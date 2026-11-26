#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
experiment_player.py - nodo ROS2 per manovre di identificazione (PoliMOVE SimChallenge).

Pubblica su /sim_in (SimIn) leggendo /sim_out (SimOut).
IMPORTANTE: fermare control_exe prima di lanciarlo (due publisher su /sim_in si disturbano).

Uso (con il workspace "sourced", cosi' custom_msgs e' visibile):

    python3 experiment_player.py --ros-args -p maneuver:=coast -p v0:=55.0 -p v_end:=35.0 -p gear:=6

Regola: i parametri float vanno scritti con il punto decimale (5.0, non 5);
'gear' e 'repeats' sono interi.

Convenzione sterzo: i parametri 'steer', 'steer_amp' sono in "convenzione utente":
positivo = curva a sinistra (omegaz > 0). Se un comando positivo produce omegaz < 0,
lancia con  -p steer_sign:=-1.0  (si applica a tutti i comandi di sterzo).

Struttura di ogni manovra:
  fase 1 APPROACH : porta la macchina a v0 con un PI di velocita' (e tiene la direzione),
                    poi attende 'settle_time' secondi a regime. (Il cerchio salta questa fase.)
  fase 2 MANEUVER : la manovra vera e propria.
  fase 3 STOPPING : sterzo a 0, frena fino a fermarsi.
  fase 4 DONE     : mantiene i freni 1 s e chiude il nodo.
Il numero di fase viene pubblicato su /experiment/phase (Int32): registralo nel bag
insieme a /sim_in e /sim_out, serve per ritagliare le finestre di dati.

Abort automatico (-> STOPPING): max_duration, |beta| > beta_abort, vx > v_abort.
"""
import math
import time

import rclpy
from rclpy.node import Node
from rclpy.signals import SignalHandlerOptions
from rclpy.qos import QoSProfile, ReliabilityPolicy, HistoryPolicy
from std_msgs.msg import Int32

from custom_msgs.msg import SimIn, SimOut

APPROACH, MANEUVER, STOPPING, DONE = 1, 2, 3, 4
MANEUVERS = ('coast', 'accel', 'brake', 'circle', 'step', 'chirp')

DEFAULTS = {
    # --- generali ---
    'maneuver': 'coast',          # coast | accel | brake | circle | step | chirp
    'topic_in': '/sim_in',
    'topic_out': '/sim_out',
    'rate_hz': 100.0,
    'max_duration': 120.0,        # [s] tempo massimo totale (sim time), poi STOPPING
    'steer_limit': 0.35,          # [rad] saturazione comando sterzo (come baseline)
    'steer_sign': 1.0,            # +1 / -1, vedi docstring
    'v0': 30.0,                   # [m/s] velocita' di approccio / di partenza
    'v_end': 10.0,                # [m/s] fine manovra (coast, brake) / (accel: soglia superiore)
    'v_max': 60.0,                # [m/s] velocita' massima della rampa (circle)
    'v_abort': 85.0,              # [m/s] abort di sicurezza
    'beta_abort': 0.5,            # [rad] abort se |beta| supera questo valore (sbandata)
    'beta_min_speed': 2.0,        # [m/s] sotto questa vx beta non e' significativo: abort disattivato
    'use_sim_stamp': False,       # False = tempo di parete (monotonic); True = header.stamp di /sim_out
    'settle_time': 3.0,           # [s] attesa a regime a fine approccio
    # --- cambio ---
    'gear': 0,                    # 0 = automatico (soglie rpm), >0 = marcia fissa DURANTE la manovra
    'upshift_rpm': 7000.0,
    'downshift_rpm': 3000.0,
    'max_gear': 6,
    # --- PI velocita' e freni ---
    'kp': 0.1,                    # [1/(m/s)] throttle per m/s di errore
    'ki': 0.05,                   # [1/(m/s*s)]
    'brake_kp': 2.0e5,            # [Pa per m/s di eccesso di velocita']
    'brake_bias_front': 0.6,
    'max_brake_pressure': 8.0e6,  # [Pa]
    'stop_brake': 1.5e6,          # [Pa] pressione totale in STOPPING
    'hh_kp': 0.5,                 # [rad/rad] mantenimento direzione (coast/accel/brake)
    # --- accel ---
    'throttle': 1.0,              # pedale (0-1) in manovra accel
    # --- brake ---
    'brake_front_pa': 0.0,        # [Pa] pressione finale anteriore (per ruota)
    'brake_rear_pa': 0.0,         # [Pa] pressione finale posteriore (per ruota)
    'ramp_time': 3.0,             # [s] rampa lineare di pressione 0 -> valore finale
    # --- circle ---
    'steer': 0.05,                # [rad] sterzo costante (convenzione utente)
    'ramp_rate': 0.2,             # [m/s^2] rampa del riferimento di velocita' (0 = velocita' costante)
    'circle_T': 30.0,             # [s] durata se ramp_rate = 0
    # --- step / chirp ---
    'steer_amp': 0.02,            # [rad] ampiezza
    't_pre': 1.0,                 # [s] prima di ogni step
    't_hold': 3.0,                # [s] durata dello step
    't_post': 2.0,                # [s] dopo ogni step
    'repeats': 3,                 # numero di step (segno alternato)
    'f0': 0.2,                    # [Hz] chirp logaritmico
    'f1': 3.0,
    'chirp_T': 40.0,              # [s]
}


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def wrap(a):
    return math.atan2(math.sin(a), math.cos(a))


class ExperimentPlayer(Node):
    def __init__(self):
        super().__init__('experiment_player')
        for k, v in DEFAULTS.items():
            self.declare_parameter(k, v)
        self.p = {k: self.get_parameter(k).value for k in DEFAULTS}
        if self.p['maneuver'] not in MANEUVERS:
            raise ValueError('maneuver deve essere uno tra %s' % (MANEUVERS,))
        if self.p['maneuver'] == 'chirp' and not (0.0 < self.p['f0'] < self.p['f1']):
            raise ValueError('chirp: serve 0 < f0 < f1')

        self.dt = 1.0 / float(self.p['rate_hz'])
        self.s = None
        self.phase = APPROACH
        self.t_start = None
        self.phase_t0 = None
        self.h_ref = 0.0
        self.i_term = 0.0
        self.settle = 0.0
        self.gear_cmd = 1
        self.last_shift = -1e9
        self.finished = False

        qos_in = QoSProfile(history=HistoryPolicy.KEEP_LAST, depth=10,
                            reliability=ReliabilityPolicy.BEST_EFFORT)
        qos_out = QoSProfile(history=HistoryPolicy.KEEP_LAST, depth=10,
                             reliability=ReliabilityPolicy.RELIABLE)
        self.create_subscription(SimOut, self.p['topic_out'], self.on_state, qos_in)
        self.pub = self.create_publisher(SimIn, self.p['topic_in'], qos_out)
        self.pub_phase = self.create_publisher(Int32, '/experiment/phase', 10)
        self.create_timer(self.dt, self.tick)
        self.get_logger().info('Manovra "%s", parametri: %s' % (self.p['maneuver'], self.p))

    # ------------------------------------------------------------------ utilita'
    def on_state(self, msg):
        self.s = msg

    def split(self, brk):
        b = self.p['brake_bias_front']
        return b * brk, (1.0 - b) * brk

    def speed_hold(self, v_ref, v):
        """PI di velocita': ritorna (throttle 0-1, pressione freno totale [Pa])."""
        p = self.p
        e = v_ref - v
        if e < -0.5:                           # troppo veloce: frena
            self.i_term = 0.0
            return 0.0, clamp(p['brake_kp'] * (-e - 0.5), 0.0, p['max_brake_pressure'])
        if abs(e) < 5.0:                       # integra solo vicino al riferimento (anti-windup)
            self.i_term = clamp(self.i_term + p['ki'] * e * self.dt, 0.0, 0.8)
        return clamp(p['kp'] * e + self.i_term, 0.0, 1.0), 0.0

    def hold_heading(self, heading):
        """Sterzo (convenzione utente) per mantenere la direzione iniziale."""
        lim = self.p['steer_limit']
        return clamp(self.p['hh_kp'] * wrap(self.h_ref - heading), -lim, lim)

    def set_phase(self, new_phase, t_sim):
        self.phase = new_phase
        self.phase_t0 = t_sim
        self.get_logger().info('-> fase %d (t=%.2f s)' % (new_phase, t_sim - self.t_start))

    def go_stop(self, t_sim, reason):
        self.get_logger().warn('STOPPING: %s' % reason)
        self.set_phase(STOPPING, t_sim)

    def update_gear(self, s, t_sim, fixed):
        p = self.p
        if fixed and p['gear'] > 0:
            self.gear_cmd = int(p['gear'])
            return
        if t_sim - self.last_shift > 0.5:
            rpm = float(s.engine_rpm)
            if rpm >= p['upshift_rpm'] and self.gear_cmd < p['max_gear']:
                self.gear_cmd += 1
                self.last_shift = t_sim
            elif rpm <= p['downshift_rpm'] and self.gear_cmd > 1:
                self.gear_cmd -= 1
                self.last_shift = t_sim

    def publish_cmd(self, steer, thr, bf, br):
        p = self.p
        m = SimIn()
        m.steering_cmd = float(clamp(p['steer_sign'] * steer, -p['steer_limit'], p['steer_limit']))
        m.accelerator_pedal_cmd = float(clamp(thr, 0.0, 1.0))
        m.gear_cmd = float(self.gear_cmd)
        bf = float(clamp(bf, 0.0, p['max_brake_pressure']))
        br = float(clamp(br, 0.0, p['max_brake_pressure']))
        m.brake_fl_cmd = bf
        m.brake_fr_cmd = bf
        m.brake_rl_cmd = br
        m.brake_rr_cmd = br
        m.push_to_pass_cmd = 0.0
        self.pub.publish(m)
        self.pub_phase.publish(Int32(data=int(self.phase)))

    def publish_safe(self, n=20):
        """Chiamata all'uscita con Ctrl+C: sterzo 0, niente gas, freni."""
        bf, br = self.split(self.p['stop_brake'])
        try:
            for _ in range(n):
                self.publish_cmd(0.0, 0.0, bf, br)
                time.sleep(0.01)
        except Exception as e:   # contesto ROS gia' chiuso: niente traceback
            print('publish_safe interrotto: %s' % e)

    # ------------------------------------------------------------------ ciclo
    def tick(self):
        s = self.s
        if s is None:
            self.get_logger().info('In attesa di %s ...' % self.p['topic_out'],
                                   throttle_duration_sec=2.0)
            return
        p = self.p
        man = p['maneuver']
        if p['use_sim_stamp']:
            t_sim = s.stamp.sec + 1e-9 * s.stamp.nanosec
        else:
            t_sim = time.monotonic()
        v = float(s.vx)
        if self.t_start is None:
            self.t_start = t_sim
            self.phase_t0 = t_sim
            self.h_ref = float(s.heading)
            self.gear_cmd = max(1, int(round(float(s.gear_data))))
            self.get_logger().info('Primo stato ricevuto: vx=%.2f, heading=%.3f' % (v, self.h_ref))
        t = t_sim - self.t_start

        # abort di sicurezza
        if self.phase in (APPROACH, MANEUVER):
            if t > p['max_duration']:
                self.go_stop(t_sim, 'max_duration raggiunta')
            elif abs(v) > p['beta_min_speed'] and abs(float(s.beta)) > p['beta_abort']:
                self.go_stop(t_sim, 'sbandata: |beta| = %.2f' % abs(float(s.beta)))
            elif v > p['v_abort']:
                self.go_stop(t_sim, 'velocita'' eccessiva: %.1f m/s' % v)

        steer, thr, bf, br = 0.0, 0.0, 0.0, 0.0
        fixed_gear = False

        # ---- fase 1: approccio
        if self.phase == APPROACH:
            if man == 'circle':
                self.set_phase(MANEUVER, t_sim)
            else:
                thr, brk = self.speed_hold(p['v0'], v)
                bf, br = self.split(brk)
                steer = self.hold_heading(float(s.heading))
                self.settle = self.settle + self.dt if abs(v - p['v0']) < 0.5 else 0.0
                if self.settle >= p['settle_time']:
                    self.set_phase(MANEUVER, t_sim)

        # ---- fase 2: manovra
        if self.phase == MANEUVER:
            tm = t_sim - self.phase_t0
            done = False
            fixed_gear = True
            if man == 'coast':
                thr, bf, br = 0.0, 0.0, 0.0
                steer = self.hold_heading(float(s.heading))
                done = v <= p['v_end']
            elif man == 'accel':
                thr = p['throttle']
                steer = self.hold_heading(float(s.heading))
                done = v >= p['v_end']
            elif man == 'brake':
                k = clamp(tm / max(p['ramp_time'], 1e-3), 0.0, 1.0)
                thr = 0.0
                bf, br = k * p['brake_front_pa'], k * p['brake_rear_pa']
                steer = self.hold_heading(float(s.heading))
                done = v <= p['v_end']
            elif man == 'circle':
                v_ref = min(p['v0'] + p['ramp_rate'] * tm, p['v_max'])
                thr, brk = self.speed_hold(v_ref, v)
                bf, br = self.split(brk)
                steer = p['steer']
                if p['ramp_rate'] > 0.0:
                    t_ramp = (p['v_max'] - p['v0']) / p['ramp_rate']
                    done = tm >= t_ramp + 3.0
                else:
                    done = tm >= p['circle_T']
            elif man == 'step':
                thr, brk = self.speed_hold(p['v0'], v)
                bf, br = self.split(brk)
                cycle = p['t_pre'] + p['t_hold'] + p['t_post']
                k = int(tm // cycle)
                if k >= p['repeats']:
                    done = True
                else:
                    loc = tm - k * cycle
                    if p['t_pre'] <= loc < p['t_pre'] + p['t_hold']:
                        steer = p['steer_amp'] * (1.0 if k % 2 == 0 else -1.0)
            elif man == 'chirp':
                thr, brk = self.speed_hold(p['v0'], v)
                bf, br = self.split(brk)
                if tm >= p['chirp_T']:
                    done = True
                else:
                    kk = math.log(p['f1'] / p['f0'])
                    phi = 2.0 * math.pi * p['f0'] * p['chirp_T'] / kk * (math.exp(tm / p['chirp_T'] * kk) - 1.0)
                    steer = p['steer_amp'] * clamp(tm / 1.0, 0.0, 1.0) * math.sin(phi)
            if done:
                self.go_stop(t_sim, 'manovra completata')

        # ---- fase 3: arresto
        if self.phase == STOPPING:
            fixed_gear = False
            steer, thr = 0.0, 0.0
            bf, br = self.split(p['stop_brake'])
            if v < 0.3:
                self.set_phase(DONE, t_sim)

        # ---- fase 4: fine
        if self.phase == DONE:
            fixed_gear = False
            steer, thr = 0.0, 0.0
            bf, br = self.split(5.0e5)
            if t_sim - self.phase_t0 > 1.0:
                self.finished = True

        self.update_gear(s, t_sim, fixed_gear)
        self.publish_cmd(steer, thr, bf, br)
        self.get_logger().info('fase %d  t=%.1f  vx=%.1f  steer=%.3f  thr=%.2f  p_f=%.2e  p_r=%.2e  gear=%d'
                               % (self.phase, t, v, steer, thr, bf, br, self.gear_cmd),
                               throttle_duration_sec=1.0)


def main():
    # NO: lascia a Python il Ctrl+C (KeyboardInterrupt) cosi' il contesto resta valido
    # e publish_safe() riesce a mandare i freni prima di uscire.
    rclpy.init(signal_handler_options=SignalHandlerOptions.NO)
    node = ExperimentPlayer()
    try:
        while rclpy.ok() and not node.finished:
            rclpy.spin_once(node, timeout_sec=0.1)
    except KeyboardInterrupt:
        node.get_logger().warn('Interrotto: comando di arresto')
        node.publish_safe()
    finally:
        node.destroy_node()
        if rclpy.ok():
            rclpy.shutdown()


if __name__ == '__main__':
    main()
