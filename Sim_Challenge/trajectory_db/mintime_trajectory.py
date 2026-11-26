#!/usr/bin/env python3
"""
Minimum-lap-time racing line  (Cars4 - PoliMOVE SimChallenge26)
=================================================================

Computes the racing line that MINIMISES LAP TIME (not curvature) inside the track
limits, respecting the physical limits of the identified vehicle model.

Method (optimal control in curvilinear coordinates, direct transcription):
  * reference = smoothed centre line of the track, abscissa s
  * states   : n  = lateral offset from the centre line [m]  (positive = left)
               xi = heading relative to the centre line [rad]
               v  = speed [m/s]
  * controls : k  = curvature of the racing line [1/m]
               ax = net longitudinal acceleration [m/s^2]
  * dynamics in s (Frenet):
       sf      = (1 - n*kc) / cos(xi)          (dt/ds = sf / v)
       dn/ds   = (1 - n*kc) * tan(xi)
       dxi/ds  = sf*k - kc
       dv/ds   = sf*ax / v
  * cost     : lap time  T = sum( sf/v * ds )  + small smoothness penalty
  * constraints (quasi-steady-state point-mass model, from the model identification):
       tyre accelerations  at = ax + c2*v^2/m   (drag helps braking, costs when accelerating)
                           ay = v^2 * k
       friction ellipse    (at/ax_lim(v))^2 + (ay/ay_lim(v))^2 <= 1
       grip vs speed       ay_lim(v) = grip*(ay0 + ay2*v^2)   (downforce)
       drive limit         at <= min(P_max/(m v), grip*(tr0 + tr2*v^2))   (power, rear traction)
       steering            |k| <= tan(steer_max)/L
       track limits        n_right + margin <= n <= n_left - margin
       closed lap          state(N) = state(0)

Usage (defaults = SimChallenge26 folders and file format):
  pip install casadi numpy scipy
  python3 mintime_trajectory.py                       # borders from trajectory_db/track_lines,
                                                      # result -> trajectory_db/traj_mintime.csv
  python3 mintime_trajectory.py --out monza_m15.csv --margin 1.5
Folders: in the ROS container borders come from ~/sim_challenge/trajectory_db/track_lines and the
result goes to ~/sim_challenge/trajectory_db. Anywhere else (VS Code, Windows) put the two border CSVs
next to this script: the result is written in the same folder. An absolute --out path is used as is.

For a new track: change --inner/--outer (and --rows if the CSV layout is different).
For a new car: edit VEHICLE below (or pass --vehicle file.json with the same keys).
Output: CSV with 3 rows (x, y, v), closed loop, starting at the first point of the track lines,
same format as optimized_racing_line_with_velocity.csv.
"""
import argparse, json, os, sys
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

# ------------------------------------------------------------------------------------------------
# Vehicle limits (identified model, see the model-identification notes)
# ------------------------------------------------------------------------------------------------
VEHICLE = dict(
    # Yas Marina car (19/11/2026): simulator constants + E1/E3/E7 runs on the Yas simulator.
    # Monza values in brackets.
    m=720.0,            # mass [kg]                                                (unchanged)
    L=3.115,            # wheelbase [m]                                            (unchanged)
    c2=0.7301,          # aerodynamic drag  F = c2 v^2 [N/(m/s)^2]                 (0.7313)
    P_max=350e3,        # effective max power at the wheels [W]: full-load torque  (260e3)
                        #   ~520 Nm above 5500 rpm, mean 368 kW in the shift band, minus shift losses
    ay0=16.3,           # lateral limit at v=0 [m/s^2]                             (13.9)
    ay2=0.0034,         # lateral limit increase per v^2 [1/m]                     (0.00264)
    tr0=9.20,           # rear-axle traction limit at v=0 [m/s^2]                  (8.76)
    tr2=0.00192,        # traction limit increase per v^2 [1/m]                    (0.0023)
                        #   same tyre/aero model as the controller (AyMax, rear axle limiting):
                        #   Magic Formula of the simulator with fs = 1 -> mu0 1.80 / 1.66, mu1 -0.10 / -0.11,
                        #   1/2 rho ClA 1.2482 / 1.0633; ay(0, 20, 40, 60 m/s) = 16.3, 17.8, 22.1, 28.5
    ax_over_ay=1.0,     # longitudinal/lateral grip ratio of the friction ellipse
    steer_max=0.32,     # max steering angle [rad]
    half_width=0.95,    # half width of the car [m] (used only to report clearance)
)


def load_xy(path, rows):
    """Reads x,y of a track border. rows='1,2' -> rows 1 and 2 of a row-wise CSV (SimChallenge);
    rows='cols' -> two columns x,y (one point per line)."""
    a = np.loadtxt(os.path.expanduser(path), delimiter=',')
    if rows == 'cols':
        return a[:, 0], a[:, 1]
    r = [int(i) for i in rows.split(',')]
    return a[r[0]], a[r[1]]


def lowpass(z, ds, lam):
    """Periodic FFT low-pass: removes wavelengths shorter than lam [m] (border noise)."""
    Z = np.fft.rfft(z)
    Z[np.fft.rfftfreq(len(z), d=ds) > 1.0 / lam] = 0
    return np.fft.irfft(Z, len(z))


def build_reference(xi, yi, xo, yo, ds, lam):
    """Smoothed centre line, left normal, curvature and available lateral space (left/right)."""
    _, j = cKDTree(np.c_[xo, yo]).query(np.c_[xi, yi])          # pair each inner point with outer
    cx, cy = (xi + xo[j]) / 2, (yi + yo[j]) / 2
    # uniform resampling of the (closed) centre line
    seg = np.hypot(np.diff(np.r_[cx, cx[0]]), np.diff(np.r_[cy, cy[0]]))
    s_raw = np.r_[0, np.cumsum(seg)]
    Ltr = s_raw[-1]
    N = int(round(Ltr / ds)); ds = Ltr / N
    s = np.arange(N) * ds
    cx = np.interp(s, s_raw, np.r_[cx, cx[0]]); cy = np.interp(s, s_raw, np.r_[cy, cy[0]])
    cx, cy = lowpass(cx, ds, lam), lowpass(cy, ds, lam)
    sp_x = CubicSpline(np.r_[s, Ltr], np.r_[cx, cx[0]], bc_type='periodic')
    sp_y = CubicSpline(np.r_[s, Ltr], np.r_[cy, cy[0]], bc_type='periodic')
    dx, dy, ddx, ddy = sp_x(s, 1), sp_y(s, 1), sp_x(s, 2), sp_y(s, 2)
    nrm = np.hypot(dx, dy)
    tx, ty = dx / nrm, dy / nrm
    nx, ny = -ty, tx                                             # left normal
    kc = (dx * ddy - dy * ddx) / nrm**3
    # signed distance of the two borders along the normal
    def side(bx, by):
        _, k = cKDTree(np.c_[bx, by]).query(np.c_[cx, cy])
        return (bx[k] - cx) * nx + (by[k] - cy) * ny
    di, do = side(xi, yi), side(xo, yo)
    left, right = (di, do) if np.median(di) > 0 else (do, di)   # works for both driving directions
    return dict(s=s, ds=ds, L=Ltr, cx=cx, cy=cy, nx=nx, ny=ny, kc=kc, left=left, right=right)


def solve(ref, veh, grip, margin, w_smooth, max_iter):
    import casadi as ca
    N, ds, kc = len(ref['s']), ref['ds'], ref['kc']
    m, c2 = veh['m'], veh['c2']
    k_max = np.tan(veh['steer_max']) / veh['L']
    lo = ref['right'] + margin
    hi = ref['left'] - margin
    if np.any(lo >= hi):
        sys.exit('margin too large: track narrower than 2*margin somewhere')

    opti = ca.Opti()
    n, xi, v = opti.variable(N), opti.variable(N), opti.variable(N)
    k, ax = opti.variable(N), opti.variable(N)

    def f(i):  # dynamics d(state)/ds at node i
        sf = (1 - n[i] * kc[i]) / ca.cos(xi[i])
        return ca.vertcat((1 - n[i] * kc[i]) * ca.tan(xi[i]), sf * k[i] - kc[i], sf * ax[i] / v[i]), sf / v[i]

    T = 0
    F = [f(i) for i in range(N)]
    for i in range(N):                          # trapezoidal collocation, closed loop (i+1 wraps)
        j = (i + 1) % N
        X_i, X_j = ca.vertcat(n[i], xi[i], v[i]), ca.vertcat(n[j], xi[j], v[j])
        opti.subject_to(X_j - X_i == ds / 2 * (F[i][0] + F[j][0]))
        T += ds / 2 * (F[i][1] + F[j][1])

    ay = v**2 * k
    at = ax + c2 * v**2 / m
    ay_lim = grip * (veh['ay0'] + veh['ay2'] * v**2)
    ax_lim = veh['ax_over_ay'] * ay_lim
    opti.subject_to((at / ax_lim)**2 + (ay / ay_lim)**2 <= 1)
    opti.subject_to(at <= veh['P_max'] / (m * v))
    opti.subject_to(at <= grip * (veh['tr0'] + veh['tr2'] * v**2))
    opti.subject_to(opti.bounded(lo, n, hi))
    opti.subject_to(opti.bounded(-k_max, k, k_max))
    opti.subject_to(opti.bounded(-0.6, xi, 0.6))
    opti.subject_to(opti.bounded(3.0, v, 100.0))

    dk = k[1:] - k[:-1]; dax = ax[1:] - ax[:-1]
    opti.minimize(T + w_smooth * (ca.sumsqr(dk) / ds * 1e3 + ca.sumsqr(dax) / ds * 1e-3))

    # initial guess: centre line at moderate speed
    opti.set_initial(n, 0); opti.set_initial(xi, 0); opti.set_initial(v, 30.0)
    opti.set_initial(k, kc); opti.set_initial(ax, 0)
    opti.solver('ipopt', {'expand': True, 'print_time': False},
                {'max_iter': max_iter, 'print_level': 5, 'tol': 1e-6, 'mu_strategy': 'adaptive'})
    try:
        sol = opti.solve()
    except RuntimeError:
        print('WARNING: IPOPT did not converge, writing the last iterate', file=sys.stderr)
        sol = opti.debug
    return sol.value(n), sol.value(v), float(sol.value(T))


def write_line(ref, n, v, out, step, lam):
    """Racing line = centre + n*normal, resampled with a periodic spline every `step` metres.
    n is low-pass filtered first (wavelengths < lam metres removed): the optimiser solution has small
    node-to-node wiggles that become curvature spikes, i.e. steering jolts, in the controller."""
    if lam > 0:
        n = lowpass(n, ref['ds'], lam)
    x = ref['cx'] + n * ref['nx']; y = ref['cy'] + n * ref['ny']
    seg = np.hypot(np.diff(np.r_[x, x[0]]), np.diff(np.r_[y, y[0]]))
    s = np.r_[0, np.cumsum(seg)]; Lr = s[-1]
    sx = CubicSpline(s, np.r_[x, x[0]], bc_type='periodic')
    sy = CubicSpline(s, np.r_[y, y[0]], bc_type='periodic')
    sv = CubicSpline(s, np.r_[v, v[0]], bc_type='periodic')
    q = np.arange(0, Lr, step)
    np.savetxt(os.path.expanduser(out), np.vstack([sx(q), sy(q), sv(q)]), delimiter=',', fmt='%.6f')
    hd = np.unwrap(np.arctan2(np.gradient(sy(q)), np.gradient(sx(q))))
    return Lr, n, float(np.max(np.abs(np.gradient(hd, q))))


# Folders: inside the ROS container borders are read from / results written to trajectory_db
# (where the launch file looks). Anywhere else (e.g. VS Code on Windows) everything happens in the
# folder of this script: put inner_track_line.csv and outer_track_line.csv next to it.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_ROS_DB = os.path.expanduser('~/sim_challenge/trajectory_db')
if os.path.isdir(_ROS_DB):
    TRAJ_DB, TRACK_DIR = _ROS_DB, os.path.join(_ROS_DB, 'track_lines')
else:
    TRAJ_DB = TRACK_DIR = SCRIPT_DIR


def main():
    p = argparse.ArgumentParser(description='Minimum-lap-time racing line')
    p.add_argument('--inner', default=os.path.join(TRACK_DIR, 'inner_track_line.csv'), help='CSV of one track border')
    p.add_argument('--outer', default=os.path.join(TRACK_DIR, 'outer_track_line.csv'), help='CSV of the other track border')
    p.add_argument('--rows', default='1,2', help="rows holding x,y ('1,2' SimChallenge) or 'cols'")
    p.add_argument('--out', default='traj_mintime.csv', help='file name (saved in trajectory_db, or next to the script) or absolute path')
    p.add_argument('--vehicle', default=None, help='optional JSON overriding VEHICLE keys')
    p.add_argument('--grip', type=float, default=0.80, help='grip scale (same meaning as the controller)')
    p.add_argument('--margin', type=float, default=1.5, help='min distance CoG-border [m]')
    p.add_argument('--ds', type=float, default=4.0, help='discretisation step [m]')
    p.add_argument('--smooth_lambda', type=float, default=25.0, help='centre-line low-pass wavelength [m]')
    p.add_argument('--w_smooth', type=float, default=1.0, help='weight of the smoothness penalty')
    p.add_argument('--max_iter', type=int, default=3000)
    p.add_argument('--line_lambda', type=float, default=30.0, help='low-pass on the line offset [m] (0 = off)')
    a = p.parse_args()
    if not os.path.isabs(os.path.expanduser(a.out)):
        os.makedirs(TRAJ_DB, exist_ok=True)
        a.out = os.path.join(TRAJ_DB, a.out)
    for f in (a.inner, a.outer):
        if not os.path.isfile(os.path.expanduser(f)):
            sys.exit(f'File not found: {f}\nPut inner_track_line.csv and outer_track_line.csv in {TRACK_DIR} '
                     f'or pass --inner / --outer')
    print(f'borders: {a.inner} | {a.outer}\noutput : {a.out}')

    veh = dict(VEHICLE)
    if a.vehicle:
        veh.update(json.load(open(os.path.expanduser(a.vehicle))))
    xi, yi = load_xy(a.inner, a.rows)
    xo, yo = load_xy(a.outer, a.rows)
    ref = build_reference(xi, yi, xo, yo, a.ds, a.smooth_lambda)
    print(f'track length {ref["L"]:.1f} m, {len(ref["s"])} nodes, '
          f'width {np.min(ref["left"]-ref["right"]):.1f}-{np.max(ref["left"]-ref["right"]):.1f} m')
    n, v, T = solve(ref, veh, a.grip, a.margin, a.w_smooth, a.max_iter)
    Lr, n, kmax = write_line(ref, n, v, a.out, 0.5, a.line_lambda)
    clear = np.minimum(ref['left'] - n, n - ref['right']) - veh['half_width']
    print(f'lap time (model) {T:.2f} s | line length {Lr:.1f} m | v min/max {v.min():.1f}/{v.max():.1f} m/s'
          f' | min wheel-border clearance {clear.min():.2f} m | max curvature {kmax:.4f} 1/m')
    print(f'written {a.out}')


if __name__ == '__main__':
    main()
