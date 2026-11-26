// ########################### Header Includes ###########################
#include "control/racing_controller.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace control {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kGeomSpan = 1.0;         // [m] half baseline for heading / curvature
constexpr double kBrakeDeadband = 300.0;  // [N] ~0.4 m/s^2 between pedal and brakes
constexpr double k_us = -0.000508; // [rad s^2/m] Gradiente sottosterzo stimato con linear regression

double WrapAngle(double a) {
  while (a > kPi) a -= 2.0 * kPi;
  while (a < -kPi) a += 2.0 * kPi;
  return a;
}

// ---------------------------------------------------------------------------------------------
// Engine map T_e(pedal, rpm) [Nm] identified in E7 (pedal 0 from E1).
// Holes of the measurement (low pedal at high rpm) are filled with the last measured value.
// ---------------------------------------------------------------------------------------------
constexpr int kNp = 7, kNn = 23;
constexpr double kPedal[kNp] = {0.0, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00};
constexpr double kRpm0 = 2500.0, kRpmStep = 250.0;  // rpm grid: 2500 .. 8000
// Yas Marina car: engine torque [Nm] from the E7 (constant pedal) and E1 (pedal 0) runs, 19/11/2026
constexpr double kTorque[kNn][kNp] = {
    {-22.4, 89.2, 94.0, 119.2, 132.2, 146.1, 153.8},  // 2500
    {-27.4, 80.3, 94.0, 120.2, 132.2, 146.1, 167.4},  // 2750
    {-32.5, 67.7, 91.2, 123.1, 142.9, 160.7, 198.6},  // 3000
    {-33.7, 60.7, 89.2, 127.0, 157.2, 186.2, 245.8},  // 3250
    {-34.5, 54.4, 89.0, 135.1, 175.5, 216.6, 308.6},  // 3500
    {-34.0, 49.8, 89.9, 144.7, 193.4, 245.8, 376.9},  // 3750
    {-33.2, 44.6, 88.4, 150.6, 206.2, 268.5, 421.9},  // 4000
    {-32.1, 41.9, 89.8, 155.1, 216.1, 284.6, 448.1},  // 4250
    {-31.1, 39.7, 94.1, 165.7, 231.9, 302.1, 468.3},  // 4500
    {-30.4, 35.7, 95.9, 180.6, 251.9, 323.9, 485.5},  // 4750
    {-31.2, 35.7, 95.1, 191.3, 267.3, 343.9, 502.7},  // 5000
    {-31.1, 35.7, 97.5, 198.7, 278.3, 359.2, 520.0},  // 5250
    {-30.6, 35.7, 102.8, 206.9, 288.2, 372.4, 537.2},  // 5500
    {-30.3, 35.7, 107.5, 216.6, 297.8, 385.2, 532.0},  // 5750
    {-30.1, 35.7, 112.2, 224.0, 304.1, 395.2, 525.9},  // 6000
    {-30.1, 35.7, 117.3, 230.6, 310.0, 403.0, 522.6},  // 6250
    {-30.1, 35.7, 118.2, 236.4, 315.2, 409.1, 520.6},  // 6500
    {-30.1, 35.7, 122.7, 240.4, 319.0, 413.3, 518.7},  // 6750
    {-30.1, 35.7, 122.7, 243.8, 322.2, 416.2, 516.8},  // 7000
    {-30.1, 35.7, 122.7, 247.9, 326.6, 420.5, 517.1},  // 7250
    {-30.1, 35.7, 122.7, 252.4, 331.3, 426.0, 518.0},  // 7500
    {-30.1, 35.7, 122.7, 256.9, 336.0, 431.1, 518.8},  // 7750
    {-30.1, 35.7, 122.7, 260.3, 338.6, 431.1, 518.8},  // 8000
};

// Torque row at a given rpm (linear interpolation along rpm)
void TorqueRow(double rpm, double out[kNp]) {
  if (!std::isfinite(rpm)) rpm = kRpm0;
  const double r = std::clamp((rpm - kRpm0) / kRpmStep, 0.0, double(kNn - 1));
  const int i0 = std::min(int(r), kNn - 2);
  const double w = r - i0;
  for (int j = 0; j < kNp; ++j) out[j] = (1.0 - w) * kTorque[i0][j] + w * kTorque[i0 + 1][j];
}

}  // namespace

// =============================================================================================
//                                       ENGINE MAP
// =============================================================================================
double RacingController::EngineTorque(double pedal, double rpm) const {
  double row[kNp];
  TorqueRow(rpm, row);
  pedal = std::clamp(pedal, 0.0, 1.0);
  for (int j = 0; j < kNp - 1; ++j) {
    if (pedal <= kPedal[j + 1]) {
      const double w = (pedal - kPedal[j]) / (kPedal[j + 1] - kPedal[j]);
      return (1.0 - w) * row[j] + w * row[j + 1];
    }
  }
  return row[kNp - 1];
}

// Inverse map: the torque is monotonic in the pedal at fixed rpm, so a piecewise linear inverse exists
double RacingController::PedalForTorque(double T, double rpm) const {
  double row[kNp];
  TorqueRow(rpm, row);
  if (T <= row[0]) return 0.0;
  for (int j = 0; j < kNp - 1; ++j) {
    if (T <= row[j + 1]) {
      const double w = (T - row[j]) / std::max(row[j + 1] - row[j], 1e-6);
      return kPedal[j] + w * (kPedal[j + 1] - kPedal[j]);
    }
  }
  return 1.0;
}

// =============================================================================================
//                                       VEHICLE LIMITS
// =============================================================================================
double RacingController::AyMax(double v) const {
  const double Fz0f = veh.m * veh.g * veh.b / veh.L, Fz0r = veh.m * veh.g * veh.a / veh.L;
  const double Fzf = Fz0f + 0.5 * veh.rho * veh.ClfA * v * v;
  const double Fzr = Fz0r + 0.5 * veh.rho * veh.ClrA * v * v;
  const double muf = std::max(veh.mu0_f * (1.0 + veh.mu1_f * (Fzf - Fz0f) / Fz0f), 0.5);
  const double mur = std::max(veh.mu0_r * (1.0 + veh.mu1_r * (Fzr - Fz0r) / Fz0r), 0.5);
  const double ay_f = muf * Fzf * veh.L / (veh.b * veh.m);  // front axle carries m ay b / L
  const double ay_r = mur * Fzr * veh.L / (veh.a * veh.m);
  return params.grip_scale * std::min(ay_f, ay_r);
}

double RacingController::AxBrakeMax(double v) const {
  const double Fz0f = veh.m * veh.g * veh.b / veh.L, Fz0r = veh.m * veh.g * veh.a / veh.L;
  const double Fzf = Fz0f + 0.5 * veh.rho * veh.ClfA * v * v;
  const double Fzr = Fz0r + 0.5 * veh.rho * veh.ClrA * v * v;
  const double muf = std::max(veh.mu0_f * (1.0 + veh.mu1_f * (Fzf - Fz0f) / Fz0f), 0.5);
  const double mur = std::max(veh.mu0_r * (1.0 + veh.mu1_r * (Fzr - Fz0r) / Fz0r), 0.5);
  const double tyre = params.grip_scale * (muf * Fzf + mur * Fzr);
  const double brakes = 2.0 * veh.K_brake * 1e5 * (params.brake_front_max + params.brake_rear_max) / veh.Rw;
  return (std::min(tyre, brakes) + veh.c2 * v * v) / veh.m;  // drag helps braking
}

double RacingController::AxDriveMax(double v) const {
  const double Fz0r = veh.m * veh.g * veh.a / veh.L;
  const double Fzr = Fz0r + 0.5 * veh.rho * veh.ClrA * v * v;
  const double mur = std::max(veh.mu0_r * (1.0 + veh.mu1_r * (Fzr - Fz0r) / Fz0r), 0.5);
  const double traction = params.traction_scale * mur * Fzr;
  double engine = 0.0;
  for (double G : veh.G) {
    const double rpm = v * G / veh.Rw * 30.0 / kPi;
    if (rpm > 7900.0) continue;
    engine = std::max(engine, EngineTorque(1.0, std::max(rpm, 2500.0)) * G / veh.Rw);
  }
  return (std::min(traction, engine) - veh.c2 * v * v) / veh.m;
}

// =============================================================================================
//                                    TRACK AND SPEED PROFILE
// =============================================================================================
std::size_t RacingController::Wrap(long i) const {
  const long n = static_cast<long>(x_.size());
  return static_cast<std::size_t>(((i % n) + n) % n);
}

double RacingController::Interp(const std::vector<double>& f, double s) const {
  s = std::fmod(s, length_);
  if (s < 0) s += length_;
  const double r = s / ds_;
  const std::size_t i0 = static_cast<std::size_t>(r) % f.size();
  const double w = r - std::floor(r);
  return (1.0 - w) * f[i0] + w * f[Wrap(static_cast<long>(i0) + 1)];
}

void RacingController::LoadTrack(const std::vector<double>& x, const std::vector<double>& y,
                                 const std::vector<double>& v) {
  if (x.size() < 10 || x.size() != y.size() || x.size() != v.size())
    throw std::runtime_error("RacingController: bad racing line");

  // 1) Arc length of the raw (closed) line
  const std::size_t n_raw = x.size();
  std::vector<double> s_raw(n_raw + 1, 0.0);
  for (std::size_t i = 1; i <= n_raw; ++i) {
    const std::size_t j = i % n_raw;
    s_raw[i] = s_raw[i - 1] + std::hypot(x[j] - x[i - 1], y[j] - y[i - 1]);
  }
  length_ = s_raw.back();

  // 2) Uniform resampling with step ds_
  const std::size_t n = static_cast<std::size_t>(std::floor(length_ / ds_));
  ds_ = length_ / static_cast<double>(n);
  x_.assign(n, 0); y_.assign(n, 0); v_csv_.assign(n, 0); s_.assign(n, 0);
  std::size_t k = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double si = i * ds_;
    while (k + 1 < s_raw.size() - 1 && s_raw[k + 1] < si) ++k;
    const double w = (si - s_raw[k]) / std::max(s_raw[k + 1] - s_raw[k], 1e-9);
    const std::size_t k1 = (k + 1) % n_raw;
    x_[i] = (1 - w) * x[k] + w * x[k1];
    y_[i] = (1 - w) * y[k] + w * y[k1];
    v_csv_[i] = (1 - w) * v[k] + w * v[k1];
    s_[i] = si;
  }

  // 3) Heading and curvature. The racing line has a few cm of noise: differentiating it over 0.5 m
  //    gives curvature spikes of 0.007 1/m on the straights (fake 140 m corners, which made the speed
  //    profile drop to ~50 m/s at the start of the lap). Heading is therefore taken over a 2*span
  //    baseline, curvature over the same baseline, then a 5 m moving average.
  const long span = std::max(1L, std::lround(kGeomSpan / ds_));
  psi_.assign(n, 0); kappa_.assign(n, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t ip = Wrap(long(i) + span), im = Wrap(long(i) - span);
    psi_[i] = std::atan2(y_[ip] - y_[im], x_[ip] - x_[im]);
  }
  std::vector<double> kraw(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t ip = Wrap(long(i) + span), im = Wrap(long(i) - span);
    kraw[i] = WrapAngle(psi_[ip] - psi_[im]) / (2.0 * span * ds_);
  }
  const int half = static_cast<int>(std::round(2.5 / ds_));
  for (std::size_t i = 0; i < n; ++i) {
    double acc = 0;
    for (int j = -half; j <= half; ++j) acc += kraw[Wrap(long(i) + j)];
    kappa_[i] = acc / (2 * half + 1);
  }

  idx_ = 0;
  localised_ = false;
  BuildSpeedProfile();
}

void RacingController::BuildSpeedProfile() {
  const std::size_t n = x_.size();
  v_prof_.assign(n, 0.0);
  std::vector<double> v_lat(n);
  // a) Cornering limit v^2 |kappa| <= ay_max(v), solved by fixed-point iteration
  for (std::size_t i = 0; i < n; ++i) {
    const double k = std::max(std::abs(kappa_[i]), 1e-5);
    double v = 30.0;
    for (int it = 0; it < 30; ++it) v = std::min(std::sqrt(AyMax(v) / k), params.v_max);
    v_lat[i] = params.use_csv_speed ? std::min({v, params.v_max, std::max(v_csv_[i], 5.0)})
                                    : std::min(v, params.v_max);
    v_prof_[i] = v_lat[i];
  }
  // Grip left for braking / traction once part of it is used laterally (friction ellipse)
  auto lon_share = [&](std::size_t i, double v) {
    const double ay = v * v * std::abs(kappa_[i]);
    const double r = std::min(ay / std::max(AyMax(v), 1e-3), 1.0);
    return std::sqrt(std::max(1.0 - r * r, 0.05));
  };
  // b) Backward pass (braking), two laps to close the loop
  for (int lap = 0; lap < 2; ++lap)
    for (long i = long(n) - 1; i >= 0; --i) {
      const std::size_t j = Wrap(i + 1);
      const double vj = v_prof_[j];
      const double a = AxBrakeMax(vj) * lon_share(j, vj);
      v_prof_[i] = std::min(v_prof_[i], std::sqrt(vj * vj + 2.0 * a * ds_));
    }
  // c) Forward pass (acceleration)
  for (int lap = 0; lap < 2; ++lap)
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t j = Wrap(long(i) + 1);
      const double vi = v_prof_[i];
      const double a = std::max(AxDriveMax(vi), 0.0) * lon_share(i, vi);
      v_prof_[j] = std::min(v_prof_[j], std::sqrt(vi * vi + 2.0 * a * ds_));
    }
}

// Projection of (x, y) on the line: arc length, signed lateral error (left positive), reference heading
void RacingController::Project(double x, double y, double& s, double& e_y, double& psi_ref) {
  const long n = static_cast<long>(x_.size());
  long best = static_cast<long>(idx_);
  double best_d = std::numeric_limits<double>::max();
  const long window = localised_ ? 200 : n;  // +-100 m local search once localised
  const long start = localised_ ? static_cast<long>(idx_) - window / 4 : 0;
  for (long k = 0; k < window; ++k) {
    const std::size_t i = Wrap(start + k);
    const double d = (x_[i] - x) * (x_[i] - x) + (y_[i] - y) * (y_[i] - y);
    if (d < best_d) { best_d = d; best = static_cast<long>(i); }
  }
  idx_ = static_cast<std::size_t>(best);
  localised_ = true;
  // Refine on the segment [i, i+1]
  const std::size_t i0 = idx_, i1 = Wrap(best + 1);
  const double tx = x_[i1] - x_[i0], ty = y_[i1] - y_[i0];
  const double seg = std::hypot(tx, ty);
  const double ux = tx / seg, uy = ty / seg;
  const double dx = x - x_[i0], dy = y - y_[i0];
  const double along = dx * ux + dy * uy;
  s = std::fmod(s_[i0] + along + length_, length_);
  e_y = -dx * uy + dy * ux;  // left of the line = positive
  psi_ref = std::atan2(uy, ux);
}

// =============================================================================================
//                                          GEARSHIFT
// =============================================================================================
int RacingController::SelectGear(const VehicleState& st, double dt) {
  // The simulator engages a new gear about 0.2 s after the command, so the controller keeps its own
  // commanded gear and takes a new decision only when the shift is complete and the cooldown is over.
  since_shift_ += dt;
  const int fb = std::clamp(st.gear, 1, 6);
  if (gear_cmd_ < 1) gear_cmd_ = fb;
  if (fb != gear_cmd_ || since_shift_ < params.shift_cooldown) return gear_cmd_;
  if (st.engine_rpm > params.upshift_rpm && gear_cmd_ < 6) {
    since_shift_ = 0.0;
    return ++gear_cmd_;
  }
  if (gear_cmd_ > 1) {
    // rpm the engine would have in the lower gear: shift down if it stays below the limiter
    const double rpm_low = st.engine_rpm * veh.G[gear_cmd_ - 2] / veh.G[gear_cmd_ - 1];
    if (rpm_low < 8000.0 - params.downshift_margin) {
      since_shift_ = 0.0;
      return --gear_cmd_;
    }
  }
  return gear_cmd_;
}

// =============================================================================================
//                                          MAIN STEP
// =============================================================================================
Command RacingController::Step(const VehicleState& st, double dt_hint) {
  Command cmd;
  const double vx = st.vx;
  const double vy = st.vy;

  // ---- 0) Simulated time step, from distance travelled (wall time != sim time)
  double dt = dt_hint;
  if (have_prev_ && vx > 3.0) {
    // Ground speed (vx, vy), not vx alone: the travelled distance includes the lateral component
    const double ds = std::hypot(st.x - prev_x_, st.y - prev_y_);
    const double v_g = std::hypot(vx, vy);
    dt = std::clamp(ds / v_g, 0.2 * dt_hint, 3.0 * dt_hint);
  }
  have_prev_ = true; prev_x_ = st.x; prev_y_ = st.y;

  // ---- 1) Localisation on the racing line
  double s, e_y, psi_ref;
  Project(st.x, st.y, s, e_y, psi_ref);
  const double e_psi = WrapAngle(st.heading - psi_ref);
  const double kappa_here = Interp(kappa_, s);

  // ---- 2) LATERAL: preview over the actuator delay, then feedforward + PID
  const double tau = params.preview_time;
  const double e_y_dot = vx * std::sin(e_psi) + st.vy * std::cos(e_psi);  // analytic derivative
  const double e_y_p = e_y + e_y_dot * tau;
  const double e_psi_p = e_psi + (st.yaw_rate - kappa_here * vx) * tau;
  const double e_y_dot_p = vx * std::sin(e_psi_p) + st.vy * std::cos(e_psi_p);
  const double s_p = s + vx * tau;

  const double kappa_ff = Interp(kappa_, s_p + vx * params.ff_lead_time);
  const double steer_ff = veh.L * kappa_ff + k_us * vx * vx * kappa_ff;

  const double v_eff = std::max(vx, 6.0);
  const double w = params.lat_omega, z = params.lat_zeta;
  const double gain = veh.L / (v_eff * v_eff);
  // Integral (anti-windup: clamp the contribution, freeze it at very low speed)
  // int_lat_ is stored already in steering units [rad]: the gain is applied before integrating, so a
  // change of speed (gain ~ 1/v^2) does not rescale the accumulated value and kick the steering
  if (vx > 5.0) {
    int_lat_ += gain * params.lat_ki * e_y_p * dt;
  }
  int_lat_ = std::clamp(int_lat_, -params.lat_int_limit, params.lat_int_limit);
  const double steer_fb = -gain * (w * w * e_y_p + 2.0 * z * w * e_y_dot_p) - int_lat_;
  cmd.steer = std::clamp(steer_ff + steer_fb, -params.steer_limit, params.steer_limit);

  // ---- 3) LONGITUDINAL: profile feedforward + PID on speed -> force -> pedal or brakes
  const double s_l = s + std::max(vx, 1.0) * params.lon_preview_time;
  const double v_ref = Interp(v_prof_, s_l);
  const double v_a = Interp(v_prof_, s_l - 5.0), v_b = Interp(v_prof_, s_l + 5.0);
  const double a_ff = (v_b * v_b - v_a * v_a) / 20.0;  // v dv/ds over 10 m (a 2 m baseline made a_ff chatter in corners)
  const double e_v = v_ref - vx;

  const double a_cmd = a_ff + params.lon_kp * e_v + params.lon_ki * int_lon_;
  const double a_hi = AxDriveMax(vx) + 2.0, a_lo = -AxBrakeMax(vx) - 2.0;
  if ((a_cmd < a_hi || e_v < 0) && (a_cmd > a_lo || e_v > 0)) int_lon_ += e_v * dt;  // conditional integration
  int_lon_ = std::clamp(int_lon_, -params.lon_int_limit / params.lon_ki, params.lon_int_limit / params.lon_ki);

  cmd.gear = SelectGear(st, dt);
  const int gi = std::clamp(st.gear, 1, 6) - 1;
  const double G = veh.G[gi];
  double F_des = veh.m * a_cmd + veh.c2 * vx * vx;         // wheel force asked for

  // Traction control: the rear axle cannot transmit more than what its friction ellipse leaves
  // after the lateral force needed for the current yaw rate (combined slip).
  {
    const double Fz0r = veh.m * veh.g * veh.a / veh.L;
    // load transfer from the acceleration the car can actually reach, not from a_cmd (at a standing
    // start a_cmd is tens of m/s^2 and inflated the rear load, so the cap never limited the torque)
    const double a_lt = std::clamp(a_cmd, 0.0, std::max(AxDriveMax(vx), 0.0));
    const double Fzr = Fz0r + 0.5 * veh.rho * veh.ClrA * vx * vx + veh.m * a_lt * veh.h / veh.L;
    const double mur = std::max(veh.mu0_r * (1.0 + veh.mu1_r * (Fzr - Fz0r) / Fz0r), 0.5);
    const double F_tot = params.traction_scale * mur * Fzr;
    const double Fy_r = veh.m * std::abs(vx * st.yaw_rate) * veh.a / veh.L;
    double F_x_max = std::sqrt(std::max(F_tot * F_tot - Fy_r * Fy_r, 0.0));
    // Safety: if the car is sliding (sideslip above the limit) cut the drive force
    const double beta = std::atan2(st.vy, std::max(vx, 1.0));
    if (std::abs(beta) > params.beta_cut) F_x_max = 0.0;
    // Measured-slip traction control: above tc_slip the drive force is reduced linearly, zero at
    // tc_slip + 0.10 (the Yas Marina engine has ~500 Nm and spins the rear wheels in the low gears)
    if (st.wheel_rl > 0.0 || st.wheel_rr > 0.0) {
      const double slip = (0.5 * (st.wheel_rl + st.wheel_rr) * veh.Rw - vx) / std::max(vx, 3.0);
      F_x_max *= std::clamp(1.0 - (slip - params.tc_slip) / 0.10, 0.0, 1.0);
    }
    F_des = std::min(F_des, F_x_max);
  }
  const double F_coast = veh.T_eb * G / veh.Rw;            // engine braking at pedal 0 (negative)

  // Dead band between pedal and brakes: the brakes engage only when the engine braking is short by
  // more than kBrakeDeadband, so small oscillations of F_des around F_coast do not toggle pedal/brake.
  if (F_des >= F_coast - kBrakeDeadband) {
    const double rpm = std::max(st.engine_rpm, 900.0);
    cmd.pedal = std::clamp(PedalForTorque(F_des * veh.Rw / G, rpm), 0.0, 1.0);
  } else {
    // Brakes: total force split with the bias, then each axle limited by its friction ellipse
    // (longitudinal force left after the lateral force the corner needs), so the wheels do not lock.
    const double F_b = F_coast - F_des;  // what the brakes must add [N]
    // lateral demand: what the car does now, or what the line asks within the next 0.5 s
    // (planned lateral acceleration v_prof^2 |kappa|: the car will be at the profile speed there)
    double ay_ahead = 0.0;
    for (double ds_a = 0.0; ds_a <= std::max(vx, 1.0) * 0.5; ds_a += 2.0) {
      const double vp = Interp(v_prof_, s + ds_a);
      ay_ahead = std::max(ay_ahead, vp * vp * std::abs(Interp(kappa_, s + ds_a)));
    }
    const double ay = std::max(std::abs(vx * st.yaw_rate), ay_ahead);
    // Patched ax_est avoiding it goes down -22 m/s^2
    const double ax_est = std::max(std::min(a_cmd, 0.0), -AxBrakeMax(vx));
    const double Fz0f = veh.m * veh.g * veh.b / veh.L, Fz0r = veh.m * veh.g * veh.a / veh.L;
    const double Fzf = Fz0f + 0.5 * veh.rho * veh.ClfA * vx * vx - veh.m * ax_est * veh.h / veh.L;
    const double Fzr = std::max(Fz0r + 0.5 * veh.rho * veh.ClrA * vx * vx + veh.m * ax_est * veh.h / veh.L, 500.0);
    const double muf = std::max(veh.mu0_f * (1.0 + veh.mu1_f * (Fzf - Fz0f) / Fz0f), 0.5);
    const double mur = std::max(veh.mu0_r * (1.0 + veh.mu1_r * (Fzr - Fz0r) / Fz0r), 0.5);
    const double Ff_tot = params.abs_scale * muf * Fzf, Fr_tot = params.abs_scale * mur * Fzr;
    const double Fyf = veh.m * ay * veh.b / veh.L, Fyr = veh.m * ay * veh.a / veh.L;
    const double Ff_max = std::max(std::sqrt(std::max(Ff_tot * Ff_tot - Fyf * Fyf, 0.0)), 0.3 * Ff_tot);
    const double Fr_max = std::max(std::sqrt(std::max(Fr_tot * Fr_tot - Fyr * Fyr, 0.0)) + F_coast, 0.0);  // engine brake uses rear grip
    // Steering priority: when the steering is close to its limit, the front tyres are needed for turning
    const double steer_use = std::abs(cmd.steer) / params.steer_limit;
    const double k_steer = std::clamp((1.0 - steer_use) / 0.4, 0.2, 1.0);
    double Ff = std::min(params.brake_bias_front * F_b, k_steer * Ff_max);
    double Fr = std::min(F_b - Ff, Fr_max);
    Ff = std::min(F_b - Fr, k_steer * Ff_max);              // give back to the front what the rear cannot take
    const double bar_per_N = veh.Rw / (2.0 * veh.K_brake) / 1e5;  // per-wheel pressure for 1 N on an axle
    cmd.brake_front = std::min(Ff * bar_per_N, params.brake_front_max);
    cmd.brake_rear = std::min(Fr * bar_per_N, params.brake_rear_max);
  }

  diag_ = Diagnostics{s, e_y, e_psi, kappa_here, v_ref, a_cmd, dt, steer_ff, steer_fb};
  return cmd;
}

}  // namespace control
