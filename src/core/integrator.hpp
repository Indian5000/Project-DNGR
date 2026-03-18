#pragma once
#include "geodesic.hpp"
#include <cmath>
#include <algorithm>

struct IntegratorConfig {
    double abs_tol = 1e-8;
    double rel_tol = 1e-6;
    double min_step = 1e-6;
    double max_step = 0.1;
};

inline double state_norm(const GeodesicState& s) {
    double max_val = 0.0;
    for (int i = 0; i < 4; ++i) {
        max_val = std::max(max_val, std::abs(s.x[i]));
        max_val = std::max(max_val, std::abs(s.p[i]));
    }
    return max_val;
}

inline GeodesicState rk45_step(
    const Metric& metric,
    const GeodesicState& s,
    double h,
    GeodesicState& error_out)
{
    static const double a21=1.0/5;
    static const double a31=3.0/40,    a32=9.0/40;
    static const double a41=44.0/45,   a42=-56.0/15,   a43=32.0/9;
    static const double a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729;
    static const double a61=9017.0/3168,  a62=-355.0/33,     a63=46732.0/5247,
                        a64=49.0/176,     a65=-5103.0/18656;
    static const double b1=35.0/384,  b3=500.0/1113, b4=125.0/192,
                        b5=-2187.0/6784, b6=11.0/84;
    static const double b1s=5179.0/57600, b3s=7571.0/16695, b4s=393.0/640,
                        b5s=-92097.0/339200, b6s=187.0/2100, b7s=1.0/40;

    GeodesicState k1 = geodesic_rhs(metric, s);
    GeodesicState k2 = geodesic_rhs(metric, {s.x+h*a21*k1.x, s.p+h*a21*k1.p});
    GeodesicState k3 = geodesic_rhs(metric, {s.x+h*(a31*k1.x+a32*k2.x),
                                              s.p+h*(a31*k1.p+a32*k2.p)});
    GeodesicState k4 = geodesic_rhs(metric, {s.x+h*(a41*k1.x+a42*k2.x+a43*k3.x),
                                              s.p+h*(a41*k1.p+a42*k2.p+a43*k3.p)});
    GeodesicState k5 = geodesic_rhs(metric, {s.x+h*(a51*k1.x+a52*k2.x+a53*k3.x+a54*k4.x),
                                              s.p+h*(a51*k1.p+a52*k2.p+a53*k3.p+a54*k4.p)});
    GeodesicState k6 = geodesic_rhs(metric, {s.x+h*(a61*k1.x+a62*k2.x+a63*k3.x+a64*k4.x+a65*k5.x),
                                              s.p+h*(a61*k1.p+a62*k2.p+a63*k3.p+a64*k4.p+a65*k5.p)});
    GeodesicState k7 = geodesic_rhs(metric, {
        s.x+h*(b1*k1.x+b3*k3.x+b4*k4.x+b5*k5.x+b6*k6.x),
        s.p+h*(b1*k1.p+b3*k3.p+b4*k4.p+b5*k5.p+b6*k6.p)});

    GeodesicState s5 = s;
    for (int i = 0; i < 4; ++i) {
        s5.x[i] += h*(b1*k1.x[i]+b3*k3.x[i]+b4*k4.x[i]+b5*k5.x[i]+b6*k6.x[i]);
        s5.p[i] += h*(b1*k1.p[i]+b3*k3.p[i]+b4*k4.p[i]+b5*k5.p[i]+b6*k6.p[i]);
    }

    GeodesicState s4 = s;
    for (int i = 0; i < 4; ++i) {
        s4.x[i] += h*(b1s*k1.x[i]+b3s*k3.x[i]+b4s*k4.x[i]+b5s*k5.x[i]+b6s*k6.x[i]+b7s*k7.x[i]);
        s4.p[i] += h*(b1s*k1.p[i]+b3s*k3.p[i]+b4s*k4.p[i]+b5s*k5.p[i]+b6s*k6.p[i]+b7s*k7.p[i]);
    }

    for (int i = 0; i < 4; ++i) {
        error_out.x[i] = s5.x[i] - s4.x[i];
        error_out.p[i] = s5.p[i] - s4.p[i];
    }
    return s5;
}

inline void integrate_geodesic_adaptive(
    const Metric& metric,
    GeodesicState& state,
    double t_end,
    IntegratorConfig cfg)
{
    double t = 0;
    double h = cfg.max_step;
    double snc = state_norm(state);

    while (t < t_end) {
        double h_try = std::min(h, t_end - t);
        GeodesicState error;
        GeodesicState next = rk45_step(metric, state, h_try, error);
        double err = state_norm(error);
        double tol = cfg.abs_tol + cfg.rel_tol * std::max(snc, 1.0);
        if (err <= tol || h_try <= cfg.min_step) {
            state = next;
            t += h_try;
            snc = state_norm(state);
        }
        double scale = 0.9 * std::pow(tol / (err + 1e-16), 0.2);
        scale = std::clamp(scale, 0.1, 5.0);
        h = std::clamp(h * scale, cfg.min_step, cfg.max_step);
    }
}