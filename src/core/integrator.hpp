#pragma once
#include "geodesic.hpp"
#include <cmath>
#include <algorithm>
#include <array> // if GeodesicState uses std::array

struct IntegratorConfig {
    double abs_tol = 1e-8;
    double rel_tol = 1e-6;
    double min_step = 1e-6;
    double max_step = 0.1;
};

// Compute infinity norm of a GeodesicState
inline double state_norm(const GeodesicState& s) {
    double max_val = 0.0;
    for(int i=0;i<4;++i){
        max_val = std::max(max_val, std::abs(s.x[i]));
        max_val = std::max(max_val, std::abs(s.u[i]));
    }
    return max_val;
}

// Single RK45 step (Dormand-Prince)
inline GeodesicState rk45_step(
    const Metric& metric,
    const GeodesicState& s,
    double h,
    GeodesicState& error_out
){
    // Dormand-Prince coefficients
    static const double a21=1.0/5;
    static const double a31=3.0/40, a32=9.0/40;
    static const double a41=44.0/45, a42=-56.0/15, a43=32.0/9;
    static const double a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729;
    static const double a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656;

    static const double b1=35.0/384, b3=500.0/1113, b4=125.0/192, b5=-2187.0/6784, b6=11.0/84;
    static const double b1s=5179.0/57600, b3s=7571.0/16695, b4s=393.0/640, b5s=-92097.0/339200, b6s=187.0/2100, b7s=1.0/40;

    // Compute K values
    GeodesicState k1 = geodesic_rhs(metric, s);
    GeodesicState k2 = geodesic_rhs(metric, GeodesicState{s.x + h*a21*k1.x, s.u + h*a21*k1.u});
    GeodesicState k3 = geodesic_rhs(metric, GeodesicState{s.x + h*(a31*k1.x + a32*k2.x), s.u + h*(a31*k1.u + a32*k2.u)});
    GeodesicState k4 = geodesic_rhs(metric, GeodesicState{s.x + h*(a41*k1.x + a42*k2.x + a43*k3.x),
                                                           s.u + h*(a41*k1.u + a42*k2.u + a43*k3.u)});
    GeodesicState k5 = geodesic_rhs(metric, GeodesicState{s.x + h*(a51*k1.x + a52*k2.x + a53*k3.x + a54*k4.x),
                                                           s.u + h*(a51*k1.u + a52*k2.u + a53*k3.u + a54*k4.u)});
    GeodesicState k6 = geodesic_rhs(metric, GeodesicState{s.x + h*(a61*k1.x + a62*k2.x + a63*k3.x + a64*k4.x + a65*k5.x),
                                                           s.u + h*(a61*k1.u + a62*k2.u + a63*k3.u + a64*k4.u + a65*k5.u)});
    // k7 for 4th-order
    GeodesicState k7 = geodesic_rhs(metric, GeodesicState{
        s.x + h*(b1*k1.x + b3*k3.x + b4*k4.x + b5*k5.x + b6*k6.x),
        s.u + h*(b1*k1.u + b3*k3.u + b4*k4.u + b5*k5.u + b6*k6.u)
    });

    // 5th-order solution
    GeodesicState s5 = s;
    for(int i=0;i<4;++i){
        s5.x[i] += h*(b1*k1.x[i] + b3*k3.x[i] + b4*k4.x[i] + b5*k5.x[i] + b6*k6.x[i]);
        s5.u[i] += h*(b1*k1.u[i] + b3*k3.u[i] + b4*k4.u[i] + b5*k5.u[i] + b6*k6.u[i]);
    }

    // 4th-order solution
    GeodesicState s4 = s;
    for(int i=0;i<4;++i){
        s4.x[i] += h*(b1s*k1.x[i] + b3s*k3.x[i] + b4s*k4.x[i] + b5s*k5.x[i] + b6s*k6.x[i] + b7s*k7.x[i]);
        s4.u[i] += h*(b1s*k1.u[i] + b3s*k3.u[i] + b4s*k4.u[i] + b5s*k5.u[i] + b6s*k6.u[i] + b7s*k7.u[i]);
    }

    // Error estimate
    for(int i=0;i<4;++i){
        error_out.x[i] = s5.x[i] - s4.x[i];
        error_out.u[i] = s5.u[i] - s4.u[i];
    }

    return s5;
}

// --------------------------
// Adaptive RK45 integration loop
// --------------------------
inline void integrate_geodesic_adaptive(
    const Metric& metric,
    GeodesicState& state,
    double t_end,
    IntegratorConfig cfg
){
    double t = 0;
    double h = cfg.max_step;
    double state_norm_cache = state_norm(state);

    while(t < t_end){
        double h_try = std::min(h, t_end - t);

        GeodesicState error;
        GeodesicState next = rk45_step(metric, state, h_try, error);

        double err = state_norm(error);
        double tol = cfg.abs_tol + cfg.rel_tol * std::max(state_norm_cache, 1.0);

        if(err <= tol || h_try <= cfg.min_step){
            // Accept step
            state = next;
            t += h_try;
            state_norm_cache = state_norm(state); // update cache
        }

        // Adaptive step size (only reduces h if necessary)
        double scale = 0.9 * std::pow(tol / (err + 1e-16), 0.2);
        scale = std::clamp(scale, 0.1, 5.0);
        h = std::clamp(h*scale, cfg.min_step, cfg.max_step);
    }
}