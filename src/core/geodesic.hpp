#pragma once
#include <array>

// State: position x^μ and 4-velocity u^μ
struct GeodesicState {
    std::array<double, 4> x;  // position (t, x, y, z)
    std::array<double, 4> u;  // 4-velocity
};

// Overload operators for arithmetic
inline std::array<double, 4> operator+(const std::array<double, 4>& a, const std::array<double, 4>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]};
}

inline std::array<double, 4> operator*(double s, const std::array<double, 4>& a) {
    return {s*a[0], s*a[1], s*a[2], s*a[3]};
}

// Abstract metric interface
struct Metric {
    virtual ~Metric() = default;
    // Compute Christoffel symbols Γ^μ_νρ at position x
    virtual double christoffel(int mu, int nu, int rho, const std::array<double,4>& x) const = 0;
};

// Geodesic equation: dx^μ/dλ = u^μ, du^μ/dλ = -Γ^μ_νρ u^ν u^ρ
inline GeodesicState geodesic_rhs(const Metric& metric, const GeodesicState& s) {
    GeodesicState ds;
    ds.x = s.u;  // dx/dλ = u
    
    // du^μ/dλ = -Γ^μ_νρ u^ν u^ρ
    for(int mu = 0; mu < 4; ++mu) {
        ds.u[mu] = 0.0;
        for(int nu = 0; nu < 4; ++nu) {
            for(int rho = 0; rho < 4; ++rho) {
                ds.u[mu] -= metric.christoffel(mu, nu, rho, s.x) * s.u[nu] * s.u[rho];
            }
        }
    }
    return ds;
}