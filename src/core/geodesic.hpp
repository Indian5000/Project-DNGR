#pragma once
#include <array>
#include <cmath>

struct GeodesicState {
    std::array<double, 4> x;  // position (t, r, θ, φ)
    std::array<double, 4> p;  // covariant momentum p_μ
};

inline std::array<double, 4> operator+(const std::array<double, 4>& a, const std::array<double, 4>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]};
}
inline std::array<double, 4> operator*(double s, const std::array<double, 4>& a) {
    return {s*a[0], s*a[1], s*a[2], s*a[3]};
}

struct Metric {
    virtual ~Metric() = default;

    // Inverse metric g^{μν} at position x
    virtual double g_inv(int mu, int nu, const std::array<double,4>& x) const = 0;

    // Numerical derivative ∂g^{μν}/∂x^α  (central differences)
    double dg_inv(int mu, int nu, int alpha, const std::array<double,4>& x) const {
        const double h = 1e-5;
        auto xp = x; xp[alpha] += h;
        auto xm = x; xm[alpha] -= h;
        return (g_inv(mu, nu, xp) - g_inv(mu, nu, xm)) / (2.0 * h);
    }
};

// Super-Hamiltonian formulation:  H = (1/2) g^{μν} p_μ p_ν
// dx^μ/dζ =  g^{μν} p_ν
// dp_μ/dζ = -(1/2) ∂g^{αβ}/∂x^μ  p_α p_β
// p_t and p_φ are exactly conserved (no coordinate dependence on t or φ)
inline GeodesicState geodesic_rhs(const Metric& metric, const GeodesicState& s) {
    GeodesicState ds;

    // dx^μ/dζ = g^{μν} p_ν
    for (int mu = 0; mu < 4; ++mu) {
        ds.x[mu] = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            ds.x[mu] += metric.g_inv(mu, nu, s.x) * s.p[nu];
    }

    // dp_μ/dζ = -(1/2) ∂g^{αβ}/∂x^μ p_α p_β
    for (int mu = 0; mu < 4; ++mu) {
        ds.p[mu] = 0.0;
        for (int alpha = 0; alpha < 4; ++alpha)
            for (int beta = 0; beta < 4; ++beta)
                ds.p[mu] -= 0.5 * metric.dg_inv(alpha, beta, mu, s.x)
                            * s.p[alpha] * s.p[beta];
    }

    return ds;
}