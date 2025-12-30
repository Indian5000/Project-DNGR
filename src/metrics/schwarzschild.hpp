#pragma once
#include "../core/geodesic.hpp"
#include <cmath>

// Schwarzschild metric in Schwarzschild coordinates (t, r, θ, φ)
// ds² = -(1-2M/r)dt² + (1-2M/r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
struct SchwarzschildMetric : public Metric {
    double M;  // Mass parameter (geometric units, G=c=1)
    
    explicit SchwarzschildMetric(double mass = 1.0) : M(mass) {}
    
    double christoffel(int mu, int nu, int rho, const std::array<double,4>& x) const override {
        // x = (t, r, θ, φ)
        double r = x[1];
        double theta = x[2];
        
        // Avoid singularities
        if (r <= 2.0 * M + 1e-6) return 0.0;
        if (std::abs(std::sin(theta)) < 1e-10) return 0.0;
        
        double f = 1.0 - 2.0 * M / r;
        double df_dr = 2.0 * M / (r * r);
        
        // Non-zero Christoffel symbols for Schwarzschild metric
        // Γ^t_{tr} = Γ^t_{rt}
        if (mu == 0 && ((nu == 0 && rho == 1) || (nu == 1 && rho == 0))) {
            return M / (r * r * f);
        }
        
        // Γ^r_{tt}
        if (mu == 1 && nu == 0 && rho == 0) {
            return M * f / (r * r);
        }
        
        // Γ^r_{rr}
        if (mu == 1 && nu == 1 && rho == 1) {
            return -M / (r * r * f);
        }
        
        // Γ^r_{θθ}
        if (mu == 1 && nu == 2 && rho == 2) {
            return -r * f;
        }
        
        // Γ^r_{φφ}
        if (mu == 1 && nu == 3 && rho == 3) {
            return -r * f * std::sin(theta) * std::sin(theta);
        }
        
        // Γ^θ_{rθ} = Γ^θ_{θr}
        if (mu == 2 && ((nu == 1 && rho == 2) || (nu == 2 && rho == 1))) {
            return 1.0 / r;
        }
        
        // Γ^θ_{φφ}
        if (mu == 2 && nu == 3 && rho == 3) {
            return -std::sin(theta) * std::cos(theta);
        }
        
        // Γ^φ_{rφ} = Γ^φ_{φr}
        if (mu == 3 && ((nu == 1 && rho == 3) || (nu == 3 && rho == 1))) {
            return 1.0 / r;
        }
        
        // Γ^φ_{θφ} = Γ^φ_{φθ}
        if (mu == 3 && ((nu == 2 && rho == 3) || (nu == 3 && rho == 2))) {
            return std::cos(theta) / std::sin(theta);
        }
        
        return 0.0;
    }
};