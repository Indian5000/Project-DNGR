#pragma once
#include "../core/geodesic.hpp"
#include <cmath>

// Kerr metric in Boyer-Lindquist coordinates (t, r, θ, φ)
// Describes a rotating black hole with mass M and spin parameter a
// ds² = -(1 - 2Mr/Σ)dt² - (4Mar sin²θ/Σ)dtdφ + (Σ/Δ)dr² + Σdθ² 
//       + ((r²+a²)² - a²Δsin²θ)/Σ sin²θ dφ²
// where Σ = r² + a²cos²θ, Δ = r² - 2Mr + a²

struct KerrMetric : public Metric {
    double M;  // Mass parameter
    double a;  // Spin parameter (a = J/M, where J is angular momentum)
    
    explicit KerrMetric(double mass = 1.0, double spin = 0.0) 
        : M(mass), a(spin) {
        // Clamp spin to physical range |a| <= M
        if (std::abs(a) > M) a = (a > 0 ? M : -M);
    }
    
    double christoffel(int mu, int nu, int rho, const std::array<double,4>& x) const override {
        // x = (t, r, θ, φ)
        double r = x[1];
        double theta = x[2];
        
        // Avoid singularities
        if (r <= 1e-6) return 0.0;
        if (std::abs(std::sin(theta)) < 1e-10) return 0.0;
        
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        double sin2 = sin_theta * sin_theta;
        double cos2 = cos_theta * cos_theta;
        
        // Helper functions
        double Sigma = r * r + a * a * cos2;
        double Delta = r * r - 2.0 * M * r + a * a;
        
        if (Sigma < 1e-10 || std::abs(Delta) < 1e-10) return 0.0;
        
        double A = (r * r + a * a) * (r * r + a * a) - a * a * Delta * sin2;
        
        // Derivatives
        double dSigma_dr = 2.0 * r;
        double dSigma_dtheta = -2.0 * a * a * cos_theta * sin_theta;
        double dDelta_dr = 2.0 * r - 2.0 * M;
        double dA_dr = 4.0 * r * (r * r + a * a) - a * a * dDelta_dr * sin2;
        double dA_dtheta = -2.0 * a * a * Delta * sin_theta * cos_theta;
        
        // Non-zero Christoffel symbols (many terms, showing main ones)
        
        // Γ^t components
        if (mu == 0) {
            // Γ^t_{tr}
            if ((nu == 0 && rho == 1) || (nu == 1 && rho == 0)) {
                return M * (r * r - a * a * cos2) / (Sigma * Sigma * Delta);
            }
            // Γ^t_{tθ}
            if ((nu == 0 && rho == 2) || (nu == 2 && rho == 0)) {
                return -2.0 * M * a * a * r * cos_theta * sin_theta / (Sigma * Sigma);
            }
            // Γ^t_{rφ}
            if ((nu == 1 && rho == 3) || (nu == 3 && rho == 1)) {
                return -a * M * sin2 * (r * r - a * a * cos2) / (Sigma * Sigma * Delta);
            }
            // Γ^t_{θφ}
            if ((nu == 2 && rho == 3) || (nu == 3 && rho == 2)) {
                return 2.0 * M * a * r * (r * r + a * a) * sin_theta * cos_theta / (Sigma * Sigma);
            }
        }
        
        // Γ^r components
        if (mu == 1) {
            // Γ^r_{tt}
            if (nu == 0 && rho == 0) {
                return M * Delta * (r * r - a * a * cos2) / (Sigma * Sigma * Sigma);
            }
            // Γ^r_{rr}
            if (nu == 1 && rho == 1) {
                return (r * (r * r + a * a) - M * (r * r - a * a * cos2)) / (Sigma * Delta);
            }
            // Γ^r_{rθ}
            if ((nu == 1 && rho == 2) || (nu == 2 && rho == 1)) {
                return -a * a * sin_theta * cos_theta / Sigma;
            }
            // Γ^r_{θθ}
            if (nu == 2 && rho == 2) {
                return -r * Delta / Sigma;
            }
            // Γ^r_{tφ}
            if ((nu == 0 && rho == 3) || (nu == 3 && rho == 0)) {
                return -a * M * sin2 * (r * r - a * a * cos2) / (Sigma * Sigma * Sigma);
            }
            // Γ^r_{φφ}
            if (nu == 3 && rho == 3) {
                return -Delta * sin2 * (A - 2.0 * M * r * Sigma) / (2.0 * Sigma * Sigma * Sigma);
            }
        }
        
        // Γ^θ components
        if (mu == 2) {
            // Γ^θ_{tt}
            if (nu == 0 && rho == 0) {
                return -2.0 * M * a * a * r * sin_theta * cos_theta / (Sigma * Sigma * Sigma);
            }
            // Γ^θ_{rr}
            if (nu == 1 && rho == 1) {
                return a * a * sin_theta * cos_theta / (Sigma * Delta);
            }
            // Γ^θ_{rθ}
            if ((nu == 1 && rho == 2) || (nu == 2 && rho == 1)) {
                return r / Sigma;
            }
            // Γ^θ_{θθ}
            if (nu == 2 && rho == 2) {
                return -a * a * sin_theta * cos_theta / Sigma;
            }
            // Γ^θ_{tφ}
            if ((nu == 0 && rho == 3) || (nu == 3 && rho == 0)) {
                return 2.0 * M * a * r * (r * r + a * a) * sin_theta * cos_theta / (Sigma * Sigma * Sigma);
            }
            // Γ^θ_{φφ}
            if (nu == 3 && rho == 3) {
                double term1 = -sin_theta * cos_theta * (A * Sigma + 2.0 * M * r * (r * r + a * a) * (r * r + a * a));
                double term2 = 2.0 * M * a * a * r * Sigma * sin_theta * cos_theta;
                return (term1 + term2) / (Sigma * Sigma * Sigma);
            }
        }
        
        // Γ^φ components
        if (mu == 3) {
            // Γ^φ_{tr}
            if ((nu == 0 && rho == 1) || (nu == 1 && rho == 0)) {
                return M * a * (a * a * cos2 - r * r) / (Sigma * Sigma * Delta);
            }
            // Γ^φ_{tθ}
            if ((nu == 0 && rho == 2) || (nu == 2 && rho == 0)) {
                return -2.0 * M * a * r * cos_theta / (Sigma * Sigma * sin_theta);
            }
            // Γ^φ_{rφ}
            if ((nu == 1 && rho == 3) || (nu == 3 && rho == 1)) {
                return (r * A - 2.0 * M * r * Sigma) / (Sigma * Delta * A);
            }
            // Γ^φ_{θφ}
            if ((nu == 2 && rho == 3) || (nu == 3 && rho == 2)) {
                return (cos_theta * A + 2.0 * M * a * a * r * sin2 * cos_theta) / (Sigma * sin_theta * A);
            }
        }
        
        return 0.0;
    }
};