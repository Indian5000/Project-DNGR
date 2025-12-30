#pragma once
#include "../core/geodesic.hpp"
#include <cmath>

// Morris-Thorne traversable wormhole in "canonical" coordinates
// ds² = -dt² + dl² + (b₀² + l²)(dθ² + sin²θ dφ²)
// where l is the "proper radial distance" and b₀ is the throat radius
// This is a simple traversable wormhole with no horizons

struct MorrisThorne : public Metric {
    double b0;  // Throat radius
    
    explicit MorrisThorne(double throat_radius = 1.0) : b0(throat_radius) {}
    
    // Helper: radial function r(l) = sqrt(b0² + l²)
    inline double r_of_l(double l) const {
        return std::sqrt(b0 * b0 + l * l);
    }
    
    // Helper: dr/dl = l / sqrt(b0² + l²)
    inline double dr_dl(double l) const {
        double r = r_of_l(l);
        return (r > 1e-10) ? l / r : 0.0;
    }
    
    double christoffel(int mu, int nu, int rho, const std::array<double,4>& x) const override {
        // x = (t, l, θ, φ) where l is proper radial distance
        double l = x[1];
        double theta = x[2];
        
        // Avoid singularities
        if (std::abs(std::sin(theta)) < 1e-10) return 0.0;
        
        double r = r_of_l(l);
        double dr_dl_val = dr_dl(l);
        
        if (r < 1e-10) return 0.0;
        
        // Non-zero Christoffel symbols for Morris-Thorne wormhole
        
        // Γ^l_{θθ}
        if (mu == 1 && nu == 2 && rho == 2) {
            return -l;
        }
        
        // Γ^l_{φφ}
        if (mu == 1 && nu == 3 && rho == 3) {
            return -l * std::sin(theta) * std::sin(theta);
        }
        
        // Γ^θ_{lθ} = Γ^θ_{θl}
        if (mu == 2 && ((nu == 1 && rho == 2) || (nu == 2 && rho == 1))) {
            return l / (b0 * b0 + l * l);
        }
        
        // Γ^θ_{φφ}
        if (mu == 2 && nu == 3 && rho == 3) {
            return -std::sin(theta) * std::cos(theta);
        }
        
        // Γ^φ_{lφ} = Γ^φ_{φl}
        if (mu == 3 && ((nu == 1 && rho == 3) || (nu == 3 && rho == 1))) {
            return l / (b0 * b0 + l * l);
        }
        
        // Γ^φ_{θφ} = Γ^φ_{φθ}
        if (mu == 3 && ((nu == 2 && rho == 3) || (nu == 3 && rho == 2))) {
            return std::cos(theta) / std::sin(theta);
        }
        
        return 0.0;
    }
};