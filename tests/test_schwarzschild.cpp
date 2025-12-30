#include "core/integrator.hpp"
#include "metrics/schwarzschild.hpp"
#include <iostream>
#include <cmath>

int main() {
    // Black hole with mass M=1 (Schwarzschild radius = 2M = 2)
    SchwarzschildMetric metric(1.0);
    
    // Light ray passing near black hole
    // Start at r=10, θ=π/2 (equatorial plane), moving tangentially
    GeodesicState state;
    double r0 = 10.0;
    double b = 6.0;  // impact parameter (closest approach distance)
    
    state.x = {0.0, r0, M_PI/2, 0.0};  // (t, r, θ, φ)
    
    // Initial velocity for light ray
    double E = 1.0;  // energy
    double L = b * E;  // angular momentum
    state.u = {E, -std::sqrt(E*E - L*L/(r0*r0) * (1.0 - 2.0/r0)), 0.0, L/(r0*r0)};
    
    IntegratorConfig cfg;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-8;
    cfg.min_step = 1e-6;
    cfg.max_step = 0.01;
    
    std::cout << "t,r,theta,phi\n";
    std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    
    // Integrate for proper time
    for(int i = 0; i < 100; ++i){
        integrate_geodesic_adaptive(metric, state, 0.5, cfg);
        std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    }
    
    return 0;
}