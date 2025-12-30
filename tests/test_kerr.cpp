#include "core/integrator.hpp"
#include "metrics/kerr.hpp"
#include <iostream>
#include <cmath>

int main() {
    // Kerr black hole: M=1, spin a=0.5 (moderate rotation)
    KerrMetric metric(1.0, 0.5);
    
    // Light ray deflection test (safe distance)
    GeodesicState state;
    double r0 = 15.0;  // Start farther out
    double b = 10.0;   // Impact parameter (closest approach)
    
    state.x = {0.0, r0, M_PI/2, 0.0};  // (t, r, θ, φ)
    
    // Initial velocity for light passing by
    double E = 1.0;
    double L = b * E;
    double f = 1.0 - 2.0 / r0;  // Schwarzschild factor at r0
    state.u = {E / f, -std::sqrt(E*E - L*L/(r0*r0) * f), 0.0, L/(r0*r0)};
    
    IntegratorConfig cfg;
    cfg.abs_tol = 1e-9;
    cfg.rel_tol = 1e-7;
    cfg.min_step = 1e-6;
    cfg.max_step = 0.02;
    
    std::cout << "t,r,theta,phi\n";
    std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    
    for(int i = 0; i < 80; ++i){
        integrate_geodesic_adaptive(metric, state, 0.5, cfg);
        
        // Stop if too close to horizon
        if (state.x[1] < 3.0) {
            std::cout << "Light too daam close to horizon,fucking program.\n";
            break;
        }
        
        std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    }
    
    return 0;
}