#include "core/integrator.hpp"
#include "metrics/wormhole.hpp"
#include <iostream>
#include <cmath>

int main() {
    // Wormhole with throat radius b0=1
    MorrisThorne metric(1.0);
    
    // Light ray passing through wormhole
    GeodesicState state;
    
    // Start on one side (l = -5), shoot toward throat
    state.x = {0.0, -5.0, M_PI/2, 0.0};  // (t, l, θ, φ)
    state.u = {1.0, 0.8, 0.0, 0.1};      // Moving radially inward with slight angular component
    
    IntegratorConfig cfg;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-8;
    cfg.min_step = 1e-6;
    cfg.max_step = 0.05;
    
    std::cout << "t,l,theta,phi\n";
    std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    
    for(int i = 0; i < 50; ++i){
        integrate_geodesic_adaptive(metric, state, 0.5, cfg);
        std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    }
    
    return 0;
}