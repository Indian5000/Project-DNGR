#include "core/integrator.hpp"
#include "metrics/flat.hpp"
#include <iostream>

int main() {
    MinkowskiMetric metric;
    
    GeodesicState state;
    state.x = {0,0,0,0};  // (t, x, y, z)
    state.u = {1,1,0,0};  // velocity along t and x
    
    IntegratorConfig cfg;
    cfg.abs_tol = 1e-8;
    cfg.rel_tol = 1e-6;
    cfg.min_step = 1e-6;
    cfg.max_step = 0.1;
    
    double dt_print = 0.1;
    int n_steps = 10;
    
    std::cout << "t,x,y,z\n";
    std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    
    // Integrate in small increments
    for(int i = 0; i < n_steps; ++i){
        integrate_geodesic_adaptive(metric, state, dt_print, cfg);
        std::cout << state.x[0] << "," << state.x[1] << "," << state.x[2] << "," << state.x[3] << "\n";
    }
    
    return 0;
}