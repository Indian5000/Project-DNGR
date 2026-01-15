#include "core/integrator.hpp"
#include "metrics/kerr.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

int main() {
    KerrMetric metric(1.0, 0.6);
    
    // Test a few specific rays
    std::cout << "Testing ray paths in Kerr spacetime:\n\n";
    
    for (double b : {0.0, 5.0, 10.0, 15.0}) {
        std::cout << "Ray with impact parameter b = " << b << "\n";
        
        double r_cam = 100.0;
        GeodesicState state;
        state.x = {0.0, r_cam, M_PI/2, 0.0};
        
        double E = 1.0;
        double L = b;
        state.u[0] = E;
        state.u[1] = -E * std::sqrt(std::max(0.0, 1.0 - L*L/(r_cam*r_cam)));
        state.u[2] = 0.0;
        state.u[3] = L / (r_cam * r_cam);
        
        IntegratorConfig cfg;
        cfg.abs_tol = 1e-7;
        cfg.rel_tol = 1e-5;
        cfg.max_step = 0.5;
        
        double r_min = r_cam;
        bool hit_horizon = false;
        
        for (int step = 0; step < 200; ++step) {
            integrate_geodesic_adaptive(metric, state, 1.0, cfg);
            
            double r = state.x[1];
            r_min = std::min(r_min, r);
            
            if (step % 50 == 0) {
                std::cout << "  Step " << step << ": r=" << r 
                         << ", theta=" << state.x[2] 
                         << ", phi=" << state.x[3] << "\n";
            }
            
            double r_horizon = 1.0 + std::sqrt(1.0 - 0.36);
            if (r < r_horizon + 0.3) {
                std::cout << "  --> HIT HORIZON at r=" << r << "\n";
                hit_horizon = true;
                break;
            }
            
            if (r > 300.0) {
                std::cout << "  --> ESCAPED to r=" << r << "\n";
                break;
            }
        }
        
        std::cout << "  Min radius: " << r_min << "\n";
        std::cout << "  Status: " << (hit_horizon ? "CAPTURED" : "ESCAPED") << "\n\n";
    }
    
    return 0;
}