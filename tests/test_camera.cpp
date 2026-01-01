#include "core/camera.hpp"
#include "core/integrator.hpp"
#include "metrics/schwarzschild.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

int main() {
    // Create black hole
    SchwarzschildMetric metric(1.0);
    
    // Create camera
    // Position: r=20, θ=π/2 (equatorial), φ=0
    // Looking toward black hole (θ=π/2, φ=0)
    Camera camera(
        {0.0, 20.0, M_PI/2, 0.0},  // position
        M_PI/2,                     // look_theta (horizontal)
        0.0,                        // look_phi (toward origin)
        M_PI/3,                     // 60 degree FOV
        64,                         // width
        64                          // height
    );
    
    SchwarzschildRayGenerator ray_gen(camera, metric);
    
    // Simple rendering: check if ray hits black hole or escapes
    std::cout << "Rendering 64x64 test image...\n";
    std::ofstream out("test_render.ppm");
    out << "P3\n" << camera.width << " " << camera.height << "\n255\n";
    
    IntegratorConfig cfg;
    cfg.abs_tol = 1e-8;
    cfg.rel_tol = 1e-6;
    cfg.max_step = 0.1;
    
    for (int y = 0; y < camera.height; ++y) {
        for (int x = 0; x < camera.width; ++x) {
            // Generate ray
            GeodesicState state = ray_gen.generate_ray(x, y);
            
            // Integrate ray
            bool hit_horizon = false;
            bool escaped = false;
            
            for (int step = 0; step < 1000 && !hit_horizon && !escaped; ++step) {
                integrate_geodesic_adaptive(metric, state, 0.5, cfg);
                
                double r = state.x[1];
                
                if (r < 2.1) {  // Near horizon
                    hit_horizon = true;
                } else if (r > 50.0) {  // Escaped to infinity
                    escaped = true;
                }
            }
            
            // Color: black if hit horizon, white if escaped
            if (hit_horizon) {
                out << "0 0 0 ";  // Black (black hole shadow)
            } else {
                out << "255 255 255 ";  // White (sky)
            }
        }
        out << "\n";
    }
    
    out.close();
    std::cout << "Rendered test_render.ppm (black hole shadow)\n";
    
    return 0;
}