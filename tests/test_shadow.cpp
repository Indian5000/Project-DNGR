#include "core/integrator.hpp"
#include "metrics/schwarzschild.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

struct RGB {
    unsigned char r, g, b;
};

int main() {
    SchwarzschildMetric metric(1.0);
    
    const int width  = 512;  // Higher resolution
    const int height = 512;
    const double fov = 15.0;  // Slightly smaller FOV for better view
    
    std::cout << "Rendering " << width << "x" << height
              << " black hole shadow...\n";
    
    std::vector<RGB> image(width * height);
    
    IntegratorConfig cfg;
    cfg.abs_tol  = 1e-8;
    cfg.rel_tol  = 1e-6;
    cfg.max_step = 0.1;
    
    const double r_cam = 100.0;
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int py = 0; py < height; ++py) {
        for (int px = 0; px < width; ++px) {
            // Screen coordinates
            double u = (2.0 * px / width  - 1.0) * fov;
            double v = (1.0 - 2.0 * py / height) * fov;
            double b = std::sqrt(u*u + v*v);
            double phi_dir = std::atan2(v, u);
            
            GeodesicState state;
            state.x = {0.0, r_cam, M_PI / 2.0, phi_dir};
            
            double E = 1.0;
            double L = b;
            
            // FIXED: Simple initial conditions for distant observer
            state.u[0] = E;
            state.u[1] = -E * std::sqrt(1.0 - L*L/(r_cam*r_cam));
            state.u[2] = 0.0;
            state.u[3] = L / (r_cam * r_cam);
            
            bool hit_horizon = false;
            
            for (int step = 0; step < 500; ++step) {
                integrate_geodesic_adaptive(metric, state, 1.0, cfg);
                
                double r = state.x[1];
                
                if (r < 2.1) {
                    hit_horizon = true;
                    break;
                }
                
                if (r > r_cam * 1.5) {
                    break;
                }
            }
            
            RGB &pix = image[py * width + px];
            if (hit_horizon) {
                pix = {0, 0, 0};  // Black hole shadow
            } else {
                pix = {255, 255, 255};  // Sky
            }
        }
        
        if (py % 64 == 0) {
            std::cout << "Row " << py << "/" << height << "\n";
        }
    }
    
    std::ofstream out("shadow.ppm");
    out << "P3\n" << width << " " << height << "\n255\n";
    for (int py = 0; py < height; ++py) {
        for (int px = 0; px < width; ++px) {
            RGB &p = image[py * width + px];
            out << int(p.r) << " " << int(p.g) << " " << int(p.b) << " ";
        }
        out << "\n";
    }
    out.close();
    
    std::cout << "Done! Saved to shadow.ppm\n";
    std::cout << "You should see a black circle (the black hole shadow)!\n";
    
    return 0;
}