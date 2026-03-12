#include "core/integrator.hpp"
#include "metrics/kerr.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <atomic>
#ifdef _OPENMP
#include <omp.h>
#endif

struct RGB {
    unsigned char r, g, b;
};

int main() {
    KerrMetric metric(1.0, 0.6);
    
    const int W = 512;
    const int H = 512;
    const double FOV = 25.0;
    const double r_cam = 100.0;
    
    std::vector<RGB> image(W * H);
    
    IntegratorConfig cfg;
    cfg.abs_tol  = 1e-8;
    cfg.rel_tol  = 1e-6;
    cfg.max_step = 0.5;
    
    std::cout << "Diagnostic: Checking what each ray hits...\n";
#ifdef _OPENMP
    std::cout << "Using " << omp_get_max_threads() << " threads\n\n";
#else
    std::cout << "WARNING: OpenMP not available!\n\n";
#endif
    
    std::atomic<int> hit_disk(0);
    std::atomic<int> hit_horizon(0);
    std::atomic<int> escaped(0);
    
    #pragma omp parallel for schedule(dynamic, 4)
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            double u = (2.0 * x / W - 1.0) * FOV;
            double v = (1.0 - 2.0 * y / H) * FOV;
            double b = std::sqrt(u*u + v*v);
            double phi_dir = std::atan2(v, u);
            
            GeodesicState s;
            s.x = {0.0, r_cam, M_PI/2.0, phi_dir};
            
            double E = 1.0;
            double L = b;
            s.u[0] = E;
            s.u[1] = -E * std::sqrt(std::max(0.0, 1.0 - L*L/(r_cam*r_cam)));
            s.u[2] = 0.0;
            s.u[3] = L / (r_cam * r_cam);
            
            int outcome = 0;
            
            for (int i = 0; i < 800; ++i) {
                integrate_geodesic_adaptive(metric, s, 0.5, cfg);
                
                double r = s.x[1];
                double theta = s.x[2];
                double z = r * std::cos(theta);
                double r_cyl = r * std::sin(theta);
                
                if (std::abs(z) < 1.0 && r_cyl > 3.0 && r_cyl < 15.0) {
                    outcome = 2;
                    break;
                }
                
                if (r < 1.8) {
                    outcome = 1;
                    break;
                }
                
                if (r > 500.0) {
                    outcome = 0;
                    break;
                }
            }
            
            RGB &pix = image[y * W + x];
            
            if (outcome == 2) {
                pix = {255, 100, 50};
                hit_disk++;
            } else if (outcome == 1) {
                pix = {0, 0, 0};
                hit_horizon++;
            } else {
                pix = {255, 255, 255};
                escaped++;
            }
        }
        
        if (y % 64 == 0) {
            #pragma omp critical
            std::cout << "Row " << y << "/" << H << "\r" << std::flush;
        }
    }
    
    std::ofstream out("diagnostic_fast.ppm");
    out << "P3\n" << W << " " << H << "\n255\n";
    for (auto &p : image) {
        out << int(p.r) << " " << int(p.g) << " " << int(p.b) << " ";
    }
    out.close();
    
    std::cout << "\n\nResults:\n";
    std::cout << "  Disk: " << hit_disk << "\n";
    std::cout << "  Horizon: " << hit_horizon << "\n";
    std::cout << "  Escaped: " << escaped << "\n";
    
    return 0;
}