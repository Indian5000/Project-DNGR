#include "core/integrator.hpp"
#include "core/starfield.hpp"
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
    static RGB from_double(double r, double g, double b) {
        return {
            (unsigned char)std::clamp(r * 255.0, 0.0, 255.0),
            (unsigned char)std::clamp(g * 255.0, 0.0, 255.0),
            (unsigned char)std::clamp(b * 255.0, 0.0, 255.0)
        };
    }
};

int main() {
    KerrMetric metric(1.0, 0.6);
    StarField stars(15000);
    
    const int W = 800;
    const int H = 800;
    const double FOV = 15.0;
    const double r_cam = 100.0;
    const double r_horizon = 1.0 + std::sqrt(1.0 - 0.36);
    
    const double disk_r_in = 3.0;
    const double disk_r_out = 12.0;
    
    std::vector<RGB> image(W * H);
    
    IntegratorConfig cfg;
    cfg.abs_tol  = 1e-8;
    cfg.rel_tol  = 1e-6;
    cfg.max_step = 0.5;
    
    std::cout << "Kerr Black Hole Renderer (a=0.6)\n";
    std::cout << "Resolution: " << W << "x" << H << "\n";
    std::cout << "Threads: " << omp_get_max_threads() << "\n";
    std::cout << "r_horizon = " << r_horizon << "\n\n";

    {
        std::cout << "=== DEBUG: center ray trace ===\n";
        GeodesicState dbg;
        dbg.x = {0.0, r_cam, 1.4, 0.0};
        dbg.u[0] = 1.0; dbg.u[1] = -1.0; dbg.u[2] = 0.0; dbg.u[3] = 0.0;
        for (int i = 0; i < 30; ++i) {
            integrate_geodesic_adaptive(metric, dbg, 0.5, cfg);
            std::cout << "  step " << i << " r=" << dbg.x[1] << " theta=" << dbg.x[2] << "\n";
            if (dbg.x[1] < r_horizon + 0.2) { std::cout << "  -> HIT HORIZON\n"; break; }
            if (dbg.x[1] > 500.0)           { std::cout << "  -> ESCAPED\n";      break; }
        }

        std::cout << "\n=== DEBUG: slight offset ray (u=1, v=0) ===\n";
        GeodesicState dbg2;
        dbg2.x = {0.0, r_cam, 1.4, 0.0};
        dbg2.u[0] = 1.0; dbg2.u[1] = -std::sqrt(1.0 - 1.0/(r_cam*r_cam));
        dbg2.u[2] = 0.0; dbg2.u[3] = 1.0;
        for (int i = 0; i < 30; ++i) {
            integrate_geodesic_adaptive(metric, dbg2, 0.5, cfg);
            std::cout << "  step " << i << " r=" << dbg2.x[1] << "\n";
            if (dbg2.x[1] < r_horizon + 0.2) { std::cout << "  -> HIT HORIZON\n"; break; }
            if (dbg2.x[1] > 500.0)           { std::cout << "  -> ESCAPED\n";      break; }
        }
        std::cout << "=== END DEBUG ===\n\n";
    }
    
    std::atomic<int> count_disk(0);
    std::atomic<int> count_horizon(0);
    std::atomic<int> count_stars(0);
    std::atomic<int> count_timeout(0);
    
    #pragma omp parallel for schedule(dynamic, 4)
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            double u = (2.0 * x / W - 1.0) * FOV;
            double v = (1.0 - 2.0 * y / H) * FOV;
            
            GeodesicState s;
            s.x = {0.0, r_cam, 1.4, 0.0};
            
            s.u[0] = 1.0;
            s.u[1] = -std::sqrt(std::max(0.0, 1.0 - (u*u + v*v)/(r_cam*r_cam)));
            s.u[2] = -v / (r_cam * r_cam);
            s.u[3] =  u / (r_cam * r_cam);
            
            bool hit_disk = false;
            bool hit_horizon = false;
            double disk_r = 0, disk_g = 0, disk_b = 0;
            double final_theta = M_PI/2.0;
            double final_phi = 0.0;
            
            double prev_theta = s.x[2];
            bool timed_out = true;
            
            for (int i = 0; i < 5000; ++i) {
                integrate_geodesic_adaptive(metric, s, 0.5, cfg);
                
                double r = s.x[1];
                double theta = s.x[2];
                
                bool crossed_equator = (prev_theta - M_PI/2.0) * (theta - M_PI/2.0) < 0;
                
                if (crossed_equator && !hit_disk) {
                    double r_cross = r;
                    double r_cyl = r_cross * std::sin(theta);
                    
                    if (r_cyl > disk_r_in && r_cyl < disk_r_out) {
                        hit_disk = true;
                        
                        double temp = (disk_r_out - r_cyl) / (disk_r_out - disk_r_in);
                        temp = std::pow(temp, 0.7);
                        
                        disk_r = 0.9 - 0.3 * temp;
                        disk_g = 0.4 * temp;
                        disk_b = 0.15 * temp;
                        timed_out = false;
                        break;
                    }
                }
                
                prev_theta = theta;
                
                if (r < r_horizon + 0.2) {
                    hit_horizon = true;
                    timed_out = false;
                    break;
                }
                
                if (r > 500.0 || (r > r_cam && s.u[1] > 0)) {
                    final_theta = theta;
                    final_phi = std::fmod(s.x[3], 2.0 * M_PI);
                    if (final_phi < 0) final_phi += 2.0 * M_PI;
                    timed_out = false;
                    break;
                }
            }
            
            if (timed_out) count_timeout++;

            RGB &pix = image[y * W + x];
            
            if (hit_disk) {
                pix = RGB::from_double(disk_r, disk_g, disk_b);
                count_disk++;
            } else if (hit_horizon) {
                pix = {0, 0, 0};
                count_horizon++;
            } else {
                auto c = stars.sample(final_theta, final_phi);
                pix = RGB::from_double(c.r, c.g, c.b);
                count_stars++;
            }
        }
        
        if (y % 100 == 0) {
            #pragma omp critical
            std::cout << "Progress: " << (100 * y / H) << "%\r" << std::flush;
        }
    }
    
    std::ofstream out("kerr_proper.ppm");
    out << "P3\n" << W << " " << H << "\n255\n";
    for (auto &p : image) {
        out << int(p.r) << " " << int(p.g) << " " << int(p.b) << " ";
    }
    out.close();
    
    std::cout << "\n\nResults:\n";
    std::cout << "  Disk hits:    " << count_disk    << "\n";
    std::cout << "  Horizon hits: " << count_horizon << "\n";
    std::cout << "  Star field:   " << count_stars   << "\n";
    std::cout << "  Timed out:    " << count_timeout << "\n";
    std::cout << "\n✨ Saved kerr_proper.ppm\n";
    
    return 0;
}