#include "core/integrator.hpp"
#include "core/starfield.hpp"
#include "metrics/kerr.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

struct RGB {
    unsigned char r, g, b;
    
    static RGB from_double(double r, double g, double b) {
        return {
            static_cast<unsigned char>(std::clamp(r * 255.0, 0.0, 255.0)),
            static_cast<unsigned char>(std::clamp(g * 255.0, 0.0, 255.0)),
            static_cast<unsigned char>(std::clamp(b * 255.0, 0.0, 255.0))
        };
    }
};

struct KerrAccretionDisk {
    double M, a;
    double inner_radius, outer_radius, thickness;
    
    KerrAccretionDisk(double mass, double spin, double r_in, double r_out, double thick)
        : M(mass), a(spin), inner_radius(r_in), outer_radius(r_out), thickness(thick) {}
    
    struct Emission {
        bool hit;
        double r, g, b;
    };
    
    Emission sample(double r, double theta) const {
        double z = r * std::cos(theta);
        if (std::abs(z) > thickness) return {false, 0, 0, 0};
        
        double r_cyl = r * std::sin(theta);
        if (r_cyl < inner_radius || r_cyl > outer_radius) return {false, 0, 0, 0};
        
        double temp_factor = (outer_radius - r_cyl) / (outer_radius - inner_radius);
        temp_factor = std::pow(temp_factor, 0.5);
        
        double red   = 0.4 + 0.6 * (1.0 - temp_factor);
        double green = 0.3 + 0.5 * temp_factor;
        double blue  = 0.2 + 0.8 * temp_factor;
        double brightness = temp_factor * 2.5;
        
        return {true, red * brightness, green * brightness, blue * brightness};
    }
};

int main() {
    double M = 1.0;
    double a = 0.6;  // Reduced from 0.9 for stability
    
    KerrMetric metric(M, a);
    StarField stars(5000);  // Reduced from 10000
    KerrAccretionDisk disk(M, a, 2.5, 12.0, 0.6);
    
    const int width  = 512;   // Reduced from 1024 for speed
    const int height = 512;
    const double fov = 20.0;  // Reduced from 25
    
    std::cout << "Kerr Black Hole Renderer\n";
    std::cout << "M=" << M << ", a=" << a << ", " << width << "x" << height << "\n\n";
    
    std::vector<RGB> image(width * height);
    
    IntegratorConfig cfg;
    cfg.abs_tol  = 1e-7;   // Relaxed from 1e-9
    cfg.rel_tol  = 1e-5;   // Relaxed from 1e-7
    cfg.min_step = 1e-5;   // Increased from 1e-6
    cfg.max_step = 0.5;    // Increased from 0.05 (10x faster!)
    
    const double r_cam = 100.0;
    
    #pragma omp parallel for schedule(dynamic, 4)
    for (int py = 0; py < height; ++py) {
        for (int px = 0; px < width; ++px) {
            double u = (2.0 * px / width  - 1.0) * fov;
            double v = (1.0 - 2.0 * py / height) * fov;
            double b = std::sqrt(u*u + v*v);
            
            GeodesicState state;
            state.x = {0.0, r_cam, M_PI / 2.0, 0.0};
            
            // Simplified initial conditions (distant observer approximation)
            double E = 1.0;
            double L = b;
            state.u[0] = E;
            state.u[1] = -E * std::sqrt(std::max(0.0, 1.0 - L*L/(r_cam*r_cam)));
            state.u[2] = 0.0;
            state.u[3] = L / (r_cam * r_cam);
            
            bool hit_horizon = false;
            bool hit_disk = false;
            double disk_r = 0, disk_g = 0, disk_b = 0;
            double final_theta = state.x[2], final_phi = state.x[3];
            
            // OPTIMIZED: Fewer, larger steps
            for (int step = 0; step < 500 && !hit_horizon && !hit_disk; ++step) {
                integrate_geodesic_adaptive(metric, state, 1.0, cfg);  // Larger interval!
                
                double r = state.x[1];
                double theta = state.x[2];
                
                // Check disk
                auto emission = disk.sample(r, theta);
                if (emission.hit) {
                    hit_disk = true;
                    disk_r = emission.r;
                    disk_g = emission.g;
                    disk_b = emission.b;
                }
                
                // Kerr horizon
                double r_horizon = M + std::sqrt(M*M - a*a);
                if (r < r_horizon + 0.3) {
                    hit_horizon = true;
                }
                
                // Escaped
                if (r > 300.0) {
                    final_theta = theta;
                    final_phi = state.x[3];
                    break;
                }
            }
            
            RGB &pix = image[py * width + px];
            
            if (hit_disk) {
                pix = RGB::from_double(disk_r, disk_g, disk_b);
            } else if (hit_horizon) {
                pix = {0, 0, 0};
            } else {
                auto color = stars.sample(final_theta, final_phi);
                pix = RGB::from_double(color.r, color.g, color.b);
            }
        }
        
        if (py % 64 == 0) {
            #pragma omp critical
            std::cout << "Row " << py << "/" << height << "\r" << std::flush;
        }
    }
    
    std::ofstream out("kerr_fast.ppm");
    out << "P3\n" << width << " " << height << "\n255\n";
    for (int py = 0; py < height; ++py) {
        for (int px = 0; px < width; ++px) {
            RGB &p = image[py * width + px];
            out << int(p.r) << " " << int(p.g) << " " << int(p.b) << " ";
        }
        out << "\n";
    }
    out.close();
    
    std::cout << "\nDone! Saved to kerr_fast.ppm\n";
    
    return 0;
}