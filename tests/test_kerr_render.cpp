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
    const double FOV       = 15.0;
    const double r_cam     = 100.0;
    const double theta_cam = 1.4;
    const double r_horizon = 1.0 + std::sqrt(1.0 - 0.36);

    const double disk_r_in  = 3.0;
    const double disk_r_out = 12.0;

    std::vector<RGB> image(W * H);

    IntegratorConfig cfg;
    cfg.abs_tol  = 1e-8;
    cfg.rel_tol  = 1e-6;
    cfg.max_step = 0.5;

        std::cout << "Kerr Black Hole Renderer (a=0.6) — Super-Hamiltonian\n";
        std::cout << "Resolution: " << W << "x" << H << "\n";
    #ifdef _OPENMP
        std::cout << "Threads: " << omp_get_max_threads() << "\n";
    #else
        std::cout << "Threads: 1 (OpenMP not enabled)\n";
    #endif
        std::cout << "r_horizon = " << r_horizon << "\n\n";

    // DEBUG: center ray — p_t and p_phi should stay constant
    {
        std::cout << "=== DEBUG: center ray ===\n";
        GeodesicState dbg;
        dbg.x = {0.0, r_cam, theta_cam, 0.0};
        dbg.p = {-1.0, -1.0, 0.0, 0.0};
        for (int i = 0; i < 30; ++i) {
            integrate_geodesic_adaptive(metric, dbg, 0.5, cfg);
            std::cout << "  step " << i
                      << " r="     << dbg.x[1]
                      << " p_t="   << dbg.p[0]
                      << " p_phi=" << dbg.p[3] << "\n";
            if (dbg.x[1] < r_horizon + 0.2) { std::cout << "  -> HIT HORIZON\n"; break; }
            if (dbg.x[1] > 500.0)           { std::cout << "  -> ESCAPED\n";      break; }
        }
        std::cout << "=== END DEBUG ===\n\n";
    }

    std::atomic<int> count_disk(0), count_horizon(0), count_stars(0), count_timeout(0);

    #pragma omp parallel for schedule(dynamic, 4)
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            double u_px = (2.0 * x / W - 1.0) * FOV;
            double v_px = (1.0 - 2.0 * y / H) * FOV;

            GeodesicState s;
            s.x = {0.0, r_cam, theta_cam, 0.0};
            s.p[0] = -1.0;
            s.p[2] = -v_px;
            s.p[3] =  u_px;

            // Exact null condition: g^{μν} p_μ p_ν = 0, solve for p_r
            {
                double gi_tt   = metric.g_inv(0, 0, s.x);
                double gi_tphi = metric.g_inv(0, 3, s.x);
                double gi_thth = metric.g_inv(2, 2, s.x);
                double gi_phph = metric.g_inv(3, 3, s.x);
                double gi_rr   = metric.g_inv(1, 1, s.x);
                double H = -gi_tt * s.p[0]*s.p[0]
                           - 2.0*gi_tphi * s.p[0]*s.p[3]
                           - gi_thth * s.p[2]*s.p[2]
                           - gi_phph * s.p[3]*s.p[3];
                s.p[1] = -std::sqrt(std::max(0.0, H / gi_rr));
            }

            bool hit_disk    = false;
            bool hit_horizon = false;
            double disk_r = 0, disk_g = 0, disk_b = 0;
            double final_theta = M_PI / 2.0;
            double final_phi   = 0.0;
            double prev_theta  = s.x[2];
            bool timed_out     = true;

            for (int i = 0; i < 5000; ++i) {
                integrate_geodesic_adaptive(metric, s, 0.5, cfg);

                double r     = s.x[1];
                double theta = s.x[2];

                bool crossed = (prev_theta - M_PI/2.0) * (theta - M_PI/2.0) < 0;

                if (crossed && !hit_disk) {
                    double r_cyl = r * std::sin(theta);
                    if (r_cyl > disk_r_in && r_cyl < disk_r_out) {
                        hit_disk = true;
                        double temp = std::pow((disk_r_out - r_cyl) /
                                               (disk_r_out - disk_r_in), 0.7);
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
                    timed_out   = false;
                    break;
                }

                if (r > 500.0 || (r > r_cam && s.p[1] > 0)) {
                    final_theta = theta;
                    final_phi   = std::fmod(s.x[3], 2.0 * M_PI);
                    if (final_phi < 0) final_phi += 2.0 * M_PI;
                    timed_out   = false;
                    break;
                }
            }

            if (timed_out) count_timeout++;

            RGB& pix = image[y * W + x];
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
    for (auto& p : image)
        out << int(p.r) << " " << int(p.g) << " " << int(p.b) << " ";
    out.close();

    std::cout << "\n\nResults:\n";
    std::cout << "  Disk hits:    " << count_disk    << "\n";
    std::cout << "  Horizon hits: " << count_horizon << "\n";
    std::cout << "  Star field:   " << count_stars   << "\n";
    std::cout << "  Timed out:    " << count_timeout << "\n";
    std::cout << "\n✨ Saved kerr_proper.ppm\n";

    return 0;
}