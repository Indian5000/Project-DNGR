#include "core/integrator.hpp"
#include "metrics/kerr.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>   // FIX: required for std::clamp

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

int main() {
    KerrMetric metric(1.0, 0.6);

    const int width  = 512;
    const int height = 512;
    const double fov = 20.0;
    const double r_cam = 100.0;

    std::cout << "DIAGNOSTIC RENDERER - Checking phi propagation\n\n";

    std::vector<RGB> image(width * height);

    IntegratorConfig cfg;
    cfg.abs_tol  = 1e-7;
    cfg.rel_tol  = 1e-5;
    cfg.max_step = 0.5;

    for (int py = 0; py < height; ++py) {
        for (int px = 0; px < width; ++px) {

            // Screen-space coordinates
            double u = (2.0 * px / width  - 1.0) * fov;
            double v = (1.0 - 2.0 * py / height) * fov;

            double impact = std::sqrt(u * u + v * v);
            double phi_dir = std::atan2(v, u);

            GeodesicState state;
            state.x = {0.0, r_cam, M_PI / 2.0, phi_dir};

            double E = 1.0;
            double L = impact;

            state.u[0] = E;
            state.u[1] = -E * std::sqrt(
                std::max(0.0, 1.0 - (L * L) / (r_cam * r_cam))
            );
            state.u[2] = 0.0;
            state.u[3] = L / (r_cam * r_cam);

            // Integrate briefly
            for (int step = 0; step < 100; ++step) {
                integrate_geodesic_adaptive(metric, state, 1.0, cfg);

                if (state.x[1] < 2.5 || state.x[1] > 300.0)
                    break;
            }

            // Color by final phi
            double final_phi = state.x[3];

            // Normalize phi from [-π, π] → [0, 1]
            double phi_norm = (final_phi + M_PI) / (2.0 * M_PI);
            phi_norm = phi_norm - std::floor(phi_norm); // wrap safely

            double r_col, g_col, b_col;

            if (phi_norm < 1.0 / 3.0) {
                r_col = 1.0 - 3.0 * phi_norm;
                g_col = 3.0 * phi_norm;
                b_col = 0.0;
            } else if (phi_norm < 2.0 / 3.0) {
                r_col = 0.0;
                g_col = 2.0 - 3.0 * phi_norm;
                b_col = 3.0 * phi_norm - 1.0;
            } else {
                r_col = 3.0 * phi_norm - 2.0;
                g_col = 0.0;
                b_col = 3.0 - 3.0 * phi_norm;
            }

            image[py * width + px] = RGB::from_double(r_col, g_col, b_col);
        }

        if (py % 64 == 0) {
            std::cout << "Row " << py << "/" << height << "\r" << std::flush;
        }
    }

    std::ofstream out("diagnostic_phi.ppm");
    out << "P3\n" << width << " " << height << "\n255\n";
    for (int py = 0; py < height; ++py) {
        for (int px = 0; px < width; ++px) {
            const RGB &p = image[py * width + px];
            out << int(p.r) << " " << int(p.g) << " " << int(p.b) << " ";
        }
        out << "\n";
    }
    out.close();

    std::cout << "\nDone! Saved to diagnostic_phi.ppm\n";
    std::cout << "\nIF YOU SEE:\n";
    std::cout << "  - Rainbow/colorful pattern → phi IS changing (good!)\n";
    std::cout << "  - Concentric rings of one color → phi NOT changing (bug in Kerr metric)\n\n";

    return 0;
}
