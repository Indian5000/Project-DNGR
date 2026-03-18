#pragma once
#include "../core/geodesic.hpp"
#include <cmath>

struct MorrisThorne : public Metric {
    double b0;
    explicit MorrisThorne(double throat_radius = 1.0) : b0(throat_radius) {}

    double g_inv(int mu, int nu, const std::array<double,4>& x) const override {
        if (mu != nu) return 0.0;
        double l     = x[1];
        double theta = x[2];
        double st    = std::sin(theta);
        if (std::abs(st) < 1e-10) return 0.0;
        double r2 = b0*b0 + l*l;
        switch (mu) {
            case 0: return -1.0;
            case 1: return  1.0;
            case 2: return  1.0 / r2;
            case 3: return  1.0 / (r2 * st * st);
            default: return 0.0;
        }
    }
};