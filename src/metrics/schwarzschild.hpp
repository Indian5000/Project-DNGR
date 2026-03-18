#pragma once
#include "../core/geodesic.hpp"
#include <cmath>

struct SchwarzschildMetric : public Metric {
    double M;
    explicit SchwarzschildMetric(double mass = 1.0) : M(mass) {}

    double g_inv(int mu, int nu, const std::array<double,4>& x) const override {
        if (mu != nu) return 0.0;
        double r     = x[1];
        double theta = x[2];
        double f     = 1.0 - 2.0 * M / r;
        if (f <= 0.0) return 0.0;
        double st = std::sin(theta);
        if (std::abs(st) < 1e-10) return 0.0;
        switch (mu) {
            case 0: return -1.0 / f;
            case 1: return  f;
            case 2: return  1.0 / (r * r);
            case 3: return  1.0 / (r * r * st * st);
            default: return 0.0;
        }
    }
};