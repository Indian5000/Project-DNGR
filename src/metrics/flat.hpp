#pragma once
#include "../core/geodesic.hpp"

// Minkowski metric in Cartesian coordinates
// g^{μν} = diag(-1, 1, 1, 1)
struct MinkowskiMetric : public Metric {
    double g_inv(int mu, int nu, const std::array<double,4>& x) const override {
        if (mu != nu) return 0.0;
        return (mu == 0) ? -1.0 : 1.0;
    }
};