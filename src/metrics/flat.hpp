#pragma once
#include "../core/geodesic.hpp"

// Minkowski metric: η_μν = diag(-1, 1, 1, 1)
// All Christoffel symbols are zero (flat spacetime)
struct MinkowskiMetric : public Metric {
    double christoffel(int mu, int nu, int rho, const std::array<double,4>& x) const override {
        return 0.0;  // All Christoffel symbols vanish in flat space
    }
};