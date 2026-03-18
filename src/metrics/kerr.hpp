#pragma once
#include "../core/geodesic.hpp"
#include <cmath>

struct KerrMetric : public Metric {
    double M, a;

    KerrMetric(double mass = 1.0, double spin = 0.0) : M(mass), a(spin) {
        if (std::abs(a) > M) a = std::copysign(M, a);
    }

    double g_inv(int mu, int nu, const std::array<double,4>& x) const override {
        double r  = x[1];
        double th = x[2];
        double ct = std::cos(th);
        double st = std::sin(th);
        if (std::abs(st) < 1e-10) return 0.0;

        double S = r*r + a*a*ct*ct;
        double D = r*r - 2.0*M*r + a*a;
        double A = (r*r+a*a)*(r*r+a*a) - a*a*D*st*st;

        if (S < 1e-10) return 0.0;

        if (mu == 0 && nu == 0) {
            if (std::abs(D) < 1e-10) return 0.0;
            return -A / (S * D);
        }
        if ((mu==0&&nu==3)||(mu==3&&nu==0)) {
            if (std::abs(D) < 1e-10) return 0.0;
            return -2.0*M*a*r / (S * D);
        }
        if (mu == 1 && nu == 1) return D / S;
        if (mu == 2 && nu == 2) return 1.0 / S;
        if (mu == 3 && nu == 3) {
            if (std::abs(D) < 1e-10) return 0.0;
            return (D - a*a*st*st) / (S * D * st*st);
        }

        return 0.0;
    }
};