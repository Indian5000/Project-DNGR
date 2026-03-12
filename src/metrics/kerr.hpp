#pragma once
#include "../core/geodesic.hpp"
#include <cmath>

struct KerrMetric : public Metric {
    double M;
    double a;

    KerrMetric(double mass=1.0, double spin=0.0)
        : M(mass), a(spin) {
        if (std::abs(a) > M) a = std::copysign(M, a);
    }

    inline double Sigma(double r, double th) const {
        return r*r + a*a*std::cos(th)*std::cos(th);
    }

    inline double Delta(double r) const {
        return r*r - 2*M*r + a*a;
    }

    double christoffel(int mu, int nu, int rho,
                       const std::array<double,4>& x) const override {

        double r  = x[1];
        double th = x[2];

        double st = std::sin(th);
        double ct = std::cos(th);

        if (st == 0.0) return 0.0;

        double S = Sigma(r, th);
        double D = Delta(r);
        double A = (r*r + a*a)*(r*r + a*a) - a*a*D*st*st;

        if (S == 0.0 || D == 0.0) return 0.0;

        // --- RADIAL ---
        if (mu==1 && nu==0 && rho==0)
            return M*(r*r - a*a*ct*ct)*D/(S*S*S);

        if (mu==1 && nu==1 && rho==1)
            return (-M*r*r + a*a*r*st*st + M*a*a*ct*ct) / (S*D);

        if (mu==1 && nu==2 && rho==2)
            return -r*D/S;

        if (mu==1 && nu==3 && rho==3)
            return -D*st*st*(r*S - M*(r*r - a*a*ct*ct))/(S*S*S);

        // --- THETA ---
        if (mu==2 && nu==1 && rho==2)
            return r/S;
        if (mu==2 && nu==2 && rho==1)
            return r/S;

        if (mu==2 && nu==3 && rho==3)
            return -st*ct*(A + 2*M*a*a*r*st*st)/(S*S);

        // --- PHI ---
        if (mu==3 && nu==0 && rho==1)
            return 2*M*a*r/(D*S);
        if (mu==3 && nu==1 && rho==0)
            return 2*M*a*r/(D*S);

        if (mu==3 && nu==1 && rho==3)
            return (r*A - 2*M*r*S)/(D*S*A);
        if (mu==3 && nu==3 && rho==1)
            return (r*A - 2*M*r*S)/(D*S*A);

        if (mu==3 && nu==2 && rho==3)
            return ct/(st*S);
        if (mu==3 && nu==3 && rho==2)
            return ct/(st*S);

        // --- TIME (FRAME DRAGGING — REQUIRED) ---
        if (mu==0 && nu==1 && rho==3)
            return 2*M*a*r/(D*S);
        if (mu==0 && nu==3 && rho==1)
            return 2*M*a*r/(D*S);

        return 0.0;
    }
};
