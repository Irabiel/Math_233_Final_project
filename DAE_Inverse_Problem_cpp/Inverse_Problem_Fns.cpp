//
// Created by irabiel on 11/13/22.
//

#include "Inverse_Problem_Fns.h"
#include "Needed_functions.h"

void Forward_solver(std::vector<double> & x, std::vector<double> & y, double m) {

    for (int i = 0; i < x.size(); i++) {

        int iter = 0;
        int max_iter = 10;
        double s = 0;
        double yst = y[i];
        while (iter < max_iter){
            s = -g(x[i], yst, m) / gy(x[i], yst, m);
            yst += s;
            if abs(s) < epislon_DAE{
                y[i] = yst;
                break;
            }
            iter += 1;
        }

        if (i != x.size() - 1){
            k1 = f(x[i], y[i], m);
            k2 = f(x[i] + .5 * dt * k1, y[i], m);
            k3 = f(x[i] + .5 * dt * k2, y[i], m);
            k4 = f(x[i] + dt * k3, y[i], m);
            x[i + 1] = x[i] + (dt / 6.) * (k1 + 2. * k2 + 2. * k3 + k4);
        }
    }
};

void Adjoint_solver(std::vector<double> & x, std::vector<double> & y, double m
                    , std::vector<double> & la, std::vector<double> & mu){
    la[la.size()-1] = 0;
    mu[mu.size()-1] = 0;
    // Backwards Euler
    for (int t = x.size() - 2 ; t >= 0; t--){
        la[t] = (1. + dt * fx(x[t + 1], y[t + 1], m)) * la[t + 1] + dt * ((x[t + 1] - d[t + 1])
                                                                         + mu[t + 1] * gx(x[t + 1], y[t + 1], m));
        mu[t] = -la[t] * (fy(x[t], y[t], m) / gy(x[t], y[t], m));
    }
};

