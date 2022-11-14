//
// Created by irabiel on 11/13/22.
//

#include "Inverse_Problem_Fns.h"
#include "Needed_functions.h"
#include "math.h"
#include <iostream>
#include <stdlib.h>     /* srand, rand */

Inv_Prb::Inv_Prb() {

}

void Inv_Prb::Forward_solver(std::vector<double> & x, std::vector<double> & y, double m) {
    double k1;
    double k2;
    double k3;
    double k4;
    for (int i = 0; i < x.size(); i++) {

        int iter = 0;
        int max_iter = 10;
        double s = 0.;
        double yst = y[i];
        while (iter < max_iter){
            s = -g(x[i], yst, m) / gy(x[i], yst, m);
            yst += s;
            if (abs(s) < epislon_DAE){
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

void Inv_Prb::Adjoint_solver(std::vector<double> & x, std::vector<double> & y, double m
                    , std::vector<double> & la, std::vector<double> & mu){
    la[la.size()-1] = 0;
    mu[mu.size()-1] = 0;
    // Backwards Euler
    for (int t = x.size() - 1 ; t >= 0; t--){
        la[t-1] = (1 + dt * fx(x[t], y[t],m) ) * la[t] + dt * ( (x[t] - d[t])
                                                                + mu[t]* gx(x[t],y[t],m) );
        mu[t-1] = - la[t-1] * (fy(x[t-1],y[t-1],m) / gy(x[t-1],y[t-1],m));
    }
};

void Inv_Prb::Grad_eval(std::vector<double> & x, std::vector<double> & y, std::vector<double> & la, std::vector<double> & mu, double m, double gnorm
        , double grad){
    grad = 0.;
    for (int t = 0; t < x.size(); t++)
        grad += la[t] * fm(x[t],y[t],m) + mu[t] * gm(x[t],y[t],m);

    grad *= dt;

    gnorm = abs(grad);
}

double Inv_Prb::cost(std::vector<double> & x, std::vector<double> & d){
    double cost = 0.;
    for (int t = 0; t < x.size(); t++)
        cost += (x[t] - d[t]) * (x[t] - d[t]);

    cost *= dt/2;
    return cost;
}

void Inv_Prb::Inc_state(std::vector<double> & x, std::vector<double> & y, double m
        , std::vector<double> & xh, std::vector<double> & yh, double mh){
    double temp_sum = 0.;
    for (int t = 0; t < x.size()-1; t++) {

        temp_sum = 1. - dt * fx(x[t + 1], y[t + 1], m) + dt * (gx(x[t + 1], y[t + 1], m) * fy(x[t + 1], y[t + 1], m)
                                                              / gy(x[t + 1], y[t + 1], m));

        xh[t + 1] = (1. / temp_sum) * (xh[t] + dt * mh * (fm(x[t + 1], y[t + 1], m)
                                                         - (gm(x[t + 1], y[t + 1], m) * fy(x[t + 1], y[t + 1], m)
                                                            / gy(x[t + 1], y[t + 1], m))));

        yh[t + 1] = -(1. / gy(x[t + 1], y[t + 1], m)) * (xh[t + 1] * gx(x[t + 1], y[t + 1], m)
                                                        + mh * gm(x[t + 1], y[t + 1], m));
    }
}

void Inv_Prb::Inc_Adj(std::vector<double> & x, std::vector<double> & y, double m
        , std::vector<double> & la, std::vector<double> & mu
        , std::vector<double> & xh, std::vector<double> & yh, double mh
        , std::vector<double> & lah, std::vector<double> & muh){
    lah[la.size()-1] = 0;
    muh[mu.size()-1] = 0;
    for (int t = x.size() - 1 ; t >= 0; t--){
        lah[t - 1] = (1. + dt * fx(x[t], y[t], m)) * lah[t] + dt * (xh[t] + muh[t] * gx(x[t], y[t], m)
                                                                   + la[t] * (xh[t] * fxx(x[t], y[t], m) +
                                                                              yh[t] * fxy(x[t], y[t], m)
                                                                              + mh * fxm(x[t], y[t], m))
                                                                   + mu[t] * (xh[t] * gxx(x[t], y[t], m) +
                                                                              yh[t] * gxy(x[t], y[t], m)
                                                                              + mh * gxm(x[t], y[t], m)));

        muh[t - 1] = -(1. / gy(x[t - 1], y[t - 1], m)) * (lah[t - 1] * fy(x[t - 1], y[t - 1], m)
                                                         + la[t - 1] * (xh[t - 1] * fyx(x[t - 1], y[t - 1], m) +
                                                                        yh[t - 1] * fyy(x[t - 1], y[t - 1], m)
                                                                        + mh * fym(x[t - 1], y[t - 1], m))
                                                         + mu[t - 1] * (xh[t - 1] * gyx(x[t - 1], y[t - 1], m) +
                                                                        yh[t - 1] * gyy(x[t - 1], y[t - 1], m)
                                                                        + mh * gym(x[t - 1], y[t - 1], m)));
    }
}

void Inv_Prb::Newton_update(std::vector<double> & x, std::vector<double> & y, double m
        , std::vector<double> & la, std::vector<double> & mu
        , std::vector<double> & xh, std::vector<double> & yh
        , std::vector<double> & lah, std::vector<double> & muh, double W, double Fn){
    W = 0.;
    Fn = 0.;
    for (int t = 0; t < x.size(); t++){
        W += w(x[t], y[t], m, la[t], mu[t]);
        Fn += Fnew(x[t], y[t], m, la[t], mu[t], xh[t], yh[t], lah[t], muh[t]);
    }

    W  *= dt;
    Fn *= dt;
}

double Inv_Prb::Newton_Linear_solve(std::vector<double> & x, std::vector<double> & y, double m
        , std::vector<double> & la, std::vector<double> & mu
        , std::vector<double> & xh, std::vector<double> & yh
        , std::vector<double> & lah, std::vector<double> & muh, double grad){
    double W, Fn;
    Newton_update(x ,y, m, la, mu, xh, yh, lah, muh, W,Fn);

    return -(Fn + grad)/W;
}

void Inv_Prb::assign_true_m(double m){
    m_true = m;
    std::vector<double> x, y;
    x.assign(N,0.);
    y.assign(N,0.);
    d.assign(N,0.);
    Forward_solver(x,y,m_true);
    Observation(x, d);
};

void Inv_Prb::Observation(std::vector<double> & x, std::vector<double> & d){
    double noise_lvl = 0.02;
    double noise_vec;

    for (int i = 0; i < x.size(); i++){
        noise_vec = noise_lvl * ( (float) rand()/RAND_MAX - 0.5) * 2.;
        d[i] = noise_vec + x[i];
    }
}

void Inv_Prb::Initialize_system(std::vector<double> & x, std::vector<double> & y
        , std::vector<double> & la, std::vector<double> & mu
        , std::vector<double> & xh, std::vector<double> & yh
        , std::vector<double> & lah, std::vector<double> & muh){
    // forward system
    x.assign(N, 0.);
    y.assign(N, 0.);
    x[0] = x0;

    //adjoint system
    la.assign(N, 0.);
    mu.assign(N, 0.);
    la[N-1] = 0.;
    mu[N-1] = 0.;

    //Incremental forward system
    xh.assign(N, 0.);
    yh.assign(N, 0.);
    xh[0] = 0.;
    yh[0] = 0.;

    //Incremental adjoint system
    lah.assign(N, 0.);
    muh.assign(N, 0.);
    lah[N-1] = 0.;
    muh[N-1] = 0.;
}

double Inv_Prb::Newton_Solver(double m) {

    int Max_new_iter = 10;
    int new_iter = 0;
    int arm_it_max = 10;
    int arm_it;
    double mh = 0.01;
    double err[Max_new_iter];
    err[0] = abs(m - m_true);

    double c = 1E-4;
    double alpha;
    bool line_search = true;
    bool suff_dec = false;

    std::vector<double> x,y,la,mu,xh,yh,lah,muh;
    Initialize_system(x,y,la,mu,xh,yh,lah,muh);

    double Jk,gnorm, grad, J;
    while (new_iter < Max_new_iter) {
        Forward_solver(x, y, m);
        Adjoint_solver(x, y, m, la, mu);
        Jk = cost(x, d);
        Grad_eval(x, y, la, mu, m,gnorm, grad);

        Inc_state(x, y, m, xh, yh, mh);
        Inc_Adj(x, y, m, la, mu, xh, yh, mh, lah, muh);
        mh = Newton_Linear_solve(x, y, m, la, mu, xh, yh, lah, muh, grad);

        //line search
        if (line_search)
            if (sign(mh * grad) > 0){
                alpha = 100.;
                mh = -grad ;
            }
            else
                alpha = 1.;

        arm_it = 0;
        while((~suff_dec) && (arm_it < arm_it_max)) {
            Forward_solver(x, y, m + alpha * mh);
            J = cost(x, d);
            if (J < Jk + c * alpha * mh * grad) {
                suff_dec = true;
                m += alpha * mh;
            }
            else {
                alpha *= 0.5;
                arm_it += 1;
            }
        }

        m += mh;

        new_iter += 1;
        err[new_iter] = abs(m - m_true);
        std::cout << new_iter << " " <<  Jk << " " <<  abs(grad) << " " << abs(mh * grad) << " " << err[new_iter] << std::endl;

        if (abs(grad) < 1E-4 || abs(mh * grad) < 1E-4) {
            std::cout <<"Termiination reason: Norm of the gradient less than tolerance" << std::endl;
            break;
        }
    }
}