//
// Created by Irabiel on 11/10/22.
//

#include "Needed_functions.h"
#include "cmath"

double f(double x, double y, double m){
    return -exp(m) * x * x + sin(y);
}

double fx(double x, double y, double m){
    return -exp(m) * x * 2.;
};
double fy(double x, double y, double m){
    return cos(y);
};
double fm(double x, double y, double m){
    return -exp(m) * x * x;
};
double fxx(double x, double y, double m){
    return -2. * exp(m);
};
double fxy(double x, double y, double m){
    return 0.;
};
double fxm(double x, double y, double m){
    return -2. * exp(m) * x;
};
double fyy(double x, double y, double m){
    return -sin(y);
};
double fym(double x, double y, double m){
    return 0.;
};
double fyx(double x, double y, double m){
    return 0.;
};
double fmx(double x, double y, double m){
    return -2. * exp(m) * x;
};
double fmy(double x, double y, double m){
    return 0.;
};
double fmm(double x, double y, double m){
    return -exp(m) * x * x;
};


double g(double x, double y, double m){
    return  exp( -( (x - 1.) * (x - 1.) )/2. ) - (y*y*y - y +1.);
};
double gx(double x, double y, double m){
    return -(x - 1.) * exp( -( (x - 1.) * (x - 1.) )/2. );
};
double gy(double x, double y, double m){
    return 1. - 3. * y * y;
};
double gm(double x, double y, double m){
    return 0.;
};
double gxx(double x, double y, double m){
    return ( ( (x - 1.) * (x - 1.)  ) - 1. ) * exp( -( (x - 1.) * (x - 1.) )/2. );
};
double gxy(double x, double y, double m){
    return 0.;
};
double gxm(double x, double y, double m){
    return 0.;
};
double gyy(double x, double y, double m){
    return -6. * y;
};
double gym(double x, double y, double m){
    return 0.;
};
double gyx(double x, double y, double m){
    return 0.;
};
double gmx(double x, double y, double m){
    return 0.;
};
double gmy(double x, double y, double m){
    return 0.;
};
double gmm(double x, double y, double m){
    return 0.;
};

//    Function required for the adjoint based newton system
double w(double x,double y, double m, double la, double mu){
    return la * fmm(x,y,m) + mu * gmm(x,y,m);
};
double Fnew(double x, double y, double m, double la, double mu, double xh, double yh, double lah, double muh){
    return  (( lah * fm(x,y,m) + muh * gm(x,y,m) + la * ( fmx(x,y,m) * xh + fmy(x,y,m) * yh ) ) + mu * ( gmx(x,y,m) * xh + gmy(x,y,m) * yh ) );
};
