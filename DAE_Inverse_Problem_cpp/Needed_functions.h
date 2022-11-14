//
// Created by Irabiel on 11/10/22.
//

#ifndef DAE_INVERSE_PROBLEM_CPP_NEEDED_FUNCTIONS_H
#define DAE_INVERSE_PROBLEM_CPP_NEEDED_FUNCTIONS_H


double f(double x, double y, double m);
double fx(double x, double y, double m);
double fy(double x, double y, double m);
double fm(double x, double y, double m);
double fxx(double x, double y, double m);
double fxy(double x, double y, double m);
double fxm(double x, double y, double m);
double fyy(double x, double y, double m);
double fym(double x, double y, double m);
double fyx(double x, double y, double m);
double fmx(double x, double y, double m);
double fmy(double x, double y, double m);
double fmm(double x, double y, double m);


double g(double x, double y, double m);
double gx(double x, double y, double m);
double gy(double x, double y, double m);
double gm(double x, double y, double m);
double gxx(double x, double y, double m);
double gxy(double x, double y, double m);
double gxm(double x, double y, double m);
double gyy(double x, double y, double m);
double gym(double x, double y, double m);
double gyx(double x, double y, double m);
double gmx(double x, double y, double m);
double gmy(double x, double y, double m);
double gmm(double x, double y, double m);

//    Function required for the adjoint based newton system
double w(double x,double y, double m, double la, double mu);
double Fnew(double x, double y, double m, double la, double mu, double xh, double yh, double lah, double muh);

// sign function
double sign(double x);

#endif //DAE_INVERSE_PROBLEM_CPP_NEEDED_FUNCTIONS_H
