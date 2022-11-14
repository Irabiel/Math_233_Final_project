#include <iostream>

#include "Inverse_Problem_Fns.h"
#include "cmath"
#include <vector>
int main() {
    Inv_Prb IP = Inv_Prb();

    // parameters
    double t_init         = 0.;
    double t_final        = 1.;
    double t_obs          = 1.;
    double dt             = .01;

    IP.assign_dt(dt);

    int n = round ((t_final - t_init) /dt);
    IP.assign_N(n);

    // DAE Newton system parameter
    double epislon_DAE = 1e-6;
    IP.assign_DAE_tol(epislon_DAE);

    //True parameter
    double mtrue = 0.693;
    IP.assign_true_m(mtrue);

    // initial guess
    double m = 4.609;
    IP.Newton_Solver(m);

    return 0;
}
