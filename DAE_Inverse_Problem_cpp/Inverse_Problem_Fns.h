//
// Created by irabiel on 11/13/22.
//

#ifndef MATH_233_FINAL_PROJECT_INVERSE_PROBLEM_FNS_H
#define MATH_233_FINAL_PROJECT_INVERSE_PROBLEM_FNS_H

#include "Needed_functions.h"
#include <vector>

// solve the foward problem
void Forward_solver(std::vector<double> & x, std::vector<double> & y, double m);
// solve the adjoint problem
void Adjoint_solver(std::vector<double> & x, std::vector<double> & y, double m
                    , std::vector<double> & la, std::vector<double> & mu);

// evaluate the gradient
void Grad_eval(std::vector<double> & x, std::vector<double> & y, double m, double gnorm
               , double grad);
// calculate the cost
double cost(std::vector<double> & x, std::vector<double> & d);

// solve the incremental state problem
void Inc_state(std::vector<double> & x, std::vector<double> & y, double m
               , std::vector<double> & xh, std::vector<double> & yh, double mh);
// solve the incremental adjoint problem
void Inc_Adj(std::vector<double> & x, std::vector<double> & y, double m
             , std::vector<double> & la, std::vector<double> & mu
             , std::vector<double> & xh, std::vector<double> & yh, double mh
             , std::vector<double> & lah, std::vector<double> & muh);

// Create the newton scaler hessian and components that don't include mh
void Newton_update(std::vector<double> & x, std::vector<double> & y, double m
                   , std::vector<double> & la, std::vector<double> & mu
                   , std::vector<double> & xh, std::vector<double> & yh
                   , std::vector<double> & lah, std::vector<double> & muh, double W, double Fn);
// solve the newton linear system
double New_linear_solve(std::vector<double> & x, std::vector<double> & y, double m
                        , std::vector<double> & la, std::vector<double> & mu
                        , std::vector<double> & xh, std::vector<double> & yh
                        , std::vector<double> & lah, std::vector<double> & muh, double grad);

// create the observation vector d
void Observation(std::vector<double> & x, std::vector<double> & d);

// Initialize the vectors need for the inverse problem
void Initialize_system(int N, std::vector<double> & x, std::vector<double> & y
                       , std::vector<double> & la, std::vector<double> & mu
                       , std::vector<double> & xh, std::vector<double> & yh
                       , std::vector<double> & lah, std::vector<double> & muh);

#endif //MATH_233_FINAL_PROJECT_INVERSE_PROBLEM_FNS_H