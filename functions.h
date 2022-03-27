#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream> 
#include <sstream>
#include <fstream>

#include "definities.h"

double f1(double t, double* U, int ind);
double df1(double t, double* U, int f_ind, int d_ind);
double fa1(double t, int ind);
void NumMultVect(double* A, double a, size_t N, double* O);
void SumVect(double* A, double* B, size_t N, double* O);

void RK4_one_step(double* y0, double* y1, Rfunc f, double tau, double t, long n);
void RK4(double* U, double tau, double t0, double T, Rfunc f, size_t n);
double* RK4_step_var(double* U, double tau, double t0, double T, Rfunc f, Afunc fa, size_t n, double tol);
double DP4_one_step(double* y0, double* y1, Rfunc f, double tau, double t, long n, long mode);
double* DormanPrince4_5(double* U, double tau0, double t0, double T, Rfunc f, Afunc fa, size_t n, double tol, double fac0, double facmin0, double facmax0);