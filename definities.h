#pragma once

typedef double(*funct) (double x);
typedef double*(*Funct) (double* x);
typedef double(*Rfunc) (double, double*, int);
typedef double(*Afunc) (double, int);
typedef double(*Df) (double, double*, int, int);

const double eps = 1e-5;

