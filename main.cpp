#include "header.h"
#include <stdlib.h> 

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "Russian");

	double tau    = atof(argv[1]);
	double fac    = strtod(argv[2], NULL);
	double facmin = strtod(argv[3], NULL);
	double facmax = strtod(argv[4], NULL);

	Rfunc f  = f1;              // функция правой части
    Afunc fa = fa1;             // аналитическое решение
    Df    df = df1;             // производная правой части

	int n = 2;                  // размерность задачи
    double* u = new double[n];  // массив с начальным приближением
	double* info = new double[3];
	clock_t time_dp, time_rk;
   
	double t0 = 0.3;           	// начальное время
	double T = 4.3;            	// конечное время
	double tol = 1e-6;

    u[0] = fa1(t0, 0);           // начальные условия
    u[1] = fa1(t0, 1);

	time_dp = clock();
	info = DormanPrince4_5(u, tau, t0, T, f, fa, n, tol, fac, facmin, facmax);
	time_dp = clock() - time_dp;

	time_rk = clock();
	info = RK4_step_var(u, tau, t0, T, f, fa, n, tol, fac, facmin, facmax);
	time_rk = clock() - time_rk;
	

	printf("time_rk = %lf\n", (double)time_rk / CLOCKS_PER_SEC);
	printf("time_dp = %lf\n", (double)time_dp / CLOCKS_PER_SEC);

	std::ofstream ZZZ("rk.dat", std::ios::app);
	ZZZ << (double)time_rk / CLOCKS_PER_SEC << "\n";
	ZZZ.close();

	ZZZ.open("dp.dat", std::ios::app);
	ZZZ << (double)time_dp / CLOCKS_PER_SEC << "\n";
	ZZZ.close();


	delete[] u;
	delete[] info;
	return 0;
}
