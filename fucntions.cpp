#include "functions.h"

//////////////////// Функции  ///////////////////////////////
double f1(double t, double* U, int ind)
{
	if (ind == 0)
		return 2 * t * U[0] * log(std::max(U[1], 1e-3));
	else
		return -2 * t * U[1] * log(std::max(U[0], 1e-3));
}


double df1(double t, double* U, int f_ind, int d_ind)
{
    double ret_value = -1;

	if (f_ind == 0)
        if (d_ind == 0)
            ret_value = 2 * t * log(std::max(U[1], 1e-3));
        else
            if(U[1] > 1e-3)
                ret_value = 2 * t * U[0] / U[1];
            else
                ret_value = 0;
    else
		if (d_ind == 0)
            if(U[0] > 1e-3)
                ret_value = -2 * t * U[1] / U[0];
            else
                ret_value = 0;
        else
            ret_value = -2 * t * log(std::max(U[0], 1e-3));
}

double fa1(double t, int ind)
{
	if (ind == 0)
		return exp(sin(t*t));
	else
		return exp(cos(t*t));
}

void NumMultVect(double* A, double a, size_t N, double* O)
{
	for (int i = 0; i < N; i++)
		O[i] = A[i] * a;
}

void SumVect(double* A, double* B, size_t N, double* O)
{
	for (int i = 0; i < N; i++)
		O[i] = A[i] + B[i];
}

//////////////////// Метод Рунге — Кутты 4-ого порядка  ///////////////////////////////
void RK4(double* U, double tau, double t0, double T, Rfunc f, size_t n)    
{

	std::cout << "RK4 launched \n\n";

    int count = 1, i, j, s;
	double  t = t0, m;
    double* y  = new double[n];
    double* y1 = new double[n];
    double* k1 = new double[n];
    double* k2 = new double[n];
    double* k3 = new double[n];
    double* k4 = new double[n];

    double* k_temp1 = new double[n];
    double* k_temp2 = new double[n];

    memcpy(y1, U, sizeof(double)*n);

	std::ofstream RK4_dat("RK4.dat");

    RK4_dat << t0 << "   ";
	for (j = 0; j < n; ++j)
		RK4_dat << y1[j] << "   ";
	RK4_dat << "\n";

	while(true) {
		if (t > T) { break; };
		RK4_dat << t << "   ";
	
		for (s = 0; s < n; ++s)
			k1[s] = f(t, y1, s);

		for (s = 0; s < n; ++s) {
            NumMultVect(k1, tau / 2., n, k_temp1);
            SumVect(y1, k_temp1, n, k_temp2);
            k2[s] = f(t + tau / 2., k_temp2, s);
        }
			

		for (s = 0; s < n; ++s) {
            NumMultVect(k2, tau / 2., n, k_temp1);   
            SumVect(y1, k_temp1, n, k_temp2);
            k3[s] = f(t + tau / 2., k_temp2, s);
        }
			

        for (s = 0; s < n; ++s) {
            NumMultVect(k3, tau, n, k_temp1);   
            SumVect(y1, k_temp1, n, k_temp2); 
            k4[s] = f(t + tau, k_temp2, s);
        }

		for (s = 0; s < n; ++s)
		{
			y[s] = y1[s] + tau * (k1[s] + 2 * k2[s] + 2 * k3[s] + k4[s]) / 6.;  
			RK4_dat << y[s] << "   ";
		} 
		
		count++;
		RK4_dat << "\n"; 
        memcpy(y1, y, sizeof(double)*n);
        t = t0 + tau * count;
	} 

	RK4_dat.close();

    delete[] y;
    delete[] y1;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] k_temp1;
    delete[] k_temp2;
}

void RK4_one_step(double* y0, double* y1, Rfunc f, double tau, double t, long n) {

    int s;

    double* k_temp1 = new double[n];
    double* k_temp2 = new double[n];
    double* k1 = new double[n];
    double* k2 = new double[n];
    double* k3 = new double[n];
    double* k4 = new double[n];
    
    for (s = 0; s < n; ++s)
        k1[s] = f(t, y0, s);

    for (s = 0; s < n; ++s) {
        NumMultVect(k1, tau / 2., n, k_temp1);
        SumVect(y0, k_temp1, n, k_temp2);
        k2[s] = f(t + tau / 2., k_temp2, s);
    }
        

    for (s = 0; s < n; ++s) {
        NumMultVect(k2, tau / 2., n, k_temp1);   
        SumVect(y0, k_temp1, n, k_temp2);
        k3[s] = f(t + tau / 2., k_temp2, s);
    }
        

    for (s = 0; s < n; ++s) {
        NumMultVect(k3, tau, n, k_temp1);   
        SumVect(y0, k_temp1, n, k_temp2); 
        k4[s] = f(t + tau, k_temp2, s);
    }

    for (s = 0; s < n; ++s)
        y1[s] = y0[s] + tau * (k1[s] + 2 * k2[s] + 2 * k3[s] + k4[s]) / 6.;  


    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] k_temp1;
    delete[] k_temp2;
    
}

double* RK4_step_var(double* U, double tau0, double t0, double T, Rfunc f, Afunc fa, size_t n, double tol, double fac0, double facmin0, double facmax0)   
{

    int i, j, s;
    double t = t0, error, sq;
    double taumax = 0.5, tau = tau0;
    double fac = fac0, facmin = facmin0, facmax = facmax0;
    double* y0 = new double[n];
    double* y01 = new double[n];
    double* y1 = new double[n];
    double* y2 = new double[n];
    double* d  = new double[n];
    double* info = new double[3]; // 0 - шаги; 1 - принятые шаги; 2 - глобальная погрешность
    
    memcpy(y0, U, sizeof(double)*n);

	std::ofstream RK4_dat("RK4.dat");
    std::ofstream ZZZ("rk.dat", std::ios::app);

    RK4_dat << t << "   ";
	for (j = 0; j < n; ++j)
		RK4_dat << y0[j] << "   ";
    RK4_dat << tau << "   ";
	RK4_dat << "\n";

    info[0] = 0; info[1] = 0; info[2] = 0;
	while(true) {
		
		if (t > T) { break; };
	
        memcpy(y01, y0, sizeof(double)*n);
		RK4_one_step(y01, y1, f, tau, t, n);     // два шага метода РК с шагом tau
        memcpy(y01, y1, sizeof(double)*n);
        RK4_one_step(y01, y1, f, tau, t + tau, n);

        RK4_one_step(y0, y2, f, tau*2, t, n);   // один шаг метода РК с шагом 2*tau

       
        for(s = 0; s < n; ++s) d[s] = 1;
        // for(s = 0; s < n; ++s) d[s] = fabs(y1[s]);

        error = fabs(y1[0] - y2[0]) / d[0];
        for(s = 1; s < n; ++s) {
            double temp = fabs(y1[s] - y2[s]) / d[s];
            if (temp > error)
                error = temp;
        }
        error = error / 15.;
        
        if(error < tol) {
            t = t + 2*tau;
            memcpy(y0, y1, sizeof(double)*n);
            RK4_dat << t << "   ";
            for (s = 0; s < n; ++s)
			    RK4_dat << y0[s] << "   ";
            RK4_dat << tau << "   ";
            RK4_dat << "\n"; 
            info[1]++;   
            facmax = facmax0;
        }
        else {
            facmax = 1;
        }

        tau = tau * std::min(facmax, std::max(facmin, fac * pow((tol / error), 1 / 5.)));
		info[0]++;
    
	} 

    info[2] = fabs(y0[0] - fa(t, 0));
    for(s = 0; s < n; ++s) {
        double temp = fabs(y0[s] - fa(t, s));
        if(info[2] < temp) info[2] = temp;
    }
    
    RK4_dat << info[0] << " " << info[0] - info[1] << " " << info[2];
    ZZZ << tau0 << " " << fac0 << " " << facmin0 << " " << facmax0 << " "  << info[0] << " " << info[0] - info[1] << " " << info[2] << " ";
    
	RK4_dat.close();
    ZZZ.close();

    delete[] info;
    delete[] d;
    delete[] y0;
    delete[] y01;
    delete[] y1;
    delete[] y2;

}

double DP4_one_step(double* y0, double* y1, Rfunc f, double tau, double t, long n, long mode) {

    int s;
    double error = 0, sq = 0;

    double* yt1 = new double[n];
    double* yt2 = new double[n];
    double* yt3 = new double[n];
    double* yt4 = new double[n];
    double* yt5 = new double[n];
    double* yt6 = new double[n];
    double* kt1 = new double[n];
    double* kt2 = new double[n];
    double* kt3 = new double[n];
    double* kt4 = new double[n];
    double* kt5 = new double[n];
    double* kt6 = new double[n];
    double* k1 = new double[n];
    double* k2 = new double[n];
    double* k3 = new double[n];
    double* k4 = new double[n];
    double* k5 = new double[n];
    double* k6 = new double[n];
    double* k7 = new double[n];
    
    for (s = 0; s < n; ++s)
        k1[s] = f(t, y0, s);

    for (s = 0; s < n; ++s) {
        NumMultVect(k1, tau / 5., n, kt1);
        SumVect(y0, kt1, n, yt1);
        k2[s] = f(t + tau * 0.2, yt1, s);
    }
        

    for (s = 0; s < n; ++s) {
        NumMultVect(k1, tau * (3 / 40.), n, kt1); 
        NumMultVect(k2, tau * (9 / 40.), n, kt2);   
        SumVect(y0, kt1, n, yt1);
        SumVect(yt1, kt2, n, yt2);
        k3[s] = f(t + tau * 0.3, yt2, s);
    }
        

    for (s = 0; s < n; ++s) {
        NumMultVect(k1, tau * (44 / 45.), n, kt1); 
        NumMultVect(k2, tau * (-56 / 15.), n, kt2);  
        NumMultVect(k3, tau * (32 / 9.), n, kt3);   
        SumVect(y0, kt1, n, yt1);
        SumVect(yt1, kt2, n, yt2);
        SumVect(yt2, kt3, n, yt3);
        k4[s] = f(t + tau * 0.8, yt3, s);
    }

    for (s = 0; s < n; ++s) {
        NumMultVect(k1, tau * (19372 / 6561.), n, kt1); 
        NumMultVect(k2, tau * (-25360 / 2187.), n, kt2);  
        NumMultVect(k3, tau * (64448 / 6561.), n, kt3); 
        NumMultVect(k4, tau * (-212 / 729.), n, kt4);   
        SumVect(y0, kt1, n, yt1);
        SumVect(yt1, kt2, n, yt2);
        SumVect(yt2, kt3, n, yt3);
        SumVect(yt3, kt4, n, yt4);
        k5[s] = f(t + tau * 8 / 9., yt4, s);
    }

    for (s = 0; s < n; ++s) {
        NumMultVect(k1, tau * (9017 / 3168.), n, kt1); 
        NumMultVect(k2, tau * (-355 / 33.), n, kt2);  
        NumMultVect(k3, tau * (46732 / 5247.), n, kt3); 
        NumMultVect(k4, tau * (49 / 176.), n, kt4);  
        NumMultVect(k5, tau * (-5103 / 18656.), n, kt5);   
        SumVect(y0, kt1, n, yt1);
        SumVect(yt1, kt2, n, yt2);
        SumVect(yt2, kt3, n, yt3);
        SumVect(yt3, kt4, n, yt4);
        SumVect(yt4, kt5, n, yt5);
        k6[s] = f(t + tau, yt5, s);
    }

    // if (mode == 1) {
        for (s = 0; s < n; ++s) {
            NumMultVect(k1, tau * (35 / 384.) , n, kt1); 
            NumMultVect(k3, tau * (500 / 1113.), n, kt3); 
            NumMultVect(k4, tau * (125 / 192.), n, kt4);  
            NumMultVect(k5, tau * (-2187 / 6784.), n, kt5);
            NumMultVect(k6, tau * (11 / 84.), n, kt6);   
            SumVect(y0, kt1, n, yt2);
            SumVect(yt2, kt3, n, yt3);
            SumVect(yt3, kt4, n, yt4);
            SumVect(yt4, kt5, n, yt5);
            SumVect(yt5, kt6, n, yt6);
            k7[s] = f(t + tau, yt6, s);
        }
        // printf("k7 = %.16f\n", k7[0]/40.*tau);
    // }

    // if (mode == 0) {
        for (s = 0; s < n; ++s)
            y1[s] = y0[s] + k1[s]* (35 / 384.)*tau+ k3[s]* (500 / 1113.)*tau+ k4[s] * (125 / 192.)*tau+ k5[s]* (-2187 / 6784.)*tau+ k6[s]* (11 / 84.)*tau;
        // printf("k0 = %f %f %f %f %f\n", k1[0], k3[0], k4[0], k5[0], k6[0]); 
        // printf("k0*coeff = %f %f %f %f %f\n", k1[0]* (35 / 384.), k3[0]* (500 / 1113.), k4[0] * (125 / 192.), k5[0]* (-2187 / 6784.), k6[0]* (11 / 84.)); 
        // double kk1 = k1[0]* (35 / 384.)*tau;
        // double kk3 = k3[0]* (500 / 1113.)*tau;
        // double kk4 = k4[0] * (125 / 192.)*tau;
        // double kk5 = k5[0]* (-2187 / 6784.)*tau;
        // double kk6 = k6[0]* (11 / 84.)*tau;
        // printf("k0*coeff*tau = %f %f %f %f %f\n", kk1, kk3, kk4, kk5, kk6); 
        // printf("k0*tau = %.16f\n", kk1 + kk3 + kk4 + kk5 + kk6);
        // printf("000 %f\n\n", y0[0] + k1[0]* (35 / 384.)*tau+ k3[0]* (500 / 1113.)*tau+ k4[0] * (125 / 192.)*tau+ k5[0]* (-2187 / 6784.)*tau+ k6[0]* (11 / 84.)*tau); 
    // }

    for(s = 0; s < n; ++s) {
        error += fabs(k7[s] - k6[s]);
        sq += k7[s]*k7[s];
    }

    // if (mode == 1) {
    //     for (s = 0; s < n; ++s)
    //         y1[s] = y0[s] +  k1[s]* (5179 / 57600.)*tau+ k3[s]* (7571 / 16695.)*tau+ k4[s]* (393 / 640.)*tau+ k5[s] * (-92097 / 339200.)*tau+ k6[s]* (187 / 2100.)*tau + k7[s] / 40.*tau;  
        // printf("k1 = %f %f %f %f %f %f\n", k1[0], k3[0], k4[0], k5[0], k6[0], k7[0]);     
        // printf("k1*coeff = %f %f %f %f %f %f\n", k1[0]* (5179 / 57600.), k3[0]* (7571 / 16695.), k4[0]* (393 / 640.), k5[0] * (-92097 / 339200.), k6[0]* (187 / 2100.) , k7[0] / 40.);
        // double kk1 = k1[0]* (5179 / 57600.)*tau;
        // double kk3 = k3[0]* (7571 / 16695.)*tau;
        // double kk4 = k4[0]* (393 / 640.)*tau;
        // double kk5 = k5[0] * (-92097 / 339200.)*tau;
        // double kk6 = k6[0]* (187 / 2100.)*tau;
        // double kk7 = k7[0] / 40.*tau;
        // printf("k1*coeff*tau = %f %f %f %f %f %f\n", kk1, kk3, kk4, kk5, kk6, kk7);  
        // printf("k1*tau = %.16f\n", kk1 + kk3 + kk4 + kk5 + kk6 + kk7);
        // printf("111 %f\n\n", y0[0] +  k1[0]* (5179 / 57600.)*tau+ k3[0]* (7571 / 16695.)*tau+ k4[0]* (393 / 640.)*tau+ k5[0] * (-92097 / 339200.)*tau+ k6[0]* (187 / 2100.)*tau + k7[0] / 40.*tau); 
    // }


    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] k5;
    delete[] k6;
    delete[] k7;
    delete[] kt1;
    delete[] kt2;
    delete[] kt3;
    delete[] kt4;
    delete[] kt5;
    delete[] kt6;
    delete[] yt1;
    delete[] yt2;
    delete[] yt3;
    delete[] yt4;
    delete[] yt5;
    delete[] yt6;

    return (error / sq);
    
    
}

double* DormanPrince4_5(double* U, double tau0, double t0, double T, Rfunc f, Afunc fa, size_t n, double tol, double fac0, double facmin0, double facmax0) {

    int i, j, s;
	double  t = t0, error, sq;
    double taumax = 0.5, tau = tau0;
    double fac = fac0, facmin = facmin0, facmax = facmax0;
    double g;
    double* y0 = new double[n];
    double* y1 = new double[n];
    double* y2 = new double[n];
    double* d  = new double[n];
    double* info = new double[3]; // 0 - шаги; 1 - принятые шаги; 2 - глобальная погрешность
    
    memcpy(y0, U, sizeof(double)*n);

	std::ofstream DP4_5("DP4_5.dat");
    std::ofstream ZZZ("dp.dat", std::ios::app);

    DP4_5 << t << "   ";
	for (j = 0; j < n; ++j)
		DP4_5 << y0[j] << "   ";
    DP4_5 << tau << "   ";
	DP4_5 << "\n";
    DP4_5.close();
    DP4_5.open("DP4_5.dat", std::ios::app);

    info[0] = 0; info[1] = 0; info[2] = 0;
	while(true) {
		
		if (t > T) { break; };
	
		error = DP4_one_step(y0, y1, f, tau, t, n, 0);   
        
        if(error < tol) {
            t = t + tau;
            memcpy(y0, y1, sizeof(double)*n);
            DP4_5 << t << "   ";
            for (s = 0; s < n; ++s)
                DP4_5 << y0[s] << "   ";
            DP4_5 << tau << "   ";
            DP4_5 << "\n"; 
            facmax = facmax0;
            g = std::min(facmax, std::max(facmin, pow(fac*(tol/error), 0.2)));
            tau = std::min(g * tau, taumax);
            info[1]++;  
        }
        else {
            facmax = 1;
            g = std::min(facmax, std::max(facmin, pow(fac*(tol/error), 0.2)));
            tau = g * tau;
        }

        info[0]++;
	} 

    info[2] = fabs(y0[0] - fa(t, 0));
    for(s = 0; s < n; ++s) {
        double temp = fabs(y0[s] - fa(t, s));
        if(info[2] < temp) info[2] = temp;
    }
    
    DP4_5 << info[0] << " " << info[0] - info[1] << " " << info[2];
    ZZZ << tau0 << " " << fac0 << " " << facmin0 << " " << facmax0 << " "  << info[0] << " " << info[0] - info[1] << " " << info[2] << " ";
    
	DP4_5.close();
    ZZZ.close();

    delete[] info;
    delete[] d;
    delete[] y0;
    delete[] y1;
    delete[] y2;

}
