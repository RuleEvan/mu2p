#include "romberg.h"
double RombergIntegrator(double (*f)(double), double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(a) + (*f)(b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}

/*double RombergSpline(double (*f)(double, gsl_spline*, gsl_interp_accel*, double), double a, double b, double p, gsl_spline *f_spline, gsl_interp_accel *acc, double tol) {
  // Numerical integrator
  int maxiter = 25;
  int maxj = 6;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(p, f_spline, acc, a) + f(p, f_spline, acc, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(p, f_spline, acc, a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("RombSpline failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}
*/

double RombergSpline(double (*f)(double, gsl_spline*, gsl_interp_accel*, double), double a, double b, double p, gsl_spline *f_spline, gsl_interp_accel *acc, double tol) {
  // Numerical integrator
/*
  size_t max_steps = 35;

  double R1[max_steps], R2[max_steps];
  double *Rp=&R1[0], *Rc=&R2[0];
  double h = (b-a);
  Rp[0] = (f(p, f_spline, acc, a) + f(p, f_spline, acc, b))*h*0.5;

  for (size_t i = 1; i < max_steps; ++i) {
	  h *= 0.5;
	  double c = 0;
	  size_t ep = 1 << (i - 1);
	  for (size_t j = 1; j <= ep; ++j) {
		  c += f(p, f_spline, acc, a + (2*j-1)*h);
	  }
	  Rc[0] = h*c + 0.5*Rp[0];

	  for (size_t j = 1; j <= i; ++j) {
		  double n_k = pow(4, j);
		  Rc[j] = (n_k*Rc[j-1] - Rp[j - 1])/(n_k - 1);
	  }
	  if (i > 1 && fabs(Rp[i-1]-Rc[i]) < tol) {
		  return Rc[i-1];
	  }

	  double *rt = Rp;
	  Rp = Rc;
	  Rc = rt;
  }
  return Rp[max_steps-1];
*/

  int maxiter = 25;
  int maxj = 6;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(p, f_spline, acc, a) + f(p, f_spline, acc, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(p, f_spline, acc, a + (2*k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;

}

double RombergSpline2FunIntegrator(double (*f)(double (*lep1) (double), double (*lep2) (double), gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , double), double (*lep1) (double), double (*lep2) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 25;
  int maxj = 6;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(lep1, lep2, gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a) + f(lep1, lep2, gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(lep1, lep2, gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}


double RombergSplineFunIntegrator(double (*f)(double (*lep) (double), gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , double), double (*lep) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 25;
  int maxj = 6;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(lep, gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a) + f(lep, gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(lep, gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}


double RombergSpline41Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 25;
  int maxj = 6;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a) + f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}


double RombergSpline42Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , double, double, double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double q1, double q2, double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 25;
  int maxj = 6;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a, q1, q2) + f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, b, q1, q2));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a + (k + k - 1)*h, q1, q2);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}


double RombergSpline4Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_spline*, gsl_spline*, gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , gsl_interp_accel* , double, double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel* acc1, gsl_interp_accel* acc2, gsl_interp_accel* acc3, gsl_interp_accel* acc4, double q, double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 25;
  int maxj = 6;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a, q) + f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, b, q));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(gmu_spline, fmu_spline, ge_spline, fe_spline, acc1, acc2, acc3, acc4, a + (k + k - 1)*h, q);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}


double RombergSpline2Integrator(double (*f)(gsl_spline*, gsl_spline*, gsl_interp_accel* , double), gsl_spline* g_spline, gsl_spline* f_spline, gsl_interp_accel* acc, double a, double b, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * (f(g_spline, f_spline, acc, a) + f(g_spline, f_spline, acc, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + f(g_spline, f_spline, acc, a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}


 double Romberg2Vars(double (*f)(double, double), double a, double b, double r, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  float h, g0, fourj, gmax, error, g1, romb;
  float g[maxj + 1];
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(a, r) + (*f)(b, r));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(a + (k + k - 1)*h, r);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  
  return romb;
}
 


double Romberg7Vars(double (*f)(double, int, double, int, double, double), double a, double b, double r, int L1, double kt, int L2, double qt, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  double h, g0, fourj, gmax, error, g1, romb;
  double *g = (double*) malloc(sizeof(double)*(maxj + 1));
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(r, L1, kt, L2, qt, a) + (*f)(r, L1, kt, L2, qt, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(r, L1, kt, L2, qt, a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  free(g);
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint 5Var failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  return romb;
}


double Romberg5Vars(double (*f)(double, int, double, double, double), double a, double b, double r, int l, double q, double m, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  double h, g0, fourj, gmax, error, g1, romb;
  double *g = (double*) malloc(sizeof(double)*(maxj + 1));
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(r, l, q, m, a) + (*f)(r, l, q, m, b));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(r, l, q, m, a + (k + k - 1)*h);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  free(g);
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint 5Var failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  return romb;
}


double Romberg3Vars(double (*f)(double, double, int), double a, double b, double p, int iv, double tol) {
  // Numerical integrator
  int maxiter = 20;
  int maxj = 5;
  double h, g0, fourj, gmax, error, g1, romb;
  double *g = (double*) malloc(sizeof(double)*(maxj + 1));
  int nint;
  int i, j, jmax, k;
  h = 0.5 * (b - a);
  gmax = h * ((*f)(a, p, iv) + (*f)(b, p, iv));
  g[0] = gmax;
  nint = 1;
  error = 1.0 * pow(10, 20);
  i = 0;
  while (1 == 1) {
    i = i + 1;
    if (i > maxiter || ((i > 9) && (fabs(error) < tol))){break;}
    g0 = 0.0;
    for (k = 1; k < nint + 1; k++) {
      g0 = g0 + (*f)(a + (k + k - 1)*h, p, iv);
    }
 
    g0 = 0.5 * g[0] + h * g0;
    h = 0.5 * h;
    nint = nint + nint;
    jmax = MIN(i, maxj);
    fourj = 1.0;
    for (j = 0; j < jmax; j++) {
      fourj = 4.0 * fourj;
      g1 = g0 + (g0 - g[j]) / (fourj - 1.0);
      g[j] = g0;
      g0 = g1;
    }

    if (fabs(g0) > tol) {
      error = 1.0 - gmax / g0;
    } else {
      error = gmax;
    }
    gmax = g0;
    g[jmax] = g0;
  }
  romb = g0;
  free(g);
  if (i > maxiter && fabs(error) > tol) {
    printf("Rombint failed to converge; integral= %g, error= %g \n", romb, error);
    
    return romb;
  }
  return romb;
}

