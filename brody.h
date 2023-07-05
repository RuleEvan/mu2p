#ifndef BRODY_H
#define BRODY_H
#include "potential.h"
double g_factor(double a, double b, double c);
double h_factor(double a, double b, double c);
double brody_mosh_zero(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int l1, int l2); 
double brody_mosh(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int n1, int l1, int n2, int l2);

// q = 0 functions
double compute_radial_matrix_element_scalar(int iv, int n1p, int l1p, int n2p, int l2p, int n1, int l1, int n2, int l2, int lamda, int s, int t);
double compute_relative_density_scalar(double r, int n1p, int l1p, int n2p, int l2p, int n1, int l1, int n2, int l2, int lamda, int s, int t);
double compute_radial_matrix_element_y2(int iv, int n1p, int l1p, int n2p, int l2p, int lambdap, int n1, int l1, int n2, int l2, int lamda, int s, int t);


// Finite q Spline functions
void compute_radial_matrix_element_GTJ(int n1p, int l1p, int n2p, int l2p, int lambdap, int n1, int l1, int n2, int l2, int lambda, int s, int t12, int J, double q, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *rm_RE, double *rm_IM);

void compute_radial_matrix_element_GTJFull(int n1p, int l1p, int n2p, int l2p, int lambdap, int n1, int l1, int n2, int l2, int lambda, int s, int t, int mu_rel, int mu_cm, int J, double q, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *rm_RE, double *rm_IM); 

#endif
