#ifndef MATRIX_ELEMENT_H
#define MATRIX_ELEMENT_H
#include "brody.h"

// Finite q Spline functions
void compute_total_matrix_element_GTJ(char* density_file, double q, int J, double alpha, double delta_e, double *m_RE, double *m_IM);

void compute_matrix_element_GTJ(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM);

void compute_total_matrix_element_FJ(char* density_file, double q, int J, double alpha, double delta_e, double *m_RE, double *m_IM);

void compute_matrix_element_FJ(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM);

void filter_spectrum(char* density_file, double J, double T);

double compute_total_matrix_element_GT0(char* density_file);
double compute_matrix_element_GT0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12); 

double compute_total_matrix_element_GT2(char* density_file);
double compute_matrix_element_GT2(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12); 


double compute_total_matrix_element_F0(char* density_file);
double compute_matrix_element_F0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12); 
void generate_bigstick_int_file(int i_model, double q);


int get_l(int n, int j);
int get_n(int index);
int get_j(int index);

double cg_fact(int l1, int l2, int k, double Ji, double Jf, int iJn);

int test_suite();

double compute_matrix_element_F1(int ina, int ija, int inb, int ijb, int L, double q);
double compute_matrix_element_F1_double(int ina, int ija, int inb, int ijb, int L1, double kt, int L2, double qt, int L);

void compute_matrix_element_FJFull(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int l, int L, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM);
void compute_matrix_element_GTJFull(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int l, int L, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM);
void compute_matrix_element_T2JFull(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int mu_rel, int mu_cm, int L, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM);


double compute_total_matrix_element_F1_double(char* density_file_i, char* density_file_f, double qt, double kt);
#endif
