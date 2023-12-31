#include "matrix_element.h"

/* Two-body matrix elements */

void compute_total_matrix_element_GTJ(char* density_file, double q, int J, double alpha, double delta_e, double *mat_RE, double *mat_IM) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline_RE = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array_RE = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
//  printf("Beginning splines\n");
  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array_RE[i] = v_NUC(J, q, r);//(v_light_nu_RE(J, q, alpha, delta_e, r);
    f_array_IM[i] = 0;//v_light_nu_IM(J, q, alpha, r);
  }
//  exit(0);

  gsl_spline_init(f_spline_RE, r_array, f_array_RE, NSPLINE + 1);
  gsl_spline_init(f_spline_IM, r_array, f_array_IM, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_GTJ(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
    *mat_RE += 4.0*M_PI*m4_RE*density;
    *mat_IM += 4.0*M_PI*m4_IM*density;

  }

  gsl_spline_free(f_spline_RE);
  gsl_spline_free(f_spline_IM);
  gsl_interp_accel_free(acc);  
  
  return;
}


void compute_matrix_element_GTJ(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4_RE = 0.0;
  double m4_IM = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
      // Reduced matrix element of tau(1)_- tau(2)_-
      m1 *= sqrt(5);
      m1 *= pow(-1.0, J)/sqrt(2.0*J + 1.0);

      m1 *= fact;
      double rm_RE = 0.0;
      double rm_IM = 0.0;
      compute_radial_matrix_element_GTJ(n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, J, q, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
      m4_RE += m1*rm_RE;
      m4_IM += m1*rm_IM;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}

  *m_RE = m4_RE;
  *m_IM = m4_IM;
 
//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, *m_RE);

  return;
}

void compute_total_matrix_element_FJ(char* density_file, double q, int J, double alpha, double delta_e, double *mat_RE, double *mat_IM) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
 
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline_RE = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array_RE = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
//  printf("Beginning splines\n");
  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array_RE[i] = v_NUC(J, q, r);//v_light_nu_RE(J, q, alpha, delta_e, r);
    f_array_IM[i] = 0;//v_light_nu_IM(J, q, alpha, r);
  }
//  exit(0);

  gsl_spline_init(f_spline_RE, r_array, f_array_RE, NSPLINE + 1);
  gsl_spline_init(f_spline_IM, r_array, f_array_IM, NSPLINE + 1);

  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_FJ(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, q, J, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
    *mat_RE += 4.0*M_PI*m4_RE*density;
    *mat_IM += 4.0*M_PI*m4_IM*density;

  }

  gsl_spline_free(f_spline_RE);
  gsl_spline_free(f_spline_IM);
  gsl_interp_accel_free(acc);  
  
  return;
}


void compute_matrix_element_FJ(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4_RE = 0.0;
  double m4_IM = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = 1.0;
      // Reduced matrix element of tau(1)_- tau(2)_-
      m1 *= sqrt(5);
      m1 *= pow(-1.0, J)/sqrt(2.0*J + 1.0);


      m1 *= fact;
      double rm_RE = 0.0;
      double rm_IM = 0.0;
      compute_radial_matrix_element_GTJ(n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t12, J, q, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
      m4_RE += m1*rm_RE;
      m4_IM += m1*rm_IM;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}

  *m_RE = m4_RE;
  *m_IM = m4_IM;
 
  return;
}



int get_l(int N, int j) {
  int l;
  if ((N - (j + 1)/2) % 2 == 0) {
    l = (j + 1)/2;
  } else {
    l = (j - 1)/2;
  }
  return l;
}

void filter_spectrum(char* density_file, double J, double T) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
 
  int i_state;
  float energy, excitation, j_tot, t_tot;

  while(fscanf(in_file, "%d %f %f %f %f\n", &i_state, &energy, &excitation, &j_tot, &t_tot) == 5) {
//    if (excitation > 25.0) {continue;} 
    if (MIN(fabs(j_tot - 0.5 - floor(j_tot - 0.5)), fabs(j_tot - 0.5 - ceil(j_tot - 0.5))) > 0.01) {continue;}
    if (MIN(fabs(t_tot - 0.5 - floor(t_tot - 0.5)), fabs(t_tot - 0.5 - ceil(t_tot - 0.5))) > 0.01) {continue;}

//        if (j_tot == 2.5 && t_tot == 2.5) {printf("%d, %d\n", i_state - 1, 1);}
    if (excitation > 25.0) {continue;}
    if (j_tot == 2.5 && t_tot == 2.5) {printf("%d, %g\n", i_state - 1, excitation);}
  }

  return;
}


double compute_total_matrix_element_GT0(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;

  double mat = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4 = compute_matrix_element_GT0(in1p, ij1p, in2p, ij2p, ij12p, it12p, in1, ij1, in2, ij2, ij12, it12); 
    mat += m4*density;
    printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12p, in1, ij1, in2, ij2, ij12, it12, m4*density);

  }

  return mat;
}

double compute_total_matrix_element_GT2(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;

  double mat = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4 = compute_matrix_element_GT2(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12); 
    mat += m4*density;
  //  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12p, in1, ij1, in2, ij2, ij12, it12, m4*density);

  }

  return mat;
}


double compute_matrix_element_GT0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int it12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  double t12 = it12/2.0;
  double t12p = it12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4 = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;
      // Reduced matrix element of tau(1)_- tau(2)_-
      m1 *= sqrt(2.0*t12 + 1.0)*pow(-1.0, 1.0 + t12)*six_j(t12,0.5,0.5,1.0,0.5,0.5)*6.0;

//      m1 *= sqrt(5);
//      m1 *= 6.0*sqrt((2.0*t12 + 1)*(2.0*t12p + 1))*nine_j(0.5, 0.5, 1, 0.5,0.5, 1, t12p, t12, 1);
      double anti_symm = 0.0;
      if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
      if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1 + l2 + lambda + s + 1);}
      if (anti_symm == 0) {continue;}
      m1 *= anti_symm;


      m1 *= fact;
      m4 += m1;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1.0/sqrt(2.0);}

//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, m4);

  return m4;
}

double compute_matrix_element_F1(int ina, int ija, int inb, int ijb, int L, double q) {

  double ja = ija/2.0;
  double jb = ijb/2.0;

  int la, lb;
   
  la = get_l(ina, ija);
  lb = get_l(inb, ijb);
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double na = (ina - la)/2.0;
  double nb = (inb - lb)/2.0;
  
  // Lambda = Lambdap and S = SP
  double fact = pow(-1, 0.5 + jb + L)*sqrt((2*ja + 1)*(2*jb + 1));
  fact *= six_j(la, ja, 0.5, jb, lb, L);
  if (fact == 0.0) {return 0.0;}
  fact *= sqrt((2*la + 1)*(2*lb + 1)*(2*L + 1)/(4*M_PI))*three_j(la, L, lb, 0, 0, 0);
  if (fact == 0.0) {return 0.0;}
  // Un-reduced matrix element of sig(1) dot sig(2)
  double m1 = fact*compute_potential_double(na, la, nb, lb, L, q, 0, 0);
  m1 *= sqrt(3.0);
  // Reduced matrix element of tau(1)_- tau(2)_-

  return m1;
}


double compute_total_matrix_element_F1_double(char* density_file_i, char* density_file_f, double qt, double kt) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file_i, "r");

  int iJi, iJf, iJn;
  int iTi, iTf, iTn;
  int iJop, iTop;

  double *mat_store; 

  int iFirst = 1;

  while (fscanf(in_file, "    %d   %d   %d   %d   %d   %d\n", &iJn, &iTn, &iJi, &iTi, &iJop, &iTop) == 6) {;
    int L2Min = (int) fabs((iJi - iJn)/2.0);
    int L2Max = (int) ((iJi + iJn)/2.0);
    int numL = L2Max - L2Min + 1;
    int Jop = iJop/2;
  //  printf("Start: %d %d %d %d %d\n", iJi, iJn, L2Min, L2Max, Jop);
    if (iFirst == 1) {mat_store = malloc(sizeof(double)*numL); iFirst = 0;}
    int toRead;
    if (Jop == 0) {
      toRead = 3;	    
    } else if (Jop == 1) {
      toRead = 7;
    } else if (Jop == 2) {
      toRead = 8;
    } else if (Jop == 3) {
      toRead = 6;
    } else if (Jop == 4) {
      toRead = 3;
    } else if (Jop == 5) {
      toRead = 1;
    } else {toRead = 0;}

    int iNa, iNb, ija, ijb;
    double density;
  
    double mat = 0;
  
    for (int i = 0; i < toRead; i++) {
      fscanf(in_file, "        %d      %d       %d        %d       %lf\n", &iNa, &ija, &iNb, &ijb, &density);
    //  printf("%d %d %d %d %g\n", iNa, ija, iNb, ijb, density);
      double matb = compute_matrix_element_F1(iNa, ija, iNb, ijb, Jop, qt);
      mat += matb*density;
    }
    int junk;
    fscanf(in_file, "       %d\n", &junk);
 //   printf("Mat: %g\n", mat);
    mat_store[Jop - L2Min] = mat;
  }

  fclose(in_file);

  in_file = fopen(density_file_f, "r");
  double total_mat = 0.0;

  while (fscanf(in_file, "    %d   %d   %d   %d   %d   %d\n", &iJf, &iTf, &iJn, &iTn, &iJop, &iTop) == 6) {;
    iJf = 5;
    int L2Min = (int) fabs((iJi - iJn)/2.0);
    int L2Max = (int) ((iJi + iJn)/2.0);
    int numL = L2Max - L2Min + 1;
    int Jop = iJop/2;
    double Jf = iJf/2.0;
 //   double Ji = Jf;
    double Ji = iJi/2.0;
    double Jn = iJn/2.0;
   // printf("Start: %d %d %d %d %d\n", iJf, iJn, L2Min, L2Max, Jop);
    int toRead;
    if (Jop == 0) {
      toRead = 3;	    
    } else if (Jop == 1) {
      toRead = 7;
    } else if (Jop == 2) {
      toRead = 8;
    } else if (Jop == 3) {
      toRead = 6;
    } else if (Jop == 4) {
      toRead = 3;
    } else if (Jop == 5) {
      toRead = 1;
    } else {toRead = 0;}

    int iNa, iNb, ija, ijb;
    double density;
  
    double mat = 0;
  
    for (int i = 0; i < toRead; i++) {
      fscanf(in_file, "        %d      %d       %d        %d       %lf\n", &iNa, &ija, &iNb, &ijb, &density);
    //  printf("%d %d %d %d %g\n", iNa, ija, iNb, ijb, density);
      for (int L2 = L2Min; L2 <= L2Max; L2++) {
        for (int L1 = abs(L2 - Jop); L1 <= L2 + Jop; L1++) {
//	  printf("Store: %d %g\n", L2 - L2Min, mat_store[L2 - L2Min]);
          if (L1 % 2 || L2 % 2) {continue;}
          double m4 = pow(4.0*M_PI, 2.5)*density*compute_matrix_element_F1_double(iNa, ija, iNb, ijb, L1, kt, L2, qt, Jop)*mat_store[L2 - L2Min]*pow(-1.0, L1/2.0 + Jop + Jf - Ji)*clebsch_gordan(L1, L2, Jop, 0, 0, 0)*cg_fact(L1, L2, Jop, Ji, Jf, iJn)*sqrt((2.0*L1 + 1)*(2.0*L2 + 1)/((2.0*Jop + 1)*4*M_PI))*sqrt(2.0*L1 + 1)/(sqrt(2.0*Jn + 1)*sqrt(2.0*Jf + 1));
	  mat += m4;
        }
      }	
    }
    int junk;
    fscanf(in_file, "       %d\n", &junk);
    total_mat += mat;
  }


  return total_mat;
}

double cg_fact(int l1, int l2, int k, double Ji, double Jf, int iJn) {
	double cg = 0.0;
        double Jn = iJn/2.0;
	double Mtot = 0.5;
	for (int m2 = -l2; m2 <= l2; m2++) {
          for (int iMn = -iJn; iMn <= iJn; iMn += 2) {
		  double Mn = iMn/2.0;
		  cg += pow(-1.0, l2 + m2)*clebsch_gordan(l1, l2, k, 0, m2, m2)*clebsch_gordan(l2, Ji, Jn, -m2, Mtot, Mn)*clebsch_gordan(k, Jn, Jf, m2, Mn, Mtot);
          }
	}
     return(cg);
}

double compute_matrix_element_F1_double(int ina, int ija, int inb, int ijb, int L1, double kt, int L2, double qt, int L) {

  double ja = ija/2.0;
  double jb = ijb/2.0;

  int la, lb;
   
  la = get_l(ina, ija);
  lb = get_l(inb, ijb);
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double na = (ina - la)/2.0;
  double nb = (inb - lb)/2.0;
  
  // Lambda = Lambdap and S = SP
  double fact = pow(-1, 0.5 + jb + L)*sqrt((2*ja + 1)*(2*jb + 1));
  fact *= six_j(la, ja, 0.5, jb, lb, L);
  if (fact == 0.0) {return 0.0;}
  fact *= sqrt((2*la + 1)*(2*lb + 1)*(2*L + 1)/(4*M_PI))*three_j(la, L, lb, 0, 0, 0);
  if (fact == 0.0) {return 0.0;}
  // Un-reduced matrix element of sig(1) dot sig(2)
  double m1 = fact*compute_potential_double(na, la, nb, lb, L1, kt, L2, qt);
  m1 *= sqrt(3.0);
  // Reduced matrix element of tau(1)_- tau(2)_-

  return m1;
}


double compute_matrix_element_GT2(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4 = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      if (s != 1) {continue;}
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0)*sqrt(2*j12p + 1.0)*pow(-1.0, lambda + 1 + j12p)*six_j(s, j12p, lambda, j12, s, 2);
      if (fact == 0.0) {continue;}
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = 2.0*sqrt(5.0);
      // Reduced matrix element of tau(1)_- tau(2)_-
      m1 *= sqrt(5);
      double anti_symm = 0.0;
      if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
      if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1 + l2 + lambda + s + 1);}
      if (anti_symm == 0) {continue;}
      m1 *= anti_symm;


      m1 *= fact;
      m4 += m1;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1.0/sqrt(2.0);}

//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, m4);

  return m4;
}


double compute_total_matrix_element_F0(char* density_file) {
  // Computes the two-body nuclear matrix element sigma_1 dot sigma_2 tau_1+ tau_2+ with arbitrary radial function specified by iv
  FILE *in_file;
  in_file = fopen(density_file, "r");
 
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  float density;

  double mat = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4 = compute_matrix_element_F0(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12); 
    mat += m4*density;

  }

  return mat;
}

double compute_matrix_element_F0(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4 = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    if ((lambda < abs(l1p - l2p)) || (lambda > l1p + l2p)) {continue;} 
    int s_max = MIN(lambda + j12, 1);
    int s_min = abs(lambda - j12);
    for (int s = s_min; s <= s_max; s++) {
      // JJ -> LS coupling factors
      double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
      fact *= sqrt((2*lambda + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
      fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
      fact *= nine_j(l1p, l2p, lambda, 0.5, 0.5, s, j1p, j2p, j12p);
      if (fact == 0.0) {continue;}
      // Un-reduce matrix element wrt total J
      fact *= sqrt(2*j12 + 1.0);
      // Un-reduced matrix element of sig(1) dot sig(2)
      double m1 = 1.0;
      // Reduced matrix element of tau(1)_- tau(2)_-
      m1 *= sqrt(5);
      double anti_symm = 0.0;
      if ((n1 == n1p) && (l1 == l1p) && (n2 == n2p) && (l2 == l2p)) {anti_symm = 1.0;}
      if ((n1 == n2p) && (l1 == l2p) && (n2 == n1p) && (l2 == l1p)) {anti_symm += pow(-1.0, t12 + l1 + l2 + lambda + s + 1);}
      if (anti_symm == 0) {continue;}
      m1 *= anti_symm;


      m1 *= fact;
      m4 += m1;
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4 *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4 *= 1.0/sqrt(2.0);}

//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, m4);

  return m4;
}


void generate_bigstick_int_file(int i_model, double q) {
int num_shells;

if (i_model == -5) {
	printf("Selecting gcn2850 model space\n");
	num_shells = 4;
} else if (i_model == -4) {
	printf("Selecting jj44b model space\n");
	num_shells = 4;
} else if (i_model == -3) {
	printf("Selecting JUN45 model space\n");
	num_shells = 4;
} else if (i_model == -2) {
	printf("Selecting Brown-Wildenthal model space\n");
	num_shells = 3;
} else if (i_model == -1) {
	printf("Selecting p + core model space\n");
	num_shells = 3;
} else if (i_model == 0) {
	printf("Selecting sd model space\n");
	num_shells = 3;
} else if (i_model == 1) {
	printf("Selecting GXPF1 model space\n");
	num_shells = 4;
} else if (i_model == 2) {
	printf("Selecting fpg + core model space\n");
	num_shells = 11;
} else if (i_model == 3) {
	printf("Selecting sd + core + extra model space\n");
	num_shells = 7;
} else if (i_model == 4) {
	printf("Selecting kb3g model space\n");
	num_shells = 4;
} else if (i_model == 5) {
	printf("Selecting o16 4hw model space\n");
	num_shells = 21;
}


int *n_shell = (int*) malloc(sizeof(int)*num_shells);
int *j_shell = (int*) malloc(sizeof(int)*num_shells);
int *l_shell = (int*) malloc(sizeof(int)*num_shells);
int *w_shell = (int*) malloc(sizeof(int)*num_shells);
int *i_core = (int*) malloc(sizeof(int)*num_shells);

if (i_model == -5) {

      n_shell[0] = 4;
      j_shell[0] = 9;
      l_shell[0] = 4;
      i_core[0] = 0;

      n_shell[1] = 3;
      j_shell[1] = 5;
      l_shell[1] = 3;
      i_core[1] = 0;

      n_shell[2] = 3;
      j_shell[2] = 3;
      l_shell[2] = 1;
      i_core[2] = 0;

      n_shell[3] = 3;
      j_shell[3] = 1;
      l_shell[3] = 1;
      i_core[3] = 0;


} else if (i_model == -4) {

      n_shell[0] = 3;
      j_shell[0] = 5;
      l_shell[0] = 3;
      i_core[0] = 0;

      n_shell[1] = 3;
      j_shell[1] = 3;
      l_shell[1] = 1;
      i_core[1] = 0;

      n_shell[2] = 3;
      j_shell[2] = 1;
      l_shell[2] = 1;
      i_core[2] = 0;

      n_shell[3] = 4;
      j_shell[3] = 9;
      l_shell[3] = 4;
      i_core[3] = 0;


} else if (i_model == -3) {

      n_shell[0] = 3;
      j_shell[0] = 3;
      l_shell[0] = 1;
      i_core[0] = 0;

      n_shell[1] = 3;
      j_shell[1] = 5;
      l_shell[1] = 3;
      i_core[1] = 0;

      n_shell[2] = 3;
      j_shell[2] = 1;
      l_shell[2] = 1;
      i_core[2] = 0;

      n_shell[3] = 4;
      j_shell[3] = 9;
      l_shell[3] = 4;
      i_core[3] = 0;


} else if (i_model == -2) {
      n_shell[0] = 2;
      j_shell[0] = 1;
      l_shell[0] = 0;
      i_core[0] = 0;

      n_shell[1] = 2;
      j_shell[1] = 3;
      l_shell[1] = 2;
      i_core[1] = 0;

      n_shell[2] = 2;
      j_shell[2] = 5;
      l_shell[2] = 2;
      i_core[2] = 0;

} else if (i_model == -1) {
      n_shell[0] = 1;
      j_shell[0] = 1;
      l_shell[0] = 1;
      i_core[0] = 0;

      n_shell[1] = 1;
      j_shell[1] = 3;
      l_shell[1] = 1;
      i_core[1] = 0;

      n_shell[2] = 0;
      j_shell[2] = 1;
      l_shell[2] = 0;
      i_core[2] = 1;

} else if (i_model == 0) {
      n_shell[0] = 2;
      j_shell[0] = 3;
      l_shell[0] = 2;

      n_shell[1] = 2;
      j_shell[1] = 5;
      l_shell[1] = 2;

      n_shell[2] = 2;
      j_shell[2] = 1;
      l_shell[2] = 0;

} else if (i_model == 1) {

  
n_shell[0] = 3;
j_shell[0] = 7;
l_shell[0] = 3;

n_shell[1] = 3;
j_shell[1] = 3;
l_shell[1] = 1;

n_shell[2] = 3;
j_shell[2] = 5;
l_shell[2] = 3;

n_shell[3] = 3;
j_shell[3] = 1;
l_shell[3] = 1;

} else if (i_model == 2) {

  
n_shell[0] = 3;
j_shell[0] = 5;
l_shell[0] = 3;
i_core[0] = 0;

n_shell[1] = 3;
j_shell[1] = 3;
l_shell[1] = 1;
i_core[1] = 0;

n_shell[2] = 3;
j_shell[2] = 1;
l_shell[2] = 1;
i_core[2] = 0;

n_shell[3] = 4;
j_shell[3] = 9;
l_shell[3] = 4;
i_core[3] = 0;

n_shell[4] = 3;
j_shell[4] = 7;
l_shell[4] = 3;
i_core[4] = 1;

n_shell[5] = 2;
j_shell[5] = 3;
l_shell[5] = 2;
i_core[5] = 1;

n_shell[6] = 2;
j_shell[6] = 1;
l_shell[6] = 0;
 i_core[6] = 1;

n_shell[7] = 2;
j_shell[7] = 5;
l_shell[7] = 2;
i_core[7] = 1;

n_shell[8] = 1;
j_shell[8] = 1;
l_shell[8] = 1;
i_core[8] = 1;

n_shell[9] = 1;
j_shell[9] = 3;
l_shell[9] = 1;
i_core[9] = 1;

n_shell[10] = 0;
j_shell[10] = 1;
l_shell[10] = 0;
i_core[10] = 1;

} else if (i_model == 3) {
      n_shell[0] = 2;
      j_shell[0] = 3;
      l_shell[0] = 2;
      i_core[0] = 0;

      n_shell[1] = 2;
      j_shell[1] = 5;
      l_shell[1] = 2;
      i_core[1] = 0;

      n_shell[2] = 2;
      j_shell[2] = 1;
      l_shell[2] = 0;
      i_core[2] = 0;

      n_shell[3] = 0;
      j_shell[3] = 1;
      l_shell[3] = 0;
      i_core[3] = 0;

      n_shell[4] = 1;
      j_shell[4] = 1;
      l_shell[4] = 1;
      i_core[4] = 0;

      n_shell[5] = 1;
      j_shell[5] = 3;
      l_shell[5] = 1;
      i_core[5] = 0;

      n_shell[6] = 3;
      j_shell[6] = 7;
      l_shell[6] = 3;
      i_core[6] = 1;
} else if (i_model == 4) {

  
n_shell[0] = 3;
j_shell[0] = 7;
l_shell[0] = 3;
i_core[0] = 0;

n_shell[1] = 3;
j_shell[1] = 3;
l_shell[1] = 1;
i_core[1] = 0;

n_shell[2] = 3;
j_shell[2] = 5;
l_shell[2] = 3;
i_core[2] = 0;

n_shell[3] = 3;
j_shell[3] = 1;
l_shell[3] = 1;
i_core[3] = 0;


} else if (i_model == 5) {

  
n_shell[0] = 0;
j_shell[0] = 1;
l_shell[0] = 0;
w_shell[0] = 1;

n_shell[1] = 1;
j_shell[1] = 1;
l_shell[1] = 1;
w_shell[1] = 2;

n_shell[2] = 1;
j_shell[2] = 3;
l_shell[2] = 1;
w_shell[2] = 2;

n_shell[3] = 2;
j_shell[3] = 1;
l_shell[3] = 0;
w_shell[3] = 3;

n_shell[4] = 2;
j_shell[4] = 3;
l_shell[4] = 2;
w_shell[4] = 3;

n_shell[5] = 2;
j_shell[5] = 5;
l_shell[5] = 2;
w_shell[5] = 3;

n_shell[6] = 3;
j_shell[6] = 1;
l_shell[6] = 1;
w_shell[6] = 4;

n_shell[7] = 3;
j_shell[7] = 3;
l_shell[7] = 1;
w_shell[7] = 4;

n_shell[8] = 3;
j_shell[8] = 5;
l_shell[8] = 3;
w_shell[8] = 4;

n_shell[9] = 3;
j_shell[9] = 7;
l_shell[9] = 3;
w_shell[9] = 4;

n_shell[10] = 4;
j_shell[10] = 1;
l_shell[10] = 0;
w_shell[10] = 5;

n_shell[11] = 4;
j_shell[11] = 3;
l_shell[11] = 2;
w_shell[11] = 5;

n_shell[12] = 4;
j_shell[12] = 5;
l_shell[12] = 2;
w_shell[12] = 5;

n_shell[13] = 4;
j_shell[13] = 7;
l_shell[13] = 4;
w_shell[13] = 5;

n_shell[14] = 4;
j_shell[14] = 9;
l_shell[14] = 4;
w_shell[14] = 5;

n_shell[15] = 5;
j_shell[15] = 1;
l_shell[15] = 1;
w_shell[15] = 6;

n_shell[16] = 5;
j_shell[16] = 3;
l_shell[16] = 1;
w_shell[16] = 6;

n_shell[17] = 5;
j_shell[17] = 5;
l_shell[17] = 3;
w_shell[17] = 6;

n_shell[18] = 5;
j_shell[18] = 7;
l_shell[18] = 3;
w_shell[18] = 6;

n_shell[19] = 5;
j_shell[19] = 9;
l_shell[19] = 5;
w_shell[19] = 6;

n_shell[20] = 5;
j_shell[20] = 11;
l_shell[20] = 5;
w_shell[20] = 6;

for (int i = 0; i < 21; i++) {
  i_core[i] = 0;
  }
} 
FILE* out_file;
out_file = fopen("kb3g_J2.int", "w");

for (int ia = 0; ia < num_shells; ia++) {
  int n1p = n_shell[ia];
  int ij1p = j_shell[ia];
  double j1p = ij1p/2.0;
  int l1p = l_shell[ia];
  for (int ib = 0; ib <= ia; ib++) {
    int n2p = n_shell[ib];
    int ij2p = j_shell[ib];
    double j2p = ij2p/2.0;
    int l2p = l_shell[ib];
    int j12p_min = abs(ij1p - ij2p)/2;
    int j12p_max = (ij1p + ij2p)/2;
    for (int ig = 0; ig < num_shells; ig++) {
      int n1 = n_shell[ig];
      int ij1 = j_shell[ig];
      double j1 = ij1/2.0;
      int l1 = l_shell[ig];
      for (int id = 0; id <= ig; id++) {
	int n2 = n_shell[id];
	int ij2 = j_shell[id];
	double j2 = ij2/2.0;
	int l2 = l_shell[id];
	if (pow(-1, l1 + l2) != pow(-1, l1p + l2p)) {continue;}
	int j12_min = abs(ij1 - ij2)/2;
	int j12_max = (ij1 + ij2)/2;
	for (int j12 = MAX(j12_min, j12p_min); j12 <= MIN(j12_max, j12p_max); j12++) {
	  for (int t12 = 0; t12 <= 1; t12++) {
	  if ((ia == ib) && ((j12 + t12) % 2 == 0)) {continue;}
	  if ((ig == id) && ((j12 + t12) % 2 == 0)) {continue;}
	  printf("%d %d %d %d\n", ia, ib, ig, id);
	  double mat =  0.0;
	  
	  if (j1p == j1 && j2p == j2) {
		   mat += 2.0*pow(-1.0, j1 + j2p + j12)*six_j(j12, j2p, j1p, 1, j1, j2)*sqrt(j1*(2*j1 + 1)*(j1 + 1))*sqrt(j2*(2*j2 + 1)*(j2 + 1));
		 //mat += pow(-1.0, 1 + t12)*six_j(t12, 0.5, 0.5, 1, 0.5, 0.5)*6.0;
	  } 

	  if (fabs(mat) > pow(10, -8)) {
	    if (ia != ig || ib != id) {mat *= 0.5;}
	    fprintf(out_file, "  %d  %d  %d  %d    %d  %d  %g\n", ia + 1, ib + 1, ig + 1, id + 1, j12, t12, mat);
	  }
	 }
	}
      }
    }
  }
}

fclose(out_file);

return;
}


void compute_matrix_element_T2JFull(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int mu_rel, int mu_cm, int L, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4_RE = 0.0;
  double m4_IM = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    for (int lambdap = abs(l1p - l2p); lambdap <= l1p + l2p; lambdap++) { 
      int s_max = MIN(lambda + j12, 1);
      s_max = MIN(s_max, lambdap + j12p);
      int s_min = abs(lambda - j12);
      s_min = MAX(s_min, abs(lambdap - j12p));
      for (int s = s_min; s <= s_max; s++) {
	if (s != 1) {continue;}
        // JJ -> LS coupling factors
        double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambdap + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        double m1 = sqrt(2.0*j12 + 1)*sqrt(2.0*j12p + 1)*sqrt(2.0*J + 1)*nine_j(lambdap, lambda, L, s, s, 2, j12p, j12, J);
	m1 *= 2.0*sqrt(5);

        // Reduced matrix element of tau(1)_- tau(2)_-
        m1 *= sqrt(5);


        m1 *= fact;
        double rm_RE = 0.0;
        double rm_IM = 0.0;
        compute_radial_matrix_element_GTJFull(n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, s, t12, mu_rel, mu_cm, L, q, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
        m4_RE += m1*rm_RE;
        m4_IM += m1*rm_IM;
      }
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}

  *m_RE = m4_RE;
  *m_IM = m4_IM;
 
//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, *m_RE);

  return;
}


void compute_matrix_element_GTJFull(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int l, int L, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4_RE = 0.0;
  double m4_IM = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    for (int lambdap = abs(l1p - l2p); lambdap <= l1p + l2p; lambdap++) { 
      int s_max = MIN(lambda + j12, 1);
      s_max = MIN(s_max, lambdap + j12p);
      int s_min = abs(lambda - j12);
      s_min = MAX(s_min, abs(lambdap - j12p));
      for (int s = s_min; s <= s_max; s++) {
        // JJ -> LS coupling factors
        double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambdap + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        double m1 = pow(-1, lambdap + s + j12 + J)*sqrt(2.0*j12 + 1)*sqrt(2.0*j12p + 1)*six_j(lambda, lambdap, J, j12p, j12, s);
	m1 *= pow(-1.0, 1.0 + s)*six_j(s,0.5,0.5,1.0,0.5,0.5)*6.0;

        // Reduced matrix element of tau(1)_- tau(2)_-
        m1 *= sqrt(5);


        m1 *= fact;
        double rm_RE = 0.0;
        double rm_IM = 0.0;
        compute_radial_matrix_element_GTJFull(n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, s, t12, l, L, J, q, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
        m4_RE += m1*rm_RE;
        m4_IM += m1*rm_IM;
      }
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}

  *m_RE = m4_RE;
  *m_IM = m4_IM;
 
//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, *m_RE);

  return;
}


void compute_matrix_element_FJFull(int in1p, int ij1p, int in2p, int ij2p, int ij12p, int in1, int ij1, int in2, int ij2, int ij12, int it12, double q, int l, int L, int J, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *m_RE, double *m_IM) {

  double j1 = ij1/2.0;
  double j2 = ij2/2.0;
  double j12 = ij12/2.0;
  double t12 = it12/2.0;
  double j1p = ij1p/2.0;
  double j2p = ij2p/2.0;
  double j12p = ij12p/2.0;

  int l1, l2, l1p, l2p;
   
  l1p = get_l(in1p, ij1p);
  l2p = get_l(in2p, ij2p);
  l1 = get_l(in1, ij1);
  l2 = get_l(in2, ij2); 
    
  // The N's listed in the input file are energy quanta, we want radial quantum numbers
  double n1 = (in1 - l1)/2.0;
  double n2 = (in2 - l2)/2.0;
  double n1p = (in1p - l1p)/2.0;
  double n2p = (in2p - l2p)/2.0;
  double m4_RE = 0.0;
  double m4_IM = 0.0;
  // Lambda = Lambdap and S = SP
  for (int lambda = abs(l1 - l2); lambda <= (l1 + l2); lambda++) {
    for (int lambdap = abs(l1p - l2p); lambdap <= l1p + l2p; lambdap++) { 
      int s_max = MIN(lambda + j12, 1);
      s_max = MIN(s_max, lambdap + j12p);
      int s_min = abs(lambda - j12);
      s_min = MAX(s_min, abs(lambdap - j12p));
      for (int s = s_min; s <= s_max; s++) {
        // JJ -> LS coupling factors
        double fact = sqrt((2*lambda + 1)*(2*s + 1)*(2*j1 + 1)*(2*j2 + 1));
        fact *= sqrt((2*lambdap + 1)*(2*s + 1)*(2*j1p + 1)*(2*j2p + 1));
        fact *= nine_j(l1, l2, lambda, 0.5, 0.5, s, j1, j2, j12);
        fact *= nine_j(l1p, l2p, lambdap, 0.5, 0.5, s, j1p, j2p, j12p);
        if (fact == 0.0) {continue;}
        double m1 = pow(-1, lambdap + s + j12 + J)*sqrt(2.0*j12 + 1)*sqrt(2.0*j12p + 1)*six_j(lambdap, j12p, s, j12, lambda, J);
        // Reduced matrix element of tau(1)_- tau(2)_-
        m1 *= sqrt(5);


        m1 *= fact;
        double rm_RE = 0.0;
        double rm_IM = 0.0;
        compute_radial_matrix_element_GTJFull(n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, s, t12, l, L, J, q, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
        m4_RE += m1*rm_RE;
        m4_IM += m1*rm_IM;
      }
    }
  }

  if ((n1 == n2) && (j1 == j2) && (l1 == l2)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}
  if ((n1p == n2p) && (j1p == j2p) && (l1p == l2p)) {m4_RE *= 1.0/sqrt(2.0); m4_IM *= 1.0/sqrt(2.0);}

  *m_RE = m4_RE;
  *m_IM = m4_IM;
 
//  printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%g\n", in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12, *m_RE);

  return;
}


int test_suite() {
  int error = 0;
  double tol = pow(10, -6);
  // Check Brody-Moshinsky brackets against values in Table I of Buck & Merchant, Nucl. Phys. A 600 (1996) 387-402
  if (fabs(+1.00000000 - brody_mosh(0, 0, 0, 0, 0, 0, 0, 0, 0)) > tol) {error = 1; printf("BM Error 1: %g\n", brody_mosh(0, 0, 0, 0, 0, 0, 0, 0, 0));} 
  if (fabs(-0.41833000 - brody_mosh(1, 0, 0, 2, 2, 0, 1, 0, 3)) > tol) {error = 1; printf("BM Error 2: %g\n", brody_mosh(1, 0, 0, 2, 2, 0, 1, 0, 3));} 
  if (fabs(+0.50000001 - brody_mosh(0, 5, 0, 1, 6, 0, 1, 0, 5)) > tol) {error = 1; printf("BM Error 3: %g\n", brody_mosh(0, 5, 0, 1, 6, 0, 1, 0, 5));} 
  if (fabs(-0.00000001 - brody_mosh(0, 3, 0, 1, 4, 0, 2, 0, 2)) > tol) {error = 1; printf("BM Error 4: %g\n", brody_mosh(0, 3, 0, 1, 4, 0, 2, 0, 2));} 
  if (fabs(+0.29880716 - brody_mosh(0, 1, 1, 3, 3, 0, 2, 0, 4)) > tol) {error = 1; printf("BM Error 5: %g\n", brody_mosh(0, 1, 1, 3, 3, 0, 2, 0, 4));} 
  if (fabs(+0.11785111 - brody_mosh(0, 2, 0, 5, 4, 0, 2, 0, 5)) > tol) {error = 1; printf("BM Error 6: %g\n", brody_mosh(0, 2, 0, 5, 4, 0, 2, 0, 5));} 
  if (fabs(-0.19394780 - brody_mosh(1, 6, 0, 3, 4, 2, 2, 1, 3)) > tol) {error = 1; printf("BM Error 7: %g\n", brody_mosh(1, 6, 0, 3, 4, 2, 2, 1, 3));} 
  if (fabs(-0.07097762 - brody_mosh(2, 5, 1, 0, 5, 2, 2, 1, 3)) > tol) {error = 1; printf("BM Error 8: %g\n", brody_mosh(2, 5, 1, 0, 5, 2, 2, 1, 3));} 
  if (fabs(-0.13328930 - brody_mosh(4, 2, 0, 2, 2, 2, 2, 1, 4)) > tol) {error = 1; printf("BM Error 9: %g\n", brody_mosh(4, 2, 0, 2, 2, 2, 2, 1, 4));} 
  if (fabs(+0.09959471 - brody_mosh(0, 4, 3, 2, 4, 2, 2, 1, 4)) > tol) {error = 1; printf("BM Error 10: %g\n", brody_mosh(0, 4, 3, 2, 4, 2, 2, 1, 4));} 
  if (error) {printf("BM error\n");} else {printf("Brody-Moshinsky Tests: Pass\n");}

  if (COR_FAC) {
    printf("Two nucleon correlation function is turned ON. Skipping integration tests.\n");}
  else {
    for (int J = 0; J <= 2; J++) {
      for (int iq = 0; iq < 10; iq++) {
        double q = iq*2.0;
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *f_spline = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
        double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
        double *f_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
        double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
        for (int i = 0; i <= NSPLINE; i++) {
          double r = RMIN + i*delta_r;
          r_array[i] = r;
          f_array[i] = v_cm_finite_q(r, J, q);
        }

        gsl_spline_init(f_spline, r_array, f_array, NSPLINE + 1);

        for (int n = 0; n <= 2; n++) {
          for (int l = 0; l <= 2; l++) {
            for (int np = 0; np <= 2; np++) {
              for (int lp = 0; lp <= 2; lp++) {
                double Mt = compute_rel_potential(np, lp, n, l, J, q, 2);
		double Mts = compute_rel_potential_spline(np, lp, n, l, f_spline, acc);
                double MtEx = HOBesselMatCM(np, lp, n, l, J, q);
                if (fabs(Mt - MtEx) > tol) {error = 1; printf("INT Error: %g %g %g\n", q, Mt, MtEx); exit(0);}
                if (fabs(Mts - MtEx) > tol) {error = 1; printf("Spline INT Error: %d %d %d %d %d %g %g %g\n", np, lp, n, l, J, q, Mts, MtEx); exit(0);}
	      }
	    }
	  }
	}
      }
    }
  }

  if (error) {printf("Integration error\n");} else {printf("Integration Tests 1: Pass\n");}

  int n1p = 0;
  int n2p = 0;
  int l1p = 0;
  int l2p = 0;
  int lambdap = 0;
  int n1 = 0;
  int n2 = 0;
  int l1 = 0;
  int l2 = 0;
  int lambda = 0;
  int s = 0;
  int t = 1;
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *f_spline_RE = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *f_array_RE = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array_RE[i] = 1.0;
    f_array_IM[i] = 1.0;
  }

  gsl_spline_init(f_spline_RE, r_array, f_array_RE, NSPLINE + 1);
  gsl_spline_init(f_spline_IM, r_array, f_array_IM, NSPLINE + 1);

  double rm_RE, rm_IM;

  compute_radial_matrix_element_GTJ(n1p, l1p, n2p, l2p, lambdap, n1, l1, n2, l2, lambda, s, t, 0, 0, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
  printf("%g, %g\n", 2*M_PI*rm_RE, 2*M_PI*rm_IM); 
  compute_radial_matrix_element_GTJ(1, l1p, n2p, l2p, lambdap, 1, l1, n2, l2, lambda, s, t, 0, 0, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
  printf("%g, %g\n", 4*M_PI*rm_RE, 4*M_PI*rm_IM);
  compute_radial_matrix_element_GTJ(1, 1, n2p, l2p, 1, 1, 1, n2, l2, 1, s, t, 0, 0, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
  printf("%g, %g\n", 4*M_PI*rm_RE, 4*M_PI*rm_IM); 
  compute_radial_matrix_element_GTJ(2, 1, 1, 1, 2, 1, 1, 2, 1, 2, s, t, 0, 0, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
  printf("%g, %g\n", 4*M_PI*rm_RE, 4*M_PI*rm_IM); 
  compute_radial_matrix_element_GTJ(2, 1, 1, 1, 2, 2, 1, 2, 1, 2, s, t, 0, 0, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
  printf("%g, %g\n", 4*M_PI*rm_RE, 4*M_PI*rm_IM); 
  compute_radial_matrix_element_GTJ(1, 1, 1, 1, 2, 2, 1, 2, 1, 2, s, t, 0, 0, f_spline_RE, f_spline_IM, acc, &rm_RE, &rm_IM);
  printf("%g, %g\n", 4*M_PI*rm_RE, 4*M_PI*rm_IM); 

  printf("\n");
  printf("Ge76 Test Results:\n");
  printf("\n");
  printf("    Result  |   Reference\n");
  printf("-------------------------\n");


  double mat_RE_tot = 0;
  double mat_IM_tot = 0;
  double mat_RE, mat_IM;

  FILE *in_file;
  in_file = fopen("density_files/ge76_se76_J0_T2.dens", "r");
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  double density;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_GTJ(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, 0, 0, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
      mat_RE = 4.0*M_PI*m4_RE;
      mat_IM = 4.0*M_PI*m4_IM;

      mat_RE_tot += mat_RE*density;
      mat_IM_tot += mat_IM*density;

      double mt0 = compute_matrix_element_GT0(in1p, ij1p, in2p, ij2p, ij12p, it12, in1, ij1, in2, ij2, ij12, it12);
      if (fabs(mt0 - mat_RE) > tol) {printf("GT Error\n"); exit(0);}
    //  printf("%g %g %g\n", mat_RE, mat_IM, mt0);
  }

  fclose(in_file);

  printf("MGT: %+.3f | %+.3f\n", fabs(mat_RE_tot), 1.278);

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    double beta = exp(-1.1*pow(B_OSC*r*sqrt(2.0), 2))*(1.0 - 0.68*pow(B_OSC*r*sqrt(2.0),2.0));

    f_array_RE[i] = v_light_limit(r)*pow(1.0 - beta, 2.0);
    f_array_IM[i] = v_light_limit(r)*pow(1.0 - beta, 2.0);
  }

  gsl_spline_init(f_spline_RE, r_array, f_array_RE, NSPLINE + 1);
  gsl_spline_init(f_spline_IM, r_array, f_array_IM, NSPLINE + 1);

  in_file = fopen("density_files/ge76_se76_J0_T2.dens", "r");
  
  mat_RE_tot = 0;
  mat_IM_tot = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_FJ(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, 0, 0, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
      mat_RE = 4.0*M_PI*m4_RE;
      mat_IM = 4.0*M_PI*m4_IM;

      mat_RE_tot += mat_RE*density;
      mat_IM_tot += mat_IM*density;


  }

  fclose(in_file);
  printf("M1 : %+.3f | %+.3f\n", mat_RE_tot, 0.129);

  in_file = fopen("density_files/ge76_se76_J0_T2.dens", "r");
  
  mat_RE_tot = 0;
  mat_IM_tot = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_GTJ(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, 0, 0, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
      mat_RE = 4.0*M_PI*m4_RE;
      mat_IM = 4.0*M_PI*m4_IM;

      mat_RE_tot += mat_RE*density;
      mat_IM_tot += mat_IM*density;

//      printf("%g %g %g\n", mat_RE, mat_IM);
  }

  fclose(in_file);
  printf("M2 : %+.3f | %+.3f\n", mat_RE_tot, -0.411);


  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    double beta = exp(-1.1*pow(B_OSC*r*sqrt(2.0), 2))*(1.0 - 0.68*pow(B_OSC*r*sqrt(2.0),2.0));

    f_array_RE[i] = v_light_limit_d(r)*pow(1.0 - beta, 2.0);
    f_array_IM[i] = v_light_limit_d(r)*pow(1.0 - beta, 2.0);
  }

  gsl_spline_init(f_spline_RE, r_array, f_array_RE, NSPLINE + 1);
  gsl_spline_init(f_spline_IM, r_array, f_array_IM, NSPLINE + 1);

  in_file = fopen("density_files/ge76_se76_J0_T2.dens", "r");
  
  mat_RE_tot = 0;
  mat_IM_tot = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_FJ(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, 0, 0, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
      mat_RE = 4.0*M_PI*m4_RE;
      mat_IM = 4.0*M_PI*m4_IM;

      mat_RE_tot += mat_RE*density;
      mat_IM_tot += mat_IM*density;

//      printf("%g %g %g\n", mat_RE, mat_IM);
  }

  fclose(in_file);
  printf("M1': %+.3f | %+.3f\n", mat_RE_tot, 0.149);

  in_file = fopen("density_files/ge76_se76_J0_T2.dens", "r");
  
  mat_RE_tot = 0;
  mat_IM_tot = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_GTJ(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, 0, 0, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
      mat_RE = 4.0*M_PI*m4_RE;
      mat_IM = 4.0*M_PI*m4_IM;

      mat_RE_tot += mat_RE*density;
      mat_IM_tot += mat_IM*density;

//      printf("%g %g %g\n", mat_RE, mat_IM);
  }

  fclose(in_file);
  printf("M2': %.3f | %.3f\n", mat_RE_tot, -0.469);

  in_file = fopen("density_files/ge76_se76_J0_T2.dens", "r");
  double mGT2p = mat_RE_tot; 
  mat_RE_tot = 0;
  mat_IM_tot = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      compute_matrix_element_T2JFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, 0, 2, 0, 2, 0, f_spline_RE, f_spline_IM, acc, &m4_RE, &m4_IM); 
      mat_RE = 4.0*M_PI*m4_RE;
      mat_IM = 4.0*M_PI*m4_IM;

      mat_RE_tot += mat_RE*density;
      mat_IM_tot += mat_IM*density;

//      printf("%g %g %g\n", mat_RE, mat_IM);
  }

  fclose(in_file);
  printf("M3': %.3f | %.3f\n", 1.0/3.0*(mGT2p) + sqrt(8.0*M_PI/15.0)*sqrt(5.0)/(2.0*sqrt(M_PI))*mat_RE_tot, -0.151);

  return error;

}

void perform_closure_analysis() {

  char density_file_i[250];
  char density_file_f[250];
  int i_state, i_zero;
  float excite;

  double Mi = 26981.5;
  double Mn = 26984.3;

  double Ebind = 0.463;

  double total_mat = 0;

  double E_avg = 10.0;
  
  double ti = 0.5;
  double tf = 2.5;
  double ji = 2.5;
  double jf = 2.5;
  double mti = -0.5;
  double mtf = -2.5;
 
  double pe_max = 92.3;
  double tn = 1.5;
  double mtn = -1.5;
  double qt = 0.0;
  //pe_max = 0.0;
  qt = pe_max;
  double alpha = M_MUON - E_avg + Mi - Mn;
  double delta_e = pe_max + E_avg - Mi + Mn;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  gsl_spline *f_spline_RE0 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM0 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE0 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM0 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE2 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM2 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE2 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM2 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE4 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM4 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE4 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM4 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE6 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM6 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE6 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM6 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE8 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM8 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE8 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM8 = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array_RE0[i] = v_light_nu_RE(0, pe_max, alpha, delta_e, r);
    f_array_IM0[i] = v_light_nu_IM(0, pe_max, alpha, r);
    f_array_RE2[i] = v_light_nu_RE(2, pe_max, alpha, delta_e, r);
    f_array_IM2[i] = v_light_nu_IM(2, pe_max, alpha, r);
    f_array_RE4[i] = v_light_nu_RE(4, pe_max, alpha, delta_e, r);
    f_array_IM4[i] = v_light_nu_IM(4, pe_max, alpha, r);
    f_array_RE6[i] = v_light_nu_RE(6, pe_max, alpha, delta_e, r);
    f_array_IM6[i] = v_light_nu_IM(6, pe_max, alpha, r);
    f_array_RE8[i] = v_light_nu_RE(8, pe_max, alpha, delta_e, r);
    f_array_IM8[i] = v_light_nu_IM(8, pe_max, alpha, r);
  }

  gsl_spline_init(f_spline_RE0, r_array, f_array_RE0, NSPLINE + 1);
  gsl_spline_init(f_spline_IM0, r_array, f_array_IM0, NSPLINE + 1);
  gsl_spline_init(f_spline_RE2, r_array, f_array_RE2, NSPLINE + 1);
  gsl_spline_init(f_spline_IM2, r_array, f_array_IM2, NSPLINE + 1);
  gsl_spline_init(f_spline_RE4, r_array, f_array_RE4, NSPLINE + 1);
  gsl_spline_init(f_spline_IM4, r_array, f_array_IM4, NSPLINE + 1);
  gsl_spline_init(f_spline_RE6, r_array, f_array_RE6, NSPLINE + 1);
  gsl_spline_init(f_spline_IM6, r_array, f_array_IM6, NSPLINE + 1);
  gsl_spline_init(f_spline_RE8, r_array, f_array_RE8, NSPLINE + 1);
  gsl_spline_init(f_spline_IM8, r_array, f_array_IM8, NSPLINE + 1);

  FILE *in_file;
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  double density;


  double mat_RE_tot = 0;
  double mat_IM_tot = 0;
  double mat_RE, mat_IM;


  printf("Performing closure sum:\n");
  printf("Computing Fermi closure sum...\n");
  double IFact0 = pow(-1, 2 + tf - ti)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1)/sqrt(2.0*jf + 1);
  
  in_file = fopen("density_files/al27_na27_usda_J0_T2_0_0.dens", "r");
 
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;
      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 0, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact0;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact0;
        }
     }
 
  }
  fclose(in_file);

  printf("MF J = 0: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);
  
  in_file = fopen("density_files/al27_na27_usda_J2_T2_0_0.dens", "r");

  double IFact2 = pow(-1, 2 + tf - ti)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1)*pow(-1, 2 + jf - ji)*clebsch_gordan(2, ji, jf, 0, 0.5, 0.5)/sqrt(2.0*jf + 1);

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 2, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact2;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact2;
        }
     }
  }
  fclose(in_file);

  printf("MF J = 2: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

  in_file = fopen("density_files/al27_na27_usda_J4_T2_0_0.dens", "r");
  double IFact4 = pow(-1, 2 + tf - ti)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1)*pow(-1, 4 + jf - ji)*clebsch_gordan(4, ji, jf, 0, 0.5, 0.5)/sqrt(2.0*jf + 1);


  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 4, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact4;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact4;
        }
     }

  }
  fclose(in_file);

  printf("MF J = 4: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);
/*
  in_file = fopen("density_files/al27_na27_usdb_J0_T2_0_0.dens", "r");

  mat_RE_tot = 0;
  mat_IM_tot = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;
      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 0, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact0;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact0;
        }
     }
 
  }
  fclose(in_file);

  double mF_RE = mat_RE_tot;
  double mF_IM = mat_IM_tot;

  printf("MGT J = 0: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

  in_file = fopen("density_files/al27_na27_usdb_J2_T2_0_0.dens", "r");

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 2, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact2;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact2;
        }
     }
  }
  fclose(in_file);

  printf("MGT J = 2: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

  in_file = fopen("density_files/al27_na27_usdb_J4_T2_0_0.dens", "r");

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 4, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact4;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact4;
        }
     }

  }
  fclose(in_file);

  printf("MGT J = 4: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);
*/

  printf("Begin Closure Check...\n");
  FILE *list_file;
  list_file = fopen("list_files/al27_mg27_na27.list", "r");

  int itn;

  double Pmat_tot = 0.0;

  while(fscanf(list_file, "%d, %f\n", &i_state, &excite) == 2) {
   if (i_state > 10000) {continue;}
   if (excite > 25.0) {continue;}
    sprintf(density_file_i, "../SPEED-DMG/DMG/examples/output/al27_mg27_usdb_0_%d.dens", i_state);
    sprintf(density_file_f, "../SPEED-DMG/DMG/examples/output/mg27_na27_usdb_%d_0.dens", i_state);
//    tn = itn/2.0;
    double kt = pe_max;// - excite;
//    kt = 0.0;
    excite = E_avg; 
    double qt = M_MUON - excite + Mi - Mn;
    
    double IFact_i = pow(-1, 1 - ti + tn)*clebsch_gordan(1, ti, tn, -1, mti, mtn)/sqrt(2.0*tn + 1);
    double IFact_f = pow(-1, 1 - tn + tf)*clebsch_gordan(1, tn, tf, -1, mtn, mtf)/sqrt(2.0*tf + 1);
    double mat = compute_total_matrix_element_F1_double(density_file_i, density_file_f, qt, kt)*IFact_i*IFact_f*qt*RNUC/HBARC/(4*M_PI*M_PI)*(-M_PI);
    total_mat += mat;
    double Pmat = 0.0;
    double deltaQ = 3.1;
    for (int iq = 0; iq < 200; iq++) {
      double qP = iq*deltaQ;
      Pmat += compute_total_matrix_element_F1_double(density_file_i, density_file_f, qP, kt)*IFact_i*IFact_f*RNUC/HBARC/(4*M_PI*M_PI)*qP/(qP - M_MUON + excite - Mi + Mn)*deltaQ;
  //    printf("Progress: %g %g\n", qP, Pmat);
    }
    Pmat_tot += Pmat;

    printf("%d %g %g %g %g %g\n", i_state, excite, mat, total_mat, Pmat, Pmat_tot);
  }
  printf("Total: %.20f\n", total_mat);

  return;
}

void compute_mu2p_transitions() {
  printf("\n");
  printf("Begin Mu2p Calculations:\n");
  printf("\n");
  printf("Computing Al27(gs) -> Na27(gs):\n");

  printf("    Light neutrino exchange: \n");
  printf("-----------------------------\n");
  printf("SM Interaction: BW84\n");

  char density_file_i[250];
  char density_file_f[250];
  int i_state, i_zero;
  float excite;

  double Mi = 26981.5;
  double Mn = 26984.3;

  double Ebind = 0.463;

  double total_mat = 0;

  double E_avg = 10.0;
  
  double ti = 0.5;
  double tf = 2.5;
  double ji = 2.5;
  double jf = 2.5;
  double mti = -0.5;
  double mtf = -2.5;
 
  double pe_max = 94.5;
  double tn = 1.5;
  double mtn = -1.5;
  double qt = 0.0;
  //pe_max = 0.0;

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  double delta_r = (RMAX - RMIN)/(1.0*NSPLINE);
  double *r_array = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  gsl_spline *f_spline_RE0 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM0 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE0 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM0 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE2 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM2 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE2 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM2 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE4 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM4 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE4 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM4 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE6 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM6 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE6 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM6 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  gsl_spline *f_spline_RE8 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  gsl_spline *f_spline_IM8 = gsl_spline_alloc(gsl_interp_cspline, NSPLINE + 1);
  double *f_array_RE8 = (double*) malloc(sizeof(double)*(NSPLINE + 1));
  double *f_array_IM8 = (double*) malloc(sizeof(double)*(NSPLINE + 1));

  FILE *in_file;
  int in1, in2, ij1, ij2, ij12, it12;
  int in1p, in2p, ij1p, ij2p, ij12p, it12p;
  double density;


  double mat_RE_tot = 0;
  double mat_IM_tot = 0;
  double mat_RE, mat_IM;


  FILE *list_file;
  list_file = fopen("al27_na27_usdb_J25.list", "r");

  char density_file[250];
  double m_sq_tot = 0.0; 
  while(fscanf(list_file, "%d, %f\n", &i_state, &excite) == 2) {

    qt = pe_max - excite;;
   
    double alpha = M_MUON - E_avg + Mi - Mn;
    double delta_e = pe_max + E_avg - Mi + Mn;

    mat_RE_tot = 0;
    mat_IM_tot = 0;


    for (int i = 0; i <= NSPLINE; i++) {
      double r = RMIN + i*delta_r;
      r_array[i] = r;
      f_array_RE0[i] = v_light_nu_RE(0, pe_max, alpha, delta_e, r);
      f_array_IM0[i] = v_light_nu_IM(0, pe_max, alpha, r);
      f_array_RE2[i] = v_light_nu_RE(2, pe_max, alpha, delta_e, r);
      f_array_IM2[i] = v_light_nu_IM(2, pe_max, alpha, r);
      f_array_RE4[i] = v_light_nu_RE(4, pe_max, alpha, delta_e, r);
      f_array_IM4[i] = v_light_nu_IM(4, pe_max, alpha, r);
      f_array_RE6[i] = v_light_nu_RE(6, pe_max, alpha, delta_e, r);
      f_array_IM6[i] = v_light_nu_IM(6, pe_max, alpha, r);
      f_array_RE8[i] = v_light_nu_RE(8, pe_max, alpha, delta_e, r);
      f_array_IM8[i] = v_light_nu_IM(8, pe_max, alpha, r);
    }

    gsl_spline_init(f_spline_RE0, r_array, f_array_RE0, NSPLINE + 1);
    gsl_spline_init(f_spline_IM0, r_array, f_array_IM0, NSPLINE + 1);
    gsl_spline_init(f_spline_RE2, r_array, f_array_RE2, NSPLINE + 1);
    gsl_spline_init(f_spline_IM2, r_array, f_array_IM2, NSPLINE + 1);
    gsl_spline_init(f_spline_RE4, r_array, f_array_RE4, NSPLINE + 1);
    gsl_spline_init(f_spline_IM4, r_array, f_array_IM4, NSPLINE + 1);
    gsl_spline_init(f_spline_RE6, r_array, f_array_RE6, NSPLINE + 1);
    gsl_spline_init(f_spline_IM6, r_array, f_array_IM6, NSPLINE + 1);
    gsl_spline_init(f_spline_RE8, r_array, f_array_RE8, NSPLINE + 1);
    gsl_spline_init(f_spline_IM8, r_array, f_array_IM8, NSPLINE + 1);

 
    double IFact0 = pow(-1, 2 + tf - ti)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1)/sqrt(2.0*jf + 1);
  
    sprintf(density_file, "../SPEED-DMG/DMG/examples/output/al27_na27_usdb_J0_T2_0_%d.dens", i_state);
    in_file = fopen(density_file, "r");
 
    while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

      if (fabs(density) < pow(10, -8)) {continue;}
        double m4_RE = 0.0;
        double m4_IM = 0.0;
        double mat_RE2 = 0.0;
        double mat_IM2 = 0.0;
        for (int l1 = 0; l1 <= 8; l1 += 2) {
          for (int l2 = 0; l2 <= 8; l2 += 2) {
            double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 0, 0, 0, 0);
            if (cg_fact_l1l2 == 0.0) {continue;}
            if (l1 == 0) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 2) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 4) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 6) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 8) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
            }
            mat_RE = 4.0*M_PI*m4_RE;
            mat_IM = 4.0*M_PI*m4_IM;
            mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact0;
            mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact0;
          }
       }
 
    }
    fclose(in_file);

    printf("MF J = 0: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

    sprintf(density_file, "../SPEED-DMG/DMG/examples/output/al27_na27_usdb_J2_T2_0_%d.dens", i_state);
    in_file = fopen(density_file, "r");

    double IFact2 = pow(-1, 2 + tf - ti)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1)*pow(-1, 2 + jf - ji)*clebsch_gordan(2, ji, jf, 0, 0.5, 0.5)/sqrt(2.0*jf + 1);

    while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

      if (fabs(density) < pow(10, -8)) {continue;}
        double m4_RE = 0.0;
        double m4_IM = 0.0;
        double mat_RE2 = 0.0;
        double mat_IM2 = 0.0;

        for (int l1 = 0; l1 <= 8; l1 += 2) {
          for (int l2 = 0; l2 <= 8; l2 += 2) {
            double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 2, 0, 0, 0);
	    if (cg_fact_l1l2 == 0.0) {continue;}
	    if (l1 == 0) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 2) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

            } else if (l1 == 4) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 6) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 8) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
            }

            mat_RE = 4.0*M_PI*m4_RE;
            mat_IM = 4.0*M_PI*m4_IM;
            mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact2;
            mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact2;
          }
       }
    }
    fclose(in_file);

    printf("MF J = 2: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

    sprintf(density_file, "../SPEED-DMG/DMG/examples/output/al27_na27_usdb_J4_T2_0_%d.dens", i_state);
    in_file = fopen(density_file, "r");

    double IFact4 = pow(-1, 2 + tf - ti)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1)*pow(-1, 4 + jf - ji)*clebsch_gordan(4, ji, jf, 0, 0.5, 0.5)/sqrt(2.0*jf + 1);


    while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

      if (fabs(density) < pow(10, -8)) {continue;}
        double m4_RE = 0.0;
        double m4_IM = 0.0;
        double mat_RE2 = 0.0;
        double mat_IM2 = 0.0;

        for (int l1 = 0; l1 <= 8; l1 += 2) {
          for (int l2 = 0; l2 <= 8; l2 += 2) {
            double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 4, 0, 0, 0);
            if (cg_fact_l1l2 == 0.0) {continue;}
            if (l1 == 0) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 2) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 4) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 6) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 8) {
              compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
            }

            mat_RE = 4.0*M_PI*m4_RE;
            mat_IM = 4.0*M_PI*m4_IM;
            mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact4;
            mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact4;
          }
       }
    }
    fclose(in_file);

    printf("MF J = 4: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

    sprintf(density_file, "../SPEED-DMG/DMG/examples/output/al27_na27_usdb_J0_T2_0_%d.dens", i_state);
    in_file = fopen(density_file, "r");

    double mF_RE = mat_RE_tot;
    double mF_IM = mat_IM_tot;

    mat_RE_tot = 0;
    mat_IM_tot = 0;

    while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

      if (fabs(density) < pow(10, -8)) {continue;}
        double m4_RE = 0.0;
        double m4_IM = 0.0;
        double mat_RE2 = 0.0;
        double mat_IM2 = 0.0;
        for (int l1 = 0; l1 <= 8; l1 += 2) {
          for (int l2 = 0; l2 <= 8; l2 += 2) {
            double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 0, 0, 0, 0);
	    if (cg_fact_l1l2 == 0.0) {continue;}
	    if (l1 == 0) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 2) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 4) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 6) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 8) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
            }

            mat_RE = 4.0*M_PI*m4_RE;
            mat_IM = 4.0*M_PI*m4_IM;
            mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact0;
            mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact0;
          }
       } 
    }
    fclose(in_file);

    printf("MGT J = 0: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

    sprintf(density_file, "../SPEED-DMG/DMG/examples/output/al27_na27_usdb_J2_T2_0_%d.dens", i_state);
    in_file = fopen(density_file, "r");

    while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

      if (fabs(density) < pow(10, -8)) {continue;}
        double m4_RE = 0.0;
        double m4_IM = 0.0;
        double mat_RE2 = 0.0;
        double mat_IM2 = 0.0;

        for (int l1 = 0; l1 <= 8; l1 += 2) {
          for (int l2 = 0; l2 <= 8; l2 += 2) {
            double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 2, 0, 0, 0);
            if (cg_fact_l1l2 == 0.0) {continue;}
            if (l1 == 0) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 2) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 4) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 6) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 8) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
            }

            mat_RE = 4.0*M_PI*m4_RE;
            mat_IM = 4.0*M_PI*m4_IM;
            mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact2;
            mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact2;
          }
       }
    }
    fclose(in_file);

  printf("MGT J = 2: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

    sprintf(density_file, "../SPEED-DMG/DMG/examples/output/al27_na27_usdb_J4_T2_0_%d.dens", i_state);
    in_file = fopen(density_file, "r");


    while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

      if (fabs(density) < pow(10, -8)) {continue;}
        double m4_RE = 0.0;
        double m4_IM = 0.0;
        double mat_RE2 = 0.0;
        double mat_IM2 = 0.0;

        for (int l1 = 0; l1 <= 8; l1 += 2) {
          for (int l2 = 0; l2 <= 8; l2 += 2) {
            double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 4, 0, 0, 0);
            if (cg_fact_l1l2 == 0.0) {continue;}
	    if (l1 == 0) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 2) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 4) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	    } else if (l1 == 6) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
            } else if (l1 == 8) {
              compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
            }

            mat_RE = 4.0*M_PI*m4_RE;
            mat_IM = 4.0*M_PI*m4_IM;
            mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact4;
            mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact4;
          }
       }
    }
    fclose(in_file);

    printf("MGT J = 4: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);
 
    double mGT_RE = mat_RE_tot;
    double mGT_IM = mat_IM_tot;

    double fVA = 1.0/1.25; // Ratio of fV/fA

    double MSq = 1.0/(2.0*ji + 1.0)*(pow(pow(fVA, 2.0)*mF_RE - mGT_RE, 2.0) + pow(pow(fVA, 2.0)*mF_IM - mGT_IM, 2.0));

    m_sq_tot += MSq;

    printf("Ex: %g MSq: %g MSqTot: %g \n", excite, MSq, m_sq_tot);
  }
  exit(0);
/*
  for (int i = 0; i <= NSPLINE; i++) {
    double r = RMIN + i*delta_r;
    r_array[i] = r;
    f_array_RE0[i] = v_NUC(0, pe_max, r);
    f_array_IM0[i] = 0.0;
    f_array_RE2[i] = v_NUC(2, pe_max, r);
    f_array_IM2[i] = 0.0;
    f_array_RE4[i] = v_NUC(4, pe_max, r);
    f_array_IM4[i] = 0.0;
    f_array_RE6[i] = v_NUC(6, pe_max, r);
    f_array_IM6[i] = 0.0;
    f_array_RE8[i] = v_NUC(8, pe_max, r);
    f_array_IM8[i] = 0.0;
  }

  gsl_spline_init(f_spline_RE0, r_array, f_array_RE0, NSPLINE + 1);
  gsl_spline_init(f_spline_IM0, r_array, f_array_IM0, NSPLINE + 1);
  gsl_spline_init(f_spline_RE2, r_array, f_array_RE2, NSPLINE + 1);
  gsl_spline_init(f_spline_IM2, r_array, f_array_IM2, NSPLINE + 1);
  gsl_spline_init(f_spline_RE4, r_array, f_array_RE4, NSPLINE + 1);
  gsl_spline_init(f_spline_IM4, r_array, f_array_IM4, NSPLINE + 1);
  gsl_spline_init(f_spline_RE6, r_array, f_array_RE6, NSPLINE + 1);
  gsl_spline_init(f_spline_IM6, r_array, f_array_IM6, NSPLINE + 1);
  gsl_spline_init(f_spline_RE8, r_array, f_array_RE8, NSPLINE + 1);
  gsl_spline_init(f_spline_IM8, r_array, f_array_IM8, NSPLINE + 1);

  mat_RE_tot = 0;
  mat_IM_tot = 0;

  in_file = fopen("density_files/ca48_ti48_kb3g_J0_T2_0_0.dens", "r");
 
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;
      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 0, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact0;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact0;
        }
     }
  }
  fclose(in_file);

  printf("MF J = 0: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);
  
  in_file = fopen("density_files/al27_na27_bw84_J2_T2_0_0.dens", "r");

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 2, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact2;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact2;
        }
     }
  }
  fclose(in_file);

  printf("MF J = 2: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

  in_file = fopen("density_files/al27_na27_bw84_J4_T2_0_0.dens", "r");

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 4, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_FJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact4;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact4;
        }
     }

  }
  fclose(in_file);

  printf("MF J = 4: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

  in_file = fopen("density_files/ca48_ti48_kb3g_J0_T2_0_0.dens", "r");

  mF_RE = mat_RE_tot;
  mF_IM = mat_IM_tot;

  mat_RE_tot = 0;
  mat_IM_tot = 0;

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;
      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 0, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 0, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact0;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact0;
        }
     }
 
  }
  fclose(in_file);

  printf("MGT J = 0: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

  in_file = fopen("density_files/al27_na27_bw84_J2_T2_0_0.dens", "r");

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 2, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 2, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact2;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact2;
        }
     }
  }
  fclose(in_file);

  printf("MGT J = 2: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);

  in_file = fopen("density_files/al27_na27_bw84_J4_T2_0_0.dens", "r");

  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%lf\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    if (fabs(density) < pow(10, -8)) {continue;}
      double m4_RE = 0.0;
      double m4_IM = 0.0;
      double mat_RE2 = 0.0;
      double mat_IM2 = 0.0;

      for (int l1 = 0; l1 <= 8; l1 += 2) {
        for (int l2 = 0; l2 <= 8; l2 += 2) {
          double cg_fact_l1l2 = sqrt(2.0*l1 + 1)*sqrt(2.0*l2 + 1)*pow(-1, (l1 + l2)/2)*clebsch_gordan(l1, l2, 4, 0, 0, 0);
	  if (cg_fact_l1l2 == 0.0) {continue;}
	  if (l1 == 0) {
            compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE0, f_spline_IM0, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 2) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE2, f_spline_IM2, acc, &m4_RE, &m4_IM); 

          } else if (l1 == 4) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE4, f_spline_IM4, acc, &m4_RE, &m4_IM); 
	  } else if (l1 == 6) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE6, f_spline_IM6, acc, &m4_RE, &m4_IM); 
          } else if (l1 == 8) {
compute_matrix_element_GTJFull(in1p, ij1p, in2p, ij2p, ij12p, in1, ij1, in2, ij2, ij12, it12, qt, l1, l2, 4, f_spline_RE8, f_spline_IM8, acc, &m4_RE, &m4_IM); 
          }

          mat_RE = 4.0*M_PI*m4_RE;
          mat_IM = 4.0*M_PI*m4_IM;
          mat_RE_tot += mat_RE*density*cg_fact_l1l2*IFact4;
          mat_IM_tot += mat_IM*density*cg_fact_l1l2*IFact4;
        }
     }

  }
  fclose(in_file);

  printf("MGT J = 4: RE: %.6f IM: %.6f\n", mat_RE_tot, mat_IM_tot);
  
  mGT_RE = mat_RE_tot;
  mGT_IM = mat_IM_tot;

  MSq = 1.0/(2.0*ji + 1.0)*(pow(pow(fVA, 2.0)*mF_RE - mGT_RE, 2.0) + pow(pow(fVA, 2.0)*mF_IM - mGT_IM, 2.0));

  printf("MSq: %g\n", MSq);
*/
  gsl_spline_free(f_spline_RE0);
  gsl_spline_free(f_spline_IM0);
  gsl_interp_accel_free(acc);  


  return;
}
