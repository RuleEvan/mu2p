#include "wfn.h"

void mu2e() {

  FILE* in_file;
  int np_mu = 359;
  int np_e = 2000;
  double j_nuc = 5.0/2.0; // nuclear spin

  // Lepton wavefunction units
  double a0e = 52917.77211; // [fm]
  double a0mu = a0e*(M_ELECTRON/M_MUON); // [fm]
  double E_hartree = 27.211386*pow(10, -6); //[MeV]

  // Parameters for numerical integration
  double x_min = 0.01; // [fm]
  double x_max = 15.0; // [fm]
  double tol = 0.00001; // [fm]

  // Arrays to hold lepton wavefunctions
  double *r_arr_mu = (double*) malloc(sizeof(double)*np_mu);
  double *g_arr_mu = (double*) malloc(sizeof(double)*np_mu);
  double *f_arr_mu = (double*) malloc(sizeof(double)*np_mu);
  
  double *r_arr_e = (double*) malloc(sizeof(double)*np_e);
  double *g_arr_e = (double*) malloc(sizeof(double)*np_e);
  double *f_arr_e = (double*) malloc(sizeof(double)*np_e);

  // Read in muon wavefunction
  in_file = fopen("isotope_data/al27/muon_wfn_al27.dat", "r");
  for (int i = 0; i < np_mu; i++) {
    double r, f, g;
    fscanf(in_file, "%lf %lf %lf\n", &r, &g, &f);
    r_arr_mu[i] = a0mu*r;
    g_arr_mu[i] = g/sqrt(a0mu);
    f_arr_mu[i] = f/sqrt(a0mu);
  }
  fclose(in_file);
  
  // Load k = -1 electron wavefunction 
  in_file = fopen("isotope_data/al27/electron_wfn_km1_al27.dat", "r");  
  for (int i = 0; i < np_e; i++) {
    double r, f, g;
    fscanf(in_file, "%lf %lf %lf\n", &r, &g, &f);
    r_arr_e[i] = a0e*r;
    g_arr_e[i] = g/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(E_hartree);
    f_arr_e[i] = f/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(E_hartree);
  }
  fclose(in_file);

  // Spline acceleration
  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc2 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc3 = gsl_interp_accel_alloc();
  gsl_interp_accel *acc4 = gsl_interp_accel_alloc();

  // Declare muon and electron wavefunction splines
  gsl_spline *g_mu = gsl_spline_alloc(gsl_interp_cspline, np_mu);
  gsl_spline *f_mu = gsl_spline_alloc(gsl_interp_cspline, np_mu);
  gsl_spline *g_e_m1 = gsl_spline_alloc(gsl_interp_cspline, np_e);
  gsl_spline *f_e_m1 = gsl_spline_alloc(gsl_interp_cspline, np_e);

  // Generate muon splines
  gsl_spline_init(g_mu, r_arr_mu, g_arr_mu, np_mu);
  gsl_spline_init(f_mu, r_arr_mu, f_arr_mu, np_mu);

  // Generate k = -1 electron splines
  gsl_spline_init(g_e_m1, r_arr_e, g_arr_e, np_e);
  gsl_spline_init(f_e_m1, r_arr_e, f_arr_e, np_e);

  // Load k = +1 electron wavefunction 
  in_file = fopen("isotope_data/al27/electron_wfn_k1_al27.dat", "r");
  for (int i = 0; i < np_e; i++) {
    double r, f, g;
    fscanf(in_file, "%lf %lf %lf\n", &r, &g, &f);
    r_arr_e[i] = a0e*r;
    g_arr_e[i] = g/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(E_hartree);
    f_arr_e[i] = f/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(E_hartree);
  }
  fclose(in_file);

  // Declare and generate k = +1 electron splines
  gsl_spline *g_e_p1 = gsl_spline_alloc(gsl_interp_cspline, np_e);
  gsl_spline *f_e_p1 = gsl_spline_alloc(gsl_interp_cspline, np_e);
  gsl_spline_init(g_e_p1, r_arr_e, g_arr_e, np_e);
  gsl_spline_init(f_e_p1, r_arr_e, f_arr_e, np_e);

  // Array of functions for MJ matrix elements
  double (**m_contact_mat) (double) = calloc(36*6, sizeof(double (*) (double)));

  m_contact_mat[0] = MJ0_contact_1s1_1s1;
  m_contact_mat[7] = MJ0_contact_1p1_1p1;
  m_contact_mat[14] = MJ0_contact_1p3_1p3;
  m_contact_mat[21] = MJ0_contact_2s1_2s1;
  m_contact_mat[28] = MJ0_contact_1d3_1d3;
  m_contact_mat[35] = MJ0_contact_1d5_1d5;

  // Compute coherent J = 0 T = 0 contact operator  
  in_file = fopen("isotope_data/ni58/density/ni58-ni58_core_1bdy_J0_T0_0_0.dens", "r");
  int j_op = 0;
  int t_op = 0;
  float density;
  double mat_w1 = 0.0;
  double mat = 0.0;
  double mat_w4 = 0.0;
  double mat_v_w1 = 0.0;
  double mat_v_w4 = 0.0;
  int in1, ij1, in1p, ij1p;
  int i1, i1p;
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &in1p, &ij1p, &in1, &ij1, &density) == 5) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j1p = ij1p/2.0;
   
    i1 = get_shell_index(in1, ij1);
    i1p = get_shell_index(in1p, ij1p);
    if (m_contact_mat[i1p + 6*i1] != NULL) {
  // Compute coherent J = 0, T = 0 contact operator
      mat_w1 += sqrt(2.0)*density*RombergSplineFunIntegrator(&lep_int_w1, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_m1, f_e_m1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);
      mat_w4 += sqrt(2.0)*density*RombergSplineFunIntegrator(&lep_int_p_w4, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);
      mat_v_w1 += sqrt(2.0)*density*RombergSplineFunIntegrator(&lep_int_v_w1, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_m1, f_e_m1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);
      mat_v_w4 += sqrt(2.0)*density*RombergSplineFunIntegrator(&lep_int_a_w4, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);


    }
  }
  fclose(in_file);
  mat_w1 *= 32.0*M_PI;///sqrt(2.0*(2*j_nuc + 1));
  mat_w4 *= 32.0*M_PI/sqrt(2.0*(2*j_nuc + 1));
  mat_v_w1 *= 32.0*M_PI/sqrt(2.0*(2*j_nuc + 1));
  mat_v_w4 *= 32.0*M_PI/sqrt(2.0*(2*j_nuc + 1));
  printf("T = 0: Is1: %g Is4: %g\n", mat_w1, mat_w4);
  printf("T = 0: Iv1: %g Iv4: %g\n", mat_v_w1, mat_v_w4);

  // Compute non-coherent J = 0, T = 1 contact operator
  in_file = fopen("al27-al27_core_1bdy_J0_T1_0_0.dens", "r");
  j_op = 0;
  t_op = 1.0;
  mat_w1 = 0.0;
  mat_w4 = 0.0;
  mat_v_w1 = 0.0;
  mat_v_w4 = 0.0;
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &in1p, &ij1p, &in1, &ij1, &density) == 5) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j1p = ij1p/2.0;
   
    i1 = get_shell_index(in1, ij1);
    i1p = get_shell_index(in1p, ij1p);
    if (m_contact_mat[i1p + 6*i1] != NULL) {
  // Compute coherent J = 0, T = 0 contact operator
      mat_w1 += sqrt(6.0)*density*RombergSplineFunIntegrator(&lep_int_w1, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_m1, f_e_m1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);
      mat_w4 += sqrt(6.0)*density*RombergSplineFunIntegrator(&lep_int_p_w4, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);
    mat_v_w1 += sqrt(6.0)*density*RombergSplineFunIntegrator(&lep_int_v_w1, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_m1, f_e_m1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);
      mat_v_w4 += sqrt(6.0)*density*RombergSplineFunIntegrator(&lep_int_a_w4, m_contact_mat[i1p + 6*i1], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, 0.0, 30, 0.00001);


    }
  }
  fclose(in_file);
  mat_w1 *= 32.0*M_PI;

//  mat_w1 *= 32.0*M_PI*clebsch_gordan(0.5, 0.5, 1, -0.5, 0.5, 0)/sqrt((2*t_op + 1)*(2*5/2+1));
  mat_w4 *= 32.0*M_PI*clebsch_gordan(0.5, 0.5, 1, -0.5, 0.5, 0)/sqrt((2*t_op + 1)*(2*5/2+1));
  mat_v_w1 *= 32.0*M_PI*clebsch_gordan(0.5, 0.5, 1, -0.5, 0.5, 0)/sqrt((2*t_op + 1)*(2*5/2+1));
  mat_v_w4 *= 32.0*M_PI*clebsch_gordan(0.5, 0.5, 1, -0.5, 0.5, 0)/sqrt((2*t_op + 1)*(2*5/2+1));

  printf("J = 0 T = 1 contact Is1: %g Is4: %g\n", mat_w1, mat_w4);
  printf("J = 0 T = 1 contact Iv1: %g Iv4: %g\n", mat_v_w1, mat_v_w4);



  /* One matrix for each J of matrix elements
     Single Particle states are indexed as
     
     State | index
     ______|______
     1s1/2 |  0
     1p1/2 |  1
     1p3/2 |  2
     2s1/2 |  3
     1d3/2 |  4
     1d5/2 |  5
 */

  int num_j = 6;
  
  double (**sigma_mat) (double) = calloc(36*num_j, sizeof(double (*) (double)));


  // J = 0 Matrix Elements

  sigma_mat[1] = sigmaJ0_1p1_1s1;
  sigma_mat[6] = sigmaJ0_1p1_1s1;

  sigma_mat[9] = sigmaJ0_2s1_1p1;
  sigma_mat[19] = sigmaJ0_2s1_1p1;

  sigma_mat[16] = sigmaJ0_1d3_1p3;
  sigma_mat[26] = sigmaJ0_1d3_1p3;


  // J = 1 Matrix Elements

  sigma_mat[0 + 36] = sigmaJ1_1s1_1s1;

  sigma_mat[7 + 36] = sigmaJ1_1p1_1p1;

  sigma_mat[8 + 36] = sigmaJ1_1p3_1p1;
  sigma_mat[13 + 36] = sigmaJ1_1p3_1p1;

  sigma_mat[14 + 36] = sigmaJ1_1p3_1p3;
 
  sigma_mat[3 + 36] = sigmaJ1_2s1_1s1;
  sigma_mat[18 + 36] = sigmaJ1_2s1_1s1;
 
  sigma_mat[4 + 36] = sigmaJ1_1d3_1s1;
  sigma_mat[24 + 36] = sigmaJ1_1d3_1s1;

  sigma_mat[21 + 36] = sigmaJ1_2s1_2s1;
 
  sigma_mat[22 + 36] = sigmaJ1_1d3_2s1;
  sigma_mat[27 + 36] = sigmaJ1_1d3_2s1;

  sigma_mat[28 + 36] = sigmaJ1_1d3_1d3;
  
  sigma_mat[34 + 36] = sigmaJ1_1d5_1d3;
  sigma_mat[29 + 36] = sigmaJ1_1d5_1d3;

  sigma_mat[35 + 36] = sigmaJ1_1d5_1d5;

// J = 2 matrix elements

  sigma_mat[12 + 2*36] = sigmaJ2_1p3_1s1;
  sigma_mat[2 + 2*36] = sigmaJ2_1p3_1s1;

  sigma_mat[15 + 2*36] = sigmaJ2_2s1_1p3;
  sigma_mat[20 + 2*36] = sigmaJ2_2s1_1p3;

  sigma_mat[25 + 2*36] = sigmaJ2_1d3_1p1;
  sigma_mat[10 + 2*36] = sigmaJ2_1d3_1p1;

  sigma_mat[16 + 2*36] = sigmaJ2_1d3_1p3;
  sigma_mat[26 + 2*36] = sigmaJ2_1d3_1p3;

  sigma_mat[31 + 2*36] = sigmaJ2_1d5_1p1;
  sigma_mat[11 + 2*36] = sigmaJ2_1d5_1p1;

  sigma_mat[17 + 2*36] = sigmaJ2_1d5_1p3;
  sigma_mat[32 + 2*36] = sigmaJ2_1d5_1p3;

// J = 3 matrix elements

  sigma_mat[14 + 3*36] = sigmaJ3_1p3_1p3;

  sigma_mat[30 + 3*36] = sigmaJ3_1d5_1s1;
  sigma_mat[5 + 3*36] = sigmaJ3_1d5_1s1;

  sigma_mat[28 + 3*36] = sigmaJ3_1d3_1d3;

  sigma_mat[33 + 3*36] = sigmaJ3_1d5_2s1;
  sigma_mat[23 + 3*36] = sigmaJ3_1d5_2s1;

  sigma_mat[34 + 3*36] = sigmaJ3_1d5_1d3;
  sigma_mat[29 + 3*36] = sigmaJ3_1d5_1d3;

  sigma_mat[35 + 3*36] = sigmaJ3_1d5_1d5;

// J = 4 matrix elements

   sigma_mat[17 + 4*36] = sigmaJ4_1d5_1p3;
   sigma_mat[32 + 4*36] = sigmaJ4_1d5_1p3;

// J = 5 matrix elements
  
  sigma_mat[35 + 5*36] = sigmaJ5_1d5_1d5;

// Compute one-pion exchange 
  in_file = fopen("al27-al27_core_1bdy_J1_T1_0_0.dens", "r");
  j_op = 1.0;
  t_op = 1.0;
  mat_w1 = 0.0;
  mat_w4 = 0.0;
  while(fscanf(in_file, "%d, %d, %d, %d, %f\n", &in1p, &ij1p, &in1, &ij1, &density) == 5) {
    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j1p = ij1p/2.0;
   
    i1 = get_shell_index(in1, ij1);
    i1p = get_shell_index(in1p, ij1p);
    if (sigma_mat[i1p + 6*i1 + 36] != NULL) {
      mat_w1 += sqrt(6.0)*density*RombergSplineFunIntegrator(&lep_int_s_w3, sigma_mat[i1p + 6*i1 + 36], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, x_min, x_max, tol);
      mat_w4 += sqrt(6.0)*density*RombergSplineFunIntegrator(&lep_int_s_w4, sigma_mat[i1p + 6*i1 + 36], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, x_min, x_max, tol);

    }
  }
  fclose(in_file);

  mat_w1 *= 45.42*G_AXIAL*32.0*M_PI*clebsch_gordan(0.5, 0.5, 1, -0.5, 0.5, 0.0)/sqrt((2*t_op + 1)*(2*j_op + 1)*(2*5/2+1));
  mat_w4 *= 45.42*G_AXIAL*32.0*M_PI*clebsch_gordan(0.5, 0.5, 1, -0.5, 0.5, 0.0)/sqrt((2*t_op + 1)*(2*j_op + 1)*(2*5/2+1));

  printf("1-Pion: J = 1 T = 0: I3: %g I4: %g\n", mat_w1, mat_w4);



// Compute two-body operators

  // Compute J = 0 T = 0 coherent operator
  in_file = fopen("al27-al27_core_J0_T0_0_0.dens", "r");

  j_op = 0;
  mat_w1 = 0.0;
  mat_w4 = 0.0;
  int in2p, ij2p, ij12p, it12p, in2, ij2, ij12, it12;
  int i2, i2p;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
   
    if (t12 != t12p) {continue;}
    i1 = get_shell_index(in1, ij1);
    i2 = get_shell_index(in2, ij2);
    i1p = get_shell_index(in1p, ij1p);
    i2p = get_shell_index(in2p, ij2p);
    // Compute J = 0 operators
    double mat_k1k2_w1 = 0.0;
    double mat_k1k2_w4 = 0.0;
    for (int k1 = 0; k1 < 6; k1++) {
      for (int k2 = 0; k2 < 6; k2++) {
        double cg_fact = clebsch_gordan(k1, k2, j_op, 0.0, 0.0, 0.0);
        if (cg_fact == 0.0) {continue;}
        if (sigma_mat[i1p + i1*6 + k1*36] != NULL && sigma_mat[i2p + i2*6 + k2*36] != NULL) {
          double phase = 1;
          if (i1 > i1p) {phase *= pow(-1.0, j1p - j1);}
          if (i2 > i2p) {phase *= pow(-1.0, j2p - j2);}
          mat_k1k2_w1 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j1, k1, j2p, j2, k2, j12p, j12, 0)*RombergSpline2FunIntegrator(&lep_int2_w1, sigma_mat[i1p + i1*6 + k1*36], sigma_mat[i2p + i2*6 + k2*36], g_mu, f_mu, g_e_m1, f_e_m1, acc1, acc2, acc3, acc4, x_min, x_max, tol);
          mat_k1k2_w4 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j1, k1, j2p, j2, k2, j12p, j12, 0)*RombergSpline2FunIntegrator(&lep_int2_p_w4, sigma_mat[i1p + i1*6 + k1*36], sigma_mat[i2p + i2*6 + k2*36], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, x_min, x_max, tol);

        }
        if (sigma_mat[i1p + i2*6 + k1*36] != NULL && sigma_mat[i2p + i1*6 + k2*36] != NULL) {
          double phase = pow(-1.0, j1 + j2 - j12 - t12);
          if (i2 > i1p) {phase *= pow(-1.0, j1p - j2);}
          if (i1 > i2p) {phase *= pow(-1.0, j2p - j1);}
          mat_k1k2_w1 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j2, k1, j2p, j1, k2, j12p, j12, 0)*RombergSpline2FunIntegrator(&lep_int2_w1, sigma_mat[i1p + i2*6 + k1*36], sigma_mat[i2p + i1*6 + k2*36], g_mu, f_mu, g_e_m1, f_e_m1, acc1, acc2, acc3, acc4, x_min, x_max, tol);
          mat_k1k2_w4 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j2, k1, j2p, j1, k2, j12p, j12, 0)*RombergSpline2FunIntegrator(&lep_int2_p_w4, sigma_mat[i1p + i2*6 + k1*36], sigma_mat[i2p + i1*6 + k2*36], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, x_min, x_max, tol);

        }
      }
    } 
    
   double mat_tot_w1 = density*mat_k1k2_w1*6.0*pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*sqrt(2*j12p + 1.0)*sqrt(2*j12 + 1);
   double mat_tot_w4 = density*mat_k1k2_w4*6.0*pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*sqrt(2*j12p + 1.0)*sqrt(2*j12 + 1);

  //  if (mat != 0.0) {printf("n1p: %d j1p: %g n2p: %d j2p: %g n1: %d j1: %g n2: %d j2: %g\n", in1p, j1p, in2p, j2p, in1, j1, in2, j2);}
   if ((in1 == in2) && (j1 == j2)) {mat_tot_w1 *= 1.0/sqrt(2.0);
                                    mat_tot_w4 *= 1.0/sqrt(2.0);}
   if ((in1p == in2p) && (j1p == j2p)) {mat_tot_w1 *= 1.0/sqrt(2.0);
                                        mat_tot_w4 *= 1.0/sqrt(2.0);}
   mat_w1 += mat_tot_w1;
   mat_w4 += mat_tot_w4;
  }
  fclose(in_file);
  mat_w1 *= 1.0/4.0*64.0*113.06*pow(G_AXIAL, 2.0)/sqrt(4.0*M_PI)*(-1)*clebsch_gordan(0.5, 0.5, 0.0, -0.5, 0.5, 0.0)/sqrt(6.0);
  mat_w4 *= 1.0/4.0*64.0*113.06*pow(G_AXIAL, 2.0)/sqrt(4.0*M_PI)*(-1)*clebsch_gordan(0.5, 0.5, 0.0, -0.5, 0.5, 0.0)/sqrt(6.0);

  printf("Two-body J = 0 T = 0: I1: %g I4: %g\n", mat_w1, mat_w4);
 
  // Load k = +1 electron wavefunction 
  in_file = fopen("dirac_e_k2.dat", "r");
  
  for (int i = 0; i < np_e; i++) {
    double r, f, g;
    fscanf(in_file, "%lf %lf %lf\n", &r, &g, &f);
    r_arr_e[i] = a0e*r;
    g_arr_e[i] = g/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(E_hartree);
    f_arr_e[i] = f/sqrt(a0e)*sqrt(2.0*ALPHA_FS)/sqrt(E_hartree);
  }
  fclose(in_file);

  gsl_spline *g_e_p2 = gsl_spline_alloc(gsl_interp_cspline, np_e);
  gsl_spline *f_e_p2 = gsl_spline_alloc(gsl_interp_cspline, np_e);
  gsl_spline_init(g_e_p2, r_arr_e, g_arr_e, np_e);
  gsl_spline_init(f_e_p2, r_arr_e, f_arr_e, np_e);
 
  // Compute J = 2 T = 0 Two-body operators
  in_file = fopen("al27-al27_core_J2_T0_0_0.dens", "r");

  j_op = 2;
  mat_w1 = 0.0;
  mat_w4 = 0.0;
  while(fscanf(in_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f\n", &in1p, &ij1p, &in2p, &ij2p, &ij12p, &it12p, &in1, &ij1, &in2, &ij2, &ij12, &it12, &density) == 13) {

    // The angular momentum are doubled in the file
    double j1 = ij1/2.0;
    double j2 = ij2/2.0;
    double j12 = ij12/2.0;
    double t12 = it12/2.0;
    double j1p = ij1p/2.0;
    double j2p = ij2p/2.0;
    double j12p = ij12p/2.0;
    double t12p = it12p/2.0;
   
    if (t12 != t12p) {continue;}
    i1 = get_shell_index(in1, ij1);
    i2 = get_shell_index(in2, ij2);
    i1p = get_shell_index(in1p, ij1p);
    i2p = get_shell_index(in2p, ij2p);

    // Compute J = 0 operators
    double mat_k1k2_w1 = 0.0;
    double mat_k1k2_w4 = 0.0;
    for (int k1 = 0; k1 < 6; k1++) {
      for (int k2 = 0; k2 < 6; k2++) {
        double cg_fact = clebsch_gordan(k1, k2, j_op, 0.0, 0.0, 0.0);
        if (cg_fact == 0.0) {continue;}
        if (sigma_mat[i1p + i1*6 + k1*36] != NULL && sigma_mat[i2p + i2*6 + k2*36] != NULL) {
          double phase = 1;
          if (i1 > i1p) {phase *= pow(-1.0, j1p - j1);}
          if (i2 > i2p) {phase *= pow(-1.0, j2p - j2);}
          mat_k1k2_w1 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j1, k1, j2p, j2, k2, j12p, j12, j_op)*sqrt(2.0*j_op + 1.0)*RombergSpline2FunIntegrator(&lep_int2_w5, sigma_mat[i1p + i1*6 + k1*36], sigma_mat[i2p + i2*6 + k2*36], g_mu, f_mu, g_e_p2, f_e_p2, acc1, acc2, acc3, acc4, x_min, x_max, tol);
          mat_k1k2_w4 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j1, k1, j2p, j2, k2, j12p, j12, j_op)*sqrt(2.0*j_op + 1.0)*RombergSpline2FunIntegrator(&lep_int2_p_w4, sigma_mat[i1p + i1*6 + k1*36], sigma_mat[i2p + i2*6 + k2*36], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, x_min, x_max, tol);

        }
        if (sigma_mat[i1p + i2*6 + k1*36] != NULL && sigma_mat[i2p + i1*6 + k2*36] != NULL) {
          double phase = pow(-1.0, j1 + j2 - j12 - t12);
          if (i2 > i1p) {phase *= pow(-1.0, j1p - j2);}
          if (i1 > i2p) {phase *= pow(-1.0, j2p - j1);}
          mat_k1k2_w1 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j2, k1, j2p, j1, k2, j12p, j12, j_op)*sqrt(2.0*j_op + 1)*RombergSpline2FunIntegrator(&lep_int2_w5, sigma_mat[i1p + i2*6 + k1*36], sigma_mat[i2p + i1*6 + k2*36], g_mu, f_mu, g_e_p2, f_e_p2, acc1, acc2, acc3, acc4, x_min, x_max, tol);
          mat_k1k2_w4 += phase*cg_fact*sqrt((2*k1 + 1)*(2*k2 + 1))*nine_j(j1p, j2, k1, j2p, j1, k2, j12p, j12, j_op)*sqrt(2.0*j_op + 1.0)*RombergSpline2FunIntegrator(&lep_int2_p_w4, sigma_mat[i1p + i2*6 + k1*36], sigma_mat[i2p + i1*6 + k2*36], g_mu, f_mu, g_e_p1, f_e_p1, acc1, acc2, acc3, acc4, x_min, x_max, tol);

        }
      }
    } 
    
   double mat_tot_w1 = density*mat_k1k2_w1*6.0*pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*sqrt(2*j12p + 1.0)*sqrt(2*j12 + 1);
   double mat_tot_w4 = density*mat_k1k2_w4*6.0*pow(-1.0, 1.0 + t12)*sqrt(2.0*t12 + 1.0)*six_j(t12, 0.5, 0.5, 1.0, 0.5, 0.5)*sqrt(2*j12p + 1.0)*sqrt(2*j12 + 1);

  //  if (mat != 0.0) {printf("n1p: %d j1p: %g n2p: %d j2p: %g n1: %d j1: %g n2: %d j2: %g\n", in1p, j1p, in2p, j2p, in1, j1, in2, j2);}
   if ((in1 == in2) && (j1 == j2)) {mat_tot_w1 *= 1.0/sqrt(2.0);
                                    mat_tot_w4 *= 1.0/sqrt(2.0);}
   if ((in1p == in2p) && (j1p == j2p)) {mat_tot_w1 *= 1.0/sqrt(2.0);
                                        mat_tot_w4 *= 1.0/sqrt(2.0);}
   mat_w1 += mat_tot_w1;
   mat_w4 += mat_tot_w4;
  }
  fclose(in_file);
  mat_w1 *= -1.0*1.0/4.0*64.0*113.06*pow(G_AXIAL, 2.0)/sqrt(4.0*M_PI)*(-1)*clebsch_gordan(0.5, 0.5, 0.0, -0.5, 0.5, 0.0)/sqrt(6.0*(2.0*j_op + 1));
  mat_w4 *= -1.0*1.0/4.0*64.0*113.06*pow(G_AXIAL, 2.0)/sqrt(4.0*M_PI)*(-1)*clebsch_gordan(0.5, 0.5, 0.0, -0.5, 0.5, 0.0)/sqrt(6.0*(2.0*j_op + 1));

  printf("Two-body J = 2 T = 0: I1: %g I4: %g\n", mat_w1, mat_w4);

  return;
}



double wfn_sq(gsl_spline* g_spline, gsl_spline* f_spline, gsl_interp_accel *acc, double r) {
  double psi2 = pow(gsl_spline_eval(g_spline, r, acc), 2.0) + pow(gsl_spline_eval(f_spline, r, acc), 2.0);

  return psi2;
}

double newlep1(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = 1.0/(3.0*sqrt(2.0))/pow(M_MUON*b_osc(27)/HBARC, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(3.0 - 2*pow(r/b_osc(27) , 2.0), 2.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double core1(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = sqrt(2.0)/pow(M_MUON*b_osc(27)/HBARC, 3.0)*exp(-pow(r/b_osc(27), 2.0))*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double core2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = 2.0*sqrt(2.0)/3.0/pow(M_MUON*b_osc(27)/HBARC, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27), 2.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double core3(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

double glep = 4.0/3.0/pow(M_MUON*b_osc(27)/HBARC, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27), 2.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}


double newlep2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 8.0/(15.0)/pow(M_MUON*b_osc(27)/HBARC, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27) , 4.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double newlep3(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 4.0/(5.0)*sqrt(2.0/3.0)/pow(M_MUON*b_osc(27)/HBARC, 3.0)*exp(-pow(r/b_osc(27), 2.0))*pow(r/b_osc(27) , 4.0)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}



double g1lep(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r, double q) {
  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*gsl_sf_bessel_j0(q*r/HBARC)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));
  return glep;
}

double g1lep2(gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r, double q1, double q2) {
  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*gsl_sf_bessel_j0(q1*r/HBARC)*gsl_sf_bessel_j0(q2*r/HBARC)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));
  return glep;
}

double lep_int_w1(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*f(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double lep_int_v_w1(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*f(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) + gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double lep_int_a_w4(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = -1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*f(r)*(gsl_spline_eval(fe_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(ge_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}


double lep_int_s_w3(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 1.0/(sqrt(6.0*M_PI))*sqrt(M_MUON)*f(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}


double lep_int_s_w4(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 1.0/(2.0*sqrt(3.0*M_PI))*sqrt(M_MUON)*f(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}


double lep_int_p_w4(double (*f) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = -1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*f(r)*(gsl_spline_eval(fe_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) + gsl_spline_eval(ge_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}


double lep_int2_w1(double (*f1) (double), double (*f2) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = 1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*f1(r)*f2(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) - gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}

double lep_int2_w5(double (*f1) (double), double (*f2) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = -1.0/sqrt(1.0*M_PI)*sqrt(M_MUON)*f1(r)*f2(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(gmu_spline, r, acc2) + gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(fmu_spline, r, acc4));

  return glep;
}


double lep_int2_p_w4(double (*f1) (double), double (*f2) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = -1.0/(2.0*sqrt(M_PI))*sqrt(M_MUON)*f1(r)*f2(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(fmu_spline, r, acc2) + gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(gmu_spline, r, acc4));

  return glep;
}


double lep_int2_p_w1(double (*f1) (double), double (*f2) (double), gsl_spline* gmu_spline, gsl_spline* fmu_spline, gsl_spline* ge_spline, gsl_spline* fe_spline, gsl_interp_accel *acc1, gsl_interp_accel *acc2, gsl_interp_accel *acc3, gsl_interp_accel *acc4, double r) {

  double glep = -1.0/(2.0*sqrt(3.0*M_PI))*sqrt(M_MUON)*f1(r)*f2(r)*(gsl_spline_eval(ge_spline, r, acc1)*gsl_spline_eval(fmu_spline, r, acc2) + gsl_spline_eval(fe_spline, r, acc3)*gsl_spline_eval(gmu_spline, r, acc4));

  return glep;
}


double M_op_sd(int ijp, int ij, int j_op, double q, double b) {
  double y = pow(b*q/(2.0*HBARC), 2.0);
  printf("b: %g y: %g\n", b, y);
  double mat = pow(4.0*M_PI, -0.5)*pow(y, (j_op - 2.0)/2.0)*exp(-y);

  if (j_op == 0) {
    if ((ij == 1) && (ijp == 1)) {
      mat *= sqrt(2.0)*(y - 4.0/3.0*y*y + 2.0/3.0*y*y*y);
    } else if ((ij == 3) && (ijp == 3)) {
      mat *= 2.0*(y - 4.0/3.0*y*y + 4.0/15.0*y*y*y);
    } else if ((ij == 5) && (ijp == 5)) {
      mat *= sqrt(6.0)*(y - 4.0/3.0*y*y + 4.0/15.0*y*y*y);
    } else {
      mat = 0.0;
    }
  } else if (j_op == 2) {
    if ((ij == 1) && (ijp == 3)) {
      mat *= sqrt(10.0)*2.0/15.0*(y - 0.5*y*y);
    } else if ((ij == 3) && (ijp == 1)) {
      mat *= -sqrt(10.0)*2.0/15.0*(y - 0.5*y*y);
    } else if ((ij == 3) && (ijp == 3)) {
      mat *= 28.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp == 5) && (ij == 3)) {
      mat *= sqrt(21.0)*4.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp == 3) && (ij == 5)) {
      mat *= -sqrt(21.0)*4.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp == 5) && (ij == 5)) {
      mat *= sqrt(21.0)*8.0/15.0*(-y + 2.0/7.0*y*y);
    } else if ((ijp = 5) && (ij == 1)) {
      mat *= sqrt(15.0)*8.0/15.0*(-y + 0.5*y*y);
    } else {
      mat = 0.0;
    }
  } else if (j_op == 4) {
    if ((ijp == 5) && (ij == 3)) {
      mat *= sqrt(14.0)*8.0/35.0*y;
    } else if ((ijp == 3) && (ij == 5)) {
      mat *= -sqrt(14.0)*8.0/35.0*y;
    } else if ((ijp == 5) && (ij == 5)) {
      mat *= sqrt(7.0)*8.0/35.0*y;
    } else {
      mat *= 0.0;
    }
  } else {
    mat *= 0.0;
  }
  return mat;
}

double MJ0_contact_1s1_1s1(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = sqrt(2.0)/pow(w, 3)*exp(-z*z);

  return mat;
}

double MJ0_contact_1p1_1p1(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = 2.0*sqrt(2.0)/3.0*pow(z, 2)/pow(w, 3)*exp(-z*z);


  return mat;
}

double MJ0_contact_1p3_1p3(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = 4.0/3.0*pow(z, 2)/pow(w, 3)*exp(-z*z);


  return mat;
}


double MJ0_contact_2s1_2s1(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = 1.0/(3.0*sqrt(2.0))/pow(w, 3)*exp(-z*z)*pow(3.0 - 2*z*z, 2.0);


  return mat;
}

double MJ0_contact_1d3_1d3(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = 8.0/(15.0)/pow(w, 3)*exp(-z*z)*pow(z, 4);


  return mat;
}

double MJ0_contact_1d5_1d5(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = 4.0/(5.0)*sqrt(2.0/3.0)/pow(w, 3)*exp(-z*z)*pow(z, 4);


  return mat;
}

double SigmaJ1_contact_2s1_2s1(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = sqrt(2.0)/3.0/pow(w, 4)*z*exp(-z*z)*(21.0 - 20*z*z + 4*pow(z, 4));


  return mat;
}

double SigmaJ1_contact_1d3_2s1(double x) {
  double z = x/b_osc(A_NUC);
  double w = b_osc(A_NUC)*M_MUON/(HBARC);

  double mat = sqrt(2.0/5.0)*2.0/3.0/pow(w, 4)*z*exp(-z*z)*(15.0 - 20*z*z + 4*pow(z, 4));


  return mat;
}

double SigmaJ1_contact_1d3_1d3(double x) {
double z = x/b_osc(A_NUC);
double w = b_osc(A_NUC)*M_MUON/(HBARC);

double mat = 16.0/(15.0*sqrt(5.0))/pow(w, 4)*pow(z, 3.0)*exp(-z*z)*(-5.0 + z*z);


return mat;
}

double SigmaJ1_contact_1d5_1d3(double x) {
double z = x/b_osc(A_NUC);
double w = b_osc(A_NUC)*M_MUON/(HBARC);

double mat = 8.0/(5.0*sqrt(5.0))/pow(w, 4)*pow(z, 3.0)*exp(-z*z)*(5.0 - 2*z*z);


return mat;
}

double SigmaJ1_contact_1d5_1d5(double x) {
double z = x/b_osc(A_NUC);
double w = b_osc(A_NUC)*M_MUON/(HBARC);

double mat = 8.0*sqrt(2.0)/(5.0*sqrt(35.0))/pow(w, 4)*pow(z, 5.0)*exp(-z*z);


return mat;
}


double MJ0_s1_s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 2.0/pow(w, 3.0)*exp(-pow(z/w, 2.0))*(w*w + 2.0*pow(w, 4) + 2.0*z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*(3.0 + 4.0*w*w + 2.0*pow(w, 4))*(exp(2.0*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= -1.0/(24.0*sqrt(2.0));

  return mat;
}

double MJ0_d3_d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/pow(w, 3.0)*exp(-pow(z/w, 2.0))*(7.0*w*w + 2.0*pow(w, 4) + 2.0*z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*(15.0 + 20.0*w*w + 4.0*pow(w, 4))*(exp(2.0*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= -1.0/(120.0);

  return mat;
}

double MJ0_d5_d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/pow(w, 3.0)*exp(-pow(z/w, 2.0))*(7.0*w*w + 2.0*pow(w, 4) + 2.0*z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*(15.0 + 20.0*w*w + 4.0*pow(w, 4))*(exp(2.0*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= -sqrt(6.0)/(240.0);

  return mat;
}

double sigmaJ1_2s1_2s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(4.0*pow(w, 6) + 2.0*pow(w, 8) - w*w*z*z + 2.0*pow(z, 4) + 3.0*pow(w, 4) + 2.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(3.0 + 4.0*w*w + 2.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(2.0));

  return mat;
}

double sigmaJ1_1d3_2s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -2.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(4.0*pow(w, 6) + 2.0*pow(w, 8) - w*w*z*z + 2.0*pow(z, 4) + 2.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(2.0*w*w + pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(12.0*sqrt(10.0));

  return mat;
}

double sigmaJ1_1d3_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(8.0*pow(w, 6) + 4.0*pow(w, 8) - 2.0*w*w*z*z + 4.0*pow(z, 4) - 15.0*pow(w, 4) + 4.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(-15.0 + 8.0*w*w + 4.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(240.0*sqrt(5.0));

  return mat;
}

double sigmaJ1_1d5_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(9.0*pow(w, 6) + 2.0*pow(w, 8) + 4.0*w*w*z*z + 2.0*pow(z, 4) + 5.0*pow(w, 4) + 2.0*pow(w, 4)*z*z);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*(5.0 + 9.0*w*w + 2.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(40.0*sqrt(5.0));

  return mat;
}

double sigmaJ1_1d5_1d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 5))*exp(-pow(z/w, 2.0))*(28.0*pow(w, 6) + 4.0*pow(w, 8) + 18.0*w*w*z*z + 4.0*pow(z, 4) + 35.0*pow(w, 4) + 4.0*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(35.0 + 28.0*w*w + 4.0*pow(w, 4))*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(80.0*sqrt(70.0));

  return mat;
}

double sigmaJ1_1s1_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*w)*exp(-pow(z/w, 2.0));
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(2.0));

  return mat;
}

double sigmaJ1_1p1_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 3))*exp(-pow(z/w, 2.0))*(2.0*pow(w, 4) + 2.0*z*z - w*w);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(2.0*w*w - 1.0)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(2.0));

  return mat;
}

double sigmaJ1_1p3_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(z*pow(w, 3))*exp(-pow(z/w, 2.0))*(pow(w, 4) + z*z + w*w);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*(w*w + 1.0)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(24.0);

  return mat;
}

double sigmaJ1_1p3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*pow(w, 3))*exp(-pow(z/w, 2.0))*(2.0*pow(w, 4) + 2.0*z*z + 5.0*w*w);
  mat += sqrt(M_PI)/(z*z)*exp(w*w)*(2.0*w*w + 5.0)*(exp(2*z)*(2*z - 1.0)*gsl_sf_erfc(w + z/w) + (2*z + 1.0)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(5.0));

  return mat;
}

double sigmaJ0_1p1_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 2.0/(pow(w, 2))*exp(-pow(z/w, 2.0));
  mat += sqrt(M_PI)/(z)*exp(w*w)*w*(exp(2*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(8.0*sqrt(3.0));

  return mat;
}

double sigmaJ2_1p3_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(pow(z*w, 2))*exp(-pow(z/w, 2.0))*(2.0*z*z + 3.0*w*w);
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(-exp(2*z)*(3.0 - 6.0*z + 4.0*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6.0*z + 4.0*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(6.0));

  return mat;
}

double sigmaJ3_1d3_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(2*w*w - 1.0) + 2*pow(w, 4)*(-5 + 8*w*w + 4*pow(w, 4))*z*z + 4*w*w*(2*w*w - 1.0)*pow(z, 4) + 8*pow(z, 6));
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(2*w*w - 1.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(160.0*sqrt(5.0));

  return mat;
}

double sigmaJ3_1d5_2s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(w*w + 2.0) + 2*pow(w, 4)*(2 + w*w)*(5 + 2*w*w)*z*z + 4*w*w*(2 + w*w)*pow(z, 4) + 4*pow(z, 6));
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(w*w + 2.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(15.0));

  return mat;
}

double sigmaJ3_1d5_1d3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(w*w + 2.0) + 2*pow(w, 4)*(2 + w*w)*(5 + 2*w*w)*z*z + 4*w*w*(2 + w*w)*pow(z, 4) + 4*pow(z, 6));
  mat += sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(w*w + 2.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(40.0*sqrt(30.0));

  return mat;
}

double sigmaJ3_1d5_1d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(w*w*pow(z*w, 3))*exp(-pow(z/w, 2.0))*(15.0*pow(w, 6)*(2*w*w + 9.0) + 2*pow(w, 4)*(45 + 28*w*w + 4*pow(w, 4))*z*z + 4*w*w*(9 + 2*w*w)*pow(z, 4) + 8*pow(z, 6));
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(2*w*w + 9.0)*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(240.0*sqrt(5.0));

  return mat;
}

double sigmaJ5_1d5_1d5(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(pow(z*w, 5))*exp(-pow(z/w, 2.0))*(945*pow(w, 8) + 210*pow(w, 6)*(3 + 2*w*w)*z*z + 4*pow(w, 4)*(63 + 28*w*w + 4*pow(w, 4))*pow(z, 4) + 8*w*w*(9 + 2*w*w)*pow(z, 6) + 16*pow(z, 8));
  mat += sqrt(M_PI)/(pow(z, 6))*exp(w*w)*pow(w, 4)*(exp(2*z)*(-945 + 1890*z - 1680*z*z + 840*z*z*z - 240*pow(z, 4) + 32*pow(z, 5))*gsl_sf_erfc(w + z/w) + (945 + 1890*z + 1680*z*z + 840*z*z*z + 240*pow(z, 4) + 32*pow(z, 5))*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(192.0*sqrt(7.0));

  return mat;
}

double sigmaJ0_2s1_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = exp(-pow(z/w, 2.0))/pow(w, 4)*(w*w - 2*pow(w, 4) - 2*z*z);
  mat -= sqrt(M_PI)/z*exp(w*w)*w*(1 + w*w)*(exp(2*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(12.0*sqrt(2.0));

  return mat;
}

double sigmaJ0_1d3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0*exp(-pow(z/w, 2.0))/pow(w, 4)*(w*w + pow(w, 4) +z*z);
  mat += sqrt(M_PI)/z*exp(w*w)*w*(5 + 2*w*w)*(exp(2*z)*gsl_sf_erfc(w + z/w) - exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(12.0*sqrt(10.0));

  return mat;
}

double sigmaJ1_2s1_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0*exp(-pow(z/w, 2.0))/(w*w*w*z)*(pow(w, 4) + z*z);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*w*w*(exp(2*z)*(2*z-1)*gsl_sf_erfc(w + z/w) + exp(-2*z)*(2*z + 1)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(16.0*sqrt(3.0));

  return mat;
}

double sigmaJ1_1d3_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0*exp(-pow(z/w, 2.0))/(w*w*w*z)*(pow(w, 4) + z*z);
  mat -= sqrt(M_PI)/(z*z)*exp(w*w)*w*w*(exp(2*z)*(2*z-1)*gsl_sf_erfc(w + z/w) + exp(-2*z)*(2*z + 1)*gsl_sf_erfc(w - z/w));

  mat *= 1.0/(8.0*sqrt(15.0));

  return mat;
}

double sigmaJ3_1d5_1s1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(pow(z*w, 3))*exp(-pow(z/w, 2.0))*(10*w*w*z*z + 4*pow(z, 4) + 15*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(10.0));

  return mat;
}


double sigmaJ2_2s1_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(3*pow(w, 4)*(1 + w*w) + 2*w*w*(1+ w*w)*z*z + 2*pow(z, 4));
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(1+ w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0);

  return mat;
}

double sigmaJ2_1d3_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(6.0*pow(w, 6) - 2*w*w*z*z + 4*pow(z, 4) - 3*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(-1+ 2*w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(48.0*sqrt(10));

  return mat;
}

double sigmaJ2_1d3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(3*pow(w, 4)*(1 + w*w) + 2*w*w*(1+ w*w)*z*z + 2*pow(z, 4));
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(1+ w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(24.0*sqrt(10));

  return mat;
}

double sigmaJ2_1d5_1p1(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(3*pow(w, 4)*(2 + w*w) + 2*w*w*(2 + w*w)*z*z + 2*pow(z, 4));
  mat -= sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(2 + w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(15.0));

  return mat;
}

double sigmaJ2_1d5_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = -4.0/(z*z*pow(w, 4))*exp(-pow(z/w, 2.0))*(6*pow(w, 6) + 14*w*w*z*z + 4*pow(z, 4) + 21*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat += sqrt(M_PI)/(z*z*z)*exp(w*w)*w*(7 + 2*w*w)*(exp(2*z)*(-3.0 + 6*z -4*z*z)*gsl_sf_erfc(w + z/w) + (3.0 + 6*z + 4*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(16.0*sqrt(210.0));

  return mat;
}

double sigmaJ3_1p3_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(pow(z*w, 3))*exp(-pow(z/w, 2.0))*(10*w*w*z*z + 4*pow(z, 4) + 15*pow(w, 4) + 4*pow(w, 4)*z*z);
  mat -= sqrt(M_PI)/(z*z*z*z)*exp(w*w)*w*w*(exp(2*z)*(-15.0 + 30*z - 24*z*z + 8*z*z*z)*gsl_sf_erfc(w + z/w) + (15 + 30*z + 24*z*z + 8*z*z*z)*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(5.0));

  return mat;
}

double sigmaJ4_1d5_1p3(double x) {
  double z = x*M_PION/(2.0*HBARC);
  double w = b_osc(A_NUC)*M_PION/(2.0*HBARC);

  double mat = 4.0/(pow(z*w, 4))*exp(-pow(z/w, 2.0))*(28*w*w*pow(z, 4) + 8*pow(z, 6) + 5*pow(w, 6)*(21 + 8*z*z) + pow(w, 4)*(70*z*z + 8*pow(z, 4)));
  mat -= sqrt(M_PI)/pow(z, 5)*exp(w*w)*w*w*w*(exp(2*z)*(-105 + 210*z - 180*z*z + 80*z*z*z - 16*pow(z, 4))*gsl_sf_erfc(w + z/w) + (105 + 210*z + 180*z*z + 80*z*z*z + 16*pow(z, 4))*exp(-2*z)*gsl_sf_erfc(w - z/w));

  mat *=1.0/(32.0*sqrt(35.0));

  return mat;
}


