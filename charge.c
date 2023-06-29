#include "charge.h"

void generate_potential_files() {
  double a0e = 52917.7;
  double a0mu = a0e*(M_ELECTRON/M_MUON);
  double E0e = 27.2114*pow(10, -6);
  int n_steps = 17;
  double r_cut;
  int z_atom; 
  FILE* OUTFILE;
  int num_grid_points = 1000;
  int log_r_min = -7;
  int log_r_max = 1;
  float delta_log_r = (float) (log_r_max - log_r_min)/num_grid_points;
  double log_r = log_r_min;
  double *a = (double*) malloc(sizeof(double)*n_steps);
  // Set coefficients for Al27
  z_atom = 13;
  r_cut = 7.0;
  a[0] = 4.3418e-2;
  a[1] = 6.0298e-2;
  a[2] = 2.8950e-3;
  a[3] = -2.3522e-2;
  a[4] = -7.9791e-3;
  a[5] = 2.3010e-3;
  a[6] = 1.0794e-3;
  a[7] = 1.2574e-4;
  a[8] = -1.3021e-4;
  a[9] = 5.6563e-5;
  a[10] = -1.8011e-5;
  a[11] = 4.2869e-6;
  a[12] = 0.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;
  // Write electron potential file
  OUTFILE = fopen("pot_e_al27.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v/E0e);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_al27.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


// Set coefficients for Si28
  z_atom = 14;
  r_cut = 8.0;
  a[0] = 0.33495e-1;
  a[1] = 0.59533e-1;
  a[2] = 0.20979e-1;
  a[3] = -0.16900e-1;
  a[4] = -0.14998e-1;
  a[5] = -0.93248e-3;
  a[6] = 0.33266e-2;
  a[7] = 0.59244e-3;
  a[8] = -0.40013e-3;
  a[9] = 0.12242e-3;
  a[10] = -0.12994e-4;
  a[11] = -0.92784e-5;
  a[12] = 0.72595e-5;
  a[13] = -0.42096e-5;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

  // Write electron potential file
  OUTFILE = fopen("pot_e_si28.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_si28.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for P31
  z_atom = 15;
  r_cut = 8.0;
  a[0] = 0.3530e-1;
  a[1] = 0.59642e-1;
  a[2] = 0.17274e-1;
  a[3] = -0.19303e-1;
  a[4] = -0.13545e-1;
  a[5] = 0.63209e-3;
  a[6] = 0.35462e-2;
  a[7] = 0.83653e-3;
  a[8] = -0.47904e-3;
  a[9] = 0.19099e-3;
  a[10] = -0.69611e-4;
  a[11] = 0.23196e-4;
  a[12] = -0.77780e-5;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_p31.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_p31.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for S32
  z_atom = 16;
  r_cut = 8.0;
  a[0] = 0.37251e-1;
  a[1] = 0.60248e-1;
  a[2] = 0.14748e-1;
  a[3] = -0.18352e-1;
  a[4] = -0.10347e-1;
  a[5] = 0.30461e-2;
  a[6] = 0.35277e-2;
  a[7] = -0.39834e-4;
  a[8] = -0.97177e-4;
  a[9] = 0.92279e-4;
  a[10] = -0.51931e-4;
  a[11] = 0.22958e-4;
  a[12] = -0.86609e-5;
  a[13] = 0.28879e-5;
  a[14] = -0.86632e-6;
  a[14] = -0.86632e-6;
  a[15] = 0.0;
  a[15] = 0.0;


// Write electron potential file
  OUTFILE = fopen("pot_e_s32.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_s32.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for C12
  z_atom = 6;
  r_cut = 8.0;
  a[0] = 0.15721e-1;
  a[1] = 0.38732e-1;
  a[2] = 0.36808e-1;
  a[3] = 0.14671e-1;
  a[4] = -0.43277e-2;
  a[5] = -0.97752e-2;
  a[6] = -0.68908e-2;
  a[7] = -0.27631e-2;
  a[8] = -0.63568e-3;
  a[9] = 0.71809e-4;
  a[10] = 0.18441e-3;
  a[11] = 0.75066e-4;
  a[12] = 0.51069e-4;
  a[13] = 0.14308e-4;
  a[14] = 0.23170e-5;
  a[15] = 0.68465e-6;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_c12.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_c12.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


// Set coefficients for O16
  z_atom = 8;
  r_cut = 8.0;
  a[0] = 0.20238e-1;
  a[1] = 0.44793e-1;
  a[2] = 0.33533e-1;
  a[3] = 0.35030e-2;
  a[4] = -0.12293e-1;
  a[5] = -0.10329e-1;
  a[6] = -0.34036e-2;
  a[7] = -0.41627e-3;
  a[8] = -0.94435e-3;
  a[9] = -0.25771e-3;
  a[10] = 0.23759e-3;
  a[11] = -0.10603e-3;
  a[12] = 0.41480e-4;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_o16.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_o16.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for Ar40
  z_atom = 18;
  r_cut = 9.0;
  a[0] = 0.30451e-1;
  a[1] = 0.55337e-1;
  a[2] = 0.20203e-1;
  a[3] = -0.16765e-1;
  a[4] = -0.13578e-1;
  a[5] = -0.43204e-4;
  a[6] = 0.91988e-3;
  a[7] = -0.41205e-3;
  a[8] = 0.11971e-3;
  a[9] = -0.19801e-4;
  a[10] = -0.43204e-5;
  a[11] = 0.61205e-5;
  a[12] = -0.37803e-5;
  a[13] = 0.18001e-5;
  a[14] = -0.77407e-6;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_ar40.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_ar40.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


// Set coefficients for Ca40
  z_atom = 20;
  r_cut = 8.0;
  a[0] = 0.44846e-1;
  a[1] = 0.61326e-1;
  a[2] = -0.16818e-2;
  a[3] = -0.26217e-1;
  a[4] = -0.29725e-2;
  a[5] = 0.85534e-2;
  a[6] = 0.35322e-2;
  a[7] = -0.48258e-3;
  a[8] = -0.39346e-3;
  a[9] = 0.20338e-3;
  a[10] = 0.25461e-4;
  a[11] = -0.17794e-4;
  a[12] = 0.67394e-5;
  a[13] = -0.21033e-5;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_ca40.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_ca40.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


  // Set coefficients for Ni58
  r_cut = 9.0;
  z_atom = 28;
  a[0] = 0.44880e-1;
  a[1] = 0.64756e-1;
  a[2] = -0.27899e-2;
  a[3] = -0.37016e-1;
  a[4] = -0.71915e-2;
  a[5] = 0.13594e-1;
  a[6] = 0.66331e-2;
  a[7] = -0.14095e-2;
  a[8] = -0.10141e-2;
  a[9] = 0.38616e-3;
  a[10] = -0.13871e-3;
  a[11] = 0.47788e-4;
  a[12] = -0.15295e-4;
  a[13] = 0.59131e-5;
  a[14] = -0.67880e-5;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_ni58.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_ni58.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


// Set coefficients for Ti48
  z_atom = 22;
  r_cut = 10.0;
  a[0] = 0.27850e-1;
  a[1] = 0.55432e-1;
  a[2] = 0.26369e-1;
  a[3] = -0.17091e-1;
  a[4] = -0.21798e-1;
  a[5] = -0.24889e-2;
  a[6] = 0.76631e-2;
  a[7] = 0.34554e-2;
  a[8] = -0.67477e-3;
  a[9] = 0.10764e-3;
  a[10] = -0.16564e-5;
  a[11] = -0.55566e-5;
  a[12] = 0.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_ti48.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_ti48.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


// Set coefficients for Cr52
  z_atom = 24;
  r_cut = 9.0;
  a[0] = 0.39287e-1;
  a[1] = 0.62477e-1;
  a[2] = 0.62482e-2;
  a[3] = -0.32885e-1;
  a[4] = -0.10648e-1;
  a[5] = 0.10520e-1;
  a[6] = 0.85478e-2;
  a[7] = -0.24003e-3;
  a[8] = -0.20499e-2;
  a[9] = -0.12001e-2;
  a[10] = -0.56649e-3;
  a[11] = 0.0;
  a[12] = 0.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_cr52.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_cr52.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for Fe56
  z_atom = 26;
  r_cut = 9.0;
  a[0] = 0.42018e-1;
  a[1] = 0.62337e-1;
  a[2] = 0.23995e-3;
  a[3] = -0.32776e-1;
  a[4] = -0.79941e-2;
  a[5] = 0.10844e-1;
  a[6] = 0.49123e-2;
  a[7] = -0.22144e-2;
  a[8] = -0.18146e-3;
  a[9] = 0.37261e-3;
  a[10] = -0.23296e-3;
  a[11] = 0.11494e-3;
  a[12] = -0.50596e-4;
  a[13] = 0.20652e-4;
  a[14] = -0.79428e-5;
  a[15] = 0.28986e-5;
  a[16] = -0.10075e-5;

// Write electron potential file
  OUTFILE = fopen("pot_e_fe56.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_fe56.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for Co59
  z_atom = 27;
  r_cut = 9.0;
  a[0] = 0.43133e-1;
  a[1] = 0.61249e-1;
  a[2] = -0.32523e-2;
  a[3] = -0.32681e-1;
  a[4] = -0.49583e-2;
  a[5] = 0.11494e-1;
  a[6] = 0.55428e-2;
  a[7] = 0.31398e-3;
  a[8] = -0.70578e-4;
  a[9] = 0.53725e-5;
  a[10] = -0.74650e-6;
  a[11] = 0.19793e-5;
  a[12] = -0.28059e-5;
  a[13] = 0.27183e-5;
  a[14] = -0.19454e-5;
  a[15] = 0.10963e-5;
  a[16] = -0.51114e-6;

// Write electron potential file
  OUTFILE = fopen("pot_e_co59.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_co59.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for Cu63
  z_atom = 29;
  r_cut = 9.0;
  a[0] = 0.45598e-1;
  a[1] = 0.60706e-1;
  a[2] = -0.78616e-2;
  a[3] = -0.31638e-1;
  a[4] = -0.14447e-2;
  a[5] = 0.10953e-1;
  a[6] = 0.42578e-2;
  a[7] = -0.24224e-3;
  a[8] = -0.30067e-3;
  a[9] = 0.23903e-3;
  a[10] = -0.12910e-3;
  a[11] = 0.60195e-4;
  a[12] = -0.25755e-4;
  a[13] = 0.10332e-4;
  a[14] = -0.39330e-5;
  a[15] = 0.14254e-5;
  a[16] = -0.49221e-6;


// Write electron potential file
  OUTFILE = fopen("pot_e_cu63.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_cu63.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


// Set coefficients for Zn64
  z_atom = 30;
  r_cut = 9.0;
  a[0] = 0.47038e-1;
  a[1] = 0.61536e-1;
  a[2] = -0.90045e-2;
  a[3] = -0.30669e-1;
  a[4] = -0.78705e-3;
  a[5] = 0.10034e-1;
  a[6] = 0.14053e-2;
  a[7] = -0.20640e-2;
  a[8] = 0.35105e-3;
  a[9] = 0.27303e-4;
  a[10] = -0.63811e-4;
  a[11] = 0.40893e-4;
  a[12] = -0.20311e-4;
  a[13] = 0.88986e-5;
  a[14] = -0.35849e-5;
  a[15] = 0.13522e-5;
  a[16] = -0.38635e-6;


// Write electron potential file
  OUTFILE = fopen("pot_e_zn64.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_zn64.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

// Set coefficients for Sr88
  z_atom = 38;
  r_cut = 9.0;
  a[0] = 0.56435e-1;
  a[1] = 0.55072e-1;
  a[2] = -0.33363e-1;
  a[3] = -0.26061e-1;
  a[4] = 0.15749e-1;
  a[5] = 0.75233e-2;
  a[6] = -0.55044e-2;
  a[7] = -0.23643e-2;
  a[8] = 0.39362e-3;
  a[9] = -0.22733e-3;
  a[10] = 0.12519e-3;
  a[11] = -0.61176e-4;
  a[12] = 0.27243e-4;
  a[13] = -0.11285e-4;
  a[14] = 0.43997e-5;
  a[15] = -0.16248e-5;
  a[16] = 0.57053e-6;


// Write electron potential file
  OUTFILE = fopen("pot_e_sr88.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_sr88.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


// Set coefficients for Zr90
  z_atom = 40;
  r_cut = 10.0;
  a[0] = 0.46188e-1;
  a[1] = 0.61795e-1;
  a[2] = -0.12315e-1;
  a[3] = -0.36915e-1;
  a[4] = 0.25175e-2;
  a[5] = 0.15234e-1;
  a[6] = -0.55146e-3;
  a[7] = -0.60631e-2;
  a[8] = -0.12198e-2;
  a[9] = 0.36200e-3;
  a[10] = -0.16466e-3;
  a[11] = 0.53305e-4;
  a[12] = -0.50873e-5;
  a[13] = -0.85658e-5;
  a[14] = 0.86095e-5;
  a[15] = 0.0;
  a[16] = 0.0;


// Write electron potential file
  OUTFILE = fopen("pot_e_zr90.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_zr90.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


  // Set coefficients for Ge74
  r_cut = 10.0;
  z_atom = 32;
  a[0] = 0.37989e-1;
  a[1] = 0.58298e-1;
  a[2] = 0.27406e-2;
  a[3] = -0.30666e-1;
  a[4] = -0.81505e-2;
  a[5] = 0.10231e-1;
  a[6] = 0.49382e-2;
  a[7] = -0.16270e-2;
  a[8] = -0.13937e-2;
  a[9] = 0.15476e-3;
  a[10] = 0.14396e-3;
  a[11] = -0.73075e-4;
  a[12] = 0.31998e-4;
  a[13] = -0.12822e-4;
  a[14] = 0.48406e-5;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_ge74.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_ge74.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Set coefficients for Mo98
  r_cut = 12.0;
  z_atom = 42;
  a[0] = 0.30483e-1;
  a[1] = 0.57207e-1;
  a[2] = 0.17888e-1;
  a[3] = -0.28388e-1;
  a[4] = -0.21778e-1;
  a[5] = 0.56780e-2;
  a[6] = 0.11236e-1;
  a[7] = 0.82176e-3;
  a[8] = -0.50390e-2;
  a[9] = -0.23877e-2;
  a[10] = 0.71492e-3;
  a[11] = 0.29839e-3;
  a[12] = -0.31408e-3;
  a[13] = 0.80177e-3;
  a[14] = 0.43682e-4;
  a[15] = -0.51394e-4;
  a[16] = 0.22293e-4;

// Write electron potential file
  OUTFILE = fopen("pot_e_mo98.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_mo98.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Set coefficients for Pd110
  r_cut = 11.0;
  z_atom = 46;
  a[0] = 0.40668e-1;
  a[1] = 0.58793e-1;
  a[2] = -0.61375e-2;
  a[3] = -0.35983e-1;
  a[4] = -0.17447e-2;
  a[5] = 0.14998e-1;
  a[6] = 0.19994e-2;
  a[7] = -0.53170e-2;
  a[8] = -0.14289e-2;
  a[9] = 0.16033e-2;
  a[10] = 0.31574e-3;
  a[11] = -0.42195e-3;
  a[12] = 0.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_pd110.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_pd110.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Set coefficients for Sm152
  r_cut = 10.5;
  z_atom = 62;
  a[0] = 0.56097e-1;
  a[1] = 0.45123e-1;
  a[2] = -0.40310e-1;
  a[3] = -0.18171e-1;
  a[4] = 0.20515e-1;
  a[5] = 0.49023e-2;
  a[6] = -0.67674e-2;
  a[7] = -0.18927e-2;
  a[8] = 0.15333e-2;
  a[9] = 0.0;
  a[10] = 0.0;
  a[11] = 0.0;
  a[12] = 0.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_sm152.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_sm152.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Set coefficients for Er166
  r_cut = 11.0;
  z_atom = 68;
  a[0] = 0.54426e-1;
  a[1] = 0.47165e-1;
  a[2] = -0.38654e-1;
  a[3] = -0.19672e-1;
  a[4] = 0.22092e-1;
  a[5] = 0.78708e-2;
  a[6] = -0.53005e-2;
  a[7] = 0.50005e-3;
  a[8] = 0.52005e-3;
  a[9] = -0.35003e-3;
  a[10] = 0.12001e-3;
  a[11] = 0.0;
  a[12] = 0.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_er166.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_er166.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


  // Set coefficients for Gd158
  r_cut = 10.5;
  z_atom = 64;
  a[0] = 0.57217e-1;
  a[1] = 0.43061e-1;
  a[2] = -0.41996e-1;
  a[3] = -0.17203e-1;
  a[4] = 0.19933e-1;
  a[5] = 0.51060e-2;
  a[6] = -0.73665e-2;
  a[7] = -0.20926e-2;
  a[8] = 0.21883e-2;
  a[9] = 0.0;
  a[10] = 0.0;
  a[11] = 0.0;
  a[12] = 0.0;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_gd158.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_gd158.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Set coefficients for Hg204
  r_cut = 12.0;
  z_atom = 80;
  a[0] = 5.0880e-2;
  a[1] = 5.0679e-2;
  a[2] = -3.9771e-2;
  a[3] = -3.1403e-2;
  a[4] = 2.8120e-2;
  a[5] = 1.0580e-2;
  a[6] = -1.6402e-2;
  a[7] = -3.1958e-3;
  a[8] = 6.9355e-3;
  a[9] = 7.0777e-4;
  a[10] = -1.4961e-4;
  a[11] = 2.4032e-4;
  a[12] = -2.9939e-4;
  a[13] = 2.0003-4;
  a[14] = -1.5870e-4;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_hg204.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_hg204.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Set coefficients for Tl205
  r_cut = 12.0;
  z_atom = 81;
  a[0] = 0.51518e-1;
  a[1] = 0.51165e-1;
  a[2] = -0.39559e-1;
  a[3] = -0.30118e-1;
  a[4] = 0.27600e-1;
  a[5] = 0.10412e-1;
  a[6] = -0.15725e-1;
  a[7] = -0.26546e-2;
  a[8] = 0.70184e-2;
  a[9] = 0.82116e-3;
  a[10] = -0.51805e-3;
  a[11] = 0.32560e-3;
  a[12] = -0.18670e-3;
  a[13] = 0.10202e-3;
  a[14] = -0.53857e-4;
  a[15] = 0.27672e-4;
  a[16] = -0.13873e-4;

// Write electron potential file
  OUTFILE = fopen("pot_e_tl205.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_tl205.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


  // Set coefficients for Pb208
  r_cut = 11.0;
  z_atom = 82;
  a[0] = 0.62732e-1;
  a[1] = 0.38542e-1;
  a[2] = -0.55105e-1;
  a[3] = -0.26990e-2;
  a[4] = 0.31016e-1;
  a[5] = -0.99486e-2;
  a[6] = -0.93012e-2;
  a[7] = 0.76653e-2;
  a[8] = 0.20885e-2;
  a[9] = -0.17840e-2;
  a[10] = 0.74876e-4;
  a[11] = 0.32278e-3;
  a[12] = -0.11353e-3;
  a[13] = 0.0;
  a[14] = 0.0;
  a[15] = 0.0;
  a[16] = 0.0;

// Write electron potential file
  OUTFILE = fopen("pot_e_pb208.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0e, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);

  // Write muon potential file
  OUTFILE = fopen("pot_mu_pb208.dat", "w");
  fprintf(OUTFILE, "%g %g\n", 0.0, 0.0);
  log_r = log_r_min;
  for (int i = 0; i < num_grid_points; i++) {
    double r = pow(10, log_r);
    double v = bessel_potential(n_steps, z_atom, a, r_cut, a0mu, r);
    fprintf(OUTFILE, "%g %g\n", r, r*v);
    log_r += delta_log_r;
  }
  fclose(OUTFILE);


  return;
}

double bessel_potential(int n_steps, int z_atom, double *a, double r_cut, double a0, double r) {
  double v = 0.0;
  if (r*a0 > r_cut) {
    v = -z_atom*ALPHA_FS*HBARC/(r*a0);
  } else {
    double bessel_fact = 0.0;
    for (int i = 0; i < n_steps; i++) {
      int n = i + 1;
      bessel_fact += a[i]*sin(n*M_PI*r*a0/r_cut)/pow(n, 3);
    }
    v = -z_atom*ALPHA_FS*HBARC/r_cut - 4*M_PI*ALPHA_FS*HBARC*pow(r_cut/M_PI, 3.0)/(r*a0)*bessel_fact;
  }
 
  return v;
}


void generate_charge_density_file(wfnData* wd, double b_osc) {

// Generates a file containing the spherically-averaged nuclear
// charge density using the supplied wavefunction data

  FILE *OUT_FILE;
  OUT_FILE = fopen("charge_density.dat", "w");
  double total = 0.0;
  int num_r_steps = 1000;
  int num_theta_steps = 1;
  double delta_theta = M_PI/(1.0*num_theta_steps);
  double r_max = 7.0; //[fm]
  double r_min = 0.001; //[fm]
  double delta_r = (r_max - r_min)/(1.0*num_r_steps);
  double *cache = (double*) malloc(sizeof(double)*wd->n_data);
  double *an = (double*) malloc(sizeof(double)*12);
  for (int i = 0; i < num_r_steps; i++) {
    double r = r_min + i*delta_r;
//    printf("r: %g total: %g\n", r, total);
  //  for (int j = 0; j < num_theta_steps; j++) {
  //    double theta = j*delta_theta;
      double theta = 0.0;
        double rho_valence = generate_proton_density(wd, r, theta, b_osc)/(4.0*M_PI);
        double rho_core =  core_wfn(2, r, b_osc)/(4.0*M_PI);
        double rho = (rho_valence + rho_core);
        for (int n = 1; n <= 12; n++) {
          an[n - 1] += r*r*delta_r*rho*gsl_sf_bessel_j0(n*M_PI*r/r_max)*2.0*pow(n*M_PI, 2.0)/pow(r_max, 3);
        }
        total += r*r*delta_r*rho*4.0*M_PI;
//      total += r*r*delta_r/pow(b_osc, 3)*harmonic_osc_square(1, 0, 0.5, 0.5, 0, r, theta, b_osc, &cache);

//      total += 2.0*M_PI*r*r*delta_r*sin(theta)*delta_theta/pow(b_osc, 3)*generate_proton_density(wd, r, theta, b_osc);
//      total += delta_r*r*r*pow(radial_osc_wfn(2, 1, r, b_osc), 2.0)/pow(b_osc, 3);
    
//    fprintf(OUT_FILE, "%g, %g, %g, %g\n", i/100.0, generate_proton_density(wd, i/100.0, 0.0, b_osc), generate_proton_density(wd, i/100.0, M_PI/4.0, b_osc), generate_proton_density(wd, i/100.0, M_PI/2.0, b_osc));
    fprintf(OUT_FILE, "%g, %g, %g, %g\n", r, rho, rho_core, rho_valence);
  }
  fclose(OUT_FILE);
  printf("total: %g\n", total);
  for (int n = 1; n <= 12; n++) {
    printf("n: %d an: %g\n", n, an[n - 1]);
  }  
  return;
}

double generate_proton_density(wfnData *wd, double r, double theta, double b_osc) {
  double rho = 0.0;

  double *cache = (double*) malloc(sizeof(double)*wd->n_shells);
  for (int i = 0; i < wd->n_shells; i++) {
    cache[i] = -100.204;;
  }

  for (int i_state = 0; i_state < wd->n_states; i_state++) {
    if (fabs(wd->bc[i_state]) < pow(10, -8)) {continue;}
    rho += pow(wd->bc[i_state], 2.0)*anti_proton_wfn(wd, i_state, r, theta, b_osc, &cache);
  }

  free(cache);
  
  return rho;
}    

double core_wfn(int n_max_quanta, double r, double b) {
  double rho = 0.0;
  for (int n = 0; n < n_max_quanta; n++) {
    for (int nr = 0; nr <= n; nr++) {
      for (int l = 0; l <= n; l++) {
        if (2*nr + l != n) {continue;}
        rho += 4.0*M_PI*(2.0*l + 1.0)/(2.0*M_PI)*pow(radial_osc_wfn(n, l, r, b), 2.0);
      }
    }
  }
  return rho;
}



double harmonic_osc_square(int n, int l, double j, double mj, int i_shell, double x, double theta, double b_osc, double **cache) {
//  printf("n: %d l: %d j: %g mj: %g shell: %d\n", n, l, j, mj, i_shell);
  double psi_sq = 0.0;
  if ((*cache)[i_shell] == -100.204) {
    if (mj + 0.5 <= l) {
//      psi_sq += pow(clebsch_gordan(l, 0.5, j, mj + 0.5, -0.5, mj)*radial_osc_wfn(n, l, x, b_osc)*gsl_sf_legendre_sphPlm(l, fabs(mj + 0.5), cos(theta)), 2.0);
      psi_sq += pow(clebsch_gordan(l, 0.5, j, mj + 0.5, -0.5, mj)*radial_osc_wfn(n, l, x, b_osc), 2.0);

    } 
  
    if (mj - 0.5 >= -l) {
   //   psi_sq += pow(clebsch_gordan(l, 0.5, j, mj - 0.5, 0.5, mj)*radial_osc_wfn(n, l, x, b_osc)*gsl_sf_legendre_sphPlm(l, fabs(mj - 0.5), cos(theta)), 2.0);
      psi_sq += pow(clebsch_gordan(l, 0.5, j, mj - 0.5, 0.5, mj)*radial_osc_wfn(n, l, x, b_osc), 2.0);

    }
    (*cache)[i_shell] = psi_sq;
  } else {
    psi_sq = (*cache)[i_shell];
  }
  return psi_sq;
}

double anti_proton_wfn(wfnData *wd, int i_state, double r, double theta, double b_osc, double **cache) {
  double w = 0.0;

  for (int i = 0; i < wd->n_proton; i++) {
    int sp_state = wd->basis[i_state]->label[i];
    w += harmonic_osc_square(wd->n_shell[sp_state], wd->l_shell[sp_state], wd->j_shell[sp_state]/2.0, wd->jz_shell[sp_state]/2.0, sp_state, r, theta, b_osc, cache);
  }
//  w *= 1.0/(wd->n_data);

  return w;
}

