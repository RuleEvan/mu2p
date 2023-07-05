#include "brody.h"
double g_factor(double a, double b, double c) {
  // Computes the G(l',l'',l) coefficient required to compute the BM coefficients
  // See Moshinsky (1959) Eq. 48
  double g = pow(-1.0, b)*sqrt(4.0*M_PI)*sqrt(gsl_sf_gamma(2.0*c + 2.0)/(gsl_sf_gamma(2.0*a + 2.0)*gsl_sf_gamma(2.0*b + 2.0)));
  return g;
}

double h_factor(double a, double b, double c) {
  // Computes the H(l',l'',l) coefficient required to compute the BM coefficients
  // See Moshinsky (1959) Eq. 54
  double h = sqrt(((2.0*a + 1.0)*(2.0*b + 1.0))/(4.0*M_PI*(2.0*c + 1.0)))*clebsch_gordan(a,b,c,0,0,0);
  return h;
}

double brody_mosh_zero(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int l1, int l2) {
  // Computes the n1 = n2 = 0 Brody-Moshinksky brackets
  // See Moshinsky (1959)
  double bm = (2.0*l_rel + 1.0)*(2.0*l_cm + 1.0)*(2.0*l1 + 1.0)*(2.0*l2 + 1.0)/(gsl_sf_gamma(l1 + 1.5)*gsl_sf_gamma(l2 + 1.5)*gsl_sf_gamma(n_rel + 1.0)*gsl_sf_gamma(n_rel + l_rel + 1.5)*gsl_sf_gamma(n_cm + 1.0)*gsl_sf_gamma(n_cm + l_cm + 1.5));
  bm = sqrt(bm);
  bm *= pow(2.0, -0.5*(l1 + l2));
  double f = 0.0;
  int l_rel_1, l_rel_2, l_cm_1, l_cm_2;
  for (l_rel_1 = 0; l_rel_1 <= l1; l_rel_1++) {
    for (l_rel_2 = 0; l_rel_2 <= l2; l_rel_2++) {
      if ((l_rel_1 + l_rel_2) != 2*n_rel + l_rel) {continue;}
      for (l_cm_1 = 0; l_cm_1 <= l1; l_cm_1++) {
        if ((l_rel_1 + l_cm_1) != l1) {continue;}
        for (l_cm_2 = 0; l_cm_2 <= l2; l_cm_2++) {
          if ((l_rel_2 + l_cm_2) != l2) {continue;}
          if ((l_cm_1 + l_cm_2) != 2*n_cm + l_cm) {continue;}
            f += pow(-1.0, l_rel_1)*g_factor(l_cm_1, l_rel_1, l1)*g_factor(l_cm_2, l_rel_2, l2)*h_factor(l_rel_1, l_rel_2, l_rel)*h_factor(l_cm_1, l_cm_2, l_cm)*nine_j(l_rel_1, l_rel_2, l_rel, l_cm_1, l_cm_2, l_cm, l1, l2, l_tot)*pow(-1.0, n_rel)*gsl_sf_gamma(0.5*(l_rel_1 + l_rel_2 - l_rel) + 1.0)*gsl_sf_gamma(0.5*(l_rel_1 + l_rel_2 + l_rel + 3.0))*pow(-1.0, n_cm)*gsl_sf_gamma(0.5*(l_cm_1 + l_cm_2 - l_cm) + 1.0)*gsl_sf_gamma(0.5*(l_cm_1 + l_cm_2 + l_cm + 3.0));
        }
      }
    }
  }
  bm *= f;

  return bm;
} 

double brody_mosh(int n_rel, int l_rel, int n_cm, int l_cm, int l_tot, int n1, int l1, int n2, int l2) {
  // Uses the r^2 recursion relation of Moshinsky to compute the Brody-Moshinsky coefficients for  // arbitrary n1, n2
  // See Moshinsky (1959)
  double bm = 0.0;
  if ((n1 == 0) && (n2 == 0)) {
    bm = brody_mosh_zero(n_rel, l_rel, n_cm, l_cm, l_tot, l1, l2);
    return bm;
  } else if (n2 == 0) {
    if (n_rel > 0) {
      bm += 0.5*sqrt(n_rel*(n_rel + l_rel + 0.5))*brody_mosh(n_rel - 1, l_rel, n_cm, l_cm, l_tot, n1 - 1, l1, n2, l2);
      if (n_cm > 0) {
        bm += sqrt(n_rel*n_cm*(l_rel + 1.0)*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm - 1, l_cm + 1, l_tot, n1 - 1, l1, n2, l2);
      }
      if (l_cm > 0) {
        bm += sqrt(n_rel*(n_cm + l_cm + 0.5)*(l_rel + 1.0)*l_cm)*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm, l_cm - 1, l_tot, n1 - 1, l1, n2, l2);
      }
    }
    if (n_cm > 0) {
      bm += 0.5*sqrt(n_cm*(n_cm + l_cm + 0.5))*brody_mosh(n_rel, l_rel, n_cm - 1, l_cm, l_tot, n1 - 1, l1, n2, l2);
      if (l_rel > 0) {
        bm += sqrt((n_rel + l_rel + 0.5)*n_cm*l_rel*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel - 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm - 1, l_cm + 1, l_tot, n1 - 1, l1, n2, l2);
      }
    }
    if ((l_rel > 0) && (l_cm > 0)) {
      bm += sqrt((n_rel + l_rel + 0.5)*(n_cm + l_cm + 0.5)*l_rel*l_cm)*pow(-1.0, l_tot + l_rel + l_cm)*six_j(l_rel, l_rel - 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm, l_cm - 1, l_tot, n1 - 1, l1, n2, l2);
    }
    bm *= 1.0/sqrt((n1)*(n1 + l1 + 0.5));
  } else { 
   if (n_rel > 0) {
      bm += 0.5*sqrt(n_rel*(n_rel + l_rel + 0.5))*brody_mosh(n_rel - 1, l_rel, n_cm, l_cm, l_tot, n1, l1, n2 - 1, l2);
      if (n_cm > 0) {
        bm -= sqrt(n_rel*n_cm*(l_rel + 1.0)*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm - 1, l_cm + 1, l_tot, n1, l1, n2 - 1, l2);
      }
      if (l_cm > 0) {
        bm -= sqrt(n_rel*(n_cm + l_cm + 0.5)*(l_rel + 1.0)*l_cm)*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel + 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel - 1, l_rel + 1, n_cm, l_cm - 1, l_tot, n1, l1, n2 - 1, l2);
      }
    }
    if (n_cm > 0) {
      bm += 0.5*sqrt(n_cm*(n_cm + l_cm + 0.5))*brody_mosh(n_rel, l_rel, n_cm - 1, l_cm, l_tot, n1, l1, n2 - 1, l2);
      if (l_rel > 0) {
        bm -= sqrt((n_rel + l_rel + 0.5)*n_cm*l_rel*(l_cm + 1.0))*pow(-1.0, l_tot + l_cm + l_rel)*six_j(l_rel, l_rel - 1, 1, l_cm + 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm - 1, l_cm + 1, l_tot, n1, l1, n2 - 1, l2);
      }
    }
    if ((l_rel > 0) && (l_cm > 0)) {
      bm -= sqrt((n_rel + l_rel + 0.5)*(n_cm + l_cm + 0.5)*l_rel*l_cm)*pow(-1.0, l_tot + l_rel + l_cm)*six_j(l_rel, l_rel - 1, 1, l_cm - 1, l_cm, l_tot)*brody_mosh(n_rel, l_rel - 1, n_cm, l_cm - 1, l_tot, n1, l1, n2 - 1, l2);
    }
    bm *= 1.0/sqrt((n2)*(n2 + l2 + 0.5));
  }
  return bm;
}

double compute_radial_matrix_element_scalar(int iv, int n1p, int l1p, int n2p, int l2p, int n1, int l1, int n2, int l2, int lambda, int s, int t) {
  // Computes the matrix element 
  // <n1p l1p n2p l2p lambdap| V(r) | n1 l1 n2 l2 lambda>
  // in terms of Moshinsky brackets and Talmi integrals
  double mat = 0.0;
  int max = 2*n1 + 2*n2 + l1 + l2;
  int maxp = 2*n1p + 2*n2p + l1p + l2p;
  for (int l_cm = 0; l_cm <= max; l_cm++) {
    for (int l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
      double sym = 1.0 + pow(-1.0, l_rel + s + 1 + t);
      if (sym == 0.0) {continue;}
      if (pow(-1.0, l_cm + l_rel) != pow(-1.0, l1 + l2)) {continue;}
      for (int n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
        int n_rel = (max - l_rel - l_cm)/2 - n_cm;
        int n_relp = n_rel + (maxp - max)/2;
        if (n_relp < 0) {continue;}
        double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
        rm *= brody_mosh(n_relp, l_rel, n_cm, l_cm, lambda, n1p, l1p, n2p, l2p);
        if (rm == 0.0) {continue;}
        rm *= sym;
        rm *= compute_potential(n_rel, n_relp, l_rel, l_rel, iv);
        mat += rm;
      }
    }
  } 

  return mat;
}

double compute_radial_matrix_element_y2(int iv, int n1p, int l1p, int n2p, int l2p, int lambdap, int n1, int l1, int n2, int l2, int lambda, int s, int t) {
  // Computes the matrix element
  // <n1p l1p n2p l2p lambdap || Y_2 V(r) ||n1 l1 n2 l2 lambda>
  // in terms of Moshinsky brackets and Talmi integrals
  int max = 2*n1 + 2*n2 + l1 + l2;
  int maxp = 2*n1p + 2*n2p + l1p + l2p;
  double mat = 0.0;
  for (int l_cm = 0; l_cm <= max; l_cm++) {
    for (int l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
      if (pow(-1.0, l_rel + l_cm) != pow(-1.0, l1 + l2)) {continue;}
      double sym = 1.0 + pow(-1.0, l_rel + s + 1 + t);
      if (sym == 0.0) {continue;}
      for (int n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
        int n_rel = (max - l_rel - l_cm)/2 - n_cm;
        for (int l_relp = l_rel - 2; l_relp <= l_rel + 2; l_relp += 2) {
          if (l_relp < 0) {continue;}
          if (pow(-1.0, l_relp + l_cm) != pow(-1.0, l1p + l2p)) {continue;}
          int n_relp = n_rel + (maxp - max + l_rel - l_relp)/2;
          if (n_relp < 0) {continue;}
          double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
          rm *= brody_mosh(n_relp, l_relp, n_cm, l_cm, lambdap, n1p, l1p, n2p, l2p);
          rm *= sqrt(5.0/(4.0*M_PI))*sqrt(2.0*l_relp + 1.0)*sqrt(2.0*l_rel + 1.0)*sqrt(2.0*lambda + 1.0)*sqrt(2.0*lambdap + 1.0)*three_j(l_relp, 2.0, l_rel, 0, 0, 0)*six_j(l_relp, l_rel, 2, lambda, lambdap, l_cm)*pow(-1.0, l_cm);
          if (rm == 0.0) {continue;}
          rm *= sym;
          rm *= compute_potential(n_rel, n_relp, l_rel, l_relp, iv);
          mat += rm;
        }
      }
    }
  }

  return mat;
}


void compute_radial_matrix_element_GTJFull(int n1p, int l1p, int n2p, int l2p, int lambdap, int n1, int l1, int n2, int l2, int lambda, int s, int t, int mu_rel, int mu_cm, int J, double q, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *rm_RE, double *rm_IM) {
  // Computes the matrix element
  // <n1p l1p n2p l2p lambdap || Y_2 V(r) ||n1 l1 n2 l2 lambda>
  // in terms of Moshinsky brackets and Talmi integrals
  int max = 2*n1 + 2*n2 + l1 + l2;
  int maxp = 2*n1p + 2*n2p + l1p + l2p;
  double mat_RE = 0.0;
  double mat_IM = 0.0;
  
  for (int l_cm = 0; l_cm <= max; l_cm++) {
    for (int l_rel = 0; l_rel <= max - l_cm; l_rel ++) {
      if (pow(-1.0, l_rel + l_cm) != pow(-1.0, l1 + l2)) {continue;}
      double sym =  1.0 + pow(-1.0, l_rel + s + 1 + t);
      if (sym == 0.0) {continue;}
      for (int n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
        int n_rel = (max - l_rel - l_cm)/2 - n_cm;
	if (n_rel < 0) {continue;}
	for (int l_cmp = 0; l_cmp <= maxp; l_cmp++) {
          for (int l_relp = 0; l_relp <= maxp - l_cmp; l_relp ++) {
            if (pow(-1.0, l_relp + l_cmp) != pow(-1.0, l1p + l2p)) {continue;}
	    for (int n_cmp = 0; n_cmp <= (maxp - l_cmp - l_relp)/2; n_cmp++) {
            int n_relp = (maxp - l_relp - l_cmp)/2 - n_cmp;
            if (n_relp < 0) {continue;}
 
 	    double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
            rm *= brody_mosh(n_relp, l_relp, n_cmp, l_cmp, lambdap, n1p, l1p, n2p, l2p);
            rm *= 1.0/(4.0*M_PI)*pow(-1.0, l_relp + l_cmp)*sqrt(2.0*J + 1.0)*sqrt(2.0*lambdap + 1.0)*sqrt(2.0*lambda + 1.0)*sqrt((2.0*l_rel + 1.0)*(2.0*l_relp + 1.0)*(2.0*l_cm + 1.0)*(2.0*l_cmp + 1.0)*(2.0*mu_rel + 1.0)*(2.0*mu_cm + 1.0))*three_j(l_relp, mu_rel, l_rel, 0.0, 0.0, 0.0)*three_j(l_cmp, mu_cm, l_cm, 0.0, 0.0, 0.0)*nine_j(l_relp, l_rel, mu_rel, l_cmp, l_cm, mu_cm, lambdap, lambda, J);
            if (rm == 0.0) {continue;}
            rm *= sym;
            double mrel_RE = compute_rel_potential_spline(n_relp, l_relp, n_rel, l_rel, f_spline_RE, acc);
            double mrel_IM = compute_rel_potential_spline(n_relp, l_relp, n_rel, l_rel, f_spline_IM, acc);

	  //  printf("np: %d lp: %d n: %d l: %d mat: %g\n", n_relp, l_relp, n_rel, l_rel, mrel_RE);
	    rm *= compute_rel_potential(n_cmp, l_cmp, n_cm, l_cm, mu_cm, q, 2);
	 //   printf("Np: %d LP: %d N: %d L: %d RCM: %g\n", n_cmp, l_cmp, n_cm, l_cm, compute_rel_potential(n_cmp, l_cmp, n_cm, l_cm, J, q, 2));
            mat_RE += rm*mrel_RE;
	    mat_IM += rm*mrel_IM;
	    }
	  }
        }
      }
    }
  }
  //printf("%d %d %d %d %d %d %d %d %d %d %d %d %g\n", n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t, mat_RE);
  *rm_RE = mat_RE;
  *rm_IM = mat_IM;

  return;
}

void compute_radial_matrix_element_GTJ(int n1p, int l1p, int n2p, int l2p, int lambdap, int n1, int l1, int n2, int l2, int lambda, int s, int t, int J, double q, gsl_spline *f_spline_RE, gsl_spline *f_spline_IM, gsl_interp_accel *acc, double *rm_RE, double *rm_IM) {
  // Computes the matrix element
  // <n1p l1p n2p l2p lambdap || V_rel(r) V_cm(R) Y_J(r) dot Y_J(R) ||n1 l1 n2 l2 lambda>
  // in terms of Moshinsky brackets and Talmi integrals
  
  // Initialize matrix element storage
  double mat_RE = 0.0;
  double mat_IM = 0.0;
 

  // Compute max oscillator quanta
  int max = 2*n1 + 2*n2 + l1 + l2;
  int maxp = 2*n1p + 2*n2p + l1p + l2p;

  // Operator is angular momentum scalar
  if (lambda != lambdap) {return;}

  // Transform to relative and CoM coordinates using Brody-Moshinsky transformation
  // Our conventions follow Moshinsky Nucl. Phys. 13 (1959) 104
  for (int l_cm = 0; l_cm <= max; l_cm++) {
    for (int l_rel = 0; l_rel <= max - l_cm; l_rel++) {
      // Parity is conserved
      if (pow(-1.0, l_rel + l_cm) != pow(-1.0, l1 + l2)) {continue;}
      // Symmetry factor arises from anti-symmetrized matrix element (exchange term)
      double sym =  1.0 + pow(-1.0, l_rel + s + 1 + t);
      if (sym == 0.0) {continue;}

      for (int n_cm = 0; n_cm <= (max - l_cm - l_rel)/2; n_cm++) {
        int n_rel = (max - l_rel - l_cm)/2 - n_cm;
	if (n_rel < 0) {continue;}
	for (int l_cmp = 0; l_cmp <= maxp; l_cmp++) {
          for (int l_relp = 0; l_relp <= maxp - l_cmp; l_relp ++) {
            if (pow(-1.0, l_relp + l_cmp) != pow(-1.0, l1p + l2p)) {continue;}
	    for (int n_cmp = 0; n_cmp <= (maxp - l_cmp - l_relp)/2; n_cmp++) {
            int n_relp = (maxp - l_relp - l_cmp)/2 - n_cmp;
            if (n_relp < 0) {continue;}
 
 	    double rm = brody_mosh(n_rel, l_rel, n_cm, l_cm, lambda, n1, l1, n2, l2);
            rm *= brody_mosh(n_relp, l_relp, n_cmp, l_cmp, lambdap, n1p, l1p, n2p, l2p);
            rm *= 1.0/(4.0*M_PI)*pow(-1.0, l_rel + l_relp + lambda)*(2.0*J + 1.0)*sqrt((2.0*l_rel + 1.0)*(2.0*l_relp + 1.0)*(2.0*l_cm + 1.0)*(2.0*l_cmp + 1.0))*three_j(l_relp, J, l_rel, 0.0, 0.0, 0.0)*three_j(l_cmp, J, l_cm, 0.0, 0.0, 0.0)*six_j(lambda, l_cmp, l_relp, J, l_rel, l_cm);
            if (rm == 0.0) {continue;}
            rm *= sym;
            double mrel_RE = compute_rel_potential_spline(n_relp, l_relp, n_rel, l_rel, f_spline_RE, acc);
            double mrel_IM = compute_rel_potential_spline(n_relp, l_relp, n_rel, l_rel, f_spline_IM, acc);

	  //  printf("np: %d lp: %d n: %d l: %d mat: %g\n", n_relp, l_relp, n_rel, l_rel, mrel_RE);
	    rm *= compute_rel_potential(n_cmp, l_cmp, n_cm, l_cm, J, q, 2);
	 //   printf("Np: %d LP: %d N: %d L: %d RCM: %g\n", n_cmp, l_cmp, n_cm, l_cm, compute_rel_potential(n_cmp, l_cmp, n_cm, l_cm, J, q, 2));
            mat_RE += rm*mrel_RE;
	    mat_IM += rm*mrel_IM;
	    }
	  }
        }
      }
    }
  }
  //printf("%d %d %d %d %d %d %d %d %d %d %d %d %g\n", n1p, l1p, n2p, l2p, lambda, n1, l1, n2, l2, lambda, s, t, mat_RE);
  *rm_RE = mat_RE;
  *rm_IM = mat_IM;

  return;
}
