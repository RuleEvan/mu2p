#include "matrix_element.h"

int main(int argc, char *argv[]) {

//  filter_spectrum("na27_usdb_all.res", 2.5, 2.5);
//exit(0);
      	compute_mu2p_transitions();

//  printf("%g\n", compute_total_matrix_element_GT0("mg26_mg26_usdb_2bdy_J0_T0_0_0.dens"));
exit(0);

	int n = 9;
	int np = 10;
	int l = 0;
	int lp = 0;
	int p = 0;

printf("B: %g\n", b_coeff(n, l, np, lp, p));	
//  generate_bigstick_int_file(4, 1);
  exit(0);


//  printf("%g\n", compute_total_matrix_element_F0("../SPEED-DMG/DMG/examples/output/na19_f19_usdb_J0_T2_0_16.dens"));
//  exit(0); 
//  filter_spectrum("ne19_usdb_clean.res", 0.5, 1.5);
//  exit(0);
//  test_suite();
  compute_mu2p_transitions();
//  perform_closure_analysis();
  exit(0);
  FILE *list_file;
  list_file = fopen("lists/al27_mg27_na27.list", "r");
//  list_file = fopen("f18.list", "r");

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
  double ji = 0.5;
  double mti = -0.5;
  double mtf = -2.5;
 
  double pe_max = 94.5;
  double tn = 1.5;
  double mtn = -1.5;
/*
  while(fscanf(list_file, "%d, %f\n", &i_state, &excite) == 2) {
//   if (i_state > 5000) {continue;}
    sprintf(density_file_i, "../SPEED-DMG/DMG/examples/output/al27_mg27_usdb_0_%d.dens", i_state);
    sprintf(density_file_f, "../SPEED-DMG/DMG/examples/output/mg27_na27_usdb_%d_0.dens", i_state);
  
    double kt = 0.001;
    excite = E_avg; 
    double qt = M_MUON - excite + Mi - Mn;
    
    double IFact_i = pow(-1, 1 - ti + tn)*clebsch_gordan(1, ti, tn, -1, mti, mtn)/sqrt(2.0*tn + 1);
    double IFact_f = pow(-1, 1 - tn + tf)*clebsch_gordan(1, tn, tf, -1, mtn, mtf)/sqrt(2.0*tf + 1);
    double mat = compute_total_matrix_element_F1_double(density_file_i, density_file_f, qt, kt)*IFact_i*IFact_f*qt*RNUC/HBARC/(4*M_PI*M_PI)*(-M_PI);

    total_mat += mat;

    printf("%d %g %g %g\n", i_state, excite, mat, total_mat);
  }
  printf("Total: %g\n", total_mat);

*/
//exit(0);

//  exit(0);
//  generate_bigstick_int_file(-2, 1.0);
//  exit(0);
//  printf("MGT: %g\n", compute_total_matrix_element_GT0("../SPEED-DMG/DMG/examples/output/al27_na27_usdb_J0_T2_0_0.dens")*clebsch_gordan(2, 0.5, 2.5, -2, -0.5, -2.5)/sqrt(2.0*2.5 + 1.0));
//  exit(0);

 
  char density_file[250];
  float m_total = 0.0;

  double IFact;

//  double pe_max = 94.5;
//  double Mi = 26981.5;
//  double Mn = 26984.3;
  double alpha = M_MUON - E_avg + Mi - Mn;

  double Stot = 0.0;

  list_file = fopen("lists/al27_na27_usdb_J25.list", "r");


  while(fscanf(list_file, "%d, %f, %f\n", &i_state, &tf, &excite) == 3) {
    sprintf(density_file, "../SPEED-DMG/DMG/examples/output/al27_na27_bw84_J0_T2_0_0.dens", i_state);
    tf = 2.5;
    if (excite > 25.0) {continue;}
    IFact = pow(-1, ti - tf)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1.0);
    double m_gt = pow(compute_total_matrix_element_GT0(density_file)*IFact, 2);
    m_total += m_gt;
//    printf("%g, %g\n", excite, m_gt);

//    printf("E: %g Bare GT: %g Bare F: %g\n", excite, compute_total_matrix_element_GT0(density_file), compute_total_matrix_element_F0(density_file));
    double pe = pe_max - excite;
    double delta_e = pe + E_avg - Mi + Mn;

    double m_RE_GT = 0.0;
    double m_IM_GT = 0.0;

    compute_total_matrix_element_GTJ(density_file, pe, 0, alpha, delta_e, &m_RE_GT, &m_IM_GT);
    compute_total_matrix_element_GTJ(density_file, pe, 2, alpha, delta_e, &m_RE_GT, &m_IM_GT);
    compute_total_matrix_element_GTJ(density_file, pe, 4, alpha, delta_e, &m_RE_GT, &m_IM_GT);


    printf("MGT: %g, %g\n", m_RE_GT, m_IM_GT);


    double m_RE_F = 0.0;
    double m_IM_F = 0.0;

    compute_total_matrix_element_FJ(density_file, pe, 0, alpha, delta_e, &m_RE_F, &m_IM_F);
    compute_total_matrix_element_FJ(density_file, pe, 2, alpha, delta_e, &m_RE_F, &m_IM_F);
    compute_total_matrix_element_FJ(density_file, pe, 4, alpha, delta_e, &m_RE_F, &m_IM_F);

    printf("MF: %g, %g\n", m_RE_F*IFact, m_IM_F*IFact);
//    printf("Re: %g Im: %g\n", (pow(1.0/1.25, 2)*m_RE_F - m_RE_GT)*IFact,  (pow(1.0/1.25, 2)*m_IM_F - m_IM_GT)*IFact);


    double MTOT = (pow(pow(1.0/1.25, 2)*m_RE_F - m_RE_GT, 2) + pow(pow(1.0/1.25, 2)*m_IM_F - m_IM_GT, 2))*pow(IFact, 2)*1.0/(2.0*ji + 1.0);

    printf("Ex: %g S: %g\n", excite, MTOT); 
    exit(0);
//    printf("%g, %g\n", excite, 4.0*1.0/(2*ji + 1.0)*m_gt); 
    printf("\n");
    Stot += pow(1.0 - excite/pe_max, 2)*MTOT;


  }

//  printf("Total: %g\n", m_total);
//  printf("Stot: %g\n", Stot);
//  printf("\n");

  exit(0);

  return 0;
}
