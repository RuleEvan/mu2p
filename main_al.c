#include "matrix_element.h"

int main(int argc, char *argv[]) {

  static int nt = 23;
  char density_file[23][50] = {"al27_na27_bw84_J0_T2_0_0.dens", "al27_na27_bw84_J0_T2_0_5.dens", "al27_na27_bw84_J0_T2_0_7.dens", "al27_na27_bw84_J0_T2_0_12.dens", "al27_na27_bw84_J0_T2_0_14.dens", "al27_na27_bw84_J0_T2_0_21.dens", "al27_na27_bw84_J0_T2_0_22.dens", "al27_na27_bw84_J0_T2_0_30.dens", "al27_na27_bw84_J0_T2_0_33.dens", "al27_na27_bw84_J0_T2_0_40.dens", "al27_na27_bw84_J0_T2_0_41.dens", "al27_na27_bw84_J0_T2_0_48.dens", "al27_na27_bw84_J0_T2_0_54.dens", "al27_na27_bw84_J0_T2_0_57.dens", "al27_na27_bw84_J0_T2_0_61.dens", "al27_na27_bw84_J0_T2_0_66.dens", "al27_na27_bw84_J0_T2_0_72.dens", "al27_na27_bw84_J0_T2_0_77.dens", "al27_na27_bw84_J0_T2_0_80.dens", "al27_na27_bw84_J0_T2_0_83.dens", "al27_na27_bw84_J0_T2_0_166.dens", "al27_na27_bw84_J0_T2_0_343.dens", "al27_na27_bw84_J0_T2_0_494.dens"};

  double Earr[23] = {0.0, 2.68685, 3.60913, 3.91678, 4.35136, 4.90839, 5.12316, 5.66790, 5.75033, 6.20264, 6.31398, 6.69324, 6.95446, 7.00781, 7.17524, 7.30895, 7.54311, 7.82888, 7.86513, 7.96331, 10.06092, 12.53328, 13.98288};
  double ji = 2.5;
  double jf = 2.5;
  double ti = 0.5;
  double tf = 2.5;
  double mti = -0.5;
  double mtf = -2.5;

  double E_avg = 10.0;
  double pe = 94.5;
  double Mi = 25126.5;
  double Mn = 25129.6;
  double alpha = M_MUON - E_avg + Mi - Mn;

  double Stot = 0.0;

  for (int i = nt - 1; i >= 0; i--) {

    double delta_e = pe - Earr[i] + E_avg - Mi + Mn;


//    printf("de: %g, alpha: %g\n", delta_e, alpha);

    double m_RE_GT = 0.0;
    double m_IM_GT = 0.0;

    compute_total_matrix_element_GTJ(density_file[i], pe, 0, alpha, delta_e, &m_RE_GT, &m_IM_GT);
    compute_total_matrix_element_GTJ(density_file[i], pe, 2, alpha, delta_e, &m_RE_GT, &m_IM_GT);
    compute_total_matrix_element_GTJ(density_file[i], pe, 4, alpha, delta_e, &m_RE_GT, &m_IM_GT);
//    printf("GT: %g, %g\n", m_RE_GT, m_IM_GT);

    double m_RE_F = 0.0;
    double m_IM_F = 0.0;

    compute_total_matrix_element_FJ(density_file[i], pe, 0, alpha, delta_e, &m_RE_F, &m_IM_F);
    compute_total_matrix_element_FJ(density_file[i], pe, 2, alpha, delta_e, &m_RE_F, &m_IM_F);
    compute_total_matrix_element_FJ(density_file[i], pe, 4, alpha, delta_e, &m_RE_F, &m_IM_F);
//    printf("MF: %g, %g\n", m_RE_F, m_IM_F);

    double IFact = pow(-1, ti - tf)*clebsch_gordan(2, ti, tf, -2, mti, mtf)/sqrt(2.0*tf + 1.0);


    double MTOT = (pow(pow(1.0/1.25, 2)*m_RE_F - m_RE_GT, 2) + pow(pow(1.0/1.25, 2)*m_IM_F - m_IM_GT, 2))*pow(IFact, 2)*1.0/(2.0*ji + 1.0)*pow(4.0*M_PI, 2);

    printf("E: %g |M|^2: %g\n", Earr[i], MTOT); 

    Stot += pow(1.0 - Earr[i]/pe, 2)*MTOT;
  }
  printf("\n");
  printf("Stot: %g\n", Stot);

  return 0;
}
