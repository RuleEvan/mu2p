#ifndef CHARGE_H
#define CHARGE_H
#include "file_io.h"
double bessel_potential(int n_steps, int z_atom, double *a, double r_cut, double a0, double r);

void generate_charge_density_file(wfnData* wd, double b_osc);

double generate_proton_density(wfnData *wd, double r, double theta, double b_osc);


double harmonic_osc_square(int n, int l, double j, double mj, int i_shell, double x, double theta, double b_osc, double **cache);

double anti_proton_wfn(wfnData *wd, int i_state, double x, double theta, double b_osc, double** cache);

double core_wfn(int n_max_quanta, double r, double b);

void generate_potential_files();

#endif

