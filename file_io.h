#ifndef FILE_IO_H
#define FILE_IO_H
#include "multipole.h"

typedef struct slater_det
{
  int* label;
} slater_det;

typedef struct wfnData
{
  int n_proton, n_neutron;
  slater_det** basis;
  long long int n_states;
  int n_shells, n_orbits, n_data;
  int n_eig;
  float b_mag;
  unsigned int n_sds_proton, n_sds_neutron;
  double *bc;
  int *n_shell, *l_shell, *j_shell, *jz_shell, *tz_shell;
  int *n_orb, *l_orb;
  float *j_orb;
  float jz;
  float *e_nuc, *j_nuc, *t_nuc;
} wfnData;

wfnData* read_binary_wfn_data(char *wfn_file, char* basis_file, int i_eig);

#endif
