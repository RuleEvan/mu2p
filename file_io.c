#include "file_io.h"

wfnData* read_binary_wfn_data(char *wfn_file, char *basis_file, int i_eig) {
  wfnData *wd = malloc(sizeof(*wd));
  FILE *in_file;
  // Read in initial wavefunction data
  printf("Opening file\n");
  in_file = fopen(wfn_file, "rb");
  int junk;
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  int vec_offset;
  fread(&vec_offset, sizeof(int), 1, in_file);
  for (int i = 0; i < 5; i++) {
    fread(&junk, sizeof(int), 1, in_file);
  }
  fread(&wd->n_proton, sizeof(int), 1, in_file);
  fread(&wd->n_neutron, sizeof(int), 1, in_file);
  wd->n_data = wd->n_proton + wd->n_neutron;
  printf("Initial state contains %d protons and %d neutrons\n", wd->n_proton, wd->n_neutron);

  fread(&junk, sizeof(int), 1, in_file);
//  fread(&junk, sizeof(int), 1, in_file);
  int n_proton_orbs, n_neutron_orbs;
  fread(&n_proton_orbs, sizeof(int), 1, in_file);
  fread(&n_neutron_orbs, sizeof(int), 1, in_file);
  wd->n_orbits = n_proton_orbs;
  printf("Model space contains %d orbitals\n", n_proton_orbs);
  wd->n_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->l_orb = (int*) malloc(sizeof(int)*wd->n_orbits);
  wd->j_orb = (float*) malloc(sizeof(float)*wd->n_orbits);
  for (int i = 0; i < n_proton_orbs; i++) {
    int j_orb, pi_orb, w_orb;
    fread(&wd->n_orb[i], sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    wd->j_orb[i] = j_orb/2.0;
    fread(&wd->l_orb[i], sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&w_orb, sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_neutron_orbs; i++) {
    int n_orb, j_orb, l_orb, pi_orb, w_orb;
    fread(&n_orb, sizeof(int), 1, in_file);
    fread(&j_orb, sizeof(int), 1, in_file);
    fread(&l_orb, sizeof(int), 1, in_file);
    fread(&pi_orb, sizeof(int), 1, in_file);
    fread(&w_orb, sizeof(int), 1, in_file);
  }
  int n_sps_proton, n_sps_neutron;
  fread(&n_sps_proton, sizeof(int), 1, in_file);
  fread(&n_sps_neutron, sizeof(int), 1, in_file);
  wd->n_shells = n_sps_proton;
  printf("Number of single particle states: %d\n", wd->n_shells);
  wd->n_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->j_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->l_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  wd->jz_shell = (int*) malloc(sizeof(int)*wd->n_shells);
  for (int i = 0; i < n_sps_proton; i++) {
    int n_shell, l_shell, j_shell, jz_shell, w_shell, pi_shell, id_shell, gid_shell;
    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }
  for (int i = 0; i < n_sps_neutron; i++) {
    int n_shell, j_shell, jz_shell, l_shell, w_shell, pi_shell, id_shell, gid_shell;

    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
    fread(&pi_shell, sizeof(int), 1, in_file);
    fread(&id_shell, sizeof(int), 1, in_file);
    fread(&gid_shell, sizeof(int), 1, in_file);
  }


  fread(&wd->jz, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&wd->n_states, sizeof(long long int), 1, in_file);
  printf("The model space has %d shells for a total basis size of %d\n", wd->n_shells, wd->n_states);
  fread(&wd->n_eig, sizeof(int), 1, in_file);
//  printf("Initial state file contains %d eigenstates\n", wd->n_eig);
  
  wd->e_nuc =  (float*) malloc(sizeof(float)*wd->n_eig);
  wd->j_nuc =  (float*) malloc(sizeof(float)*wd->n_eig);
  wd->t_nuc =  (float*) malloc(sizeof(float)*wd->n_eig);
  wd->bc = malloc(sizeof(double)*wd->n_states*wd->n_eig);

  printf("Reading in initial state wavefunction coefficients\n");
  for (int i = 0; i < wd->n_eig; i++) {
    fread(&junk, sizeof(int), 1, in_file);
    int vec_index;
    float j_nuc, t_nuc;
    fread(&vec_index, sizeof(int), 1, in_file);
    fread(&wd->e_nuc[i], sizeof(float), 1, in_file);
    fread(&j_nuc, sizeof(float), 1, in_file);
    fread(&t_nuc, sizeof(float), 1, in_file);
    t_nuc = -0.5 + 0.5*sqrt(1 + 4.0*t_nuc);
    int ij_nuc = round(2*j_nuc);
    int it_nuc = round(2*t_nuc);
    wd->j_nuc[i] = ij_nuc/2.0;
    wd->t_nuc[i] = it_nuc/2.0;
    if (i == i_eig) {
      for (int j = 0; j < wd->n_states; j++) {
        float bc_float;
        fread(&bc_float, sizeof(float), 1, in_file);
        wd->bc[j] = bc_float;
      }
    } else {
      fseek(in_file, sizeof(float)*wd->n_states, SEEK_CUR);
    }
  }
  fclose(in_file);
  printf("Done.\n");
  printf("Reading in basis\n");
  in_file = fopen(basis_file, "rb");
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&junk, sizeof(int), 1, in_file);
  fread(&vec_offset, sizeof(int), 1, in_file);
  fseek(in_file, vec_offset + 16, SEEK_SET);
  for (int i = 0; i < wd->n_shells; i++) {
    int w_shell;
    fread(&junk, sizeof(int), 1, in_file);
    fread(&wd->n_shell[i], sizeof(int), 1, in_file);
    fread(&wd->l_shell[i], sizeof(int), 1, in_file);
    fread(&wd->j_shell[i], sizeof(int), 1, in_file);
    fread(&wd->jz_shell[i], sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
 //   printf("%d, %d, %d, %d, %d\n", i+1, wd->n_shell[i], wd->l_shell[i], wd->j_shell[i], wd->jz_shell[i]);
  }
  for (int i = 0; i < wd->n_shells; i++) {
    int n_shell, j_shell, l_shell, jz_shell, w_shell;
    fread(&junk, sizeof(int), 1, in_file);
    fread(&n_shell, sizeof(int), 1, in_file);
    fread(&l_shell, sizeof(int), 1, in_file);
    fread(&j_shell, sizeof(int), 1, in_file);
    fread(&jz_shell, sizeof(int), 1, in_file);
    fread(&w_shell, sizeof(int), 1, in_file);
  }
  wd->basis = (slater_det**) calloc(wd->n_states, sizeof(slater_det*));
//  int* im_set = (int*) malloc(sizeof(int)*wd->n_data); 
//  int* compress = (int*) malloc(sizeof(int)*wd->n_data);
  for (int i = 0; i < wd->n_states; i++) {
    wd->basis[i] = malloc(sizeof(slater_det*));
//    wd->basis[i]->m_val = malloc(sizeof(float)*wd->n_data);
//    wd->basis[i]->n_val = malloc(sizeof(int)*wd->n_data);
//    wd->basis[i]->j_val = malloc(sizeof(float)*wd->n_data);
//    wd->basis[i]->l_val = malloc(sizeof(int)*wd->n_data);
    wd->basis[i]->label = malloc(sizeof(int)*wd->n_data);
    for (int j = 0; j < wd->n_data; j++) {
      int i_state;
      fread(&i_state, sizeof(int), 1, in_file);
    //  im_set[j] = wd->jz_shell[i_state - 1] + 2*wd->j_orb[0];
//      wd->basis[i]->m_val[j] = wd->jz_shell[i_state - 1]/2.0;
 //     wd->basis[i]->n_val[j] = wd->n_shell[i_state - 1];
//      wd->basis[i]->j_val[j] = wd->j_shell[i_state - 1];
 //     wd->basis[i]->l_val[j] = wd->l_shell[i_state - 1];
      wd->basis[i]->label[j] = i_state - 1;
    }
  /*  perm_compress(im_set, &compress, wd->n_data);
    int phase = perm_sign(compress, wd->n_data);
    order_perm(&im_set, wd->n_data);
    wd->bc[i] *= phase;
    for (int j = 0; j < wd->n_data; j++) {
      wd->basis[i]->m_val[j] = im_set[j]/2.0 - wd->j_orb[0];
    }*/
  }
  printf("Done.\n");
  fclose(in_file);
  
  return wd;
}
