#ifndef MEPHIT_FORTFEM_H
#define MEPHIT_FORTFEM_H

#include <complex.h>

int mephit_fortfem_init(
  int npoint,
  const double *node_R,
  const double *node_Z,
  int ntri,
  const int *tri_node,
  int nedge,
  const int *edge_node);

int mephit_fortfem_l2(
  int nedge,
  const double complex *mephit_dofs,
  double *norm);

int mephit_fortfem_prepare(int toroidal_mode);

int mephit_fortfem_solve(
  int nedge,
  int npoint,
  const double complex *current_dofs,
  double complex *magnetic_dofs,
  double complex *potential_R,
  double complex *potential_Z);

void mephit_fortfem_deinit(void);

#endif
