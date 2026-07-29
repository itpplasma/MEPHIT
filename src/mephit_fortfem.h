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

void mephit_fortfem_deinit(void);

#endif
