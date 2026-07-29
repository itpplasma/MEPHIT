#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "mephit_fortfem.h"

int main(void)
{
  const double node_R[4] = {0.0, 1.0, 1.0, 0.0};
  const double node_Z[4] = {0.0, 0.0, 1.0, 1.0};
  const int tri_node[6] = {1, 2, 3, 1, 4, 3};
  const int edge_node[10] = {
    2, 1, 2, 3, 4, 3, 1, 4, 3, 1
  };
  const double complex scale = 1.0 + 2.0 * I;
  double complex dofs[5];
  double norm;
  int edge;
  int status;

  for (edge = 0; edge < 5; ++edge) {
    const int a = edge_node[2 * edge] - 1;
    const int b = edge_node[2 * edge + 1] - 1;
    const double dR = node_R[b] - node_R[a];
    const double dZ = node_Z[b] - node_Z[a];
    const double qR_average = 1.0 + node_R[a] + node_R[b];
    const double qZ_average = 3.0 + node_Z[a] + node_Z[b];

    dofs[edge] = scale * (dZ * qR_average - dR * qZ_average);
  }

  status = mephit_fortfem_init(
    4, node_R, node_Z, 2, tri_node, 5, edge_node);
  if (status != 0) {
    fprintf(stderr, "adapter initialization failed: %d\n", status);
    return 1;
  }
  status = mephit_fortfem_l2(5, dofs, &norm);
  mephit_fortfem_deinit();
  if (status != 0) {
    fprintf(stderr, "adapter L2 evaluation failed: %d\n", status);
    return 2;
  }
  if (fabs(norm - sqrt(5.0 * 62.0 / 3.0)) > 1.0e-13) {
    fprintf(stderr, "adapter L2 mismatch: %.17g\n", norm);
    return 3;
  }
  return 0;
}
