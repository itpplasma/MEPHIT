#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "mephit_fortfem.h"

int main(void)
{
  const double node_R[4] = {2.0, 3.0, 3.0, 2.0};
  const double node_Z[4] = {0.0, 0.0, 1.0, 1.0};
  const int tri_node[6] = {1, 2, 3, 1, 4, 3};
  const int edge_node[10] = {
    2, 1, 2, 3, 4, 3, 1, 4, 3, 1
  };
  const double complex scale = 1.0 + 2.0 * I;
  double complex AnR[4];
  double complex AnZ[4];
  double complex Bn[5];
  double complex dofs[5];
  const double complex zero_dofs[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double magnitude;
  double norm;
  double phase_error;
  int edge;
  int status;

  for (edge = 0; edge < 5; ++edge) {
    const int a = edge_node[2 * edge] - 1;
    const int b = edge_node[2 * edge + 1] - 1;
    const double dR = node_R[b] - node_R[a];
    const double dZ = node_Z[b] - node_Z[a];
    const double qR_average = node_R[a] + node_R[b] - 3.0;
    const double qZ_average = 3.0 + node_Z[a] + node_Z[b];

    dofs[edge] = scale * (dZ * qR_average - dR * qZ_average);
  }

  status = mephit_fortfem_init(
    4, node_R, node_Z, 2, tri_node, 5, edge_node);
  if (status != 0) {
    fprintf(stderr, "adapter initialization failed: %d\n", status);
    return 1;
  }
  status = mephit_fortfem_prepare(2);
  if (status != 0) {
    fprintf(stderr, "adapter preparation failed: %d\n", status);
    return 2;
  }
  status = mephit_fortfem_solve(5, 4, zero_dofs, Bn, AnR, AnZ);
  if (status != 0) {
    fprintf(stderr, "adapter solve failed: %d\n", status);
    return 3;
  }
  if (cabs(Bn[0]) + cabs(Bn[1]) + cabs(Bn[2]) + cabs(Bn[3]) +
      cabs(Bn[4]) + cabs(AnR[0]) + cabs(AnR[1]) + cabs(AnR[2]) +
      cabs(AnR[3]) + cabs(AnZ[0]) + cabs(AnZ[1]) + cabs(AnZ[2]) +
      cabs(AnZ[3]) > 1.0e-14) {
    fprintf(stderr, "zero current produced a nonzero field\n");
    return 4;
  }
  status = mephit_fortfem_solve(5, 4, dofs, Bn, AnR, AnZ);
  if (status != 0) {
    fprintf(stderr, "nonzero adapter solve failed: %d\n", status);
    return 5;
  }
  magnitude = 0.0;
  phase_error = 0.0;
  for (edge = 0; edge < 5; ++edge) {
    magnitude += cabs(Bn[edge]);
    phase_error += fabs(cimag(Bn[edge]) + 0.5 * creal(Bn[edge]));
  }
  for (edge = 0; edge < 4; ++edge) {
    magnitude += cabs(AnR[edge]) + cabs(AnZ[edge]);
    phase_error += fabs(cimag(AnR[edge]) - 2.0 * creal(AnR[edge]));
    phase_error += fabs(cimag(AnZ[edge]) - 2.0 * creal(AnZ[edge]));
  }
  if (!(magnitude > 1.0e-20) ||
      phase_error > 1.0e-12 * magnitude) {
    fprintf(stderr,
            "real operator violated complex phase: %.17g %.17g\n",
            magnitude, phase_error);
    return 6;
  }
  status = mephit_fortfem_l2(5, dofs, &norm);
  mephit_fortfem_deinit();
  if (status != 0) {
    fprintf(stderr, "adapter L2 evaluation failed: %d\n", status);
    return 7;
  }
  if (fabs(norm - sqrt(5.0 * 62.0 / 3.0)) > 1.0e-13) {
    fprintf(stderr, "adapter L2 mismatch: %.17g\n", norm);
    return 8;
  }
  return 0;
}
