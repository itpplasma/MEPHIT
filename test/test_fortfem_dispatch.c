#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "mephit_fem.h"

int main(void)
{
  const double node_R[4] = {2.0, 3.0, 3.0, 2.0};
  const double node_Z[4] = {0.0, 0.0, 1.0, 1.0};
  const int tri_node[6] = {1, 2, 3, 1, 4, 3};
  const int edge_node[10] = {
    2, 1, 2, 3, 4, 3, 1, 4, 3, 1
  };
  const complex_double scale = 1.0 + 2.0 * I;
  complex_double AnR[4];
  complex_double AnZ[4];
  complex_double Bn[5];
  complex_double current[5];
  double magnitude = 0.0;
  double phase_error = 0.0;
  int edge;

  for (edge = 0; edge < 5; ++edge) {
    const int a = edge_node[2 * edge] - 1;
    const int b = edge_node[2 * edge + 1] - 1;
    const double dR = node_R[b] - node_R[a];
    const double dZ = node_Z[b] - node_Z[a];
    const double qR_average = node_R[a] + node_R[b] - 3.0;
    const double qZ_average = 3.0 + node_Z[a] + node_Z[b];

    current[edge] =
      scale * (dZ * qR_average - dR * qZ_average);
  }

  FEM_init_fortfem_mesh(
    4, node_R, node_Z, 2, tri_node, 5, edge_node);
  FEM_init(2, 5, 4, 0);
  FEM_compute_magfn(5, 4, current, Bn, AnR, AnZ);
  for (edge = 0; edge < 5; ++edge) {
    magnitude += cabs(Bn[edge]);
    phase_error += fabs(cimag(Bn[edge]) + 0.5 * creal(Bn[edge]));
  }
  for (edge = 0; edge < 4; ++edge) {
    magnitude += cabs(AnR[edge]) + cabs(AnZ[edge]);
    phase_error += fabs(cimag(AnR[edge]) - 2.0 * creal(AnR[edge]));
    phase_error += fabs(cimag(AnZ[edge]) - 2.0 * creal(AnZ[edge]));
  }
  FEM_deinit();

  if (!(magnitude > 1.0e-20) ||
      phase_error > 1.0e-12 * magnitude) {
    fprintf(stderr,
            "FortFEM dispatch violated complex phase: %.17g %.17g\n",
            magnitude, phase_error);
    return 1;
  }
  return 0;
}
