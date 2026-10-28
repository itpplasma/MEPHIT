#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mephit_fortfem.h"

enum { divisions = 40 };

static void free_arrays(
  double *node_R,
  double *node_Z,
  int *tri_node,
  int *edge_node,
  double complex *current,
  double complex *magnetic,
  double complex *potential_R,
  double complex *potential_Z)
{
  free(node_R);
  free(node_Z);
  free(tri_node);
  free(edge_node);
  free(current);
  free(magnetic);
  free(potential_R);
  free(potential_Z);
}

int main(void)
{
  const int points_per_side = divisions + 1;
  const int npoint = points_per_side * points_per_side;
  const int ntri = 2 * divisions * divisions;
  const int nedge =
    2 * divisions * points_per_side + divisions * divisions;
  const double complex scale = 1.0 + 2.0 * I;
  double complex *current = calloc((size_t) nedge, sizeof(*current));
  double complex *magnetic = calloc((size_t) nedge, sizeof(*magnetic));
  double complex *potential_R =
    calloc((size_t) npoint, sizeof(*potential_R));
  double complex *potential_Z =
    calloc((size_t) npoint, sizeof(*potential_Z));
  double *node_R = malloc((size_t) npoint * sizeof(*node_R));
  double *node_Z = malloc((size_t) npoint * sizeof(*node_Z));
  int *tri_node = malloc(3U * (size_t) ntri * sizeof(*tri_node));
  int *edge_node = malloc(2U * (size_t) nedge * sizeof(*edge_node));
  double magnitude;
  double norm;
  double phase_error;
  int edge;
  int i;
  int j;
  int status;
  int triangle;

  if (node_R == NULL || node_Z == NULL || tri_node == NULL ||
      edge_node == NULL || current == NULL || magnetic == NULL ||
      potential_R == NULL || potential_Z == NULL) {
    fprintf(stderr, "production-mesh allocation failed\n");
    free_arrays(
      node_R, node_Z, tri_node, edge_node, current, magnetic,
      potential_R, potential_Z);
    return 1;
  }

  for (j = 0; j < points_per_side; ++j) {
    for (i = 0; i < points_per_side; ++i) {
      const int point = j * points_per_side + i;

      node_R[point] = 2.0 + (double) i / (double) divisions;
      node_Z[point] = (double) j / (double) divisions;
    }
  }

  triangle = 0;
  for (j = 0; j < divisions; ++j) {
    for (i = 0; i < divisions; ++i) {
      const int lower_left = j * points_per_side + i + 1;
      const int lower_right = lower_left + 1;
      const int upper_left = lower_left + points_per_side;
      const int upper_right = upper_left + 1;

      tri_node[3 * triangle] = lower_left;
      tri_node[3 * triangle + 1] = lower_right;
      tri_node[3 * triangle + 2] = upper_right;
      ++triangle;
      tri_node[3 * triangle] = lower_left;
      tri_node[3 * triangle + 1] = upper_right;
      tri_node[3 * triangle + 2] = upper_left;
      ++triangle;
    }
  }

  edge = 0;
  for (j = 0; j < points_per_side; ++j) {
    for (i = 0; i < divisions; ++i) {
      edge_node[2 * edge] = j * points_per_side + i + 1;
      edge_node[2 * edge + 1] = edge_node[2 * edge] + 1;
      ++edge;
    }
  }
  for (j = 0; j < divisions; ++j) {
    for (i = 0; i < points_per_side; ++i) {
      edge_node[2 * edge] = j * points_per_side + i + 1;
      edge_node[2 * edge + 1] =
        edge_node[2 * edge] + points_per_side;
      ++edge;
    }
  }
  for (j = 0; j < divisions; ++j) {
    for (i = 0; i < divisions; ++i) {
      edge_node[2 * edge] = j * points_per_side + i + 1;
      edge_node[2 * edge + 1] =
        edge_node[2 * edge] + points_per_side + 1;
      ++edge;
    }
  }
  if (edge != nedge || triangle != ntri) {
    fprintf(stderr, "production-mesh topology count mismatch\n");
    free_arrays(
      node_R, node_Z, tri_node, edge_node, current, magnetic,
      potential_R, potential_Z);
    return 2;
  }

  for (edge = 0; edge < nedge; ++edge) {
    const int a = edge_node[2 * edge] - 1;
    const int b = edge_node[2 * edge + 1] - 1;
    const double dR = node_R[b] - node_R[a];
    const double dZ = node_Z[b] - node_Z[a];
    const double qR_average = node_R[a] + node_R[b] - 3.0;
    const double qZ_average = 3.0 + node_Z[a] + node_Z[b];

    current[edge] =
      scale * (dZ * qR_average - dR * qZ_average);
  }

  status = mephit_fortfem_init(
    npoint, node_R, node_Z, ntri, tri_node, nedge, edge_node);
  if (status != 0) {
    fprintf(stderr, "production-mesh initialization failed: %d\n", status);
    free_arrays(
      node_R, node_Z, tri_node, edge_node, current, magnetic,
      potential_R, potential_Z);
    return 3;
  }

  status = mephit_fortfem_l2(nedge, current, &norm);
  if (status != 0 || fabs(norm - sqrt(5.0 * 62.0 / 3.0)) > 2.0e-12) {
    fprintf(stderr, "production-mesh RT norm mismatch: %.17g\n", norm);
    mephit_fortfem_deinit();
    free_arrays(
      node_R, node_Z, tri_node, edge_node, current, magnetic,
      potential_R, potential_Z);
    return 4;
  }

  status = mephit_fortfem_prepare(2);
  if (status != 0) {
    fprintf(stderr, "production-mesh preparation failed: %d\n", status);
    mephit_fortfem_deinit();
    free_arrays(
      node_R, node_Z, tri_node, edge_node, current, magnetic,
      potential_R, potential_Z);
    return 5;
  }

  magnitude = 0.0;
  status = mephit_fortfem_solve(
    nedge, npoint, magnetic, current, potential_R, potential_Z);
  for (edge = 0; edge < nedge; ++edge) {
    magnitude += cabs(current[edge]);
  }
  for (i = 0; i < npoint; ++i) {
    magnitude += cabs(potential_R[i]) + cabs(potential_Z[i]);
  }
  if (status != 0 || magnitude > 1.0e-14) {
    fprintf(stderr,
            "production-mesh zero solve failed: %d %.17g\n",
            status, magnitude);
    mephit_fortfem_deinit();
    free_arrays(
      node_R, node_Z, tri_node, edge_node, current, magnetic,
      potential_R, potential_Z);
    return 6;
  }

  for (edge = 0; edge < nedge; ++edge) {
    const int a = edge_node[2 * edge] - 1;
    const int b = edge_node[2 * edge + 1] - 1;
    const double dR = node_R[b] - node_R[a];
    const double dZ = node_Z[b] - node_Z[a];
    const double qR_average = node_R[a] + node_R[b] - 3.0;
    const double qZ_average = 3.0 + node_Z[a] + node_Z[b];

    current[edge] =
      scale * (dZ * qR_average - dR * qZ_average);
  }
  status = mephit_fortfem_solve(
    nedge, npoint, current, magnetic, potential_R, potential_Z);
  magnitude = 0.0;
  phase_error = 0.0;
  for (edge = 0; edge < nedge; ++edge) {
    magnitude += cabs(magnetic[edge]);
    phase_error +=
      fabs(cimag(magnetic[edge]) + 0.5 * creal(magnetic[edge]));
  }
  for (i = 0; i < npoint; ++i) {
    magnitude += cabs(potential_R[i]) + cabs(potential_Z[i]);
    phase_error +=
      fabs(cimag(potential_R[i]) - 2.0 * creal(potential_R[i]));
    phase_error +=
      fabs(cimag(potential_Z[i]) - 2.0 * creal(potential_Z[i]));
  }
  mephit_fortfem_deinit();
  free_arrays(
    node_R, node_Z, tri_node, edge_node, current, magnetic,
    potential_R, potential_Z);
  if (status != 0 || !(magnitude > 1.0e-20) ||
      phase_error > 2.0e-12 * magnitude) {
    fprintf(stderr,
            "production-mesh complex phase failed: %d %.17g %.17g\n",
            status, magnitude, phase_error);
    return 7;
  }
  return 0;
}
