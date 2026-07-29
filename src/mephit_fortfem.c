#include <complex.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <fortfem.h>

#include "mephit_fortfem.h"

static int mesh_handle;
static int edge_count;
static int *fort_dof;
static int *edge_sign;
static fortfem_complex *fort_dofs;

void mephit_fortfem_deinit(void)
{
  int status;

  if (mesh_handle != 0) {
    fortfem_triangle_mesh_free(mesh_handle, &status);
  }
  free(fort_dof);
  free(edge_sign);
  free(fort_dofs);
  mesh_handle = 0;
  edge_count = 0;
  fort_dof = NULL;
  edge_sign = NULL;
  fort_dofs = NULL;
}

int mephit_fortfem_init(
  const int npoint,
  const double *node_R,
  const double *node_Z,
  const int ntri,
  const int *tri_node,
  const int nedge,
  const int *edge_node)
{
  double *vertices;
  int *triangles;
  int *edges;
  int *triangle_edge_dofs;
  int *triangle_edge_signs;
  int edge;
  int fort_edge;
  int found;
  int nedge_returned;
  int status;
  int triangle;

  mephit_fortfem_deinit();
  if (npoint < 3 || ntri < 1 || nedge < 1 || node_R == NULL ||
      node_Z == NULL || tri_node == NULL || edge_node == NULL) {
    return -1;
  }

  vertices = malloc(2U * (size_t) npoint * sizeof(*vertices));
  triangles = malloc(3U * (size_t) ntri * sizeof(*triangles));
  edges = malloc(2U * (size_t) nedge * sizeof(*edges));
  triangle_edge_dofs =
    malloc(3U * (size_t) ntri * sizeof(*triangle_edge_dofs));
  triangle_edge_signs =
    malloc(3U * (size_t) ntri * sizeof(*triangle_edge_signs));
  fort_dof = malloc((size_t) nedge * sizeof(*fort_dof));
  edge_sign = malloc((size_t) nedge * sizeof(*edge_sign));
  fort_dofs = malloc((size_t) nedge * sizeof(*fort_dofs));
  if (vertices == NULL || triangles == NULL || edges == NULL ||
      triangle_edge_dofs == NULL || triangle_edge_signs == NULL ||
      fort_dof == NULL || edge_sign == NULL || fort_dofs == NULL) {
    free(vertices);
    free(triangles);
    free(edges);
    free(triangle_edge_dofs);
    free(triangle_edge_signs);
    mephit_fortfem_deinit();
    return -2;
  }

  for (edge = 0; edge < npoint; ++edge) {
    vertices[2 * edge] = node_R[edge];
    vertices[2 * edge + 1] = node_Z[edge];
  }
  for (triangle = 0; triangle < ntri; ++triangle) {
    int a = tri_node[3 * triangle] - 1;
    int b = tri_node[3 * triangle + 1] - 1;
    int c = tri_node[3 * triangle + 2] - 1;
    double coordinate_scale;
    double signed_area_twice;

    if (a < 0 || a >= npoint || b < 0 || b >= npoint ||
        c < 0 || c >= npoint) {
      status = -1;
      goto cleanup;
    }
    signed_area_twice =
      (node_R[b] - node_R[a]) * (node_Z[c] - node_Z[a]) -
      (node_Z[b] - node_Z[a]) * (node_R[c] - node_R[a]);
    coordinate_scale = fmax(
      fmax(fabs(node_R[b] - node_R[a]), fabs(node_Z[b] - node_Z[a])),
      fmax(fabs(node_R[c] - node_R[a]), fabs(node_Z[c] - node_Z[a])));
    if (fabs(signed_area_twice) <=
        DBL_EPSILON * coordinate_scale * coordinate_scale) {
      status = -1;
      goto cleanup;
    }
    if (signed_area_twice < 0.0) {
      const int temporary = b;
      b = c;
      c = temporary;
    }
    triangles[3 * triangle] = a;
    triangles[3 * triangle + 1] = b;
    triangles[3 * triangle + 2] = c;
  }

  fortfem_triangle_mesh_create(
    npoint, vertices, ntri, triangles, &mesh_handle, &nedge_returned, &status);
  if (status != 0 || nedge_returned != nedge) {
    if (status == 0) {
      status = -3;
    }
    goto cleanup;
  }
  fortfem_triangle_mesh_edges(
    mesh_handle, nedge, &nedge_returned, edges,
    triangle_edge_dofs, triangle_edge_signs, &status);
  if (status != 0) {
    goto cleanup;
  }

  for (edge = 0; edge < nedge; ++edge) {
    const int start = edge_node[2 * edge] - 1;
    const int end = edge_node[2 * edge + 1] - 1;

    found = 0;
    for (fort_edge = 0; fort_edge < nedge; ++fort_edge) {
      if (start == edges[2 * fort_edge] &&
          end == edges[2 * fort_edge + 1]) {
        fort_dof[edge] = fort_edge;
        edge_sign[edge] = 1;
        found = 1;
        break;
      }
      if (start == edges[2 * fort_edge + 1] &&
          end == edges[2 * fort_edge]) {
        fort_dof[edge] = fort_edge;
        edge_sign[edge] = -1;
        found = 1;
        break;
      }
    }
    if (!found) {
      status = -4;
      goto cleanup;
    }
  }

  edge_count = nedge;
  status = 0;

cleanup:
  free(vertices);
  free(triangles);
  free(edges);
  free(triangle_edge_dofs);
  free(triangle_edge_signs);
  if (status != 0) {
    mephit_fortfem_deinit();
  }
  return status;
}

int mephit_fortfem_l2(
  const int nedge,
  const double complex *mephit_dofs,
  double *norm)
{
  int edge;
  int status;

  if (mesh_handle == 0 || nedge != edge_count ||
      mephit_dofs == NULL || norm == NULL) {
    return -1;
  }
  for (edge = 0; edge < nedge; ++edge) {
    fort_dofs[fort_dof[edge]] =
      (double) edge_sign[edge] * mephit_dofs[edge];
  }
  fortfem_rt0_l2_norm(mesh_handle, nedge, fort_dofs, norm, &status);
  return status;
}
