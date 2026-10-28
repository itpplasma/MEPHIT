#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <fortfem.h>

#include "mephit_fortfem.h"

static int mesh_handle;
static int extended_mesh_handle;
static int factor_handle;
static int point_count;
static int triangle_count;
static int edge_count;
static int extended_edge_count;
static int interior_count;
static int toroidal_mode;
static double *core_vertices;
static int *core_triangles;
static int *core_edge_nodes;
static int *fort_dof;
static int *edge_sign;
static int *extended_dof;
static int *extended_sign;
static int *mixed_col_ptr;
static int *mixed_row_ind;
static double *mixed_values;
static fortfem_complex *fort_dofs;
static fortfem_complex *extended_rt_dofs;
static fortfem_complex *load;
static fortfem_complex *potential_dofs;
static fortfem_complex *interior_rhs;
static fortfem_complex *interior_solution;

static void release_prepared(void)
{
  int status;

  if (factor_handle != 0) {
    fortfem_factor_free(factor_handle, &status);
  }
  if (extended_mesh_handle != 0) {
    fortfem_triangle_mesh_free(extended_mesh_handle, &status);
  }
  free(extended_dof);
  free(extended_sign);
  free(mixed_col_ptr);
  free(mixed_row_ind);
  free(mixed_values);
  free(extended_rt_dofs);
  free(load);
  free(potential_dofs);
  free(interior_rhs);
  free(interior_solution);
  extended_mesh_handle = 0;
  factor_handle = 0;
  extended_edge_count = 0;
  interior_count = 0;
  toroidal_mode = 0;
  extended_dof = NULL;
  extended_sign = NULL;
  mixed_col_ptr = NULL;
  mixed_row_ind = NULL;
  mixed_values = NULL;
  extended_rt_dofs = NULL;
  load = NULL;
  potential_dofs = NULL;
  interior_rhs = NULL;
  interior_solution = NULL;
}

void mephit_fortfem_deinit(void)
{
  int status;

  release_prepared();
  if (mesh_handle != 0) {
    fortfem_triangle_mesh_free(mesh_handle, &status);
  }
  free(core_vertices);
  free(core_triangles);
  free(core_edge_nodes);
  free(fort_dof);
  free(edge_sign);
  free(fort_dofs);
  mesh_handle = 0;
  point_count = 0;
  triangle_count = 0;
  edge_count = 0;
  core_vertices = NULL;
  core_triangles = NULL;
  core_edge_nodes = NULL;
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

  core_vertices =
    malloc(2U * (size_t) npoint * sizeof(*core_vertices));
  core_triangles =
    malloc(3U * (size_t) ntri * sizeof(*core_triangles));
  core_edge_nodes =
    malloc(2U * (size_t) nedge * sizeof(*core_edge_nodes));
  edges = malloc(2U * (size_t) nedge * sizeof(*edges));
  triangle_edge_dofs =
    malloc(3U * (size_t) ntri * sizeof(*triangle_edge_dofs));
  triangle_edge_signs =
    malloc(3U * (size_t) ntri * sizeof(*triangle_edge_signs));
  fort_dof = malloc((size_t) nedge * sizeof(*fort_dof));
  edge_sign = malloc((size_t) nedge * sizeof(*edge_sign));
  fort_dofs = malloc((size_t) nedge * sizeof(*fort_dofs));
  if (core_vertices == NULL || core_triangles == NULL ||
      core_edge_nodes == NULL || edges == NULL ||
      triangle_edge_dofs == NULL || triangle_edge_signs == NULL ||
      fort_dof == NULL || edge_sign == NULL || fort_dofs == NULL) {
    free(edges);
    free(triangle_edge_dofs);
    free(triangle_edge_signs);
    mephit_fortfem_deinit();
    return -2;
  }

  for (edge = 0; edge < npoint; ++edge) {
    core_vertices[2 * edge] = node_R[edge];
    core_vertices[2 * edge + 1] = node_Z[edge];
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
    core_triangles[3 * triangle] = a;
    core_triangles[3 * triangle + 1] = b;
    core_triangles[3 * triangle + 2] = c;
  }

  fortfem_triangle_mesh_create(
    npoint, core_vertices, ntri, core_triangles,
    &mesh_handle, &nedge_returned, &status);
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

    core_edge_nodes[2 * edge] = start;
    core_edge_nodes[2 * edge + 1] = end;
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

  point_count = npoint;
  triangle_count = ntri;
  edge_count = nedge;
  status = 0;

cleanup:
  free(edges);
  free(triangle_edge_dofs);
  free(triangle_edge_signs);
  if (status != 0) {
    mephit_fortfem_deinit();
  }
  return status;
}

int mephit_fortfem_prepare(const int mode)
{
  const long long boundary_count =
    2LL * (long long) edge_count - 3LL * (long long) triangle_count;
  double *outer_vertices = NULL;
  double *operator_values = NULL;
  double radial_min;
  double radial_max;
  double vertical_min;
  double vertical_max;
  fortfem_complex *interior_values = NULL;
  int *edges = NULL;
  int *operator_col_ptr = NULL;
  int *operator_row_ind = NULL;
  int *triangle_edge_dofs = NULL;
  int *triangle_edge_signs = NULL;
  int edge;
  int extended_triangle_count;
  int extended_vertex_count;
  int found;
  int interior_nnz;
  int matrix_dofs;
  int matrix_nnz;
  int nnz_capacity;
  int outer_count;
  int status;

  release_prepared();
  if (mesh_handle == 0 || boundary_count < 3 ||
      boundary_count > INT_MAX ||
      triangle_count > INT_MAX / 9) {
    return -1;
  }
  outer_count =
    (int) floor(0.125 * (double) boundary_count + 0.5);
  if (outer_count < 8) {
    outer_count = 8;
  }
  outer_vertices =
    malloc(2U * (size_t) outer_count * sizeof(*outer_vertices));
  if (outer_vertices == NULL) {
    return -2;
  }

  radial_min = core_vertices[0];
  radial_max = core_vertices[0];
  vertical_min = core_vertices[1];
  vertical_max = core_vertices[1];
  for (edge = 1; edge < point_count; ++edge) {
    radial_min = fmin(radial_min, core_vertices[2 * edge]);
    radial_max = fmax(radial_max, core_vertices[2 * edge]);
    vertical_min = fmin(vertical_min, core_vertices[2 * edge + 1]);
    vertical_max = fmax(vertical_max, core_vertices[2 * edge + 1]);
  }
  for (edge = 0; edge < outer_count; ++edge) {
    const double angle =
      2.0 * acos(-1.0) * (double) edge / (double) outer_count;
    const double radial_midpoint = 0.5 * (radial_min + radial_max);
    const double vertical_midpoint = 0.5 * (vertical_min + vertical_max);
    const double radial_radius = radial_max - radial_min;
    const double vertical_radius = vertical_max - vertical_min;

    outer_vertices[2 * edge] =
      radial_midpoint + radial_radius * cos(angle);
    outer_vertices[2 * edge + 1] =
      vertical_midpoint + vertical_radius * sin(angle);
  }

  fortfem_triangle_mesh_extend(
    mesh_handle, outer_count, outer_vertices, &extended_mesh_handle,
    &extended_vertex_count, &extended_triangle_count,
    &extended_edge_count, &interior_count, &status);
  free(outer_vertices);
  outer_vertices = NULL;
  if (status != 0) {
    goto cleanup;
  }
  if (extended_edge_count < edge_count ||
      interior_count < edge_count ||
      extended_triangle_count > INT_MAX / 9) {
    status = -3;
    goto cleanup;
  }

  edges =
    malloc(2U * (size_t) extended_edge_count * sizeof(*edges));
  triangle_edge_dofs = malloc(
    3U * (size_t) extended_triangle_count * sizeof(*triangle_edge_dofs));
  triangle_edge_signs = malloc(
    3U * (size_t) extended_triangle_count * sizeof(*triangle_edge_signs));
  extended_dof = malloc((size_t) edge_count * sizeof(*extended_dof));
  extended_sign = malloc((size_t) edge_count * sizeof(*extended_sign));
  if (edges == NULL || triangle_edge_dofs == NULL ||
      triangle_edge_signs == NULL || extended_dof == NULL ||
      extended_sign == NULL) {
    status = -2;
    goto cleanup;
  }
  fortfem_triangle_mesh_edges(
    extended_mesh_handle, extended_edge_count, &matrix_dofs, edges,
    triangle_edge_dofs, triangle_edge_signs, &status);
  if (status != 0 || matrix_dofs != extended_edge_count) {
    if (status == 0) {
      status = -3;
    }
    goto cleanup;
  }
  for (edge = 0; edge < edge_count; ++edge) {
    int fort_edge;
    const int start = core_edge_nodes[2 * edge];
    const int end = core_edge_nodes[2 * edge + 1];

    found = 0;
    for (fort_edge = 0; fort_edge < extended_edge_count; ++fort_edge) {
      if (start == edges[2 * fort_edge] &&
          end == edges[2 * fort_edge + 1]) {
        extended_dof[edge] = fort_edge;
        extended_sign[edge] = 1;
        found = 1;
        break;
      }
      if (start == edges[2 * fort_edge + 1] &&
          end == edges[2 * fort_edge]) {
        extended_dof[edge] = fort_edge;
        extended_sign[edge] = -1;
        found = 1;
        break;
      }
    }
    if (!found) {
      status = -4;
      goto cleanup;
    }
  }

  nnz_capacity = 9 * extended_triangle_count;
  operator_col_ptr = malloc(
    ((size_t) extended_edge_count + 1U) * sizeof(*operator_col_ptr));
  operator_row_ind =
    malloc((size_t) nnz_capacity * sizeof(*operator_row_ind));
  operator_values =
    malloc((size_t) nnz_capacity * sizeof(*operator_values));
  if (operator_col_ptr == NULL || operator_row_ind == NULL ||
      operator_values == NULL) {
    status = -2;
    goto cleanup;
  }
  toroidal_mode = mode != 0 ? mode : 2;
  fortfem_nedelec_axisymmetric_fourier_csc(
    extended_mesh_handle, toroidal_mode, 7, nnz_capacity,
    &matrix_dofs, &matrix_nnz, operator_col_ptr,
    operator_row_ind, operator_values, &status);
  if (status != 0 || matrix_dofs != extended_edge_count) {
    if (status == 0) {
      status = -3;
    }
    goto cleanup;
  }

  mixed_col_ptr = malloc(
    ((size_t) extended_edge_count + 1U) * sizeof(*mixed_col_ptr));
  mixed_row_ind =
    malloc((size_t) nnz_capacity * sizeof(*mixed_row_ind));
  mixed_values =
    malloc((size_t) nnz_capacity * sizeof(*mixed_values));
  if (mixed_col_ptr == NULL || mixed_row_ind == NULL ||
      mixed_values == NULL) {
    status = -2;
    goto cleanup;
  }
  fortfem_nedelec_rt0_mass_csc(
    extended_mesh_handle, 2, nnz_capacity, &matrix_dofs, &matrix_nnz,
    mixed_col_ptr, mixed_row_ind, mixed_values, &status);
  if (status != 0 || matrix_dofs != extended_edge_count) {
    if (status == 0) {
      status = -3;
    }
    goto cleanup;
  }

  {
    int *interior_col_ptr = malloc(
      ((size_t) interior_count + 1U) * sizeof(*interior_col_ptr));
    int *interior_row_ind =
      malloc((size_t) nnz_capacity * sizeof(*interior_row_ind));

    interior_values =
      malloc((size_t) nnz_capacity * sizeof(*interior_values));
    if (interior_col_ptr == NULL || interior_row_ind == NULL ||
        interior_values == NULL) {
      free(interior_col_ptr);
      free(interior_row_ind);
      status = -2;
      goto cleanup;
    }
    interior_nnz = 0;
    for (edge = 0; edge < interior_count; ++edge) {
      int entry;

      interior_col_ptr[edge] = interior_nnz;
      for (entry = operator_col_ptr[edge];
           entry < operator_col_ptr[edge + 1]; ++entry) {
        if (operator_row_ind[entry] < interior_count) {
          interior_row_ind[interior_nnz] = operator_row_ind[entry];
          interior_values[interior_nnz] = operator_values[entry];
          ++interior_nnz;
        }
      }
    }
    interior_col_ptr[interior_count] = interior_nnz;
    fortfem_complex_factor_csc(
      interior_count, interior_nnz, interior_col_ptr,
      interior_row_ind, interior_values, &factor_handle, &status);
    free(interior_col_ptr);
    free(interior_row_ind);
    if (status != 0) {
      goto cleanup;
    }
  }

  extended_rt_dofs = calloc(
    (size_t) extended_edge_count, sizeof(*extended_rt_dofs));
  load = malloc((size_t) extended_edge_count * sizeof(*load));
  potential_dofs = calloc(
    (size_t) extended_edge_count, sizeof(*potential_dofs));
  interior_rhs = malloc((size_t) interior_count * sizeof(*interior_rhs));
  interior_solution =
    malloc((size_t) interior_count * sizeof(*interior_solution));
  if (extended_rt_dofs == NULL || load == NULL ||
      potential_dofs == NULL || interior_rhs == NULL ||
      interior_solution == NULL) {
    status = -2;
    goto cleanup;
  }

  free(edges);
  free(triangle_edge_dofs);
  free(triangle_edge_signs);
  free(operator_col_ptr);
  free(operator_row_ind);
  free(operator_values);
  free(interior_values);
  return 0;

cleanup:
  free(outer_vertices);
  free(edges);
  free(triangle_edge_dofs);
  free(triangle_edge_signs);
  free(operator_col_ptr);
  free(operator_row_ind);
  free(operator_values);
  free(interior_values);
  release_prepared();
  return status;
}

int mephit_fortfem_solve(
  const int nedge,
  const int npoint,
  const double complex *current_dofs,
  double complex *magnetic_dofs,
  double complex *potential_R,
  double complex *potential_Z)
{
  const double rhs_scale = 4.0 * acos(-1.0) / 29979245800.0;
  int edge;
  int point;
  int status;

  if (factor_handle == 0 || nedge != edge_count ||
      npoint != point_count || current_dofs == NULL ||
      magnetic_dofs == NULL || potential_R == NULL ||
      potential_Z == NULL) {
    return -1;
  }
  for (edge = 0; edge < extended_edge_count; ++edge) {
    extended_rt_dofs[edge] = 0.0;
    load[edge] = 0.0;
    potential_dofs[edge] = 0.0;
  }
  for (edge = 0; edge < edge_count; ++edge) {
    extended_rt_dofs[extended_dof[edge]] =
      (double) extended_sign[edge] * current_dofs[edge];
  }
  for (edge = 0; edge < extended_edge_count; ++edge) {
    int entry;

    for (entry = mixed_col_ptr[edge];
         entry < mixed_col_ptr[edge + 1]; ++entry) {
      load[mixed_row_ind[entry]] +=
        mixed_values[entry] * extended_rt_dofs[edge];
    }
  }
  for (edge = 0; edge < interior_count; ++edge) {
    interior_rhs[edge] = rhs_scale * load[edge];
  }
  fortfem_complex_solve(
    factor_handle, interior_count, interior_rhs,
    interior_solution, &status);
  if (status != 0) {
    return status;
  }
  for (edge = 0; edge < interior_count; ++edge) {
    potential_dofs[edge] = interior_solution[edge];
  }
  for (edge = 0; edge < edge_count; ++edge) {
    magnetic_dofs[edge] =
      (double) extended_sign[edge] * I * (double) toroidal_mode *
      potential_dofs[extended_dof[edge]];
  }
  for (point = 0; point < point_count; ++point) {
    fortfem_complex curl;
    fortfem_complex value[2];

    fortfem_nedelec_evaluate_vertex(
      extended_mesh_handle, point, extended_edge_count,
      potential_dofs, value, &curl, &status);
    if (status != 0) {
      return status;
    }
    potential_R[point] = value[0];
    potential_Z[point] = value[1];
  }
  return 0;
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
