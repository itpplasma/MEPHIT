#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include "triangle.h"
#include "mephit_util.h"
#include "mephit_fem.h"

void gauss_legendre_unit_interval(int order, double *points, double *weights)
{
  gsl_integration_glfixed_table *table;
  size_t i, n;

  n = (size_t) order;
  table = gsl_integration_glfixed_table_alloc(n);
  for (i = 0; i < n; ++i) {
    gsl_integration_glfixed_point(0.0, 1.0, i, &points[i], &weights[i], table);
  }
  gsl_integration_glfixed_table_free(table);
}

void FEM_triangulate_external(const int npt_inner,
                              const int npt_outer,
                              const double *bdry_R,
                              const double *bdry_Z,
                              const double R_mid,
                              const double Z_mid,
                              const char *fname)
{
  int k;
  FILE *fid;
  struct triangulateio in, out, vorout;

  // initialize all fields to zero or NULL
  memset(&in, 0, sizeof(struct triangulateio));
  memset(&out, 0, sizeof(struct triangulateio));
  memset(&vorout, 0, sizeof(struct triangulateio));

  in.numberofpoints = npt_inner + npt_outer;
  in.numberofsegments = npt_inner + npt_outer;
  in.numberofholes = 1;
  in.pointlist = (REAL *) calloc(2 * (size_t) in.numberofpoints, sizeof(REAL));
  in.pointmarkerlist = (int *) calloc((size_t) in.numberofpoints, sizeof(int));
  in.segmentlist = (int *) calloc(2 * (size_t) in.numberofsegments, sizeof(int));
  in.holelist = (REAL *) calloc(2 * (size_t) in.numberofholes, sizeof(REAL));
  for (k = 0; k < in.numberofpoints; ++k) {
    in.pointlist[2 * k] = bdry_R[k];
    in.pointlist[2 * k + 1] = bdry_Z[k];
  }
  for (k = 0; k < npt_inner; ++k) {
    in.segmentlist[2 * k] = k;
    in.segmentlist[2 * k + 1] = (k + 1) % npt_inner;
  }
  for (k = 0; k < npt_outer; ++k) {
    in.segmentlist[2 * (npt_inner + k)] = npt_inner + k;
    in.segmentlist[2 * (npt_inner + k) + 1] = npt_inner + (k + 1) % npt_outer;
  }
  in.holelist[0] = R_mid;
  in.holelist[1] = Z_mid;

  // triangulate options:
  // B - omit boundary markers in output
  // e - generate edge list (to be used later)
  // j - clean point list
  // n - generate neighbor list (to be used later)
  // p - triangulate from given polygon boundary
  // q - minimum angle of 20 degrees
  // Y - don't modify boundary edges
  // z - use zero indexing
  triangulate("BejnpqYz", &in, &out, &vorout);

  fid = fopen(fname, "w");
  fprintf(fid, "MFEM mesh v1.0\ndimension\n2\n");
  fprintf(fid, "elements\n%i\n", out.numberoftriangles);
  for (k = 0; k < out.numberoftriangles; ++k) {
    // <element attribute> <geometry type> <vertex indices ...>
    // attribute: 2 for surrounding volume, geometry type: 2 for triangle
    fprintf(fid, "2 2 %i %i %i\n",
            out.trianglelist[3 * k],
            out.trianglelist[3 * k + 1],
            out.trianglelist[3 * k + 2]);
  }
  fprintf(fid, "boundary\n%i\n", out.numberofsegments);
  for (k = 0; k < out.numberofsegments; ++k) {
    // <boundary element attribute> <geometry type> <vertex indices ...>
    // attribute: 2 for surrounding volume, geometry type: 1 for segment
    fprintf(fid, "2 1 %i %i\n",
            out.segmentlist[2 * k],
            out.segmentlist[2 * k + 1]);
  }
  fprintf(fid, "vertices\n%i\n2\n", out.numberofpoints);
  for (k = 0; k < out.numberofpoints; ++k) {
    fprintf(fid, "%.16e %.16e\n",
            out.pointlist[2 * k],
            out.pointlist[2 * k + 1]);
  }
  fclose(fid);

  free(in.pointlist);
  free(in.pointmarkerlist);
  free(in.segmentlist);
  free(in.holelist);  // same as out.holelist
  trifree(out.pointlist);
  trifree(out.trianglelist);
  trifree(out.neighborlist);
  trifree(out.segmentlist);
  trifree(out.edgelist);
}
