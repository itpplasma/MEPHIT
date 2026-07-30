#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include "triangle.h"
#include "mephit_util.h"
#include "mephit_fem.h"

char shared_namedpipe[path_max];

void send_long0_to_FreeFem(const char *namedpipe, const long int *long0)
{
  int fd;
  ssize_t bytes_written, total_bytes_written, bytes_expected;

  if (!(namedpipe && strlen(namedpipe))) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for namedpipe");
  }
  if (!long0) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL value for long0");
  }
  fd = open(namedpipe, O_WRONLY);
  if (fd < 0) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to open pipe %s for writing", namedpipe);
  }
  bytes_expected = (ssize_t) sizeof(long int);
  total_bytes_written = (ssize_t) 0;
  do {
    bytes_written = write(fd, (char *) long0 + total_bytes_written,
                          (size_t) (bytes_expected - total_bytes_written));
    if (bytes_written < (ssize_t) 0) {
      errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO,
                "Failed to write to pipe %s", namedpipe);
    }
    total_bytes_written += bytes_written;
  } while (total_bytes_written < bytes_expected);
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close write end of pipe %s", namedpipe);
  }
}

void receive_long0_from_FreeFem(const char *namedpipe, long int *long0)
{
  int fd;
  ssize_t bytes_read, total_bytes_read, bytes_expected;

  if (!(namedpipe && strlen(namedpipe))) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for namedpipe");
  }
  if (!long0) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL value for long0");
  }
  fd = open(namedpipe, O_RDONLY);
  if (fd < 0) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to open pipe %s for reading", namedpipe);
  }
  bytes_expected = (ssize_t) sizeof(long int);
  total_bytes_read = (ssize_t) 0;
  do {
    bytes_read = read(fd, (char *) long0 + total_bytes_read,
                      (size_t) (bytes_expected - total_bytes_read));
    if (bytes_read < (ssize_t) 0) {
      errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO,
                "Failed to read from pipe %s", namedpipe);
    }
    total_bytes_read += bytes_read;
  } while (total_bytes_read < bytes_expected);
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close read end of pipe %s", namedpipe);
  }
}

void send_double1_to_FreeFem(const char *namedpipe, const int size, const double *double1)
{
  long int announced_size = size;
  int fd;
  ssize_t bytes_written, total_bytes_written, bytes_expected;

  if (!(namedpipe && strlen(namedpipe))) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for namedpipe");
  }
  if (!double1) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL value for double1");
  }
  fd = open(namedpipe, O_WRONLY);
  if (fd < 0) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to open pipe %s for writing", namedpipe);
  }
  bytes_expected = (ssize_t) sizeof(long int);
  total_bytes_written = (ssize_t) 0;
  do {
    bytes_written = write(fd, (char *) &announced_size + total_bytes_written,
                          (size_t) (bytes_expected - total_bytes_written));
    if (bytes_written < (ssize_t) 0) {
      errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO,
                "Failed to write to pipe %s", namedpipe);
    }
    total_bytes_written += bytes_written;
  } while (total_bytes_written < bytes_expected);
  bytes_expected = (ssize_t) ((unsigned) size * sizeof(double));
  total_bytes_written = (ssize_t) 0;
  do {
    bytes_written = write(fd, (char *) double1 + total_bytes_written,
                          (size_t) (bytes_expected - total_bytes_written));
    if (bytes_written < (ssize_t) 0) {
      errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO,
                "Failed to write to pipe %s", namedpipe);
    }
    total_bytes_written += bytes_written;
  } while (total_bytes_written < bytes_expected);
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close write end of pipe %s", namedpipe);
  }
}

void receive_double1_from_FreeFem(const char *namedpipe, const int size, double *double1)
{
  long int announced_size;
  int fd;
  ssize_t bytes_read, total_bytes_read, bytes_expected;

  if (!(namedpipe && strlen(namedpipe))) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for namedpipe");
  }
  if (!double1) {
    errno_msg(_exit, __FILE__, __LINE__, EINVAL, "NULL value for double1");
  }
  fd = open(namedpipe, O_RDONLY);
  if (fd < 0) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to open pipe %s for reading", namedpipe);
  }
  bytes_expected = (ssize_t) sizeof(long int);
  total_bytes_read = (ssize_t) 0;
  do {
    bytes_read = read(fd, (char *) &announced_size + total_bytes_read,
                      (size_t) (bytes_expected - total_bytes_read));
    if (bytes_read < (ssize_t) 0) {
      errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO,
                "Failed to read from pipe %s", namedpipe);
    }
    total_bytes_read += bytes_read;
  } while (total_bytes_read < bytes_expected);
  if (announced_size != size) {
    errno_msg(_exit, __FILE__, __LINE__, EIO,
              "Pipe %s contains %li double precision values, expected %i",
              namedpipe, announced_size, size);
  }
  bytes_expected = (ssize_t) ((unsigned) size * sizeof(double));
  total_bytes_read = (ssize_t) 0;
  do {
    bytes_read = read(fd, (char *) double1 + total_bytes_read,
                      (size_t) (bytes_expected - total_bytes_read));
    if (bytes_read < (ssize_t) 0) {
      errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO,
                "Failed to read from pipe %s", namedpipe);
    }
    total_bytes_read += bytes_read;
  } while (total_bytes_read < bytes_expected);
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close read end of pipe %s", namedpipe);
  }
}

void FEM_init(const int tormode, const int nedge, const int npoint, const int runmode)
{
  long int long_tormode = tormode ? tormode : 2;
  long int long_nedge = nedge;
  long int long_npoint = npoint;
  long int long_runmode = runmode ? runmode : (1 << 0 | 1 << 1 | 1 << 2);

  send_long0_to_FreeFem(shared_namedpipe, &long_tormode);
  receive_long0_from_FreeFem(shared_namedpipe, &long_tormode);
  send_long0_to_FreeFem(shared_namedpipe, &long_nedge);
  receive_long0_from_FreeFem(shared_namedpipe, &long_nedge);
  send_long0_to_FreeFem(shared_namedpipe, &long_npoint);
  receive_long0_from_FreeFem(shared_namedpipe, &long_npoint);
  send_long0_to_FreeFem(shared_namedpipe, &long_runmode);
  receive_long0_from_FreeFem(shared_namedpipe, &long_runmode);
}

void FEM_extend_mesh(void)
{
  long int flag = 0L;

  send_long0_to_FreeFem(shared_namedpipe, &flag);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
}

void FEM_compute_magfn(const int nedge,
                       const int npoint,
                       const complex_double *Jn,
                       complex_double *Bn,
                       complex_double *AnR,
                       complex_double *AnZ)
{
  long int flag = -1L;
  int size = 2 * nedge;

  send_long0_to_FreeFem(shared_namedpipe, &flag);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
  send_double1_to_FreeFem(shared_namedpipe, size, (double *) Jn);
  receive_double1_from_FreeFem(shared_namedpipe, size, (double *) Bn);
  size = 2 * npoint;
  receive_double1_from_FreeFem(shared_namedpipe, size, (double *) AnR);
  receive_double1_from_FreeFem(shared_namedpipe, size, (double *) AnZ);
}

void FEM_compute_L2int(const int nedge, const complex_double *elem, double *L2int)
{
  long int flag = -2L;
  int size = 2 * nedge;

  send_long0_to_FreeFem(shared_namedpipe, &flag);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
  send_double1_to_FreeFem(shared_namedpipe, size, (double *) elem);
  receive_double1_from_FreeFem(shared_namedpipe, 1, L2int);
}

void FEM_deinit(void)
{
  long int flag = -3L;

  send_long0_to_FreeFem(shared_namedpipe, &flag);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
}

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
  fprintf(fid, "%i %i %i\n",
          out.numberofpoints,
          out.numberoftriangles,
          out.numberofsegments);
  for (k = 0; k < out.numberofpoints; ++k) {
    fprintf(fid, "%.16e %.16e 0\n",
            out.pointlist[2 * k],
            out.pointlist[2 * k + 1]);
  }
  for (k = 0; k < out.numberoftriangles; ++k) {
    fprintf(fid, "%i %i %i 1\n",
            out.trianglelist[3 * k] + 1,
            out.trianglelist[3 * k + 1] + 1,
            out.trianglelist[3 * k + 2] + 1);
  }
  for (k = 0; k < out.numberofsegments; ++k) {
    fprintf(fid, "%i %i 2\n",
            out.segmentlist[2 * k] + 1,
            out.segmentlist[2 * k + 1] + 1);
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
