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
