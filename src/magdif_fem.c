#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "magdif_util.h"
#include "magdif_fem.h"

char shared_namedpipe[path_max];

void send_long0_to_FreeFem(const char *namedpipe, const long int *long0)
{
  int fd;
  ssize_t bytes_written;

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
  bytes_written = write(fd, (void *) long0, sizeof(long int));
  if (bytes_written < (ssize_t) sizeof(long int)) {
    errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO, "Failed to write to pipe %s", namedpipe);
  }
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close write end of pipe %s", namedpipe);
  }
}

void receive_long0_from_FreeFem(const char *namedpipe, long int *long0)
{
  int fd;
  ssize_t bytes_read;

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
  bytes_read = read(fd, (void *) long0, sizeof(long int));
  if (bytes_read < (ssize_t) sizeof(long int)) {
    errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO, "Failed to read from pipe %s", namedpipe);
  }
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close read end of pipe %s", namedpipe);
  }
}

void send_double1_to_FreeFem(const char *namedpipe, const int size, const double *double1)
{
  long int dum = size;
  int fd;
  ssize_t bytes_written;

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
  bytes_written = write(fd, (void *) &dum, sizeof(long int));
  if (bytes_written < (ssize_t) sizeof(long int)) {
    errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO, "Failed to write to pipe %s", namedpipe);
  }
  bytes_written = write(fd, (void *) double1, size * sizeof(double));
  if (bytes_written < (ssize_t) (size * sizeof(double))) {
    errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO, "Failed to write to pipe %s", namedpipe);
  }
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close write end of pipe %s", namedpipe);
  }
}

void receive_double1_from_FreeFem(const char *namedpipe, const int size, double *double1)
{
  long int dum;
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
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to open pipe %s for writing", namedpipe);
  }
  bytes_read = read(fd, (void *) &dum, sizeof(long int));
  if (bytes_read < (ssize_t) sizeof(long int)) {
    errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO, "Failed to read from pipe %s", namedpipe);
  }
  if (dum < size) {
    errno_msg(_exit, __FILE__, __LINE__, EIO, "Pipe %s only contains %li double precision values, "
            "expected %i.\n", namedpipe, dum, size);
  }
  bytes_expected = (ssize_t) (size * sizeof(double));
  total_bytes_read = (ssize_t) 0;
  do {
    bytes_read = read(fd, (void *) double1 + total_bytes_read, bytes_expected - total_bytes_read);
    total_bytes_read += bytes_read;
  } while (bytes_read != 0 && total_bytes_read < bytes_expected);
  if (total_bytes_read < bytes_expected) {
    errno_msg(_exit, __FILE__, __LINE__, errno ? errno : EIO, "Received %li bytes from pipe %s, "
              "expected %li.\n", total_bytes_read, namedpipe, bytes_expected);
  }
  if (close(fd)) {
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close read end of pipe %s", namedpipe);
  }
}

void FEM_init(const int tormode, const int runmode)
{
  long int n = tormode, long_runmode = runmode;
  long int flag;
 
  send_long0_to_FreeFem(shared_namedpipe, &n);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
  send_long0_to_FreeFem(shared_namedpipe, &long_runmode);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
}

void FEM_extend_mesh(void)
{
  long int flag = 0L;

  send_long0_to_FreeFem(shared_namedpipe, &flag);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
}

void FEM_compute_Bn(const int *restrict shape, const complex_double *restrict Jn, complex_double *restrict Bn)
{
  long int flag = -1L;
  int size = shape ? 2 * shape[0] * shape[1] : 0;

  send_long0_to_FreeFem(shared_namedpipe, &flag);
  receive_long0_from_FreeFem(shared_namedpipe, &flag);
  send_double1_to_FreeFem(shared_namedpipe, size, (double *) Jn);
  receive_double1_from_FreeFem(shared_namedpipe, size, (double *) Bn);
}

void FEM_compute_L2int(const int *restrict shape, const complex_double *restrict elem, double *restrict L2int)
{
  long int flag = -2L;
  int size = shape ? 2 * shape[0] * shape[1] : 0;

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
