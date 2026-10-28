#include <signal.h>
#include <stdio.h>
#include <string.h>
#include "mephit_fem.h"
#include "mephit_util.h"
extern char shared_namedpipe[path_max];

void receive_long0_from_FreeFem(const char *, long int *);
void receive_double1_from_FreeFem(const char *, int, double *);
void send_double1_to_FreeFem(const char *, int, const double *);

int main(int argc, char **argv)
{
  long int value;
  static double data[262144] = {0.0};
  const double expected[2] = {1.25, -2.5};
  if (argc != 3) return 2;
  signal(SIGPIPE, SIG_IGN);
  if (!strcmp(argv[1], "protocol")) {
    const complex_double input[2] = {1.25 + 2.0 * I, -3.0 + 0.5 * I};
    complex_double b[2], ar[2], az[2];
    double integral;
    snprintf(shared_namedpipe, path_max, "%s", argv[2]);
    FEM_init(2, 2, 2, 2);
    for (int k = 0; k < 100; ++k) {
      FEM_compute_magfn(2, 2, input, b, ar, az);
      if (memcmp(input, b, sizeof(input)) || memcmp(input, ar, sizeof(input)) ||
          memcmp(input, az, sizeof(input))) return 3;
      FEM_compute_L2int(2, input, &integral);
      if (integral != 7.25) return 3;
    }
    FEM_deinit();
    return 0;
  }
  if (!strcmp(argv[1], "scalar")) {
    receive_long0_from_FreeFem(argv[2], &value);
    return value == 42 ? 0 : 3;
  }
  if (!strcmp(argv[1], "vector")) {
    receive_double1_from_FreeFem(argv[2], 2, data);
    return memcmp(data, expected, sizeof(expected)) == 0 ? 0 : 3;
  }
  send_double1_to_FreeFem(argv[2], 262144, data);
  return 0;
}
