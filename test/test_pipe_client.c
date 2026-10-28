#include <signal.h>
#include <string.h>
#include "mephit_fem.h"

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
