#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "mephit_util.h"

char shared_namedpipe[path_max];

void mephit_run(const int runmode, const char *config, const char *suffix)
{
  (void) runmode;
  (void) suffix;
  if (!strcmp(config, "fail")) exit(23);
  if (!strcmp(config, "wait")) {
    for (;;) pause();
  }
}
