#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mephit_util.h"

void errno_msg(void (*exit_func)(int), const char *file, int line, int errnum, const char *msg_fmt, ...)
{
  int bytes_written = 0;
  size_t size = 0;
  char *msg = NULL;
  va_list ap;

  va_start(ap, msg_fmt);
  bytes_written = vsnprintf(msg, size, msg_fmt, ap);
  va_end(ap);
  if (bytes_written > 0) {
    size = (size_t) bytes_written + 1;
    msg = malloc(size);
    if (msg) {
      va_start(ap, msg_fmt);
      bytes_written = vsnprintf(msg, size, msg_fmt, ap);
      va_end(ap);
      if (bytes_written <= 0) {
        free(msg);
        msg = NULL;
      }
    }
  }
  if (msg) {
    fprintf(stderr, "%s:%i:%s: %s.\n", file ? file : "", line, msg, strerror(errnum));
  } else {
    fprintf(stderr, "%s:%i:Error: %s.\n", file ? file : "", line, strerror(errnum));
  }
  if (exit_func) {
    exit_func(1);
  }
}
