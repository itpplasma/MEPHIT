#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mephit_util.h"

void timestamp(char *buffer) {
  time_t t = time(NULL);
  struct tm *tm = localtime(&t);
  snprintf(buffer, 72, "%04i-%02i-%02i %02i:%02i:%02i",
           tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday,
           tm->tm_hour, tm->tm_min, tm->tm_sec);
}

void errno_msg(void (*exit_func)(int), const char *file, int line, int errnum, const char *msg_fmt, ...)
{
  int bytes_written = 0;
  size_t size = 0;
  char *msg = NULL;
  char now[72] = "1990-06-11 20:47:00";
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
  timestamp(now);
  fprintf(stderr, "[%s] %s:%i:%s: %s.\n",
          now, file ? file : "", line, msg ? msg : "Error", strerror(errnum));
  if (exit_func) {
    exit_func(1);
  }
}
