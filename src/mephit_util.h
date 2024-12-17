#ifndef MEPHIT_UTIL_H
#define MEPHIT_UTIL_H

#ifdef __cplusplus
#include <cstdarg>
#else
#include <stdarg.h>
#endif

#define path_max 1024

#ifdef __cplusplus
extern "C" {
#endif

void timestamp(char *buffer);
void errno_msg(void (*exit_func)(int), const char *file, int line, int errnum, const char *msg_fmt, ...);
void gsl_errno_msg(const char *reason, const char *file, int line, int gsl_errno);

#ifdef __cplusplus
}
#endif

#endif  // MEPHIT_UTIL_H
