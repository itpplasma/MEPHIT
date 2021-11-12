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

void errno_msg(void (*exit_func)(int), const char *file, int line, int errnum, const char *msg_fmt, ...);

#ifdef __cplusplus
}
#endif

#endif  // MEPHIT_UTIL_H
