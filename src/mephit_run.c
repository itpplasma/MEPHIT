#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <signal.h>
#include <gsl/gsl_errno.h>
#include "mephit_util.h"

static volatile sig_atomic_t caught_signal = 0;

static void catch_signal(int signum)
{
  if (!caught_signal) {
    caught_signal = (sig_atomic_t) signum;
  }
}

typedef struct child_process {
  pid_t pid;
  int exited;
  int status;
} child_process_t;
static child_process_t fem;
static child_process_t mephit;

void wait_for_exit(child_process_t *child, const char *name)
{
  siginfo_t siginfo;

  if (!child || child->exited) {
    return;
  }
  if (waitid(P_PID, (id_t) child->pid, &siginfo, WEXITED | WNOHANG)) {
    child->exited = 1;
    child->status = -1;
    fprintf(stderr, "Failed to wait for process %s to exit: %s.\n", name ? name : "", strerror(errno));
    return;
  }
  if (siginfo.si_signo == SIGCHLD && siginfo.si_pid == child->pid) {
    child->pid = (pid_t) 0;
    child->exited = 1;
    child->status = siginfo.si_status;
    switch (siginfo.si_code) {
    case CLD_EXITED:
      if (siginfo.si_status) {
        fprintf(stderr, "Process %s exited with error code %i.\n", name ? name : "", siginfo.si_status);
      }
      break;
    case CLD_KILLED:
    case CLD_DUMPED:
      fprintf(stderr, "Process %s terminated abnormally (signal %s).\n", name ? name : "", strsignal(siginfo.si_status));
      break;
    default:
      fprintf(stderr, "Process %s received signal %s.\n", name ? name : "", strsignal(siginfo.si_status));
    }
  }
}

extern void mephit_run(const int runmode, const char *config_file, const char *suffix);
extern char shared_namedpipe[path_max];

int main(int argc, char *argv[])
{
  char *config = NULL, *suffix = NULL, *tmpdir = NULL, *scriptpath = NULL;
  int argi = 0, runmode = 0, bytes_written, errno_save;
  struct sigaction infanticide, prev_sigint, prev_sigterm;

  // first argument
  if (argc > ++argi) {
    if (!strlen(argv[argi])) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "Empty string for numeric runmode");
    }
    runmode = (int) strtol(argv[argi], NULL, 0);
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected numeric runmode as first argument");
  }
  // second argument
  if (argc > ++argi) {
    if (!strlen(argv[argi])) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "Empty string for config file");
    }
    config = argv[argi];
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected path to config file as second argument");
  }
  // third argument
  if (argc > ++argi) {
    suffix = argv[argi];
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected file basename suffix as third argument");
  }
  // fourth argument
  if (argc > ++argi) {
    if (!strlen(argv[argi])) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "Empty string for temporary directory");
    }
    tmpdir = argv[argi];
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected path to temporary directory as fourth argument");
  }
  // fifth argument
  if (argc > ++argi) {
    if (!strlen(argv[argi])) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "Empty string for FreeFem script file");
    }
    scriptpath = argv[argi];
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected path to FreeFem script file as fifth argument");
  }
  /* TODO: if exclusive flock on mephit.h5 fails, exit with error; else, release acquired lock */
  bytes_written = snprintf(shared_namedpipe, path_max, "%s/MEPHIT_0x%.8x.dat", tmpdir, getpid());
  if (bytes_written < 0) {
    errno_msg(exit, __FILE__, __LINE__, errno, "Failed to generate FIFO name");
  } else if (bytes_written >= path_max) {
    errno_msg(exit, __FILE__, __LINE__, ENAMETOOLONG, "Failed to generate FIFO name");
  }
  if (mkfifo(shared_namedpipe, 0600)) {
    errno_msg(exit, __FILE__, __LINE__, errno, "Failed to create FIFO %s", shared_namedpipe);
  }
  fem.pid = fork();
  if (fem.pid == (pid_t) 0) {
    execl(scriptpath, scriptpath, shared_namedpipe, suffix, (char *) NULL);
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to execute FreeFem");
  } else if (fem.pid == (pid_t) -1) {
    errno_msg(exit, __FILE__, __LINE__, errno, "Failed to create child process for FreeFem");
  }
 mephit_fork: mephit.pid = fork();
  if (mephit.pid == (pid_t) 0) {
    gsl_set_error_handler(gsl_errno_msg);
    mephit_run(runmode, config, suffix);
    exit(0);
  } else if (mephit.pid == (pid_t) -1) {
    errno_save = errno;
    if (kill(fem.pid, SIGTERM)) {
      errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to kill process FreeFem");
    }
    errno_msg(exit, __FILE__, __LINE__, errno_save, "Failed to create child process for MEPHIT");
  }
  infanticide.sa_handler = catch_signal;
  infanticide.sa_flags = SA_RESTART;
  sigemptyset(&infanticide.sa_mask);
  if (sigaction(SIGINT, &infanticide, &prev_sigint)) {
    errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to register signal handler for SIGINT");
  }
  if (sigaction(SIGTERM, &infanticide, &prev_sigterm)) {
    errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to register signal handler for SIGTERM");
  }
  do {
    if (caught_signal) {
      if (fem.pid) {
        if (kill(fem.pid, caught_signal)) {
          errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to signal process FreeFem");
        }
      }
      if (mephit.pid) {
        if (kill(mephit.pid, caught_signal)) {
          errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to signal process MEPHIT");
        }
      }
    }
    wait_for_exit(&fem, "FreeFem");
    wait_for_exit(&mephit, "MEPHIT");
    if (fem.exited && fem.status && !mephit.exited) {
      if (kill(mephit.pid, SIGTERM)) {
        errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to terminate process MEPHIT");
      }
    }
    if (mephit.exited && mephit.status && !fem.exited) {
      if (kill(fem.pid, SIGTERM)) {
        errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to terminate process FreeFem");
      }
    }
  } while (!(fem.exited && mephit.exited));
  if (sigaction(SIGINT, &prev_sigint, NULL)) {
    errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to unregister signal handler for SIGINT");
  }
  if (sigaction(SIGTERM, &prev_sigterm, NULL)) {
    errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to unregister signal handler for SIGTERM");
  }
  if (unlink(shared_namedpipe)) {
    errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to delete FIFO");
  }
  shared_namedpipe[0] = '\0';

  if (caught_signal) {
    return 128 + (int) caught_signal;
  }
  if (fem.status || mephit.status) {
    return 1;
  }
  return 0;
}
