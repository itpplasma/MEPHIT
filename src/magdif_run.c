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
#include "magdif_util.h"

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

int main(int argc, char *argv[])
{
  const char ff[] = "FreeFem++-mpi";
  const char ff_log[] = "freefem.out";
  char *binpath = NULL, *config = NULL, *tmpdir = NULL, *scriptpath = NULL;
  char fifo[path_max] = "";
  int bytes_written, log_fd, errno_save;
  struct sigaction infanticide, prev_sigint, prev_sigterm;

  if (argc > 1) {
    binpath = argv[1];
    if (!(binpath && strlen(binpath))) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for Fortran binary path");
    }
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected path to Fortran binary as first argument");
  }
  if (argc > 2) {
    config = argv[2];
    if (!(config && strlen(config))) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for config file");
    }
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected path to config file as second argument");
  }
  if (argc > 3) {
    tmpdir = argv[3];
    if (!(tmpdir && strlen(tmpdir))) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for temporary directory");
    }
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected path to temporary directory as third argument");
  }
  if (argc > 4) {
    scriptpath = argv[4];
    if (!(scriptpath && strlen(scriptpath))) {
      errno_msg(exit, __FILE__, __LINE__, EINVAL, "NULL or empty value for FreeFem script file");
    }
  } else {
    errno_msg(exit, __FILE__, __LINE__, EINVAL, "Expected path to FreeFem script file as fourth argument");
  }
  /* TODO: if exclusive flock on magdif.h5 fails, exit with error; else, release acquired lock */
  bytes_written = snprintf(fifo, path_max, "%s/MEPHIT_0x%.8x.dat", tmpdir, getpid());
  if (bytes_written < 0) {
    errno_msg(exit, __FILE__, __LINE__, errno, "Failed to generate FIFO name");
  } else if (bytes_written >= path_max) {
    errno_msg(exit, __FILE__, __LINE__, ENAMETOOLONG, "Failed to generate FIFO name");
  }
  if (mkfifo(fifo, 0600)) {
    errno_msg(exit, __FILE__, __LINE__, errno, "Failed to create FIFO %s", fifo);
  }
  fem.pid = fork();
  if (fem.pid == (pid_t) 0) {
    log_fd = open(ff_log, O_WRONLY | O_CREAT, 0644);
    if (log_fd < 0) {
      errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to create FreeFem logfile");
    }
    if (0 > dup2(log_fd, STDOUT_FILENO)) {
      errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to redirect stdout");
    }
    if (0 > dup2(log_fd, STDERR_FILENO)) {
      errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to redirect stderr");
    }
    if (close(log_fd)) {
      errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to close duplicate file descriptor");
    }
    execlp(ff, ff, "-nw", scriptpath, "-P", fifo, (char *) NULL);
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to execute FreeFem");
  } else if (fem.pid == (pid_t) -1) {
    errno_msg(exit, __FILE__, __LINE__, errno, "Failed to create child process for FreeFem");
  }
  mephit.pid = fork();
  if (mephit.pid == (pid_t) 0) {
    execl(binpath, binpath, config, fifo, (char *) NULL);
    errno_msg(_exit, __FILE__, __LINE__, errno, "Failed to execute MEPHIT");
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
  if (unlink(fifo)) {
    errno_msg(NULL, __FILE__, __LINE__, errno, "Failed to delete FIFO");
  }
  fifo[0] = '\0';

  if (caught_signal) {
    return 128 + (int) caught_signal;
  }
  if (fem.status || mephit.status) {
    return 1;
  }
  return 0;
}
