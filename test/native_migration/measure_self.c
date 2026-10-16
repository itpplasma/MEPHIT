/* Task-local observation only: own post-exec Linux VmHWM and self CPU at exit.
 * No cumulative child RSS is presented as aggregate simultaneous memory. */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/resource.h>
#include <time.h>
static struct timespec started;
__attribute__((constructor)) static void begin(void) {
  clock_gettime(CLOCK_MONOTONIC,&started);
}
__attribute__((destructor)) static void finish(void) {
  const char *directory=getenv("MEPHIT_METRIC_DIR");
  if (!directory) return;
  char exe[4096]={0};
  ssize_t length=readlink("/proc/self/exe",exe,sizeof(exe)-1);
  if(length<0) return;
  const char *name=strrchr(exe,'/'); name=name?name+1:exe;
  if(strcmp(name,"mephit_run.x") && strcmp(name,"FreeFem++") &&
     strcmp(name,"FreeFem++-mpi") && strcmp(name,"mephit_sparse_probe")) return;
  long hwm=-1; char line[512];
  FILE *status=fopen("/proc/self/status","r");
  if(status) {while(fgets(line,sizeof(line),status)) {
    if(sscanf(line,"VmHWM: %ld kB",&hwm)==1) break;
  } fclose(status);}
  int runmode=-1;
  if(!strcmp(name,"mephit_run.x")) {
    FILE *args=fopen("/proc/self/cmdline","r");
    if(args) { char command[8192]; size_t n=fread(command,1,sizeof(command)-1,args);
      command[n]=0; size_t first=strnlen(command,n);
      if(first+1<n) runmode=atoi(command+first+1);
      fclose(args);
    }
  }
  struct rusage usage; getrusage(RUSAGE_SELF,&usage);
  struct timespec ended; clock_gettime(CLOCK_MONOTONIC,&ended);
  double wall=(ended.tv_sec-started.tv_sec)+1e-9*(ended.tv_nsec-started.tv_nsec);
  double user=usage.ru_utime.tv_sec+1e-6*usage.ru_utime.tv_usec;
  double system=usage.ru_stime.tv_sec+1e-6*usage.ru_stime.tv_usec;
  char path[4096]; snprintf(path,sizeof(path),"%s/%d.json",directory,getpid());
  FILE *out=fopen(path,"w"); if(!out) return;
  fprintf(out,"{\"pid\":%d,\"ppid\":%d,\"runmode\":%d,\"executable\":\"%s\",\"wall_s\":%.9f,\"user_s\":%.9f,\"system_s\":%.9f,\"cpu_s\":%.9f,\"own_postexec_vmhwm_kib\":%ld}\n",getpid(),getppid(),runmode,exe,wall,user,system,user+system,hwm);
  fclose(out);
}
