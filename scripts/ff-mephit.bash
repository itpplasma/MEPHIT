#!/bin/bash
scriptname=$0
scriptdir=$(dirname "$0")
scriptdir=$(realpath "$scriptdir")
shared_namedpipe=$1

# check whether FreeFem++ runs with MPI and
# use subshell and change error code to suppress error message
(FreeFem++-mpi -check_plugin MUMPS || false) > /dev/null 2>&1
err_MPI=$?

# when connected via SSH or not working in an X environmnent,
# suppress graphics to avoid setting an error code
graphical=1
if [ -n "$SSH_CLIENT" ]; then
    echo "$scriptname: Environment variable SSH_CLIENT is set and not empty, disabling FreeFem++ graphics."
    graphical=0
fi
if [ -z "$DISPLAY" ]; then
    echo "$scriptname: Environment variable DISPLAY is empty or not set, disabling FreeFem++ graphics."
    graphical=0
fi

ff_cmd=FreeFem++
if [ $err_MPI -eq 0 ]; then
    ff_cmd=$ff_cmd-mpi
fi
ff_args=( "-fglut" "fglut_dump" )
if [ $graphical -ne 0 ]; then
    ff_args=( "-wg" "-wait" "${ff_args[@]}" )
else
    ff_args=( "-nw" "-nowait" "${ff_args[@]}" )
fi
ff_script=$scriptdir/maxwell_daemon.edp
ff_script_args=( "-P" "$shared_namedpipe" )

exec $ff_cmd "${ff_args[@]}" "$ff_script" "${ff_script_args[@]}"
exit $?
