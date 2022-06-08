#!/bin/bash
scriptdir=$(realpath $(dirname $0))
shared_namedpipe="$1"

# check whether FreeFem++ runs with MPI and
# use subshell and change error code to suppress error message
(FreeFem++-mpi -check_plugin MUMPS || false) > /dev/null 2>&1
err_MPI=$?

# when connected via SSH or not working in an X environmnent,
# suppress graphics to avoid setting an error code
graphical=1
if [ -n "$SSH_CLIENT" ]; then
    echo "Environment variable SSH_CLIENT is set and not empty, disabling FreeFem++ graphics."
    graphical=0
fi
if [ -z "$DISPLAY" ]; then
    echo "Environment variable DISPLAY is empty or not set, disabling FreeFem++ graphics."
    graphical=0
fi

ff_cmd=FreeFem++
if [ $err_MPI -eq 0 ]; then
    ff_cmd="$ff_cmd-mpi"
fi
ff_args=-nw
if [ $graphical -ne 0 ]; then
    ff_args="$ff_args -wg -wait"
fi
ff_script="$scriptdir/maxwell_daemon.edp"
ff_script_args="-P $shared_namedpipe"

exec $ff_cmd $ff_args $ff_script $ff_script_args
exit $?
