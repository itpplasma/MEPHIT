#!/bin/bash
scriptname=${0##*/}
scriptdir=$(dirname "$0")
scriptdir=$(realpath "$scriptdir")
shared_namedpipe=$1
suffix=$2

# check whether FreeFem++ runs with MPI and
# use subshell to suppress error message and change error code
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
ff_args=( "-ne" "-fglut" "fglut_dump$suffix" )
if [ $graphical -ne 0 ]; then
    ff_args=( "-wg" "-wait" "${ff_args[@]}" )
else
    ff_args=( "-nw" "-nowait" "${ff_args[@]}" )
fi
ff_script=$scriptdir/maxwell_daemon.edp
# MUMPS is optional; prefer 64-bit factors when falling back to UMFPACK.
if ! ($ff_cmd -nw -nowait -check_plugin MUMPS) > /dev/null 2>&1; then
    # The default 32-bit UMFPACK backend can report out_of_memory when its
    # factor indices overflow, even when the Slurm memory allowance is ample.
    fallback=""
    if ($ff_cmd -nw -nowait -check_plugin UMFPACK64) > /dev/null 2>&1; then
        fallback='load "UMFPACK64"'
        echo "$scriptname: MUMPS unavailable; using UMFPACK64 (64-bit factors)."
    else
        echo "$scriptname: MUMPS and UMFPACK64 unavailable; using the default sparse solver (large factors may exceed its index range)."
    fi
    ff_script=$(mktemp "${TMPDIR:-/tmp}/maxwell_daemon_XXXXXX.edp")
    sed -e "s/^load \"MUMPS\"/$fallback/" "$scriptdir/maxwell_daemon.edp" > "$ff_script"
fi
ff_script_args=( "-P" "$shared_namedpipe" "-S" "$suffix" )

# Keep daemon diagnostics independently readable, including binary/NUL messages.
# exec preserves the supervised PID; tee also keeps diagnostics in mephit.log.
exec $ff_cmd "${ff_args[@]}" "$ff_script" "${ff_script_args[@]}" \
    > >(tee "freefem$suffix.log") 2>&1
exit $?
