#!/bin/bash
scriptname=${0##*/}
scriptdir=$(dirname "$0")
scriptdir=$(realpath "$scriptdir")
shared_namedpipe=$1
suffix=$2

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

ff_cmd="@FREEFEM_EXECUTABLE@"
ff_args=( "-ne" "-fglut" "fglut_dump$suffix" )
if [ $graphical -ne 0 ]; then
    ff_args=( "-wg" "-wait" "${ff_args[@]}" )
else
    ff_args=( "-nw" "-nowait" "${ff_args[@]}" )
fi
ff_script=$scriptdir/maxwell_daemon.edp
# The real Maxwell factors are always project-owned native FortSparse.
# Missing/incompatible plugin is a hard load error; no external fallback.
export FF_LOADPATH="$scriptdir;${FF_LOADPATH:-}"
ff_script_args=( "-P" "$shared_namedpipe" "-S" "$suffix" )

# Keep daemon diagnostics independently readable, including binary/NUL messages.
# exec preserves the supervised PID; tee also keeps diagnostics in mephit.log.
exec $ff_cmd "${ff_args[@]}" "$ff_script" "${ff_script_args[@]}" \
    > >(tee "freefem$suffix.log") 2>&1
exit $?
