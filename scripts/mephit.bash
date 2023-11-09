#!/bin/bash

if [[ "$OSTYPE" == "darwin"* ]]; then
    PATH=/opt/local/libexec/gnubin:$PATH
fi

# if first argument is absolute path, return it unchanged
# if first argument is relative path, prefix it with second argument and return absolute path
# if second argument is empty, use current path
absolutize() {
    case "$1" in
        /*)
            echo "$1"
            ;;
        *)
            if [ -z "$2" ]; then
                realpath -m "$1"
            else
                realpath -m "$2/$1"
            fi
            ;;
    esac
}

replace_first_in_line() {
    sed -Ee "$2 s/^ *'?([^' ]*)'?/$3/" -i "$1"
}


mephit_init() {
    config=
    geqdsk=
    type=
    convexwall=
    TEMP=$(getopt -o 'c:g:t:' --long 'config:,g-eqdsk:,type:' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-c'|'--config')
                config=$2
                shift 2
                continue
                ;;
            '-g'|'--g-eqdsk')
                geqdsk=$2
                shift 2
                continue
                ;;
            '-t'|'--type')
                type=$2
                shift 2
                continue
                ;;
            '--')
                shift
                break
                ;;
            *)
                echo "$scriptname: unrecognized option '$1'" >&2
                exit 1
                ;;
        esac
    done

    if [ -z "$config" ]; then
        echo "$scriptname: no config file given" >&2
        anyerr+=1
    else
        config=$(absolutize "$config")
        if [ ! -r "$config" ]; then
            echo "$scriptname: cannot read config file '$config'" >&2
            anyerr+=1
        fi
    fi
    if [ -z "$geqdsk" ]; then
        echo "$scriptname: no GEQDSK file given" >&2
        anyerr+=1
    else
        geqdsk=$(absolutize "$geqdsk")
        if [ ! -r "$geqdsk" ]; then
            echo "$scriptname: cannot read GEQDSK file '$geqdsk'" >&2
            anyerr+=1
        fi
    fi
    if [ -z "$type" ]; then
        echo "$scriptname: no type given" >&2
        anyerr+=1
    else
        case "$type" in
            asdex|kilca)
                convexwall=$datadir/convexwall_$type.dat
                ;;
            *)
                echo "$scriptname: unrecognized type '$type'" >&2
                anyerr+=1
                ;;
        esac
    fi
    if [ $anyerr -gt 0 ]; then
        return
    fi

    for workdir; do
        if [ -d "$workdir" ]; then
            echo "$workdir already exists, contained files may be overwritten."
        else
            echo "Creating $workdir."
            mkdir -p "$workdir"
        fi

        echo "Copying files to $workdir."
        cp -t "$workdir" \
           "$datadir/field_divB0.inp" \
           "$datadir/preload_for_SYNCH.inp" \
           "$convexwall" \
           "$geqdsk"
        if [ "$type" != "kilca" ]; then
            ln -s "$datadir/AUG_B_coils.h5" "$workdir/AUG_B_coils.h5"
        fi
        cp "$config" "$workdir/mephit.in"
        replace_first_in_line "$workdir/field_divB0.inp" 7 "'${geqdsk##*/}'"     # gfile
        replace_first_in_line "$workdir/field_divB0.inp" 9 "'${convexwall##*/}'" # convex
    done
}

mephit_run() {
    meshing=0
    preconditioner=0
    iterations=0
    debug=0
    TEMP=$(getopt -o 'mpi' --long 'meshing,preconditioner,iterations,debug,memcheck,test' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-m'|'--meshing')
                meshing=1
                shift
                continue
                ;;
            '-p'|'--preconditioner')
                preconditioner=1
                shift
                continue
                ;;
            '-i'|'--iterations')
                iterations=1
                shift
                continue
                ;;
            '--debug')
                debug=1
                shift
                continue
                ;;
            '--memcheck')
                debug=2
                shift
                continue
                ;;
            '--test')
                debug=3
                shift
                continue
                ;;
            '--')
                shift
                break
                ;;
            *)
                echo "$scriptname: unrecognized option '$1'" >&2
                exit 1
                ;;
        esac
    done
    if [ $meshing -eq 0 ] && [ $preconditioner -eq 0 ] && [ $iterations -eq 0 ]; then
        # if no options are set, default to go through all phases
        meshing=1
        preconditioner=1
        iterations=1
    fi
    runmode=$(( iterations << 2 | preconditioner << 1 | meshing << 0 ))

    # default to current directory if none is given on the command line
    arglist=( "$(pwd)" )
    if [ $# -gt 0 ]; then
        arglist=( "$@" )
    fi
    for arg in "${arglist[@]}"; do
        if [ -d "$arg" ]; then
            workdir="$arg"
            pushd "$workdir"
            configlist=( mephit*.in )
        elif [ -e "$arg" ]; then
            workdir="${arg%/*}"
            configlist=( "${arg##*/}" )
            pushd "$workdir"
        else
            echo "$scriptname: skipping nonexistent directory/file '$arg'."
            anyerr+=1
            continue
        fi
        if [ -f "field_divB0_unprocessed.inp" ]; then
            # backwards compatibility for directories
            # initialized with previous version of this script
            mv -b -f field_divB0_unprocessed.inp field_divB0.inp
        fi
        export GFORTRAN_ERROR_BACKTRACE=1
        for config in "${configlist[@]}"; do
            log="${config%.*}.log"
            rm -f "$log"
            suffix="${config#mephit}"
            suffix="${suffix%.*}"
            lasterr=0
            case "$debug" in
                '3')
                    "$bindir/mephit_test.x" "$config" "$suffix" 2>&1 | tee -a "$log"
                    ;;
                '2')
                    valgrind -v \
                             --leak-check=full \
                             --show-leak-kinds=all \
                             --track-origins=yes \
                             --log-file="valgrind_%p_%n.log" \
                             "$bindir/mephit_run.x" \
                             $runmode \
                             "$config" \
                             "$suffix" \
                             "$tmpdir" \
                             "$scriptdir/ff-mephit.bash" \
                             2>&1 | tee -a "$log"
                    lasterr=$?
                    ;;
                '1')
                    gdb -x "$scriptdir/mephit.gdb" --args \
                        "$bindir/mephit_run.x" \
                        $runmode \
                        "$config" \
                        "$suffix" \
                        "$tmpdir" \
                        "$scriptdir/ff-mephit.bash"
                    lasterr=$?
                    ;;
                *)
                    "$bindir/mephit_run.x" \
                        $runmode \
                        "$config" \
                        "$suffix" \
                        "$tmpdir" \
                        "$scriptdir/ff-mephit.bash" \
                        2>&1 | tee -a "$log"
                    lasterr=$?
                    ;;
            esac
            if [ "$lasterr" -ne 0 ]; then
                echo "$scriptname: MEPHIT exited with code $lasterr while running '$arg'." | tee -a "$log" >&2
                anyerr+=1
                continue
            fi
        done
        popd
    done
}

mephit_plot() {
    config=mephit.in
    data=mephit.h5
    log=mephit.log

    # default to current directory if none is given on the command line
    workdirs=( "$(pwd)" )
    if [ $# -gt 0 ]; then
        workdirs=( "$@" )
    fi
    for workdir in "${workdirs[@]}"; do
        if [ ! -d "$workdir" ]; then
            echo "$scriptname: skipping nonexistent directory '$workdir'."
            anyerr+=1
            continue
        fi
        pushd "$workdir"
        python3 "$scriptdir/magdifplot.py" "$(pwd)" "$data"
        popd
    done
}

mephit_clean() {
    # default to current directory if none is given on the command line
    workdirs=( "$(pwd)" )
    if [ $# -gt 0 ]; then
        workdirs=( "$@" )
    fi
    for workdir in "${workdirs[@]}"; do
        if [ ! -d "$workdir" ]; then
            echo "$scriptname: skipping nonexistent directory '$workdir'."
            anyerr+=1
            continue
        fi
        pushd "$workdir"
        # files from mephit_run
        rm -f fort.* mephit.h5 mephit.log core_plasma.msh outer.msh maxwell.msh
        # files from mephit_plot
        rm -f plot*.pdf convergence.pdf Bmn*.pdf currmn*.pdf
        popd
    done
}

mephit_help() {
    echo For now, please look at the README.md for how to run this script...
}

# context
scriptdir=$(dirname "$0")
scriptdir=$(realpath "$scriptdir")
bindir=$(realpath -m "$scriptdir/../bin")
datadir=$(realpath -m "$scriptdir/../data")
if [ -d "@tmpdir@" ]; then
    tmpdir=@tmpdir@
else
    tmpdir=/tmp
fi

set -m
set -o pipefail
scriptname=$0
anyerr=0
case "$1" in
    init|convert|run|plot|clean)
        mode=$1
        shift
        mephit_$mode "$@"
        ;;
    'help'|'--help'|'-h'|'-?')
        shift
        mephit_help "$@"
        ;;
    *)
        echo "$scriptname: unrecognized mode '$1'"
        mephit_help "$@"
        anyerr+=1
        ;;
esac
exit $anyerr
