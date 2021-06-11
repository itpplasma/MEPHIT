#!/bin/bash

if [[ "$OSTYPE" == "darwin"* ]]; then
    PATH=/opt/local/libexec/gnubin:$PATH
fi

# if first argument is absolute path, return it unchanged
# if first argument is relative path, prefix it with second argument and return absolute path
# if second argument is empty, use current path
absolutize() {
    case $1 in
        /*)
            echo $1
            ;;
        *)
            if [ -z "$2" ]; then
                echo $(realpath -m $1)
            else
                echo $(realpath -m $2/$1)
            fi
            ;;
    esac
}

replace_first_in_line() {
    sed -Ee "$2 s/^ *'?([^' ]*)'?/$3/" -i $1
}

read_first_in_line() {
    sed -Ene "$2 s/^ *'?([^' ]*)'? +.*/\1/ p" $1
}

nml_read_integer() {
    sed -Ene "s/ *$2 *= *([-0-9]+).*/\1/i p" $1
}


magdif_init() {
    config=
    geqdsk=
    convexwall=
    vacfield=
    TEMP=$(getopt -o 'c:g:v:w:' --long 'config:,g-eqdsk:,vacuum-field:,convex-wall:' -n "$scriptname" -- "$@")
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
            '-w'|'--convex-wall')
                convexwall=$2
                shift 2
                continue
                ;;
            '-v'|'--vacuum-field')
                vacfield=$2
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
           $(absolutize "$config") \
           $(absolutize "$geqdsk") \
           $(absolutize "$convexwall")
        if [ -n "$vacfield" ]; then
            ln -s $(absolutize "$vacfield") "$workdir/${vacfield##*/}"
        fi

        mv "$workdir/${config##*/}" "$workdir/magdif.inp"
        replace_first_in_line "$workdir/field_divB0.inp" 7 "'${geqdsk##*/}'"     # gfile
        replace_first_in_line "$workdir/field_divB0.inp" 8 "'${vacfield##*/}'"   # pfile
        replace_first_in_line "$workdir/field_divB0.inp" 9 "'${convexwall##*/}'" # convex
        cp "$workdir/field_divB0.inp" "$workdir/field_divB0_unprocessed.inp"
    done
}

magdif_convert() {
    in_type=$1
    out_type=$2
    in_dir=$3
    out_dir=$4
    if [ -z "$in_dir" ]; then
        echo "$scriptname: expected directory at third position after convert"
        anyerr=1
        return
    fi
    if [ ! -d "$in_dir" ]; then
        echo "$scriptname: directory $in_dir does not exist"
        anyerr=1
        return
    fi
    in_dir=$(absolutize "$in_dir")
    if [ -n "$out_dir" ]; then
        if [ ! -d "$out_dir" ]; then
            mkdir -p "$out_dir"
        fi
        out_dir=$(absolutize "$out_dir")
    else
        out_dir=$in_dir
    fi
    "$bindir/vacfield.x" "$in_type" "$out_type" "$in_dir" "$out_dir"
    lasterr=$?
    if [ $lasterr -ne 0 ]; then
        echo "$scriptname: error $lasterr during coil geometry / vacuum field conversion"
        anyerr=$lasterr
    fi
}

magdif_prepare() {
    config=magdif.inp
    log=magdif.log
    # maybe implement mesh file options when magdif_mesh_mod.f90 is complete
    mesh=inputformaxwell.msh
    extended=inputformaxwell_ext.msh

    for workdir; do
        pushd "$workdir"
        rm -f "$log"
        cp field_divB0_unprocessed.inp field_divB0.inp
        unprocessed=$(read_first_in_line field_divB0_unprocessed.inp 7)
        replace_first_in_line field_divB0.inp 7 "'\1_processed'"  # gfile
        "$bindir/magdif_mesher.x" "$config" "$unprocessed" 2>&1 | tee "$log"
        lasterr=$?
        if [ $lasterr -ne 0 ]; then
            echo "$scriptname: error $lasterr during mesh generation in $workdir" | tee -a "$log" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        if [ -n "$SSH_CLIENT" -o -z "$DISPLAY" ]; then
            # when connected via SSH or not working in an X environmnent,
            # suppress graphics to avoid setting an error code
            FreeFem++ -nw "$scriptdir/extmesh.edp" 2>&1 | tee -a "$log"
        else
            FreeFem++ "$scriptdir/extmesh.edp" 2>&1 | tee -a "$log"
        fi
        lasterr=$?
        if [ $lasterr -ne 0 ]; then
            echo "$scriptname: error $lasterr during mesh generation in $workdir" | tee -a "$log" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        popd
    done
}

magdif_run() {
    config=magdif.inp
    log=magdif.log
    batchmode=false
    TEMP=$(getopt -o 'b' --long 'batchmode' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-b'|'--batchmode')
                batchmode=true
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

    for workdir; do
        pushd "$workdir"
        "$bindir/magdif_run.x" \
            "$config" \
            "$tmpdir" \
            "$scriptdir/maxwell_daemon.edp" \
            2>&1 | tee -a "$log"
        lasterr=$?
        if [ "$lasterr" -ne 0 ]; then
            echo "$scriptname: magdif_run.x exited with code $lasterr during run in $workdir" | tee -a "$log" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        popd
    done
}

magdif_plot() {
    config=magdif.inp
    data=magdif.h5
    log=magdif.log

    for workdir; do
        pushd "$workdir"
        python3 "$scriptdir/magdifplot.py" $(absolutize .) "$config" "$data"
        popd
    done
}

magdif_clean() {
    for workdir; do
        pushd "$workdir"
        # files from magdif_prepare
        rm -f fort.* magdif.h5 inputformaxwell_ext.msh inputformaxwell.msh box_size_axis.dat btor_rbig.dat flux_functions.dat twodim_functions.dat optpolres.dat j0_gs.dat j0_amp.dat cmp_prof.dat check_q_step.dat check_q_cont.dat cmp_vac.dat cmp_RT0.dat
        # files from magdif_run
        rm -f fort.* magdif.h5 freefem.out magdif.log Bmn*.dat currmn*.dat currn_par*.dat
        # files from magdif_plot
        rm -f plot*.pdf convergence.pdf Bmn*.pdf currmn*.pdf
        popd
    done
}

magdif_help() {
    echo For now, please look at the README.md for how to run this script...
}

# context
bindir=$(realpath $(dirname $0))
scriptdir=$(realpath -m "$bindir/../../scripts")
datadir=$(realpath -m "$bindir/../../data")
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
    init|convert|prepare|run|plot|clean)
        mode=$1
        shift
        magdif_$mode "$@"
        ;;
    'help'|'--help'|'-h'|'-?')
        shift
        magdif_help "$@"
        ;;
    *)
        echo "$scriptname: unrecognized mode '$1'"
        exit 1
        ;;
esac
exit $anyerr
