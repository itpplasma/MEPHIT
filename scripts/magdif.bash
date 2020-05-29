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

magdif_prepare() {
    config=magdif.inp
    log=magdif.log
    # maybe implement mesh file options when magdif_mesh_mod.f90 is complete
    mesh=inputformaxwell.msh
    extended=inputformaxwell_ext.msh

    for workdir; do
        pushd "$workdir"
        cp field_divB0_unprocessed.inp field_divB0.inp
        unprocessed=$(read_first_in_line field_divB0_unprocessed.inp 7)
        replace_first_in_line field_divB0.inp 7 "'\1_processed'"  # gfile
        replace_first_in_line field_divB0.inp 1 0  # ipert
        replace_first_in_line field_divB0.inp 2 1  # iequil
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
        kilca_scale_factor=$(nml_read_integer "$config" kilca_scale_factor)
        if [ -z "$kilca_scale_factor" ]; then
            kilca_scale_factor=0
        fi
        if [ "$kilca_scale_factor" -eq 0 ]; then
            replace_first_in_line field_divB0.inp 1 1  # ipert
            # replace_first_in_line field_divB0.inp 2 0  # iequil
            ## uncomment to use at most half of available RAM and go to swap instead
            ## MemFree=$(sed -ne '/MemFree/ s/MemFree: *\([0-9]*\) kB/\1/gp' /proc/meminfo)
            ## systemd-run --scope --user -p MemoryHigh=$((MemFree / 2048))M \
            "$bindir/vacfield.x" "$config" 2>&1 | tee -a "$log"
            lasterr=$?
            replace_first_in_line field_divB0.inp 1 0  # ipert
            # replace_first_in_line field_divB0.inp 2 1  # iequil
            if [ $lasterr -ne 0 ]; then
                echo "$scriptname: error $lasterr during mesh generation in $workdir" | tee -a "$log" >&2
                popd
                anyerr=$lasterr
                continue
            fi
        fi
        popd
    done
}

magdif_run() {
    config=magdif.inp
    log=magdif.log

    for workdir; do
        pushd "$workdir"
        rm -f magdif.log convergence.dat
        "$bindir/magdif_test.x" "$config" "$scriptdir" 2>&1 | tee -a "$log"
        lasterr=$?
        if [ $lasterr -ne 0 ]; then
            echo "$scriptname: error $lasterr during run in $workdir" | tee -a "$log" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        popd
    done
}

magdif_plot() {
    config=magdif.inp
    log=magdif.log
    # implement mesh file options when geomint_mesh.f90 is complete
    mesh=inputformaxwell.msh
    extended=inputformaxwell_ext.msh

    for workdir; do
        pushd "$workdir"
        python3 "$scriptdir/magdifplot.py" . "$config" "$mesh" 2>&1 | tee -a "$log"
        popd
    done
}

magdif_clean() {
    for workdir; do
        pushd "$workdir"
        # files from magdif_prepare
        rm -f fort.* Bn_flux.dat inputformaxwell_ext.msh inputformaxwell.msh points.fmt mesh.dat mesh_new.asc box_size_axis.dat btor_rbig.dat flux_functions.dat twodim_functions.dat
        # files from magdif_run
        rm -f fort.* plot*.dat fluxvar.dat j0phi.dat j0_gs.dat j0_amp.dat cmp_prof.dat Bn*.dat eigvec*.dat presn*.dat currn*.dat freefem.out convergence.dat rect*.dat *.log Bmn*.dat currmn*.dat
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

set -o pipefail
scriptname=$0
anyerr=0
case "$1" in
    init|prepare|run|plot|clean)
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
