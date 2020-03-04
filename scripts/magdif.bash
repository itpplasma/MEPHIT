#!/bin/bash

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

extract_psi() {
    sed -e 's/ \+/\t/g; s/^\t//g' $1 | cut -f 9-10 > $2
}

extract_theta() {
    sed -e 's/ \+/\t/g; s/^\t//g' $1 | cut -f 11-12 > $2
}


magdif_init() {
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
           "$datadir/flt.inp" \
           "$datadir/17151.wall_full" \
           "$datadir/kin2d.inp" \
           $(absolutize "$config") \
           $(absolutize "$geqdsk") \
           $(absolutize "$convexwall")
        if [ -n $vacfield ]; then
            ln -s $(absolutize "$vacfield") "$workdir/${vacfield##*/}"
        fi

        replace_first_in_line "$workdir/field_divB0.inp" 7 "'${geqdsk##*/}'"     # gfile
        replace_first_in_line "$workdir/field_divB0.inp" 8 "'${vacfield##*/}'"   # pfile
        replace_first_in_line "$workdir/field_divB0.inp" 9 "'${convexwall##*/}'" # convex
        cp "$workdir/field_divB0.inp" "$workdir/field_divB0_unprocessed.inp"
    done
}

magdif_prepare() {
    TEMP=$(getopt -o 'c:' --long 'config:' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-c'|'--config')
                config=$2
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

    # implement mesh file options when geomint_mesh.f90 is complete
    mesh=inputformaxwell.msh
    extended=inputformaxwell_ext.msh

    for workdir; do
        pushd "$workdir"
        cp field_divB0_unprocessed.inp field_divB0.inp
        unprocessed=$(read_first_in_line field_divB0_unprocessed.inp 7)
        replace_first_in_line field_divB0.inp 7 "'\1_processed'"  # gfile
        replace_first_in_line field_divB0.inp 1 0  # ipert
        replace_first_in_line field_divB0.inp 2 1  # iequil
        kilca_scale_factor=$(nml_read_integer "$config" kilca_scale_factor)
        if [ -z $kilca_scale_factor ]; then
            kilca_scale_factor=0
        fi
        if [ $kilca_scale_factor -eq 0 ]; then
            "$bindir/standardise_equilibrium.x" "$unprocessed" "${unprocessed}_processed"
            "$bindir/axis.x" && \
                "$bindir/tri_mesh.x" && \
                "$bindir/readcarre_m.x" && \
                "$bindir/writeout.x" && \
                python3 "$scriptdir/extmesh.py" ${mesh%.*}
            lasterr=$?
            if [ $lasterr -ne 0 ]; then
                echo "$scriptname: error $lasterr during mesh generation in $workdir" >&2
                popd
                anyerr=$lasterr
                continue
            fi
            replace_first_in_line field_divB0.inp 1 1  # ipert
            # replace_first_in_line field_divB0.inp 2 0  # iequil
            "$bindir/vacfield.x"
            lasterr=$?
            replace_first_in_line field_divB0.inp 1 0  # ipert
            # replace_first_in_line field_divB0.inp 2 1  # iequil
            if [ $lasterr -ne 0 ]; then
                echo "$scriptname: error $lasterr during mesh generation in $workdir" >&2
                popd
                anyerr=$lasterr
                continue
            fi
        else
            "$bindir/magdif_mesher.x" "$config" "$unprocessed" && \
                FreeFem++ "$scriptdir/extmesh.edp"
            lasterr=$?
            if [ $lasterr -ne 0 ]; then
                echo "$scriptname: error $lasterr during mesh generation in $workdir" >&2
                popd
                anyerr=$lasterr
                continue
            fi
        fi
        popd
    done
}

magdif_run() {
    TEMP=$(getopt -o 'c:' --long 'config:' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-c'|'--config')
                config=$2
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

    tempfile=temp.dat
    for workdir; do
        pushd "$workdir"
        rm -f "$config.log" convergence.dat
        "$bindir/magdif_test.x" "$config" "$scriptdir"
        lasterr=$?
        if [ $lasterr -ne 0 ]; then
            echo "$scriptname: error $lasterr during run in $workdir" >&2
            popd
            anyerr=$lasterr
            continue
        fi
        kilca_scale_factor=$(nml_read_integer "$config" kilca_scale_factor)
        if [ -z $kilca_scale_factor ]; then
            kilca_scale_factor=0
        fi
        if [ $kilca_scale_factor -eq 0 ]; then
            infile=Bn.dat
            infile=plot_$infile
            outfile=Bmn_psi.dat
            extract_psi "$infile" $tempfile
            "$bindir/result_spectrum.x" $tempfile > $outfile
            rm $tempfile
            infile=Bn_vac.dat
            infile=plot_$infile
            outfile=Bmn_vac_psi.dat
            extract_psi "$infile" $tempfile
            "$bindir/result_spectrum.x" $tempfile > $outfile
            rm $tempfile
            infile=currn.dat
            infile=plot_${infile%.*}_000.${infile##*.}
            outfile=currmn_000_theta.dat
            extract_theta "$infile" $tempfile
            "$bindir/result_spectrum.x" $tempfile > $outfile
            rm $tempfile
        fi
        popd
    done
}

magdif_plot() {
    TEMP=$(getopt -o 'c:' --long 'config:' -n "$scriptname" -- "$@")
    eval set -- "$TEMP"
    unset TEMP
    while true; do
        case "$1" in
            '-c'|'--config')
                config=$2
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

    # implement mesh file options when geomint_mesh.f90 is complete
    mesh=inputformaxwell.msh
    extended=inputformaxwell_ext.msh

    for workdir; do
        pushd "$workdir"
        python3 "$scriptdir/magdifplot.py" . "$config" "$mesh"
        popd
    done
}

magdif_clean() {
    for workdir; do
        pushd "$workdir"
        # files from magdif_prepare
        rm -f fort.* bmod_psipol_vac.dat Bn_flux.dat BR_im.dat BR_re.dat BZ_im.dat BZ_re.dat hpsi_vac.dat p0.dat inputformaxwell_ext.msh mesh.png bmod_psipol.dat inputformaxwell.msh boundary_sqz.dat mesh_all_sqz.dat points_asc_sqz.dat mesh_all.dat points_asc.dat source.dat boundary.dat thermostat.dat triangles.dat wall_loads.dat points.dat wsrr_rec_fix.dat neighbour.fmt boundary_qq.fmt 17151.aug_qq.ele triangles.fmt points_m.fmt index_sol.pnt points.fmt index_pfr.pnt index_core.pnt RZpsirho.dat Opoint.dat
        # files from magdif_run
        rm -f fort.* plot*.dat fluxvar.dat j0phi.dat j0_gs.dat j0_amp.dat cmp_prof.dat Bn*.dat eigvec*.dat presn*.dat currn*.dat freefem.out convergence.dat rect*.dat *.log Bmn*.dat currmn*.dat
        # files from magdif_plot
        rm -f plot*.pdf convergence.pdf
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
