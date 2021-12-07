#!/bin/bash
############################################################################################################
set -e
############################################################################
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
############################################################################################################
############################################################################################################
if [ -z "$1" ]; then
    echo ">> PROCESSING POTENTIAL..."

    cd data
    mkfifo pot.plot

    if [ -f potential.tgz ]; then
        tar -O -zxf potential.tgz | awk '{print $4, $5}' > pot.plot &
    else
        cat potential_*.dat | awk '{print $4, $5}' > pot.plot &
    fi

    gnuplot <<-EOF
	set term png
	set output 'pot.png'
	unset key
	set title 'Potential'
	set xlabel 'r/A'
	set ylabel 'V (GeV)'
	plot 'pot.plot' w d
	EOF
    rm pot.plot

    [[ -f potential_1.dat ]] && tar -zcf potential.tgz --remove-files potential_*.dat

    ############################################################################################################
    echo ">> PROCESSING DECAY TABLE..."
    gnuplot <<-EOF
	set term png
	set output 'decay.png'
	unset key
	set title 'Evolution of the energies'
	set xlabel 't (GeV^{-1})'
	set ylabel 'E (GeV)'
	plot 'decay.dat' u 2:5 w lp, 'ground_state.out' u 2:5, 'first_excited_state.out' u 2:5, 'second_excited_state.out' u 2:5
	EOF
    ############################################################################################################
    "${DIR}/postprocess.sh" wavefunction
    ############################################################################################################
    CASES=$(find snapshot -name "wavefunction_*_0_1.dat" -o -name "wavefunction_*_0.tgz" | sort -V)
    set +e
    for DAT in $CASES
    do
        DAT="${DAT%_0_1.dat}"
        echo "${DAT%_0.tgz}"
    done | xargs -P8 -n1 "${DIR}/postprocess.sh"
############################################################################################################
else
    NAME=$1
    for i in 0 1 2
    do
        TMP1="${NAME}_wf_${i}.plot"

        if [ -f "${NAME}_${i}_1.dat" ]; then
            cat "${NAME}_${i}"_*.dat | awk '{print $4, $5, $5*$5+$6*$6}' > "${TMP1}"
        elif [ -f "${NAME}_${i}_1.dat.gz" ]; then
            zcat "${NAME}_${i}"_*.dat | awk '{print $4, $5, $5*$5+$6*$6}' > "${TMP1}"
        elif [ -f "${NAME}_${i}.tgz" ]; then
            tar -O -zxf "${NAME}_${i}.tgz" | awk '{print $4, $5, $5*$5+$6*$6}' > "${TMP1}"
        else
            continue
        fi

        echo ">> PROCESSING CASE ${NAME}, WAVEFUNCTION n=${i}..."

        set +e
        gnuplot <<-EOF
		set term png
		unset key
		set xlabel 'r/A'
		set title '|wf_${i}|^2'
		set ylabel '|wf_${i}|^2'
		set output '${NAME}_wf2_${i}.png'
		plot '${TMP1}' u 1:3 w d
		set ylabel 'Re (wf_${i})'
		set title 'Re (wf_${i})'
		set output '${NAME}_re_wf_${i}.png'
		plot '${TMP1}' u 1:2 w d
		set logscale y
		set title '|wf_${i}|^2'
		set ylabel '|wf_${i}|^2'
		set output '${NAME}_wf2_${i}_log.png'
		plot '${TMP1}' u 1:3 w d
		set ylabel 'Re (wf_${i})'
		set title 'Re (wf_${i})'
		set output '${NAME}_re_wf_${i}_log.png'
		plot '${TMP1}' u 1:(abs(\$2)) w d
		EOF

        set -e

        rm "${TMP1}"
        [[ -f "${NAME}_${i}.tgz" ]] || tar -zcf "${NAME}_${i}.tgz" --remove-files "${NAME}_${i}"_*.dat
    done
fi
############################################################################################################
