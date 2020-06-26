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
	cd data
	############################################################################################################
	"${DIR}/postprocess_sym.sh" wavefunction
	############################################################################################################
	CASES=$(find snapshot -name "wavefunction_*_0_all_norm.dat.gz" | sort -V)
	set +e
	for DAT in $CASES
	do
		echo "${DAT%_0_all_norm.dat.gz}"
	done | xargs -P8 -n1 "${DIR}/postprocess_sym.sh"
############################################################################################################
else
	NAME=$1
	echo " >> NAME: ${NAME} << "
	for i in 0 1 2
	do
		for j in "all" "m" "pk" "p"
		do

			CASE="${i}_${j}_norm"
			TMP1="${NAME}_${CASE}.dat"
			FILE=${TMP1}.gz

			echo " >> TESTING FILE ${FILE} << "

			[[ -f "${FILE}" ]] || continue
			[[ -f "${TMP1}" ]] || gunzip -k "${FILE}"
			[[ -s "${TMP1}" ]] || (rm -f "${TMP1}" && continue)

			echo ">> PROCESSING WAVEFUNCTION n=${i} OF ${NAME}..."

			if [ "$j" == "pk" ]; then
				echo " *Special pk plotting...."
				OUT_re="4:(\$6*(\$4/\$5))"
				OUT_log_re="4:(abs(\$6*(\$4/\$5)))"
				OUT_wf2="4:((\$6*\$6+\$7*\$7)*(\$4/\$5)*(\$4/\$5))"
			else
				OUT_re="4:5"
				OUT_log_re="4:(abs(\$5))"
				OUT_wf2="4:(\$5*\$5+\$6*\$6)"
			fi

			set +e

			gnuplot <<-EOF
			set term png
			unset key
			set xlabel 'r/A'
			set title '|wf_${i}|^2'
			set ylabel '|wf_${i}|^2'
			set output '${NAME}_wf2_${i}_${j}.png'
			plot '${TMP1}' u $OUT_wf2 w d
			set ylabel 'Re (wf_${i})'
			set title 'Re (wf_${i})'
			set output '${NAME}_re_wf_${i}_${j}.png'
			plot '${TMP1}' u $OUT_re w d
			set logscale y
			set title '|wf_${i}|^2'
			set ylabel '|wf_${i}|^2'
			set output '${NAME}_wf2_${i}_${j}_log.png'
			plot '${TMP1}' u $OUT_wf2 w d
			set ylabel 'Re (wf_${i})'
			set title 'Re (wf_${i})'
			set output '${NAME}_re_wf_${i}_${j}_log.png'
			plot '${TMP1}' u $OUT_log_re w d
			EOF

			[[ $? -ne 0 ]] && echo ">>> CASE ${NAME}_${CASE}.dat.gz has failed!!! <<<"

			set -e

			rm "${TMP1}"
		done
	done
fi
############################################################################################################
