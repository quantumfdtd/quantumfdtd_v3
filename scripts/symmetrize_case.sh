#!/bin/bash

############################################################################
suffix="p m pk all"
rsuffix="p m all"
osuffix="pk"
############################################################################
CAS=$1
[[ -z "$CAS" ]] && CAS="wavefunction"
############################################################################
set -e
############################################################################
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" 
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
############################################################################
echoerr() { echo "$@" 1>&2; exit 1;}
############################################################################

for i in 0 1 2
do
	[[ -f "${CAS}_${i}_all_norm.dat.gz" ]] && continue # Do not rerun for computed cases!!

	echo "Procesing case ${CAS}_${i}_%d.dat.gz"...

	if [ -f "${CAS}_${i}_all.dat.gz" ];
	then
		for j in ${suffix}
		do
			obj="${CAS}_${i}_${j}.dat.gz"
			if [ ! -f "${obj}" ];
			then
				echoerr "Error file ${obj} does not exist!!!" 
			else
				gunzip "${obj}"
			fi
		done
	else
		[[ -f "${CAS}_${i}.tgz" ]] && tar -zxf ${CAS}_${i}.tgz
		"${DIR}/symmetrize_wf.py" "${CAS}_${i}_%d.dat"

		for j in $suffix
		do
			obj="${CAS}_${i}_${j}.dat"
			[[ -f "${obj}" ]] || echoerr "Error file ${obj} has not been created!!!"
		done
	fi
		
	for j in $rsuffix
	do
		in="${CAS}_${i}_${j}.dat"
		out="${CAS}_${i}_${j}_norm.dat"
		echo "Normalizing wavefunction ${in}: "
		if [[ ! -s "${in}" ]]; then
		       echo "Null wavefunction ${in}!!"
		       continue
	        fi
          
		"${DIR}/normalize.py" "${in}" "${out}"

		[[ -f "${out}" ]] || echoerr "Error, normalized wavefunction ${out} has not been created!!!"
		[[ -s "${out}" ]] || echoerr "Error, normalized wavefunction ${out} is empty!!!"
		echo "Normalized wavefunction ${out}!!"
	done

	for j in $osuffix
	do
		in="${CAS}_${i}_${j}.dat"
		out="${CAS}_${i}_${j}_norm.dat"
		echo "Normalizing wavefunction ${in}"
		if [[ ! -s "${in}" ]]; then
		       echo "Null wavefunction ${in}!!"
		       continue
	        fi
          
		"${DIR}/normalize.py" 5 "${in}" "${out}"

		[[ -f "${out}" ]] || echoerr "Error, normalized wavefunction ${out} has not been created!!!"
		[[ -s "${out}" ]] || echoerr "Error, normalized wavefunction ${out} is empty!!!"
		echo "Normalized wavefunction ${out}!!"
	done

	for j in $suffix
	do
		[[ -s "${CAS}_${i}_${j}.dat" ]]      && gzip "${CAS}_${i}_${j}.dat"
		[[ -s "${CAS}_${i}_${j}_norm.dat" ]] && gzip "${CAS}_${i}_${j}_norm.dat"
	done
	
	[[ -f "${CAS}_${i}.tgz" ]] || tar -zcf ${CAS}_${i}.tgz --remove-files ${CAS}_${i}_*.dat
done
