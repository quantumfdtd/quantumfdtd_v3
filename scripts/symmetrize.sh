#!/bin/bash

############################################################################
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

############################################################################

if [ -z "$1" ]; then
    cd data || exit
    [[ -f potential_1.dat ]] && tar -zcf potential.tgz --remove-files potential_*.dat
    "${DIR}/symmetrize_case.sh" wavefunction

    if [ -d snapshot ]; then
        echo "Post processing snapshots..."
        for cas in $(find snapshot -name "wavefunction_[0-9]*_0_1.dat" -o -name "wavefunction_[0-9]*_1.tgz" | sort -V)
        do
            OUT="${cas%_1.tgz}"
            echo "${OUT%_0_1.dat}"
        done | xargs -P8 -n1 "${DIR}/symmetrize.sh"
    fi
else
    echo "Processing case $1"
    "${DIR}/symmetrize_case.sh" "${1}" > "${1}.out" 2> "${1}.err"
    gzip "${1}.out"
    gzip "${1}.err"
fi

exit
