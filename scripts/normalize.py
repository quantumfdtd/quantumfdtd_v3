#!/usr/bin/env python3

import sys
import os
import numpy as np
import pandas as pd
import re


def process_wf(wf, out, ncols):

    norm = np.sqrt((wf['re']**2 + wf['im']**2).sum())
    wf['re'] /= norm
    wf['im'] /= norm

    for row in wf.itertuples(index=True):
        if ncols == 4:
            index, x, y, z, r, wre, wim = row
            out.write(f'{x}\t {y}\t {z}\t {r}\t {wre}\t {wim}\n')
        elif ncols == 5:
            index, x, y, z, r, pk, wre, wim = row
            out.write(f'{x}\t {y}\t {z}\t {r}\t {pk}\t {wre}\t {wim}\n')
        else:
            sys.stderr.write("ERROR!!!Unexpected number of columns...\n")
            sys.exit(1)

    return


if __name__ == "__main__":
    ncols = 4
    if ((len(sys.argv) < 3) or (len(sys.argv) > 4)):
        sys.stderr.write("CORRECT USAGE: ../scripts/symmetrize_wf.py [PREVCOL] wavefunction.dat out.dat\n")
        sys.stderr.write("[PREVCOL]: Number of columns before Re(wf) and Im(wf). 4 by deault.\n")
        sys.exit(1)
    elif len(sys.argv) == 4:
        ncols = int(sys.argv[1])
        wf = sys.argv[2]
        out = open(sys.argv[3], 'w')
    else:
        wf = sys.argv[1]
        out = open(sys.argv[2], 'w')

    if ncols == 4:
        headers = ["x", "y", "z", "r", "re", "im"]
        typ = {"x": np.int32, "y": np.int32, "z": np.int32, "r": np.float64, "re": np.float64, "im": np.float64}
    elif ncols == 5:
        headers = ["x", "y", "z", "r", "pk", "re", "im"]
        typ = {"x": np.int32, "y": np.int32, "z": np.int32, "r": np.float64, "pk": np.float64,
               "re": np.float64, "im": np.float64}
    else:
        sys.stderr.write("ERROR!!!Unexpected number of columns...\n")
        sys.exit(1)

    wfs = pd.read_csv(wf, dtype=typ, delim_whitespace=True, header=None, names=headers)
    process_wf(wfs, out, ncols)

    del wfs
    out.close()

    print('Finished evaluation...')
