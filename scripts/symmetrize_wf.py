#!/usr/bin/env python3

import sys
import os
import numpy as np
import pandas as pd
import re
import quantumfdtd

# ###################################################### #
# IF \int |W|^2<WRITELIM, THE WF IS NOT WRITTEN TO HDD
WRITELIM = 1.e-5
# ###################################################### #


def compute_r_fft(NUM, x, y, z):
    ix, iy, iz = [((NUM + 1 - i) % NUM) + 1 for i in [x, y, z]]
    return (iz - 1) + NUM * ((iy - 1) + NUM * (ix - 1))


def compute_rp_rm_fft(NUM, x, y, z, i):
    ix, iy, iz = [((NUM + 1 - i) % NUM) + 1 for i in [x, y, z]]

    if i == 0:
        rp = (iz - 1) + NUM * ((iy - 1) + NUM * (x - 1))
        rm = (z - 1) + NUM * ((y - 1) + NUM * (ix - 1))
        pk = x
    elif i == 1:
        rp = (iz - 1) + NUM * ((y - 1) + NUM * (ix - 1))
        rm = (z - 1) + NUM * ((iy - 1) + NUM * (x - 1))
        pk = y
    elif i == 2:
        rp = (z - 1) + NUM * ((iy - 1) + NUM * (ix - 1))
        rm = (iz - 1) + NUM * ((y - 1) + NUM * (x - 1))
        pk = z
    else:
        sys.stderr.write(f'INTERNAL ERROR!!!\nBAD VALUE OF i IN compute_rp_rm_fft\n')
        sys.exit(1)

    pk = np.where(pk <= NUM / 2, pk - 1, pk - NUM - 1)

    return [rp, rm, pk]


def compute_rp_rm_CN2(NUM, x, y, z, i):
    ix, iy, iz = [NUM + 1 - i for i in [x, y, z]]

    if i == 0:
        rp = (iz - 1) + NUM * ((iy - 1) + NUM * (x - 1))
        rm = (z - 1) + NUM * ((y - 1) + NUM * (ix - 1))
        pk = x
    elif i == 1:
        rp = (iz - 1) + NUM * ((y - 1) + NUM * (ix - 1))
        rm = (z - 1) + NUM * ((iy - 1) + NUM * (x - 1))
        pk = y
    elif i == 2:
        rp = (z - 1) + NUM * ((iy - 1) + NUM * (ix - 1))
        rm = (iz - 1) + NUM * ((y - 1) + NUM * (x - 1))
        pk = z
    else:
        sys.stderr.write(f'INTERNAL ERROR!!!\nBAD VALUE OF i IN compute_rp_rm_fft\n')
        sys.exit(1)

    pk = pk - (NUM + 1) / 2.

    return [rp, rm, pk]


def process_wf(centered_on_lattice, wf, NUM, out_p, out_m, out_pk, out_all, ii=0):
    int_p = 0.
    int_m = 0.
    int_pk = 0.
    int_all = 0.

    disp = 1 if centered_on_lattice else 2

    # Add index column and compute needed index values for the reordering
    wf['rowIndex'] = wf.index
    wf['rindex'] = NUM ** 3 - wf['rowIndex'] - 1 if centered_on_lattice else compute_r_fft(NUM, wf['x'], wf['y'], wf['z'])
    wf['rpindex'], wf['rmindex'], wf['pk'] = compute_rp_rm_CN2(NUM, wf['x'], wf['y'], wf['z'], ii) if centered_on_lattice else compute_rp_rm_fft(NUM, wf['x'], wf['y'], wf['z'], ii)

    # Get the values from the reordered wf array
    wf['iwre'] = wf.iloc[wf['rindex']][['re']].reset_index(drop=True)
    wf['iwim'] = wf.iloc[wf['rindex']][['im']].reset_index(drop=True)
    wf['ppwre'] = wf.iloc[wf['rpindex']][['re']].reset_index(drop=True)
    wf['ppwim'] = wf.iloc[wf['rpindex']][['im']].reset_index(drop=True)
    wf['pmwre'] = wf.iloc[wf['rmindex']][['re']].reset_index(drop=True)
    wf['pmwim'] = wf.iloc[wf['rmindex']][['im']].reset_index(drop=True)

    # At this point we don't need the columns rindex, rpindex, and rmindex any more
    wf.drop(columns=['rindex', 'rpindex', 'rmindex'], inplace=True)

    # Compute (intermediate) results
    wf['pwre'] = 0.5 * (wf['re'] + wf['iwre'])
    wf['mwre'] = 0.5 * (wf['re'] - wf['iwre'])
    wf['pwim'] = 0.5 * (wf['im'] + wf['iwim'])
    wf['mwim'] = 0.5 * (wf['im'] - wf['iwim'])
    wf['mp_wre'] = 0.25 * ((wf['re'] + wf['ppwre']) - (wf['pmwre'] + wf['iwre']))
    wf['mp_wim'] = 0.25 * ((wf['im'] + wf['ppwim']) - (wf['pmwim'] + wf['iwim']))
    wf['int_p'] = wf['pwre'] * wf['pwre'] + wf['pwim'] * wf['pwim']
    wf['int_m'] = wf['mwre'] * wf['mwre'] + wf['mwim'] * wf['mwim']
    wf['int_pk'] = wf['mp_wre'] * wf['mp_wre'] + wf['mp_wim'] * wf['mp_wim']
    wf['int_all'] = wf['re'] * wf['re'] + wf['im'] * wf['im']

    # Compute the return values of the function
    int_p = sum(wf['int_p'])
    int_m = sum(wf['int_m'])
    int_pk = sum(wf['int_pk'])
    int_all = sum(wf['int_all'])

    write_p = int_p > WRITELIM
    write_m = int_m > WRITELIM
    write_pk = int_pk > WRITELIM
    write_all = int_all > WRITELIM

    # At this point we don't need the columns rowIndex, iwre, iwim, ppwre, ppwim, pmwre, pmwim,
    #                                         int_p, int_m, int_pk, and int_all any more
    wf.drop(columns=['rowIndex', 'iwre', 'iwim', 'ppwre', 'ppwim', 'pmwre', 'pmwim',
                     'int_p', 'int_m', 'int_pk', 'int_all'], inplace=True)

    # Write results to disk
    for row in wf.itertuples(index=False):

        x, y, z, r, re, im, pk, pwre, mwre, pwim, mwim, mp_wre, mp_wim = row

        if write_p:
            out_p.write(f'{x}\t {y}\t {z}\t {r}\t {pwre}\t {pwim}\n')
        if write_m:
            out_m.write(f'{x}\t {y}\t {z}\t {r}\t {mwre}\t {mwim}\n')
        if write_pk:
            out_pk.write(f'{x}\t {y}\t {z}\t {r}\t {pk}\t {mp_wre}\t {mp_wim}\n')
        if write_all:
            out_all.write(f'{x}\t {y}\t {z}\t {r}\t {re}\t {im}\n')

    return (int_p, int_m, int_pk, int_all)


if __name__ == "__main__":
    NUM = -1
    centered_on_lattice = None

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        sys.stderr.write("ERROR!!!\nCORRECT USAGE: ../scripts/symmetrize_wf.py [params.txt] wavefunction_i_%d.dat\n\n")
        sys.stderr.write("%d points to the split index. I.e., for wf0 (ground state), it is usually wavefunction_0_%d.dat\n")
        sys.stderr.write("If not specified, the config file is stored in ../input/params.txt\n\n")
        sys.exit(1)

    p_config = sys.argv[1] if len(sys.argv) == 3 else "../input/params.txt"
    wf = sys.argv[-1]

    if not os.path.exists(p_config):
        sys.stderr.write(f'Error: config file {p_config} does not exist!!\n')
        sys.exit(1)

    conf = quantumfdtd.quantumfdtd('..', p_config)
    NUM = float(conf.conf['NUM'])
    centered_on_lattice = conf.conf['center_on_lattice']
    ii = int(conf.conf.get('INITCONDAXIS', 0))
    del conf

    if centered_on_lattice:
        print("\nOrigin of symmetry: (N/2,N/2,N/2)\n")
    else:
        print("\nOrigin of symmetry: (0,0,0)\n")

    folder = os.path.dirname(wf)
    if folder == '':
        folder = '.'

    regexpr = re.sub('%d', '[0-9]*', os.path.basename(wf))
    regexpr = '^' + regexpr + '$'
    wfs = [folder + '/' + f for f in os.listdir(folder) if re.search(regexpr, f)]
    NUMX = len(wfs)

    if NUMX == 0:
        sys.stderr.write('Error: no valid wavefunction found in this folder!!! Check the provided name...\n')
        sys.exit(1)

    out_p = open(re.sub('%d', 'p', wf), 'w')
    out_m = open(re.sub('%d', 'm', wf), 'w')
    out_pk = open(re.sub('%d', 'pk', wf), 'w')
    out_all = open(re.sub('%d', 'all', wf), 'w')

    headers = ["x", "y", "z", "r", "re", "im"]
    typ = {"x": np.int32, "y": np.int32, "z": np.int32, "r": np.float64, "re": np.float64, "im": np.float64}

    wfs = pd.concat([pd.read_csv(fil_wf, dtype=typ, delim_whitespace=True, header=None, names=headers)
                     for fil_wf in [re.sub('%d', str(iwf + 1), wf) for iwf in range(NUMX)]],
                    ignore_index=True)

    int_p, int_m, int_pk, int_all = process_wf(centered_on_lattice, wfs, NUM, out_p, out_m, out_pk, out_all, ii)

    del wfs

    out_p.close()
    out_m.close()
    out_pk.close()
    out_all.close()

    print(f'Finished evaluation of {wf}...')
    print(f'> {wf} : <wf|wf>, <wf|P+|wf>, <wf|P-|wf>: <wf|Pi|wf>')
    print(f'> {wf} : {int_all}, {int_p}, {int_m}, {int_pk}')

    out_proj = open(re.sub('%d', 'proj', wf), 'w')
    out_proj.write(f'all\t Pp\t Pm\t Ppk\n')
    out_proj.write(f'{int_all}\t {int_p}\t {int_m}\t {int_pk}\n')
    out_proj.close()
