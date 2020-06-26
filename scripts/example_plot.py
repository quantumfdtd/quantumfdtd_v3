#!/usr/bin/env python3

import sys

import numpy as np
import pandas as pd

import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt

# IMPORTANT NOTE!!!
# THIS PATH SHOULD POINT TO THE FOLDER WHERE YOU
# INSTALL quantumfdtd.py
sys.path.append('./quantumfdtd_v3/scripts')
# THE quantumfdtd LIBRARY SHOULD BE LOADED AFTER ADDING ITS PATH
import quantumfdtd as qu

matplotlib.use("Agg")

rcParams.update({'figure.autolayout': True})
rcParams.update({'figure.dpi': 200})
rcParams.update({'legend.fontsize': 8, 'legend.handlelength': 2})

# #################################################################################

# INPORTANT NOTE!!!
# PLEASE, FEEL FREE TO ADDAPT THE FOLLOWING VARIABLES TO YOUR

# ACTUAL NEEDS
# ATTENTION: FOR SPECIFIC PARAMETERS, WE NEED KIN0 "BEFORE" MIKE'S IN ORDER TO KEEP
# THE SAME AXIS RANGES ON THE PAPER
CASES = ['KIN0', 'KIN0_MIKE_PARITY', 'KIN0_MIKE_PARITY_ANTISYM', 'KIN1', 'KIN2',
         'KIN3', 'KIN2_HEAVYMASS', 'KIN3_HEAVYMASS', 'KIN0_128', 'KIN2_128']
NO_THS = ['KIN3', 'KIN3_HEAVYMASS']
REDUCED_XRANGE = ['KIN2_HEAVYMASS', 'KIN3_HEAVYMASS']

STATES = [0, 1, 2]
figsize = (4., 2.8)

# #################################################################################


def wf_1s(rA, MASS, A):
    #
    r = rA * A
    wf = np.exp(-MASS * r)
    #
    norm = np.sqrt(np.abs(wf * wf).sum())
    #
    avg_rA = (np.abs(wf * wf) * rA / (norm * norm)).sum()
    #
    rA = np.linspace(min(rA), max(rA), 50)
    r = rA * A
    wf = np.exp(-MASS * r) / norm
    return (rA, wf, avg_rA)


def wf_2s(rA, MASS, A):
    #
    r = rA * A
    wf = (2. - MASS * r) * np.exp(-MASS * r / 2.)
    #
    norm = np.sqrt(np.abs(wf * wf).sum())
    #
    avg_rA = (np.abs(wf * wf) * rA / (norm * norm)).sum()
    #
    rA = np.linspace(min(rA), max(rA), 50)
    r = rA * A
    wf = (2. - MASS * r) * np.exp(-MASS * r / 2.) / norm
    return (rA, wf, avg_rA)


def wf_2p_cos(rA, pk, MASS, A):
    #
    r = rA * A
    cos = pk / rA
    wf = MASS * r * cos * np.exp(-MASS * r / 2.)
    #
    norm = np.sqrt(np.abs(wf * wf).sum())
    #
    avg_rA = (np.abs(wf * wf) * rA / (norm * norm)).sum()
    #
    rA = np.linspace(min(rA), max(rA), 50)
    r = rA * A
    wf = MASS * r * np.exp(-MASS * r / 2.) / norm
    return (rA, wf, avg_rA)


if __name__ == "__main__":

    for CASE in CASES:

        print(f'Processing case {CASE}')
        quant = qu.quantumfdtd(CASE)
        quant.load_proj_components()

        SNAPS = quant.get_snaps()

        # PLOT EVOLUTION OF EIGEN-ENERGIES
        fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=figsize)

        quant.e_wf.plot(ax=axes, x='time', y=['re_e_wf0', 're_e_wf1', 're_e_wf2'],
                        label=[f'${{\\rm E}}_0$', f'${{\\rm E}}_1$', f'${{\\rm E}}_2$'])

        # SPECIFIC CODE
        if CASE is 'KIN0':
            limits = axes.get_ylim()
            plt.axhline(y=0., color='k', linestyle='dotted')
        elif CASE in ['KIN0_MIKE_PARITY', 'KIN0_MIKE_PARITY_ANTISYM']:
            axes.set_ylim(limits)
            plt.axhline(y=0., color='k', linestyle='dotted')

        plt.xlabel(f'$\\tau$ (GeV$^{{-1}}$)')
        axes.yaxis.set_label_coords(-0.06, 1.04)
        plt.ylabel(f'$E$ (GeV)', rotation='horizontal')

        plt.grid(True)

        plt.savefig(f'plots/energ_{CASE}.png')
        plt.close('all')

        # PLOT POTENTIAL
        fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=figsize)

        potential = quant.load_potential()
        potential.plot(ax=axes, x='r', y='re_v', legend=None)

        plt.xlabel(f'$r$ (GeV$^{{-1}}$)')
        axes.yaxis.set_label_coords(-0.06, 1.04)
        plt.ylabel(f'$V(r)$ (GeV)', rotation='horizontal')

        plt.grid(True)

        plt.savefig(f'plots/potential_{CASE}.png')
        plt.close('all')

        # PLOT WAVE-FUNCTIONS
        for st in STATES:
            fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=figsize)

            quant.e_wf.plot(ax=axes, x='time', y=[f'Pp_proj_{st}', f'Pm_proj_{st}', f'Ppk_proj_{st}'],
                            label=[f'Positive parity', f'Negative parity', f'Pk'], logy=True)

            plt.axhline(y=1., color='k', linestyle='dotted')

            plt.xlabel(f'$\\tau$ (GeV$^{{-1}}$)', )
            plt.ylabel(f'Parity projection of $\\Psi_{st}$')
            plt.grid(True)

            plt.savefig(f'plots/proj_{CASE}_{st}.png')
            plt.close('all')

        for snap in SNAPS:
            fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=figsize)

            print(f'>> PRINTING CASE {CASE}, snap={snap}')

            if CASE in REDUCED_XRANGE:
                axes.set_xlim([0., 10.])

            Itab0, axis_name0 = quant.load_wf(snap, 0, 'p', 0.3)
            Itab0pk, axis_name0p = quant.load_wf(snap, 0, 'pk', 0.3)
            Itab1p, axis_name1p = quant.load_wf(snap, 1, 'p', 0.3)
            Itab1pk, axis_name1pk = quant.load_wf(snap, 1, 'pk', 0.3)
            Itab2p, axis_name2p = quant.load_wf(snap, 2, 'p', 0.3)
            Itab2pk, axis_name2pk = quant.load_wf(snap, 2, 'pk', 0.3)

            raw_r = None
            raw_pk = None

            if Itab2p is not None:
                tab2p = quant.fix_positive_wf(Itab2p)
                tab2p.plot(ax=axes, x='r', y='re_wf_wgt', color='g', label=f'$P^+\\Psi_{2}$', alpha=1.)
                print(f'CASE={CASE}; n=2; proj=p; <r/A>={quant.get_avg_r(Itab2p)/quant.conf["A"]}')
                raw_r = Itab2p['r']
                del Itab2p
                del tab2p

            if Itab2pk is not None:
                tab2pk = quant.fix_positive_wf(Itab2pk)
                tab2pk.plot(ax=axes, x='r', y='re_wf_wgt', color='c', label=f'$P^-_{{\\vec p_k}}\\Psi_{2}$', alpha=1.)
                print(f'CASE={CASE}; n=2; proj=pk; <r/A>={quant.get_avg_r(Itab2pk)/quant.conf["A"]}')
                raw_pk = Itab2pk['pk']
                del Itab2pk
                del tab2pk

            if Itab1p is not None:
                tab1p = quant.fix_positive_wf(Itab1p)
                tab1p.plot(ax=axes, x='r', y='re_wf_wgt', color='r', label=f'$P^+\\Psi_{1}$', alpha=0.7)
                print(f'CASE={CASE}; n=1; proj=p; <r/A>={quant.get_avg_r(Itab1p)/quant.conf["A"]}')
                raw_r = Itab1p['r']
                del Itab1p
                del tab1p

            if Itab1pk is not None:
                tab1pk = quant.fix_positive_wf(Itab1pk)
                tab1pk.plot(ax=axes, x='r', y='re_wf_wgt', color='b', label=f'$P^-_{{\\vec p_k}}\\Psi_{1}$', alpha=0.7)
                print(f'CASE={CASE}; n=1; proj=pk; <r/A>={quant.get_avg_r(Itab1pk)/quant.conf["A"]}')
                raw_pk = Itab1pk['pk']
                del Itab1pk
                del tab1pk

            if Itab0 is not None:
                tab0 = quant.fix_positive_wf(Itab0)
                tab0.plot(ax=axes, x='r', y='re_wf_wgt', color='k', label=f'$P^+\\Psi_{0}$', alpha=0.9)
                print(f'CASE={CASE}; n=0; proj=p; <r/A>={quant.get_avg_r(Itab0)/quant.conf["A"]}')
                raw_r = Itab0['r']
                del Itab0
                del tab0

            if Itab0pk is not None:
                tab0pk = quant.fix_positive_wf(Itab0pk)
                tab0pk.plot(ax=axes, x='r', y='re_wf_wgt', color='m', label=f'$P^-_{{\\vec p_k}}\\Psi_{0}$', alpha=0.9)
                print(f'CASE={CASE}; n=0; proj=pk; <r/A>={quant.get_avg_r(Itab0pk)/quant.conf["A"]}')
                del Itab0pk
                del tab0pk

            if (CASE not in NO_THS) and (raw_r is not None):

                if raw_pk is not None:
                    r, th_2p, avg_rA_2p = wf_2p_cos(raw_r, raw_pk, quant.conf['MASS'], quant.conf['A'])
                    axes.plot(r, th_2p, linestyle=':', label='$\\Psi^{{\\rm th}}_{{2p}}$', color='k')
                    print(f'AVG(r*|PSI(2p)|^2): {avg_rA_2p}')
                    del th_2p
                    del r
                    del raw_pk

                r, th_2s, avg_rA_2s = wf_2s(raw_r, quant.conf['MASS'], quant.conf['A'])
                axes.plot(r, th_2s, linestyle='-.', label=f'$\\Psi^{{\\rm th}}_{{2s}}$', color='k')
                print(f'AVG(r*|PSI(2s)|^2): {avg_rA_2s}')
                del th_2s
                del r

                r, th_1s, avg_rA_1s = wf_1s(raw_r, quant.conf['MASS'], quant.conf['A'])
                axes.plot(r, th_1s, linestyle='--', label=f'$\\Psi^{{\\rm th}}_{{1s}}$', color='k')
                print(f'AVG(r*|PSI(1s)|^2): {avg_rA_1s}')
                del th_1s
                del r

                del raw_r

            plt.axhline(y=0., color='k', linestyle='dotted')
            plt.grid(True)

            tau = float(quant.e_wf.loc[quant.e_wf['step'] == snap]['time'].head(1))

            axes.xaxis.set_label_coords(1.05, -0.025)
            plt.xlabel(f'$r/A$')

            axes.yaxis.set_label_coords(-0.06, 1.04)
            plt.ylabel(f'${{\\rm Re}}\\Psi$', rotation='horizontal')

            handles, labels = axes.get_legend_handles_labels()

            # reverse the order
            axes.legend(handles[::-1], labels[::-1], loc='upper right', bbox_to_anchor=(1.16, 1.18))

            plt.title(f'$step={snap}, \\tau={tau}$')
            plt.savefig(f'plots/figure_{CASE}_{snap}.png')
            plt.close('all')

            del tau
            del handles
            del labels
            del axes
            del fig

        print()
