#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import tarfile


# #################################################################################


class quantumfdtd:
    def __init__(self):
        self.case = None
        self.conf = None
        self.e_dec = None
        self.e_wf = None

    def __init__(self, case, config=None):
        self.case = case
        self.load_case(self.case, config)

    def load_case(self, case, config=None):
        self.case = case
        self.e_dec = None
        self.e_wf = None

        name = config if config is not None else f'{case}/input/params.txt'

        print(f'Loading config file {name}')
        self.conf = pd.read_csv(name, comment='/', delim_whitespace=True, header=None).set_index(0).T.to_dict('index')[1]
        self.conf['center_on_lattice'] = int(self.conf['POTENTIAL']) < 100

        self.load_energies()

    def load_energies(self):
        fil_dec_typ = {"step": np.int32, "time": np.float64, "re_e_b": np.float64, "im_e_b": np.float64,
                       "re_e": np.float64, "im_e": np.float64}

        fil_dec_names = ["step", "time", "re_e_b", "im_e_b", "re_e", "im_e"]

        fil_e_wf_typ = {"step": np.int32, "time": np.float64, "re_e_b": np.float64, "im_e_b": np.float64,
                        "re_e": np.float64, "im_e": np.float64, "xAvg": np.float64, "yAvg": np.float64,
                        "zAvg": np.float64, "T": np.float64, "XI": np.float64, "eB": np.float64, "Kx": np.float64}

        fil_e_wf_names = ["step", "time", "re_e_b", "im_e_b", "re_e", "im_e",
                          "xAvg", "yAvg", "zAvg", "T", "XI", "eB", "Kx"]

        self.e_dec = pd.read_csv(f'{self.case}/data/decay.dat', delim_whitespace=True, header=None, dtype=fil_dec_typ,
                                 names=fil_dec_names)[['step', 'time', 're_e', 'im_e']]

        fil_e_wf0 = pd.read_csv(f'{self.case}/data/ground_state.out', delim_whitespace=True, header=None,
                                dtype=fil_e_wf_typ, names=fil_e_wf_names)[['step', 'time', 're_e', 'im_e']]
        fil_e_wf0.rename({'re_e': 're_e_wf0', 'im_e': 'im_e_wf0'}, axis='columns', inplace=True)

        fil_e_wf1 = pd.read_csv(f'{self.case}/data/first_excited_state.out', delim_whitespace=True, header=None,
                                dtype=fil_e_wf_typ, names=fil_e_wf_names)[['re_e', 'im_e']].add_suffix('_wf1')

        fil_e_wf2 = pd.read_csv(f'{self.case}/data/second_excited_state.out', delim_whitespace=True, header=None,
                                dtype=fil_e_wf_typ, names=fil_e_wf_names)[['re_e', 'im_e']].add_suffix('_wf2')

        self.e_wf = fil_e_wf0.join([fil_e_wf1, fil_e_wf2])

        return [self.e_dec, self.e_wf]

    def load_potential(self):
        pot_names = ["x", "y", "z", "r", "re_v", "im_v"]

        pot_typ = {"x": np.int32, "y": np.int32, "z": np.int32, "r": np.float64,
                   "re_v": np.float64, "im_v": np.float64}

        pot_file = f'{self.case}/data/potential.tgz'
        if (not tarfile.is_tarfile(pot_file)):
            raise Exception(f'Error: file {pot_file} does not exist or is not a tarball')

        tar = tarfile.open(pot_file, 'r:gz')
        csv_files = [tar.extractfile(f) for f in tar.getmembers() if f.name.endswith('.dat')]

        potential = pd.concat([pd.read_csv(fil, delim_whitespace=True, header=None,
                                           dtype=pot_typ, names=pot_names) for fil in csv_files])

        return potential

    def load_proj_components(self):
        if self.e_wf is None:
            self.load_energies()

        for i in range(3):
            self.e_wf[f'case_{i}'] = \
                self.e_wf['step'].apply(lambda st: f'{self.case}/data/snapshot/wavefunction_{st}_{i}_proj.dat'
                                        if st > 0 else f'{self.case}/data/wavefunction_{i}_proj.dat')

            self.e_wf[f'proj_{i}'] = self.e_wf[f'case_{i}'].apply(lambda f: pd.read_csv(f, delim_whitespace=True))
            self.e_wf.drop(columns=f'case_{i}', inplace=True)

            concat = pd.concat(self.e_wf[f'proj_{i}'].values, axis=0, ignore_index=True)
            self.e_wf.drop(columns=f'proj_{i}', inplace=True)

            self.e_wf = self.e_wf.join(concat.add_suffix(f'_proj_{i}'))
            del concat

        return self.e_wf

    def get_snaps(self):
        return self.e_wf.step.unique()

    def load_wf(self, snap, state, figure, min_sep_edge=-1., normalization=True):
        NUM = int(self.conf['NUM'])
        center_on_lattice = self.conf['center_on_lattice']

        if snap < 0:
            name = f'{self.case}/data/wavefunction_{state}_{figure}_norm.dat.gz'
        else:
            name = f'{self.case}/data/snapshot/wavefunction_{snap}_{state}_{figure}_norm.dat.gz'

        if not os.path.exists(name) or os.path.getsize(name) == 0:
            return [None, None]

        if figure != 'pk':
            names = ['x', 'y', 'z', 'r', 're_wf', 'im_wf']
            dtype = {"x": np.int32, "y": np.int32, "z": np.int32, "r": np.float64,
                     "re_wf": np.float64, "im_wf": np.float64}
        else:
            names = ['x', 'y', 'z', 'r', 'pk', 're_wf', 'im_wf']
            dtype = {"x": np.int32, "y": np.int32, "z": np.int32, "r": np.float64,
                     "pk": np.float64, "re_wf": np.float64, "im_wf": np.float64}

        print(f'Reading file {name}')

        tab = pd.read_csv(name, delim_whitespace=True, header=None, dtype=dtype, names=names)

        if ((min_sep_edge < 1.0) and (min_sep_edge > 0.0)):
            lim = int(NUM * (1.0 - min_sep_edge))
            if center_on_lattice:
                cen = NUM * 0.5
                lin = int(cen - lim * 0.5)
                lsu = int(cen + lim * 0.5)
                tab.drop(tab[(tab.x < lin) | (tab.x > lsu) | (tab.y < lin) | (tab.y > lsu)
                             | (tab.z < lin) | (tab.z > lsu)].index, inplace=True)
            else:
                tab.drop(tab[(tab.x > lim) | (tab.y > lim) | (tab.z > lim)].index, inplace=True)

        print(f'Processing file {name}')

        if figure != 'pk':
            tab['re_wf_wgt'] = tab['re_wf']
            tab['im_wf_wgt'] = tab['im_wf']
            axis_name = f'$\\Psi_{state}$'
        else:
            tab['re_wf_wgt'] = tab['re_wf'] * tab['r'] / tab['pk']
            tab['im_wf_wgt'] = tab['im_wf'] * tab['r'] / tab['pk']
            axis_name = f'$(\\cos^{{-1}}\\theta)\\,\\Psi_{state}$'

        if normalization:
            norm = 1. / np.sqrt((tab['re_wf'] ** 2 + tab['im_wf'] ** 2).sum())
            for col in ['re_wf', 'im_wf', 're_wf_wgt', 'im_wf_wgt']:
                tab[col] = norm * tab[col]

        return[tab, axis_name]

    def fix_positive_wf(self, wf, lim=0.05):
        ret = wf
        compar = ret[ret['r'] < max(1., lim * self.conf['NUM']) * self.conf['A']][['r', 're_wf_wgt']]
        if compar['re_wf_wgt'].sum() < 0.:
            ret['re_wf_wgt'] = -1 * ret['re_wf_wgt']
            ret['im_wf_wgt'] = -1 * ret['im_wf_wgt']
        return ret

    def get_avg_r(self, wf):
        num = self.conf['A'] * ((wf['r'] * (wf['re_wf'] ** 2 + wf['im_wf'] ** 2)).sum())
        den = (wf['re_wf'] ** 2 + wf['im_wf'] ** 2).sum()
        return num / den
