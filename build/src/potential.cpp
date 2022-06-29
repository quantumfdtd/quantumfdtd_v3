/*

 potential.cpp

 Copyright (c) Michael Strickland

 GNU General Public License (GPLv3)
 See detailed text in license directory

 */

#include <cmath>

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "grid.h"
#include "potential.h"
#include "externalv.h"
#include "latextv.h"

// global variable useful for subroutines
double dx, dy, dz, raw_r2, r;

dcomp potential(int sx, int sy, int sz)
{
    double temp, iV, rV;
    double res, err;
    double m12, B, wc, e, rho;

    // note that this has no effect on loading of external potentials
    get_pos(sx, sy, sz, &raw_r2, &r, &dx, &dy, &dz);
    rho = A * sqrt(dx * dx + dy * dy);

    switch (POTENTIAL % 100)
    {
    case 0:
        // none
        return 0.;
        break;
    case 1:
        // cubic well
        if ((dx > -NUM / 4 && dx <= NUM / 4) && (dy > -NUM / 4 && dy <= NUM / 4) && (dz > -NUM / 4 && dz <= NUM / 4))
            return -10.0;
        else
            return 0.0;
        break;
    case 2:
        // coulomb
        if (r < A)
            r = A;
        return -1. / r + 1. / A;
    case 3:
        // elliptical coulomb
        dz *= 2;
        r = A * sqrt(dx * dx + dy * dy + dz * dz);
        if (r < A)
            r = A;
        return -1. / r + 1. / A;
        break;
    case 4:
        // 3d harmonic oscillator
        return r * r / 2;
        break;
    case 5:
        // Complex 3d harmonic oscillator
        return dcomp(1., 1.) * dcomp(r * r / 2, 0.);
        break;
    case 6:
        // cornell potential
        // units here are GeV for energy/momentum and GeV^(-1) for distance
        if (r > 5.5745)
            r = 5.5745;
        else if (r < A)
            r = A;
        return -0.385 / r + SIGMA * r + 4 * MASS;
        break;
    // file based
    case 90:
        // read external potential, as a table separated by spaces
        // i j k Re(V) Im(V)
        return read_external_cartes_v(sx, sy, sz);
        break;
    case 91:
        // read external potential, as a table separated by spaces
        // R^2 Re(V) Im(V)
        return read_external_radial_v(raw_r2);
        break;
    case 92:
        // read external potential, in lattice format
        //  The columns are

        // dt/a t/a r_{enc} r_{imp}/a aV_{eff}, d aV_{eff}
        //
        // r_{enc} is an encoding of the distance. r_{imp} is the improved distance in lattice units.
        // For you, the former is the one you want to use, i.e.
        // r_{enc}=30201 means 3 steps in one direction, 2 steps in another and 1 step in a third.
        // You can decode it like that. The encoding is meaningful, so we should keep it if your code can be adapted to deal with it.
        //
        // You can ignore any results with dt/a>1.
        // Then you can take for each r_{enc} an average of the results for t/a>4 as an estimate of the static energy. That should be sufficient for our purposes.
        return read_latext_v(sx, sy, sz);
    default:
        return 0.;
        break;
    }
}

// returns value of potential which should be subtracted when computing binding energies
dcomp potentialSub(int sx, int sy, int sz)
{
    double iV, rV;

    // note that this has no effect on loading of external potentials
    get_pos(sx, sy, sz, &raw_r2, &r, &dx, &dy, &dz);

    if (r < A)
        r = A;

    switch (POTENTIAL % 100)
    {
    case 0:
    case 1:
        return 0.;
        break;
    case 2:
    case 3:
        return 1. / A;
        break;
    case 4:
    case 5:
        return 0.;
        break;
    case 6:
        r = 5.5745;
        return -0.385 / r + SIGMA * r + 4 * MASS;
        break;
    // file based
    case 90:
    case 91:
    case 92:
        return 0.;
        break;
    default:
        return 0.;
        break;
    }
    return 0.;
}

// Initialize potential, important for external potential file
void initialize_potential()
{
    switch (POTENTIAL % 100)
    {
    case 90:
        load_external_cartes_v(EXTPOT);
        break;
    case 91:
        load_external_radial_v(EXTPOT);
        break;
    case 92:
        load_latext_v(EXTPOT);
        break;
    default:
        break;
    }
}

// Destroy potential, important for external potential file
void destroy_potential()
{
    if (POTENTIAL % 100 == 90)
        destroy_external_cartes_v();
    if (POTENTIAL % 100 == 91)
        destroy_external_v_eval();
}
