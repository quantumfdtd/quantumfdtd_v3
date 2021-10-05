/*

   grid.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <cmath>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "latextv.h"
#include "grid.h"
#include "outputroutines.h"
#include "initialconditions.h"
#include "potential.h"

// FDTD 3d Schroedinger Eq Solver

// this holds the current values of the wavefunction
dcomp ***w;

// this holds the updated values of the wavefunction
dcomp ***W;

// this holds the current values of the 1st excited wavefunction
dcomp ***w1;

// this holds the current values of the 2nd excited wavefunction
dcomp ***w2;

// this holds the snapshots of the wavefunction
dcomp ****wstore;

// this holds the potential array
dcomp ***v;

// this holds the kinetic term array
dcomp ***t_kin;

// this holds the potential sub array
dcomp ***vsub;

// these hold the alpha and beta arrays which are used during updates
dcomp ***a, ***b;

// temp for point swap
dcomp ***tmp;

// allocate memory arrays
void allocateMemory()
{

    // cout << "==> Allocating memory\n";
    // indices x,y,z

    w = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        w[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            w[sx][sy] = new dcomp[NUM + 2];

    W = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        W[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            W[sx][sy] = new dcomp[NUM + 2];

    w1 = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        w1[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            w1[sx][sy] = new dcomp[NUM + 2];

    w2 = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        w2[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            w2[sx][sy] = new dcomp[NUM + 2];

    v = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        v[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            v[sx][sy] = new dcomp[NUM + 2];

    t_kin = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        t_kin[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            t_kin[sx][sy] = new dcomp[NUM + 2];

    if (POTSUBARRAY())
    {
        vsub = new dcomp **[NUMX + 2];
        for (int sx = 0; sx < NUMX + 2; sx++)
            vsub[sx] = new dcomp *[NUM + 2];
        for (int sx = 0; sx < NUMX + 2; sx++)
            for (int sy = 0; sy < NUM + 2; sy++)
                vsub[sx][sy] = new dcomp[NUM + 2];
    }

    a = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        a[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            a[sx][sy] = new dcomp[NUM + 2];

    b = new dcomp **[NUMX + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        b[sx] = new dcomp *[NUM + 2];
    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            b[sx][sy] = new dcomp[NUM + 2];

    int snaps = 2; // currently only use two snapshots in order to reduce memory footprint
    wstore = new dcomp ***[snaps];
    for (int n = 0; n < snaps; n++)
        wstore[n] = new dcomp **[NUMX + 2];
    for (int n = 0; n < snaps; n++)
        for (int sx = 0; sx < NUMX + 2; sx++)
            wstore[n][sx] = new dcomp *[NUM + 2];
    for (int n = 0; n < snaps; n++)
        for (int sx = 0; sx < NUMX + 2; sx++)
            for (int sy = 0; sy < NUM + 2; sy++)
                wstore[n][sx][sy] = new dcomp[NUM + 2];

    return;
}

// "copies" updated arrays using pointer swap
void copyDown()
{
    tmp = w;
    w = W;
    W = tmp;
}

// initializes the potentials from potential.cpp
void loadPotentialArrays()
{
    int sx, sy, sz;

    initialize_potential();

    for (sx = 0; sx <= NUMX + 1; sx++)
        for (sy = 0; sy <= NUM + 1; sy++)
            for (sz = 0; sz <= NUM + 1; sz++)
            {
                v[sx][sy][sz] = potential(sx, sy, sz);
                b[sx][sy][sz] = 1. / (1. + EPS * v[sx][sy][sz] / ((dcomp)2.));
                a[sx][sy][sz] = (1. - EPS * v[sx][sy][sz] / ((dcomp)2.)) * b[sx][sy][sz];
            }

    destroy_potential();

    if (POTSUBARRAY())
    {

        for (sx = 0; sx <= NUMX + 1; sx++)
            for (sy = 0; sy <= NUM + 1; sy++)
                for (sz = 0; sz <= NUM + 1; sz++)
                {
                    vsub[sx][sy][sz] = potentialSub(sx, sy, sz);
                }
    }
}

// compute energy of a wave function
dcomp wfncEnergy(dcomp ***wfnc)
{
    dcomp res = 0.;

    // First, compute kinetic term
    wfncKinetic(wfnc);

    //R: ADDED CODE FOR USING FLAG RELATIVISTIC AND SWITCHING TO RELATIVISTIC KINETIC TERM IF NEEDED: 20/03/2018
    for (int sx = 1; sx <= NUMX; sx++)
        for (int sy = 1; sy <= NUM; sy++)
            for (int sz = 1; sz <= NUM; sz++)
            {
                res += v[sx][sy][sz] * conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz] +
                       conj(wfnc[sx][sy][sz]) * t_kin[sx][sy][sz];
            }
    return res;
}

// a convenience method for getting energy of the current wave function
dcomp computeEnergy()
{
    return wfncEnergy(w);
}

// compute norm squared
dcomp wfncNorm2(dcomp ***wfnc)
{
    dcomp norm = 0;
    for (int sx = 1; sx <= NUMX; sx++)
        for (int sy = 1; sy <= NUM; sy++)
            for (int sz = 1; sz <= NUM; sz++)
            {
                norm += conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz];
            }
    return norm;
}

// compute expectation value of vinfinity
dcomp vInfinityExpectationValue(dcomp ***wfnc)
{
    dcomp expectation = 0;

    if (POTSUBARRAY())
    {
        for (int sx = 1; sx <= NUMX; sx++)
            for (int sy = 1; sy <= NUM; sy++)
                for (int sz = 1; sz <= NUM; sz++)
                {
                    expectation += conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz] * vsub[sx][sy][sz];
                }
    }
    else
    {
        for (int sx = 1; sx <= NUMX; sx++)
            for (int sy = 1; sy <= NUM; sy++)
                for (int sz = 1; sz <= NUM; sz++)
                {
                    expectation += conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz] * potentialSub(sx, sy, sz);
                }
    }
    return expectation;
}

// compute expectation value of r^2
dcomp r2ExpectationValue(dcomp ***wfnc)
{
    dcomp expectation = 0;
    for (int sx = 1; sx <= NUMX; sx++)
        for (int sy = 1; sy <= NUM; sy++)
            for (int sz = 1; sz <= NUM; sz++)
            {
                expectation += conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz] * distsq(sx, sy, sz);
            }
    return expectation;
}

// compute expectation value of x
dcomp xExpectationValue(dcomp ***wfnc)
{
    dcomp expectation = 0;
    for (int sx = 1; sx <= NUMX; sx++)
        for (int sy = 1; sy <= NUM; sy++)
            for (int sz = 1; sz <= NUM; sz++)
            {
                expectation += conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz] * distx(sx);
            }
    return expectation;
}

// compute expectation value of y
dcomp yExpectationValue(dcomp ***wfnc)
{
    dcomp expectation = 0;
    for (int sx = 1; sx <= NUMX; sx++)
        for (int sy = 1; sy <= NUM; sy++)
            for (int sz = 1; sz <= NUM; sz++)
            {
                expectation += conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz] * disty(sy);
            }
    return expectation;
}

// compute expectation value of z
dcomp zExpectationValue(dcomp ***wfnc)
{
    dcomp expectation = 0;
    for (int sx = 1; sx <= NUMX; sx++)
        for (int sy = 1; sy <= NUM; sy++)
            for (int sz = 1; sz <= NUM; sz++)
            {
                expectation += conj(wfnc[sx][sy][sz]) * wfnc[sx][sy][sz] * distz(sz);
            }
    return expectation;
}

inline dcomp updateRule(int sx, int sy, int sz, double step)
{
    return w[sx][sy][sz] * a[sx][sy][sz] - b[sx][sy][sz] * step * t_kin[sx][sy][sz];
}

// load updated left and right boundaries into W,
// CALL wfncKinetic(w) FIRST!!!
void updateBoundaries(double step)
{

    for (int sx = 1; sx <= NUMX; sx += NUMX - 1)
        for (int sy = 0; sy < NUM + 2; sy++)
            for (int sz = 0; sz < NUM + 2; sz++)
            {
                if (sy >= 1 && sy <= NUM && sz >= 1 && sz <= NUM)
                    W[sx][sy][sz] = updateRule(sx, sy, sz, step);
                else
                    W[sx][sy][sz] = w[sx][sy][sz];
            }
}

// update the grid; note you should always call updatedBondaries before calling this routine
// CALL wfncKinetic(w) FIRST!!!
void updateInterior(double step)
{

    for (int sx = 0; sx < NUMX + 2; sx++)
        for (int sy = 0; sy < NUM + 2; sy++)
            for (int sz = 0; sz < NUM + 2; sz++)
            {
                // no need to update the two slices which were already loaded by updatedBoundaries
                if (sx == 1 || sx == NUMX)
                    continue;
                if (sx >= 1 && sx <= NUMX && sy >= 1 && sy <= NUM && sz >= 1 && sz <= NUM)
                    W[sx][sy][sz] = updateRule(sx, sy, sz, step);
                else
                    W[sx][sy][sz] = w[sx][sy][sz];
            }
}

void recordSnapshot(dcomp ***wfnc, int snap)
{
    for (int sx = 0; sx <= NUMX + 1; sx++)
        for (int sy = 0; sy <= NUM + 1; sy++)
            for (int sz = 0; sz <= NUM + 1; sz++)
                wstore[snap][sx][sy][sz] = wfnc[sx][sy][sz];
}

void normalizeWavefunction(dcomp ***wfnc)
{
    dcomp norm = sqrt(normalizationCollect);
    for (int sx = 0; sx <= NUMX + 1; sx++)
        for (int sy = 0; sy <= NUM + 1; sy++)
            for (int sz = 0; sz <= NUM + 1; sz++)
                wfnc[sx][sy][sz] /= norm;
}
