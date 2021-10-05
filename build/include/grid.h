/*

   grid.h

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __grid_h__
#define __grid_h__

// this holds the current values of the wavefunction
extern dcomp ***w;

// this holds the updated values of the wavefunction
extern dcomp ***W;

// this holds the current values of the 1st excited wavefunction
extern dcomp ***w1;

// this holds the current values of the 2nd excited wavefunction
extern dcomp ***w2;

// this holds the snapshots of the wavefunction
extern dcomp ****wstore;

// this is a temporary pointer used during array swap in copyDown
extern dcomp ***tmp;

// this holds the potential array
extern dcomp ***v;

// this holds the kinetic term after a call to wfncKinetic
extern dcomp ***t_kin;

// these hold the alpha and beta arrays which are used during updates
extern dcomp ***a, ***b;

// number of snapshots to capture in memory
extern int nsnaps;

// return energy of the passed wavefnc
dcomp wfncEnergy(dcomp ***wfnc);

// return norm squared of the passed wavefnc
dcomp wfncNorm2(dcomp ***wfnc);

// return expectation value of vInfinity
dcomp vInfinityExpectationValue(dcomp ***wfnc);

// return expectation value of r^2
dcomp r2ExpectationValue(dcomp ***wfnc);

// return expectation value of x, y, and z
dcomp xExpectationValue(dcomp ***wfnc);
dcomp yExpectationValue(dcomp ***wfnc);
dcomp zExpectationValue(dcomp ***wfnc);

// return energy of the current wavefnc
dcomp computeEnergy();

// other methods
void allocateMemory();
void updateBoundaries(double eps);
void updateInterior(double eps);
void copyDown();
void loadPotentialArrays();
void normalizeWavefunction(dcomp ***wfnc);
void recordSnapshot(dcomp ***wfnc, int step);

#endif /* __grid_h__ */
