/*

   outputroutines.h

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __outputroutines_h__
#define __outputroutines_h__

#include <iostream>
#include <fstream>

const int dwidth = 25;

void outputMeasurements(double time, int step);
void outputSnapshot(dcomp ***wfnc, char* label);
void outputWavefunction(dcomp ***wfnc, char* label);
void outputSummaryData(string label, double time, int step);
void outputPotential(char* label);
void dumpPotential();
void print_line();

// output files
extern std::fstream energy_out;

#endif /* __outputroutines_h__ */
