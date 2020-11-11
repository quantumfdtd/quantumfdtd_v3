/*

   ouputroutines.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <complex>

using namespace std;

#include "mpisolve.h"
#include "grid.h"
#include "outputroutines.h"
#include "potential.h"
#include "latextv.h"

void safe_open(fstream &str, const char *f_name, ios_base::openmode mode)
{
  str.open(f_name, mode);
  if(!str){
    cerr << " >> FATAL ERROR WHEN OPENING FILE: " << f_name << "!!!" << endl;
    exit(1);
  }
}

void outputMeasurements(const double time, const int step) {

	dcomp ener = energyCollect/normalizationCollect;
	dcomp binding = ener - vInfinityCollect/normalizationCollect;
	dcomp rRMS2 = rRMS2Collect/normalizationCollect;  
        //dcomp yAvg = yAvgCollect/normalizationCollect;  

	// output to screen

	cout.precision(8);
	cout.width(dwidth); cout << time;
	cout.width(dwidth); cout << setprecision (7) << ener;
	cout.width(dwidth); cout << setprecision (7) << binding;
	cout.width(dwidth); cout << setprecision (7) << sqrt(real(rRMS2));  
	//cout.width(dwidth); cout << setprecision (7) << real(yAvg);  
	cout << endl;

	if (SAVEDECAY>0){
  		char fname[255];
		std::fstream outfile;

  		sprintf(fname,"%s/decay.dat",DATAFOLD);
		safe_open(outfile, fname, std::ios_base::app);
		if (outfile.is_open()){
			outfile.precision(10);
  			outfile << step << '\t'
  				<< time << '\t'
				<< ener.real() << '\t'
				<< ener.imag() << '\t' 
				<< binding.real() << '\t'
				<< binding.imag() << endl;
			outfile.close();
		}
	}
}

void outputSummaryData(string label, const double time, const int step) {


      dcomp ener = energyCollect/normalizationCollect;
      dcomp binding = ener - vInfinityCollect/normalizationCollect;
      dcomp rRMS2 = rRMS2Collect/normalizationCollect;  
      dcomp xAvg = xAvgCollect/normalizationCollect;  
      dcomp yAvg = yAvgCollect/normalizationCollect;  
      dcomp zAvg = zAvgCollect/normalizationCollect;  

      print_line();
      cout << "==> " << label << " Energy : " << ener << endl;
      cout << "==> " << label << " Binding Energy : " << binding << endl;
      cout << "==> " << label << " r_RMS : " << sqrt(real(rRMS2)) << endl;  
      cout << "==> " << label << " <x> : " << xAvg << endl;  
      cout << "==> " << label << " <y> : " << yAvg << endl;  
      cout << "==> " << label << " <z> : " << zAvg << endl;  
      cout << "==> " << label << " L/r_RMS : " << float(NUM)/sqrt(real(rRMS2)) << endl;  
      print_line();

      // convert label to lower case and convert spaces to underscores for file name
      const int length = label.length();
      for(int i=0; i < length; ++i)
      {
        label[i] = std::tolower(label[i]);
	if (label[i]==' ') label[i]='_';
      }

      // output summary data to output file for later use
      fstream out;
      char fname[255];
      sprintf(fname,"%s/%s.out",DATAFOLD,label.c_str());
      safe_open(out, fname, std::ios_base::app);
      out.precision(10);
      out << step << "\t";
      out << time << "\t";
      out << real(ener) << "\t";
      out << imag(ener) << "\t";
      out << real(binding) << "\t";
      out << imag(binding) << "\t";
      out << real(xAvg) << "\t";
      out << real(yAvg) << "\t";
      out << real(zAvg) << "\t";
      out << NUM << "\t";
      out << A << "\t";
      out << MASS << "\t";
      out << SIGMA << endl;
      out.close();
}

void outputSnapshot(dcomp ***wfnc, char* label) {

  int x;
  static int h=NUM/2;
  static int hx=NUMX/2;
  dcomp data;

  fstream out;
  char fname[255];
  sprintf(fname,"%s/snapshot/wavefunction_%s.dat",DATAFOLD,label);
  safe_open(out, fname, ios::out);
  out.precision(10);

  if (DUMPSNAPS%2 == 1){ // output slices suitable for 2d viewing
    for (int s=0;s<=NUMX+1;s++) {
      x=(nodeID-1)*NUMX + s;	  
      out << x << "\t";
      data = 0.5*(wfnc[s][h][h]+wfnc[s][h+1][h+1]);
      out << scientific << real(data) << "\t";
      out << scientific << imag(data);
      out << endl;
    }
    out << "&&" << endl;
    for (int s=0;s<=NUM+1;s++) {
      out << s << "\t";
      data = 0.5*(wfnc[hx][s][h]+wfnc[hx+1][s][h+1]);
      out << scientific << real(data) << "\t";
      out << scientific << imag(data);
      out << endl;
    }
    out << "&&" << endl;
    for (int s=0;s<=NUM+1;s++) {
      out << s << "\t";
      data = 0.5*(wfnc[hx][h][s]+wfnc[hx+1][h+1][s]);
      out << scientific << real(data) << "\t";
      out << scientific << imag(data);
      out << endl;
    }
  } else { // dump wavefunction
    for (int sx=1;sx<=NUMX;sx++) {
      x=(nodeID-1)*NUMX + sx;
      for (int sy=1;sy<=NUM;sy++) {
        for (int sz=1; sz<=NUM;sz++) {
                  out << x  << "\t";
                  out << sy << "\t";
                  out << sz << "\t";
          	out << get_pos(sx, sy, sz)/A << "\t"; // R IN LATTICE UNITS!!!
                  out << real(wfnc[sx][sy][sz]) << "\t";
                  out << imag(wfnc[sx][sy][sz]);
                  out << endl;
    }}}
  }
  out.close();
  
  return;
}

void outputWavefunction(dcomp ***wfnc, char* label) {

  int x;
  fstream out;
  char fname[255];

  // output full 3d wfnc
  sprintf(fname,"%s/wavefunction_%s.dat",DATAFOLD,label);

  cout << "==> Dumping wave function to " << fname << endl;

  safe_open(out, fname, ios::out);
  out.precision(10);
  for (int sx=1;sx<=NUMX;sx++) {
    x=(nodeID-1)*NUMX + sx;
    for (int sy=1;sy<=NUM;sy++) {
      for (int sz=1; sz<=NUM;sz++) {
                out << x  << "\t";
                out << sy << "\t";
                out << sz << "\t";
		out << get_pos(sx, sy, sz)/A << "\t";
                out << real(wfnc[sx][sy][sz]) << "\t";
                out << imag(wfnc[sx][sy][sz]);
                out << endl;
  }}}
  out.close();

  return;

}

// output v 3d
void outputPotential(char* label) {

  int x;
  fstream out;
  char fname[255];

  // output full 3d wfnc
  sprintf(fname,"%s/potential_%s.dat",DATAFOLD,label);

  cout << "==> Dumping potential to " << fname << endl;

  safe_open(out, fname, ios::out);
  out.precision(8);
  for (int sx=1;sx<=NUMX;sx++) {
    x=(nodeID-1)*NUMX + sx;
    for (int sy=1;sy<=NUM;sy++) {
      for (int sz=1; sz<=NUM;sz++) {
                out << x  << "\t";
                out << sy << "\t";
                out << sz << "\t";
		out << get_pos(sx, sy, sz)/A << "\t"; // R IN LATTICE UNITS!!!
                out << real(v[sx][sy][sz]) << "\t";
                out << imag(v[sx][sy][sz]);
                out << endl;
  }}}
  out.close();

  return;

}

// output v along principal axes
void dumpPotential() {
  char fname[255];
  int h=NUM/2;
  fstream out;

  sprintf(fname,"%s/potential.dat",DATAFOLD);

  safe_open(out, fname, ios::out);
  for (int s=0;s<=NUM+1;s++) {
                out << s << "\t";
                out << v[s][h][h] << "\t";
                out << v[h][s][h] << "\t";
                out << v[h][h][s] << "\t";
                out << endl;
  }
  out.close();
  return;

}

void print_line() {
        for (int i=0;i<4*dwidth;i++) cout << "-"; cout << endl;
        return;
}

