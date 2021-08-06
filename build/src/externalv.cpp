/*

   externalv.cpp

   Copyright (c) Rafael L. Delgado 

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cerrno>
#include <complex>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

#include "mpisolve.h"
#include "externalv.h"

void load_external_file(const char *filename);

// This holds the external potential v vs. r2

struct external_v_int_t{
  gsl_interp_accel *acc_re, *acc_im;
  gsl_spline *spline_re, *spline_im;
  
  double *vec_r2;
  double *vec_re, *vec_im;

  int length;
} ev_int = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};


dcomp read_external_radial_v(double raw_r2)
{
  double buff_re, buff_im;

  buff_re = gsl_spline_eval(ev_int.spline_re, raw_r2, ev_int.acc_re);
  buff_im = gsl_spline_eval(ev_int.spline_im, raw_r2, ev_int.acc_im);

  return dcomp(buff_re, buff_im);
}

void load_external_radial_v(const char *filename)
{
  if(!filename) filename="potential.txt";
  load_external_file(filename);
}

void load_external_file(const char *filename)
{
  ifstream file(filename);

  cerr << "CHECKING FILE..." << endl;
  if(!file){
    cerr << "EXTERNAL POTENTIAL DOES NOT EXIST!!!" << endl;
    exit(-1);
  }


  int maxline = 128; // maximum line length used in the buffer for reading
  char buffer[maxline];
  double buff_r2, buff_re, buff_im;

  int counter = 0;

  ev_int.length = 0;

  while(!file.eof()){
    file.getline(buffer, maxline, '\n');
    if (sscanf(buffer, "%*e %*e %*e") != EOF){
      ev_int.length++;
    }
  }
  file.clear();
  file.seekg(0, ios_base::beg);

  ev_int.vec_r2 = new double[ev_int.length];
  ev_int.vec_re = new double[ev_int.length];
  ev_int.vec_im = new double[ev_int.length]; 

  ev_int.acc_re = gsl_interp_accel_alloc();
  ev_int.acc_im = gsl_interp_accel_alloc();
  
  ev_int.spline_re = gsl_spline_alloc(gsl_interp_cspline, ev_int.length);
  ev_int.spline_im = gsl_spline_alloc(gsl_interp_cspline, ev_int.length);
  
  while(!file.eof()){
    file.getline(buffer, maxline, '\n');
    if (sscanf(buffer, "%le %le %le", &buff_r2, &buff_re, &buff_im) != EOF){
      if (counter >= ev_int.length){
	cerr << "ERROR: UNCORRECTLY PREDICTED EXTERNAL POTENTIAL FILE LENGTH!!!" << endl;
	exit(1);
      }
      ev_int.vec_r2[counter] = buff_r2;
      ev_int.vec_re[counter] = buff_re;
      ev_int.vec_im[counter] = buff_im;
      // cerr << counter << " >> r2, re = " << ev_int.vec_r2[counter] << '\t' << ev_int.vec_re[counter] << endl;
      counter++;
    }
  }

  gsl_spline_init(ev_int.spline_re, ev_int.vec_r2, ev_int.vec_re, ev_int.length);
  gsl_spline_init(ev_int.spline_im, ev_int.vec_r2, ev_int.vec_im, ev_int.length);

  return;
}

void destroy_external_v_eval()
{
  gsl_spline_free(ev_int.spline_re);
  gsl_spline_free(ev_int.spline_im);

  gsl_interp_accel_free(ev_int.acc_re);
  gsl_interp_accel_free(ev_int.acc_im);

  delete ev_int.vec_r2;
  delete ev_int.vec_re;
  delete ev_int.vec_im;
  
  ev_int.spline_re = NULL;
  ev_int.spline_im = NULL;
  
  ev_int.acc_re = NULL;
  ev_int.acc_im = NULL;

  ev_int.vec_r2 = NULL;
  ev_int.vec_re = NULL;
  ev_int.vec_im = NULL;

  return;
}
