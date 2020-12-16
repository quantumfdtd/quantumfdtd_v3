/*

   latextv.cpp

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
#include <vector>
#include <gsl/gsl_multifit.h>

using namespace std;

#include "mpisolve.h"
#include "latextv.h"

const int Nrenc = 100; // Length of the external potential position, for reading renc
const int maxline = 128; // maximum line length used in the buffer for reading

struct t_latext_v_data{
  bool full_adj;

  double s_aVeff[Nrenc][Nrenc][Nrenc];
  int n[Nrenc][Nrenc][Nrenc];

  double aVeff_inf;
  double adj_sigma;

  double adjf_a, adjf_b, adjf_sigma;
} ltv_data = {false};

struct t_cartes_v_data{
  bool full_adj;
  int Nmax;

  dcomp ***s_Veff;
  int ***n;

  dcomp Veff_inf;
  dcomp adj_sigma;

  dcomp adjf_a, adjf_b, adjf_sigma;
} l_cartes_v_data = {false};


void get_pos_c000(const int&, const int&, const int&, int*, int*, int*);
void reorder(const int&, const int&, const int&, int*, int*, int*);

bool origin_center_lattice(void)
{
  return POTENTIAL < 100;
}

void get_pos(const int& sx, const int& sy, const int& sz, double *raw_r2, double* Ar, double* dx, double* dy, double* dz)
{
  double iraw_r2, idx, idy, idz;

  if (raw_r2==NULL) raw_r2 = &iraw_r2;
  if (dx==NULL) dx = &idx;
  if (dy==NULL) dy = &idy;
  if (dz==NULL) dz = &idz;

  if ( origin_center_lattice() ){
    // coordinate system is centered in simulation volume 
    *dx = ((double) sx) - ((double)NUMX+1.)/2. + ( ((double)nodeID) - (((double)numNodes)+1.)/2. )*NUMX;
    *dy = ((double) sy) - ((double)NUM+1.)/2.;
    *dz = ((double) sz) - ((double)NUM+1.)/2.;
    *raw_r2 = *dx*(*dx) + *dy*(*dy) + *dz*(*dz);
  } else if ( POTENTIAL < 200 ){
    int nidx, nidy, nidz;
    get_pos_c000(sx, sy, sz, &nidx, &nidy, &nidz);
    *dx = (double) nidx;
    *dy = (double) nidy;
    *dz = (double) nidz;

    *raw_r2 = *dx*(*dx) + *dy*(*dy) + *dz*(*dz);
  } else{
    cerr << "ERROR: INVALID POTENTIAL!!!" << endl;
    exit(1);
  }

  if (Ar!=NULL) *Ar = A*sqrt(*raw_r2);
  return;
}

// determines x distance to center of simulation volume in lattice units
double distx(int sx) 
{ 
	double rA, dx;
	get_pos(sx, 0, 0, NULL, NULL, &dx, NULL, NULL);
	return dx;
}

// determines y distance to center of simulation volume in lattice units
double disty(int sy) 
{ 
	double rA, dy;
	get_pos(0, sy, 0, NULL, NULL, NULL, &dy, NULL);
	return dy;
}

// determines z distance to center of simulation volume in lattice units
double distz(int sz) 
{ 
	double rA, dz;
	get_pos(0, 0, sz, NULL, NULL, NULL, NULL, &dz);
	return dz;
}

// determines square of distance to center of simulation volume in lattice units
double distsq(int sx, int sy, int sz)
{
	double raw_r2;
	get_pos(sx, sy, sz, &raw_r2, NULL, NULL, NULL, NULL);
	return raw_r2;
}

// determines distance to center of simulation volume in physical units
double get_pos(const int& sx, const int& sy, const int& sz)
{
	double rA, dx, dy, dz;
	get_pos(sx, sy, sz, NULL, &rA, NULL, NULL, NULL);
	return rA;
}

void get_pos_c000(const int& sx, const int& sy, const int& sz, int* kx, int* ky, int* kz)
{
  int rsx = sx + NUMX*(nodeID-1);

  if (rsx == 0) *kx = 1;
  else if (rsx == NUM+1) *kx = 0;
  else if (2*rsx<NUM+1) *kx = rsx - 1;
  else *kx = -(NUM - rsx + 1);

  if (sy == 0) *ky = 1;
  else if (sy == NUM+1) *ky = 0;
  else if (2*sy<NUM+1) *ky = sy - 1;
  else *ky = -(NUM - sy + 1);

  if (sz == 0) *kz = 1;
  else if (sz == NUM+1) *kz = 0;
  else if (2*sz<NUM+1) *kz = sz - 1;
  else *kz = -(NUM - sz + 1);
}

inline void reorder(const int& kx, const int& ky, const int& kz, int* k1, int* k2, int* k3)
{
  if( (kx>=ky) && (kx>=kz) ){
    *k1=kx;
    if(ky>=kz){
      *k2=ky;
      *k3=kz;
    }else{
      *k2=kz;
      *k3=ky;
    }
  }else if( (ky>=kx) && (ky>=kz) ){
    *k1=ky;
    if(kx>=kz){
      *k2=kx;
      *k3=kz;
    }else{
      *k2=kz;
      *k3=kx;
    }
  }else{
    *k1=kz;
    if(kx>=ky){
      *k2=kx;
      *k3=ky;
    }else{
      *k2=ky;
      *k3=kx;
    }
  }
}

int get_r2_c000(int sx, int sy, int sz)
{
  int kx, ky, kz;

  get_pos_c000(sx, sy, sz, &kx, &ky, &kz);
  return kx*kx + ky*ky + kz*kz;
}

dcomp read_latext_v(int sx, int sy, int sz)
{
  int kx, ky, kz, k1, k2, k3, n;
  double aVeff;
  bool found=false;

  if (origin_center_lattice()){
    kx = abs((2*sx + 2*NUMX*(nodeID-1) - (NUM+1) )/2);
    ky = abs((2*sy - (NUM+1) )/2);
    kz = abs((2*sz - (NUM+1) )/2);
  }else{
    get_pos_c000(sx, sy, sz, &kx, &ky, &kz);
  }
      
  if (kx<0) kx=-kx;
  if (ky<0) ky=-ky;
  if (kz<0) kz=-kz;

  if( kx<Nrenc && ky<Nrenc && kz<Nrenc ){
   reorder(kx, ky, kz, &k1, &k2, &k3);
   n = ltv_data.n[k1][k2][k3];
    if (n!=0){
       aVeff = ltv_data.s_aVeff[k1][k2][k3]/((double)n);
       found=true;
    }
  }

  if (!found){
    double r = sqrt(kx*kx + ky*ky + kz*kz);

    if (POTFLATR>0 && r>POTFLATR) r = POTFLATR; 

    aVeff = ltv_data.full_adj ?
        ltv_data.adjf_a + ltv_data.adjf_b/r + ltv_data.adjf_sigma*r :
        -0.385/(A*r) + ltv_data.adj_sigma*(A*r); 
  }
  
  return aVeff/A;
}

void charge_latext_v(const char *filename)
{
  int read_lines=0, discarded_lines=0;

  std::vector<int> a_dyn_r2;
  std::vector<double> a_dyn_V;

  if(!filename) filename = "effpot.dat";
  
  ifstream file(filename);
  if(!file){
    cerr << "EXTERNAL POTENTIAL DOES NOT EXIST!!!" << endl;
    exit(-1);
  }

  if (nodeID==1) cout << "Loading potential lattice-format file " << filename << endl;

  char buffer[maxline];
  double buff_r2, buff_re, buff_im;
  int maxk2=0, maxk2_prev=0;

  ltv_data.aVeff_inf = 0.;
  ltv_data.full_adj = false;

  for (int k1=0; k1<Nrenc; k1++)
    for (int k2=0; k2<Nrenc; k2++)
      for (int k3=0; k3<Nrenc; k3++){
	ltv_data.s_aVeff[k1][k2][k3] = 0.;
	ltv_data.n[k1][k2][k3] = 0;
      }

  while(!file.eof()){
    int dta, ta, renc;
    double rimp, aVeff, daVeff;
    
    file.getline(buffer, maxline, '\n');
    if (sscanf(buffer, "%d %d %d %le %le %le",
	       &dta, &ta, &renc, &rimp, &aVeff, &daVeff) != EOF){
      
      ++read_lines;

      if (dta>1 || ta < 5){
	++discarded_lines;
	continue;
      }

      int k1, k2, k3, tmp;
      tmp = renc;
      
      k3 = tmp % Nrenc;
      tmp /= Nrenc;
      k2 = tmp % Nrenc;
      tmp /= Nrenc;
      k1 = tmp % Nrenc;

      ltv_data.s_aVeff[k1][k2][k3] += aVeff;
      ++ltv_data.n[k1][k2][k3];

      //if (nodeID==1) printf("k = %d ; %d ; %d : Veff = %le n(ki) = %d \n", k1,k2,k3,aVeff,ltv_data.n[k1][k2][k3]);

      maxk2 = k1*k1+k2*k2+k3*k3;

      if (POTCRITR*POTCRITR < maxk2){
        a_dyn_r2.push_back(maxk2);
        a_dyn_V.push_back(aVeff);
      }
      if (maxk2>=maxk2_prev){
	maxk2_prev=maxk2;
	ltv_data.aVeff_inf = ltv_data.s_aVeff[k1][k2][k3]/((double)ltv_data.n[k1][k2][k3]);
      }
    }
  }

  int n = a_dyn_r2.size();
  {
    gsl_matrix *a_X, *a_cov;
    gsl_vector *a_y, *a_c;
    double a_chisq;
    gsl_multifit_linear_workspace *fit = NULL;

    if (n>4){
     //if (nodeID==1) cout << "DEBUG MULTIFIT npoints=" << n << endl;
     a_y = gsl_vector_alloc(n);
     a_c = gsl_vector_alloc(3);
     a_X = gsl_matrix_alloc(n,3);
     a_cov = gsl_matrix_alloc(3,3);

     for (int i=0; i<n; i++){
       double r = sqrt(a_dyn_r2[i]);
       double V = a_dyn_V[i];
       //if (nodeID==1) cout << "DEBUG multifit: r=" << r << "\t V=" << V << endl;
       gsl_vector_set(a_y, i, V);
       gsl_matrix_set(a_X, i, 0, 1.);
       gsl_matrix_set(a_X, i, 1, 1./r);
       gsl_matrix_set(a_X, i, 2, r);
     }

     fit = gsl_multifit_linear_alloc(n,3);
     gsl_multifit_linear(a_X, a_y, a_c, a_cov, &a_chisq, fit);
     gsl_multifit_linear_free(fit);

     ltv_data.adjf_a = gsl_vector_get(a_c, 0);
     ltv_data.adjf_b = gsl_vector_get(a_c, 1);
     ltv_data.adjf_sigma = gsl_vector_get(a_c, 2);

     gsl_vector_free(a_y);
     gsl_vector_free(a_c);
     gsl_matrix_free(a_X);
     gsl_matrix_free(a_cov);

     ltv_data.full_adj=true;
    }
  }
  {
    double r = A*sqrt(maxk2_prev);
    ltv_data.adj_sigma = (ltv_data.aVeff_inf+0.385/r)/r;
  }
  if (nodeID==1){
    cout << "File read (node 1)" << endl;
    cout << "Read lines: " << read_lines << "    Discarded lines: " << discarded_lines << endl;
    cout << "maxk2: " << maxk2 << "    aVeff_inf: " << ltv_data.aVeff_inf << endl;
    cout << "POTCRITR: " << POTCRITR << "  N CONSIDERED POINTS: " << n << endl;
    if (ltv_data.full_adj)
      cout << "Used full_adj with a=" << ltv_data.adjf_a << " b=" << ltv_data.adjf_b << " sigma=" << ltv_data.adjf_sigma << endl;
  }
}

dcomp read_external_cartes_v(int sx, int sy, int sz)
{
  int kx, ky, kz, n;
  dcomp Veff;
  bool found=false;

  if (origin_center_lattice()){
    kx = abs((2*sx + 2*NUMX*(nodeID-1) - (NUM+1) )/2);
    ky = abs((2*sy - (NUM+1) )/2);
    kz = abs((2*sz - (NUM+1) )/2);
  }else{
    get_pos_c000(sx, sy, sz, &kx, &ky, &kz);
  }

  double r = sqrt(kx*kx + ky*ky + kz*kz);

  if( abs(kx)<l_cartes_v_data.Nmax && abs(ky)<l_cartes_v_data.Nmax && abs(kz)<l_cartes_v_data.Nmax ){

   kx+=l_cartes_v_data.Nmax;  
   ky+=l_cartes_v_data.Nmax;  
   kz+=l_cartes_v_data.Nmax;  

   n = l_cartes_v_data.n[kx][ky][kz];

   if (n>0){
      Veff = l_cartes_v_data.s_Veff[kx][ky][kz]/((double)n);
      cout << "READING POTENTIAL Veff=" << Veff << endl << "n=" << n << endl << "kx ky kz" << kx << ' ' << ky << ' ' << kz << endl;
      found=true;
   }
  }

  if (!found){
    if (POTFLATR>0 && r>POTFLATR) r = POTFLATR; 
    Veff = l_cartes_v_data.full_adj ?
        l_cartes_v_data.adjf_a + l_cartes_v_data.adjf_b/r + l_cartes_v_data.adjf_sigma*r :
        -0.385/(A*(r>0?r:0.5)) + l_cartes_v_data.adj_sigma*(A*r); 
  }
  
  return Veff;
}

void charge_external_cartes_v(const char *filename)
{
  int read_lines=0;

  std::vector<int> a_dyn_r2;
  std::vector<dcomp> a_dyn_V;

  if(!filename) filename = "effpot.dat";
  
  ifstream file(filename);
  if(!file){
    cerr << "EXTERNAL POTENTIAL DOES NOT EXIST!!!" << endl;
    exit(-1);
  }

  if (nodeID==1) cout << "Loading potential lattice-format file " << filename << endl;

  char buffer[maxline];
  double buff_r2, buff_re, buff_im;
  int maxk2=0, maxk2_prev=0;

  int k1,k2,k3;
  double ReV, ImV;

  l_cartes_v_data.Veff_inf = 0.;
  l_cartes_v_data.full_adj = false;
  l_cartes_v_data.Nmax = 0;
  int Nmax = l_cartes_v_data.Nmax;

  while (!file.eof()){
    file.getline(buffer, maxline, '\n');
    if (sscanf(buffer, "%d %d %d %le %le",
	       &k1, &k2, &k3, &ReV, &ImV) != EOF){

       ++read_lines;
       if (abs(k1)>Nmax) Nmax = abs(k1);
       if (abs(k2)>Nmax) Nmax = abs(k2);
       if (abs(k3)>Nmax) Nmax = abs(k3);
    }
  }
  
  l_cartes_v_data.Nmax = Nmax;

  file.clear();
  file.seekg(0, ios::beg);

  l_cartes_v_data.n = new int**[2*Nmax+2];
  for (int sx=0;sx<2*Nmax+2;sx++) l_cartes_v_data.n[sx]=new int*[2*Nmax+2];
  for (int sx=0;sx<2*Nmax+2;sx++) for (int sy=0;sy<2*Nmax+2;sy++) l_cartes_v_data.n[sx][sy] = new int[2*Nmax+2];

  l_cartes_v_data.s_Veff = new dcomp**[2*Nmax+2];
  for (int sx=0;sx<2*Nmax+2;sx++) l_cartes_v_data.s_Veff[sx]=new dcomp*[2*Nmax+2];
  for (int sx=0;sx<2*Nmax+2;sx++) for (int sy=0;sy<2*Nmax+2;sy++) l_cartes_v_data.s_Veff[sx][sy] = new dcomp[2*Nmax+2];

  for (int k1=2*Nmax+1; k1>=0; --k1)
    for (int k2=2*Nmax+1; k2>=0; --k2)
      for (int k3=2*Nmax+1; k3>=0; --k3){
	l_cartes_v_data.s_Veff[k1][k2][k3] = 0.;
	l_cartes_v_data.n[k1][k2][k3] = 0;
      }

  while(!file.eof()){
    file.getline(buffer, maxline, '\n');
    if (sscanf(buffer, "%d %d %d %le %le",
	       &k1, &k2, &k3, &ReV, &ImV) != EOF){
      
      maxk2 = k1*k1+k2*k2+k3*k3;

      //if (nodeID==1) printf("k = %d ; %d ; %d : Veff = (%le,%le), n(ki) = %d \n", k1,k2,k3,ReV,ImV,l_cartes_v_data.n[k1][k2][k3]);

      k1+=Nmax;
      k2+=Nmax;
      k3+=Nmax;

      l_cartes_v_data.s_Veff[k1][k2][k3] += dcomp(ReV,ImV);
      ++l_cartes_v_data.n[k1][k2][k3];

      if (POTCRITR*POTCRITR < maxk2 && maxk2>0){
        a_dyn_r2.push_back(maxk2);
        a_dyn_V.push_back(dcomp(ReV,ImV));
      }
      if (maxk2>=maxk2_prev){
	maxk2_prev=maxk2;
	l_cartes_v_data.Veff_inf = l_cartes_v_data.s_Veff[k1][k2][k3]/((double)l_cartes_v_data.n[k1][k2][k3]);
      }
    }
  }

  int n = a_dyn_r2.size();
  {
    gsl_matrix *a_X, *a_cov;
    gsl_vector *a_y, *a_c;
    double a_chisq;
    gsl_multifit_linear_workspace *fit = NULL;

    if (n>4){
     //if (nodeID==1) cout << "DEBUG MULTIFIT npoints=" << n << endl;
     a_y = gsl_vector_alloc(n);
     a_y = gsl_vector_alloc(n);
     a_c = gsl_vector_alloc(3);
     a_X = gsl_matrix_alloc(n,3);
     a_cov = gsl_matrix_alloc(3,3);

     for (int i=0; i<n; i++){
       double r = sqrt(a_dyn_r2[i]);
       double V = a_dyn_V[i].real();
       //if (nodeID==1) cout << "DEBUG multifit: r=" << r << "\t V=" << V << endl;
       gsl_vector_set(a_y, i, V);
       gsl_matrix_set(a_X, i, 0, 1.);
       gsl_matrix_set(a_X, i, 1, r>0 ? 1./r : 0.5*0.5*0.5);
       gsl_matrix_set(a_X, i, 2, r);
     }

     fit = gsl_multifit_linear_alloc(n,3);
     gsl_multifit_linear(a_X, a_y, a_c, a_cov, &a_chisq, fit);
     gsl_multifit_linear_free(fit);

     l_cartes_v_data.adjf_a = gsl_vector_get(a_c, 0);
     l_cartes_v_data.adjf_b = gsl_vector_get(a_c, 1);
     l_cartes_v_data.adjf_sigma = gsl_vector_get(a_c, 2);

     for (int i=0; i<n; i++){
       double r = sqrt(a_dyn_r2[i]);
       double V = a_dyn_V[i].imag();
       //if (nodeID==1) cout << "DEBUG multifit: r=" << r << "\t V=" << V << endl;
       gsl_vector_set(a_y, i, V);
     }

     fit = gsl_multifit_linear_alloc(n,3);
     gsl_multifit_linear(a_X, a_y, a_c, a_cov, &a_chisq, fit);
     gsl_multifit_linear_free(fit);

     l_cartes_v_data.adjf_a += dcomp(0,gsl_vector_get(a_c, 0));
     l_cartes_v_data.adjf_b += dcomp(0,gsl_vector_get(a_c, 1));
     l_cartes_v_data.adjf_sigma += dcomp(0,gsl_vector_get(a_c, 2));

     gsl_vector_free(a_y);
     gsl_vector_free(a_c);
     gsl_matrix_free(a_X);
     gsl_matrix_free(a_cov);

     l_cartes_v_data.full_adj=true;
    }
  }
  {
    double r = A*sqrt(maxk2_prev);
    l_cartes_v_data.adj_sigma = (l_cartes_v_data.Veff_inf+0.385/r)/r;
  }
  if (nodeID==1){
    cout << "File read (node 1)" << endl;
    cout << "Read lines: " << read_lines << endl;
    cout << "maxk2: " << maxk2 << "    Veff_inf: " << l_cartes_v_data.Veff_inf << endl;
    cout << "POTCRITR: " << POTCRITR << "  N CONSIDERED POINTS: " << n << endl;
    if (l_cartes_v_data.full_adj){
      cout << "Used full_adj with a=" << l_cartes_v_data.adjf_a << " b=" << l_cartes_v_data.adjf_b << " sigma=" << l_cartes_v_data.adjf_sigma << endl;
      cout << "Used full_adj with a=" << l_cartes_v_data.adjf_a << " b=" << l_cartes_v_data.adjf_b << " sigma=" << l_cartes_v_data.adjf_sigma << endl;
    }
  }
}

void destroy_external_cartes_v(void)
{
  int Nmax = l_cartes_v_data.Nmax;

  for (int sx=0;sx<2*Nmax+2;sx++) for (int sy=0;sy<2*Nmax+2;sy++) delete l_cartes_v_data.n[sx][sy];
  for (int sx=0;sx<2*Nmax+2;sx++) delete l_cartes_v_data.n[sx];
  delete l_cartes_v_data.n;

  for (int sx=0;sx<2*Nmax+2;sx++) for (int sy=0;sy<2*Nmax+2;sy++) delete l_cartes_v_data.s_Veff[sx][sy];
  for (int sx=0;sx<2*Nmax+2;sx++) delete l_cartes_v_data.s_Veff[sx];
  delete l_cartes_v_data.s_Veff;
}
