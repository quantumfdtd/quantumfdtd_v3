/*

   potential.h

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __potential_h__
#define __potential_h__

dcomp potential(int sx, int sy, int sz);
dcomp potentialSub(int sx, int sy, int sz);
double alphas(double mu);
double alphasp(double mu);
double mu(double t, double fac);
double mup(double t, double fac);


double phir(double z);
double psi1(double z);
double psi2(double z);

void initialize_potential(void);
void destroy_potential(void);

void initialize_potential_sub(void);
void destroy_potential_sub(void);

inline bool POTSUBARRAY(void){ return POTENTIAL == 36;}

#endif /* __potential_h__ */
