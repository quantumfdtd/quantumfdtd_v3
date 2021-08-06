/*

   latextv.h

   Copyright (c) Rafael L. Delgado

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __latextv_h__
#define __latextv_h__


double distx(int sx);
double disty(int sy);
double distz(int sz);
double distsq(int sx, int sy, int sz);
double get_pos(const int& sx, const int& sy, const int& sz);

bool origin_center_lattice(void);
void get_pos(const int& sx, const int& sy, const int& sz, double *raw_r2, double* Ar, double* dx, double* dy, double* dz);

bool center_000(void);

dcomp read_latext_v(int sx, int sy, int sz);
dcomp read_latext_v_c000(int sx, int sy, int sz);
void load_latext_v(const char *filename);

dcomp read_external_cartes_v(int sx, int sy, int sz);
void load_external_cartes_v(const char *filename);
void destroy_external_cartes_v(void);

#endif /* __latextv_h__ */
