/*

   externalv.h

   Copyright (c) Rafael L. Delgado

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#ifndef __externalv_h__
#define __externalv_h__

dcomp read_external_v(double raw_r2);
dcomp read_external_v_sub(double raw_r2);
void charge_external_v(const char *filename);
void charge_external_v_sub(const char *filename);
void destroy_external_v_eval();

#endif /* __externalv_h__ */
