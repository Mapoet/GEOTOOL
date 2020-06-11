//C source of Mapoet Niphy 
/*
Author    :  Mapoet Niphy
Date      :  2020
Institude :  SHAO
Name      :  
Project   :  

  Created by Mapoet Niphy on 2020/06/10.
  Copyright © 2020年 Mapoet Niphy. All rights reserved.

*/
#ifndef __XFORM__H__
#define __XFORM__H__
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define sign(x,y) ((x)*(y)>0?1:-1)
/* constants -----------------------------------------------------------------*/
#ifndef M_PI
#define M_PI		3.1415926535897932				/* M_PI */
#endif
#define DEG2RAD		(M_PI/180.0)					/* deg->rad */
#define RAD2DEG		(180.0/M_PI)					/* rad->deg */
#define SEC2RAD		4.8481368110953598E-6			/* "->rad */

#define	M_C			299792458.0						/* speed of light(m/sec) */
#define GME			3.986004415E+14					/* geogravity(JGM-3) */
#define GMS			1.32712440017987E+20			/* solar-gravity(DE405) */
#define GMM			4.902801E+12					/* lunar-gravity(DE405) */
/*xform-----------------------------------------------------------------------*/
extern void geod2ecef(const double *gpos, double *epos, double *E);
extern void ecef2geod(const double *epos, double *gpos, double *E);
extern void ecsf2satf(const double *state, double *E);
extern int ele2state(const double *ele, double *state, double *dsde);
extern void state2ele(const double *state, double *ele);
extern void sphfunc_azel(double az, double el, int nmax, double *fc, double *fs);
extern void sphfunc_latlon(double lat, double lon, int nmax, double *fc, double *fs);
extern void cart2pol(double*in, double*lat, double*lon, double*radial);
extern void pol2cart(double lat, double lon, double radial, double*out);
/*
type :gei geo gse gsm mag sm 2 others
time: times[0]=yyyyddd,times[2]=sod*1000 ms.
*/
extern void xform(char* type, int*time, double *in, double *out);
#endif