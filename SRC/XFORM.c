
//C source of Mapoet Niphy 
/*
Author    :  Mapoet Niphy
Date      :  2020
Institude :  SHAO
Name      :  XFORM
Project   :  GEOTOOL

  Created by Mapoet Niphy on 2020/06/10.
  Copyright © 2020年 Mapoet Niphy. All rights reserved.

*/

#include "XFORM.h"
#include <math.h>
#ifdef UNIX
#define FILESEP     "/"								/* file path separator */
#else
#define FILESEP     "\\"
#endif

/* basic functions/operations ------------------------------------------------*/
#define NI(x)		(sizeof(x)/sizeof(*x))			/* array item counts */
#define MOD(x,y)	((x)-(y)*floor((x)/(y)))		/* modulo */
#define ABS(x)		((x)<0?-(x):(x))				/* absolute value */
#define MAX(x,y)	((x)>(y)?(x):(y))				/* maximum value */
#define MIN(x,y)	((x)<(y)?(x):(y))				/* minimum value */
#define SWAP(x,y)	do {double _t; _t=y; y=x; x=_t;} while (0) /* swap values */

#define COPY(X,n,m,Y) memcpy(Y,(void *)X,(n)*(m)*sizeof(double))

/* memory allocation/deallocation --------------------------------------------*/
#define VEC(n)		((double *)malloc((n)*    sizeof(double)))
#define MAT(n,m)	((double *)malloc((n)*(m)*sizeof(double)))
#define ZVEC(n)		((double *)calloc((n),    sizeof(double)))
#define ZMAT(n,m)	((double *)calloc((n)*(m),sizeof(double)))
#define FreeMat(mat) free(mat)

/* dot product of 2d/3d vectors ----------------------------------------------*/
#define DOT2(x,y)	((x)[0]*(y)[0]+(x)[1]*(y)[1])
#define DOT(x,y)	((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2])

/* cross product of 3d vectors -----------------------------------------------*/
#define CROSS(x,y,z) do { \
	(z)[0]=(x)[1]*(y)[2]-(x)[2]*(y)[1]; \
	(z)[1]=(x)[2]*(y)[0]-(x)[0]*(y)[2]; \
	(z)[2]=(x)[0]*(y)[1]-(x)[1]*(y)[0]; \
} while (0)

/* norm of 2d/3d vector ------------------------------------------------------*/
#define NORM2(x)	(sqrt(DOT2(x,x)))
#define NORM(x)		(sqrt(DOT(x,x)))

/* normalize 3d vector -------------------------------------------------------*/
#define NORMV(x) do { \
	double _xx=NORM(x); (x)[0]/=_xx; (x)[1]/=_xx; (x)[2]/=_xx; \
} while (0)

/* zero 3d matrix ------------------------------------------------------------*/
#define ZERO(X) do { \
	(X)[0]=(X)[1]=(X)[2]=(X)[3]=(X)[4]=(X)[5]=(X)[6]=(X)[7]=(X)[8]=0.0; \
} while (0)

/* unit 3d matrix ------------------------------------------------------------*/
#define EYE(X) do { \
	(X)[1]=(X)[2]=(X)[3]=(X)[5]=(X)[6]=(X)[7]=0.0;(X)[0]=(X)[4]=(X)[8]=1.0; \
} while (0)

/* transpose 3d matrix -------------------------------------------------------*/
#define Tr(X,Y) do { \
	(Y)[0]=(X)[0];(Y)[1]=(X)[3];(Y)[2]=(X)[6];(Y)[3]=(X)[1];(Y)[4]=(X)[4]; \
	(Y)[5]=(X)[7];(Y)[6]=(X)[2];(Y)[7]=(X)[5];(Y)[8]=(X)[8]; \
} while (0)

/* product of 3d matrix and vector -------------------------------------------*/
#define Mv(X,y,z) do { \
	(z)[0]=(X)[0]*(y)[0]+(X)[3]*(y)[1]+(X)[6]*(y)[2]; \
	(z)[1]=(X)[1]*(y)[0]+(X)[4]*(y)[1]+(X)[7]*(y)[2]; \
	(z)[2]=(X)[2]*(y)[0]+(X)[5]*(y)[1]+(X)[8]*(y)[2]; \
} while (0)

/* product of 3d matrixes ----------------------------------------------------*/
#define MM(X,Y,Z) do { \
	Mv(X,Y,Z); Mv(X,Y+3,Z+3); Mv(X,Y+6,Z+6); \
} while (0)

/* coordinate rotation matrix ------------------------------------------------*/
#define Rx(t,X) do { \
	(X)[0] = 1.0; (X)[1] = (X)[2] = (X)[3] = (X)[6] = 0.0; \
	(X)[4] = (X)[8] = cos(t); (X)[7] = sin(t); (X)[5] = -(X)[7]; \
} while (0)

#define Ry(t,X) do { \
	(X)[4] = 1.0; (X)[1] = (X)[3] = (X)[5] = (X)[7] = 0.0; \
	(X)[0] = (X)[8] = cos(t); (X)[2] = sin(t); (X)[6] = -(X)[2]; \
} while (0)

#define Rz(t,X) do { \
	(X)[8] = 1.0; (X)[2] = (X)[5] = (X)[6] = (X)[7] = 0.0; \
	(X)[0] = (X)[4] = cos(t); (X)[3] = sin(t); (X)[1] = -(X)[3]; \
} while (0)


extern void ecsf2satf(const double *state, double *E)
{
	double crt[3],alt[3],nrad,ncrt,nalt;
	
	CROSS(state,state+3,crt);
	CROSS(crt,state,alt);
	nrad=NORM(state);
	nalt=NORM(alt);
	ncrt=NORM(crt);
	E[0]=state[0]/nrad;
	E[3]=state[1]/nrad;
	E[6]=state[2]/nrad;
	E[1]=alt[0]/nalt;
	E[4]=alt[1]/nalt;
	E[7]=alt[2]/nalt;
	E[2]=crt[0]/ncrt;
	E[5]=crt[1]/ncrt;
	E[8]=crt[2]/ncrt;
}

extern void geod2ecef(const double *gpos, double *epos, double *E)
{
	static const double Re=6378137.0;			/* WGS84 */
	static const double f=1.0/298.257223563;	/* WGS84 */
	double e2,sinp,cosp,sinl,cosl,n;

	e2=f*(2.0-f);
	sinp=sin(gpos[0]*DEG2RAD); cosp=cos(gpos[0]*DEG2RAD);
	sinl=sin(gpos[1]*DEG2RAD); cosl=cos(gpos[1]*DEG2RAD);
	n=Re/sqrt(1.0-e2*sinp*sinp);
	epos[0]=(n+gpos[2])*cosp*cosl;
	epos[1]=(n+gpos[2])*cosp*sinl;
	epos[2]=((1.0-e2)*n+gpos[2])*sinp;

	/* ecef to local tangental coordinate transformation matrix */
	if (E!=0) {
	    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
	    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
	    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
	}
}
extern void ecef2geod(const double *epos, double *gpos, double *E)
{
	static const double Re=6378137.0;			/* WGS84 */
	static const double f=1.0/298.257223563;	/* WGS84 */
	double e2,r2,z,zz,n,sinp,cosp,sinl,cosl,lat;

	e2=f*(2.0-f);
	r2=epos[0]*epos[0]+epos[1]*epos[1];
	z=epos[2]; zz=0; n=Re;
	while (ABS(z-zz)>=1E-4) {
	    zz=z;
	    sinp=z/sqrt(r2+z*z);
	    n=Re/sqrt(1.0-e2*sinp*sinp);
	    z=epos[2]+n*e2*sinp;
	}
	if (r2<1E-12) lat=90.0; else lat=atan(z/sqrt(r2))*180.0/M_PI;
	gpos[0]=lat;
	gpos[1]=atan2(epos[1],epos[0])*180.0/M_PI;
	gpos[2]=sqrt(r2+z*z)-n;
	
	/* ecef to local tangental coordinate transformation matrix */
	if (E!=NULL) {
		sinp=sin(gpos[0]*DEG2RAD); cosp=cos(gpos[0]*DEG2RAD);
		sinl=sin(gpos[1]*DEG2RAD); cosl=cos(gpos[1]*DEG2RAD);
		E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
		E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
		E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
	}
}

double t0(int*time, int*iyr, int*iday, double*ut){
	*iyr = time[0] / 1000;
	*iday = time[0] - *iyr * 1000;
	*ut = time[1] / 3600000.0;
	double fracday = *ut / 24.0;
	double rmjd = 45.0 + (*iyr - 1859)*365.0 + (*iyr - 1861) / 4 + 1.0
		+ *iday - 1.0 + fracday;
	return (rmjd - 51544.5) / 36525.0;
}
void t1(int*time, double*in, double*out, int iverse){
	int iyr, iday;
	double ut,temp = t0(time, &iyr, &iday, &ut),theta=(100.461+36000.770*temp+15.04107*ut)*DEG2RAD,
		X[9];
	Rz(iverse*theta, X);
	//matmul("NN", 3, 3, 1,1.0,X,in,0.0,out);
	Mv(X,in,out);
}
void t2(int*time, double*in, double*out, int iverse){
	int iyr, iday;
	double ut, tt0 = t0(time, &iyr, &iday, &ut), epsion = (23.439 - 0.013*tt0)*DEG2RAD,
		m=(357.528+35999.05*tt0+0.04107*ut)*DEG2RAD,
		cgamma=280.46+36000.772*tt0+0.04107*ut,
		lambdas=(cgamma+(1.915-0.0048*tt0)*sin(m)+0.02*sin(2.0*m))*DEG2RAD,
		X[9],temp[3];
	if (iverse == 1){
		Rx(epsion, X);
		//matmul("NN", 3, 3, 1, 1.0, X, in, 0.0, temp);
		Mv(X,in,temp);
		Rz(lambdas, X);
		//matmul("NN", 3, 3, 1, 1.0, X, temp, 0.0, out);
		Mv(X,temp,out);
	}
	else{
		Rz(iverse*lambdas, X);
		//matmul("NN", 3, 3, 1, 1.0, X, in, 0.0, temp);
		Mv(X,in,temp);
		Rx(iverse*epsion, X);
		//matmul("NN", 3, 3, 1, 1.0, X, temp, 0.0, out);
		Mv(X,temp,out);
	}
}
extern void  pol2cart(double lat, double lon, double radial, double*out){
	double coslat = cos(lat);
	out[0] = radial*coslat*cos(lon);
	out[1] = radial*coslat*sin(lon);
	out[2] = radial*sin(lat);
}
extern void  cart2pol(double*in, double*lat, double*lon, double*radial){
	double row = sqrt(in[0] * in[0] + in[1] * in[1]);
	*radial = sqrt(row*row + in[2] * in[2]);
	*lat = atan2(in[2], row);
	*lon = atan2(in[1], in[0]);
	*lon = *lon + (1.0 -sign(1,*lon))*M_PI;
}
void get_q_c(int*time, double*q_c){
	double q_g[3], temp[3];
	int iyr = time[0] / 1000;
	int iday = time[0] - iyr * 1000;
	double ut = time[1] / 3600000.0;
	double fracday = ut / 24.0;
	double rmjd = 45.0 + (iyr - 1859)*365.0 + (iyr - 1861) / 4 + 1.0
		+ iday - 1.0 + fracday;
	double factor = (rmjd - 46066.0) / 365.25;
	double phi = (78.8 + 4.283e-2*factor)*DEG2RAD;
	double lamda = (289.1 - 1.413e-2*factor)*DEG2RAD;
	pol2cart(phi, lamda, 1.0, q_g);
	t1(time, q_g, temp, -1);
	t2(time, temp, q_c, 1);
}
void t3(int*time, double*in, double*out, int iverse){
	double q_c[3],psi,X[9];
	get_q_c(time, q_c);
	if (q_c[2] = 0.0)
		psi = -sign(1.5707963, q_c[1]);
	else
		psi = -atan(q_c[1] / q_c[2]);
	Rx(iverse*psi, X);
	//matmul("NN", 3, 3, 1, 1.0, X, in, 0.0, out);
	Mv(X,in,out);
}
void t4(int*time, double*in, double*out, int iverse){
	double q_c[3], mu,X[9];
	get_q_c(time, q_c);
	mu = -atan(q_c[0] / sqrt(q_c[1] * q_c[1] + q_c[2] * q_c[2]));
	Rz(iverse*mu, X);
	//matmul("NN", 3, 3, 1, 1.0, X, in, 0.0, out);
	Mv(X,in,out);
}
void t5(int*time, double*in, double*out, int iverse){
	double temp[3],X[9];
	int iyr = time[0] / 1000;
	int iday = time[0] - iyr * 1000;
	double ut = time[1] / 3600000.0;
	double fracday = ut / 24.0;
	double rmjd = 45.0 + (iyr - 1859)*365.0 + (iyr - 1861) / 4 + 1.0
		+ iday - 1.0 + fracday;
	double factor = (rmjd - 46066.0) / 365.25;
	double phi = (78.8 + 4.283e-2*factor)*DEG2RAD;
	double lamda = (289.1 - 1.413e-2*factor)*DEG2RAD;
	if (iverse == 1){
		Rz(lamda, X);
		//matmul("NN", 3, 3, 1, 1.0, X, in, 0.0, temp);
		Mv(X,in,temp);
		Ry(phi-M_PI/2., X);
		//matmul("NN", 3, 3, 1, 1.0, X, temp, 0.0, out);
		Mv(X,temp,out);
	}
	else{
		Ry(iverse*(phi - M_PI / 2.), X);
		//matmul("NN", 3, 3, 1, 1.0, X, in, 0.0, temp);
		Mv(X,in,temp);
		Rz(iverse*lamda, X);
		//matmul("NN", 3, 3, 1, 1.0, X, temp, 0.0, out);
		Mv(X,temp,out);
	}
}
extern void xform(char* type, int*time, double *in, double *out){
	double temp1[3],temp2[3],temp3[3],temp4[3];
	if (strcmp(type, "gei2geo") == 0){
		t1(time, in, out,1);
		return;
	}
	else if (strcmp(type, "gei2gse") == 0){
		t2(time, in, out, 1);
		return;
	}
	else if (strcmp(type, "gei2gsm") == 0){
		t2(time, in, temp1, 1);
		t3(time, temp1,out, 1);
		return;
	}
	else if (strcmp(type, "gei2mag") == 0){
		t1(time, in, temp1, 1);
		t5(time, temp1, out, 1);
		return;
	}
	else if (strcmp(type, "gei2sm") == 0){
		t2(time, in, temp1, 1);
		t3(time, temp1, temp2, 1);
		t4(time, temp2, out, 1);
		return;
	}
	else if (strcmp(type, "geo2gei") == 0){
		t1(time, in, out, -1);
		return;
	}
	else if (strcmp(type, "geo2gse") == 0){
		t1(time, in, temp1, -1);
		t2(time, temp1, out, 1);
		return;
	}
	else if (strcmp(type, "geo2gsm") == 0){
		t1(time, in, temp1, -1);
		t2(time, temp1, temp2, 1);
		t3(time, temp2, out, 1);
		return;
	}
	else if (strcmp(type, "geo2mag") == 0){
		t5(time, in, out, 1);
		return;
	}
	else if (strcmp(type, "geo2sm") == 0){
		t1(time, in, temp1, -1);
		t2(time, temp1, temp2, 1);
		t3(time, temp2, temp3, 1);
		t4(time, temp3, out, 1);
		return;
	}
	else if (strcmp(type, "gse2gei") == 0){
		t2(time, in, out, -1);
		return;
	}
	else if (strcmp(type, "gse2geo") == 0){
		t2(time, in, temp1, -1);
		t1(time, temp1, out, 1);
		return;
	}
	else if (strcmp(type, "gse2gsm") == 0){
		t3(time, in, out, 1);
		return;
	}
	else if (strcmp(type, "gse2mag") == 0){
		t2(time, in, temp1, -1);
		t1(time, temp1, temp2, 1);
		t5(time, temp2, out, 1);
		return;
	}
	else if (strcmp(type, "gse2sm") == 0){
		t3(time, in, temp1,  1);
		t4(time, temp1, out, 1);
		return;
	}
	else if (strcmp(type, "gse2gei") == 0){
		t2(time, in, out, -1);
		return;
	}
	else if (strcmp(type, "gsm2geo") == 0){
		t3(time, in, temp1, -1);
		t2(time, temp1, temp2, -1);
		t1(time, temp2, out, 1);
		return;
	}
	else if (strcmp(type, "gsm2gse") == 0){
		t3(time, in, out, -1);
		return;
	}
	else if (strcmp(type, "gsm2mag") == 0){
		t3(time, in, temp1, -1);
		t2(time, temp1, temp2, -1);
		t1(time, temp2, temp3, 1);
		t5(time, temp3, out, 1);
		return;
	}
	else if (strcmp(type, "mag2gei") == 0){
		t5(time, in, temp1, -1);
		t1(time, temp1, out, -1);
		return;
	}
	else if (strcmp(type, "mag2geo") == 0){
		t5(time, in, out, -1);
		return;
	}
	else if (strcmp(type, "mag2gse") == 0){
		t5(time, in, temp1, -1);
		t1(time, temp1, temp2, -1);
		t2(time, temp2,out, 1);
		return;
	}
	else if (strcmp(type, "mag2gsm") == 0){
		t5(time, in, temp1, -1);
		t1(time, temp1, temp2, -1);
		t2(time, temp2, temp3, 1);
		t3(time, temp3, out, 1);
		return;
	}
	else if (strcmp(type, "mag2sm") == 0){
		t5(time, in, temp1, -1);
		t1(time, temp1, temp2, -1);
		t2(time, temp2, temp3, 1);
		t3(time, temp3, temp4, 1);
		t4(time, temp4, out, 1);
		return;
	}
}