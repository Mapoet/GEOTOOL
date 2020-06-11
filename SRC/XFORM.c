
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

/* matrix(3x2)xvector(2x1) ---------------------------------------------------*/
#define Mv2(X,y,z) \
{ \
	(z)[0]=(X)[0]*(y)[0]+(X)[3]*(y)[1]; \
	(z)[1]=(X)[1]*(y)[0]+(X)[4]*(y)[1]; \
	(z)[2]=(X)[2]*(y)[0]+(X)[5]*(y)[1]; \
}

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

/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite orbit element to position/velocity
% [func]   : transform satellite orbit elements to state(position/velocity)
%            calculate partial derivatives of state by orbit elements
% [argin]  : ele   = orbit elements(m,deg) [a,e,i,OMG,omg,M]
% [argout] : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
%            (dsde) = partial derivatives of state by orbit elements(6x6)
% [note]   : error if a<=0,e<=0,1<=e
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (ç«, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
*-----------------------------------------------------------------------------*/
/* transform satellite orbit elements to state -------------------------------*/
extern int ele2state(const double *ele, double *state, double *dsde)
{
	int k;
	double a,e,i,OMG,omg,M,E,f,rr,sf,n,r[3],v[3],PQW[9],R1[9],R2[9],R3[9];
	double ar2,ar3,drda[2],drde[2],drdM[2],dvda[2],dvde[2],dvdM[2];
	double dPQdi[6],dPQdO[6],dPQdo[6];

	a=ele[0]; e=ele[1]; i=ele[2]*M_PI/180.0; OMG=ele[3]*M_PI/180.0;
	omg=ele[4]*M_PI/180.0; M=ele[5]*M_PI/180.0;
	if (a<=0.0||e<=0.0||1.0<=e) return 0;
	
	/* solve Kepler equation */
	for (k=1,E=M;k<20;k++) {
	    f=E-e*sin(E)-M; E=E-f/(1.0-e*cos(E)); if (ABS(f)<=1E-14) break;
	}
	rr=a*(1.0-e*cos(E)); f=1.0-e*e; sf=sqrt(f); n=sqrt(GME/(a*a*a));
	r[0]=a*(cos(E)-e);
	r[1]=a*sf*sin(E);
	v[0]=-n*a*a/rr*sin(E);
	v[1]=n*a*a/rr*sf*cos(E);
	r[2]=v[2]=0.0;
	Rz(-OMG,R1); Rx(-i,R2); MM(R1,R2,R3);
	Rz(-omg,R1); MM(R3,R1,PQW);
	Mv(PQW,r,state);
	Mv(PQW,v,state+3);
	
	/* partial derivatives of state by orbit elements */
	if (dsde!=NULL) {
		ar2=a*a/(rr*rr); ar3=ar2*a/rr;
		drda[0]=r[0]/a;
		drda[1]=r[1]/a;
		drde[0]=-a-r[1]*r[1]/(rr*f);
		drde[1]=r[0]*r[1]/(rr*f);
		drdM[0]=v[0]/n;
		drdM[1]=v[1]/n;
		dvda[0]=-v[0]/(2*a);
		dvda[1]=-v[1]/(2*a);
		dvde[0]=v[0]*ar2*(2*r[0]/a+e/f*r[1]*r[1]/(a*a));
		dvde[1]=n/sf*ar2*(r[0]*r[0]/rr-r[1]*r[1]/(a*f));
		dvdM[0]=-n*ar3*r[0];
		dvdM[1]=-n*ar3*r[1];
		dPQdO[0]=-PQW[1];
		dPQdO[1]=PQW[0];
		dPQdO[2]=dPQdO[5]=0.0;
		dPQdO[3]=-PQW[4];
		dPQdO[4]=PQW[3];
		for (k=0;k<3;k++) {
			dPQdi[k]=sin(omg)*PQW[6+k];
			dPQdi[3+k]=cos(omg)*PQW[6+k];
			dPQdo[k]=PQW[3+k];
			dPQdo[3+k]=-PQW[k];
		}
		Mv2(PQW,drda,dsde);		/* ds/da */
		Mv2(PQW,dvda,dsde+3);
		Mv2(PQW,drde,dsde+6);	/* ds/de */
		Mv2(PQW,dvde,dsde+9);
		Mv2(dPQdi,r,dsde+12);	/* ds/di */
		Mv2(dPQdi,v,dsde+15);
		Mv2(dPQdO,r,dsde+18);	/* ds/d¶ */
		Mv2(dPQdO,v,dsde+21);
		Mv2(dPQdo,r,dsde+24);	/* ds/dÖ */
		Mv2(dPQdo,v,dsde+27);
		Mv2(PQW,drdM,dsde+30);	/* ds/dM */
		Mv2(PQW,dvdM,dsde+33);
	}
	return 1;
}
/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite state to orbit element
% [func]   : convert satellite position/velocity to orbit element
% [argin]  : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
% [argout] : ele   = orbit element(m,deg) [a,e,i,OMG,omg,M]
% [note]   : tangental orbit element
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
%----------------------------------------------------------------------------*/
/* satellite state to orbit element -----------------------------------------*/
extern void state2ele(const double *state, double *ele)
{
	int j;
	double r[3],v[3],h[3],w[3],hh,rr,aa,a,e,i,OMG,omg,M,E;
	
	r[0]=state[0]; r[1]=state[1]; r[2]=state[2];
	v[0]=state[3]; v[1]=state[4]; v[2]=state[5];
	rr=NORM(r);
	CROSS(r,v,h);
	aa=2.0/rr-DOT(v,v)/GME;
	if (NORM(h)<1E-12||aa<=0) {
		for (j=0;j<6;j++) ele[j]=NAN; /* no rotating orbit */
		return;
	}
	a=1.0/aa;
	hh=NORM(h);
	w[0]=h[0]/hh; w[1]=h[1]/hh; w[2]=h[2]/hh;
	i=atan(NORM2(w)/w[2]);
	OMG=atan2(h[0],-h[1]);
	e=sqrt(1.0-DOT(h,h)/(GME*a));
	E=atan2(DOT(r,v)/sqrt(GME*a),1.0-rr/a);
	M=E-e*sin(E);
	omg=atan2(r[2],-r[0]*w[1]+r[1]*w[0])-atan2(sqrt(1-e*e)*sin(E),cos(E)-e);
	if (i<0.0) i=i+M_PI;
	if (OMG<0.0) OMG=OMG+2.0*M_PI;
	if (omg<0.0) omg=omg+2.0*M_PI;
	if (M<0.0) M=M+2.0*M_PI;
	ele[0]=a;
	ele[1]=e;
	ele[2]=i*180.0/M_PI;
	ele[3]=OMG*180.0/M_PI;
	ele[4]=omg*180.0/M_PI;
	ele[5]=M*180.0/M_PI;
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