
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

#define NM(n,m)		((n)+(m)*(nmax+1))		/* pnm,dpnm index */
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
/*
!-------------------------------------------------------------------------
!
! Writes TLE-records to an output-file.
! The TLE_Data-fields nDO20, nDDO60 and BSTAR are ignored.
!
! Input: toFile = name of outputfile
!    NumWrite = number of records to write
!    TLE_Data = holds 'NumWrite' TLE information records to be written.
!    Replace = optional; Replace existing outputfile? (default: append).
! Output: IOSTAT = optional; returns I/O status, 0 if no error occured.
!
SUBROUTINE Write_TLE_Rec(toFile, NumWrite, TLE_Data, &
              IOSTAT, Replace)
   CHARACTER (LEN=*), INTENT(IN) :: toFile
   INTEGER, INTENT(IN) :: NumWrite
   TYPE (TLE_Set), INTENT(IN) :: TLE_Data(MAXSat)
   INTEGER, OPTIONAL, INTENT(OUT) :: IOSTAT
   LOGICAL, OPTIONAL, INTENT(IN) :: Replace

   LOGICAL :: Opened, Existing
   INTEGER :: unit, ios, i
   CHARACTER (LEN=10) :: Status, Position

   Status='old'
   Position='append'
   IF (PRESENT(Replace)) THEN
      IF (Replace) Status = 'replace'
   ENDIF
   INQUIRE(FILE=TRIM(toFile), EXIST=Existing, OPENED=Opened)
   IF (.NOT. Existing) Status = 'replace'
   IF (.NOT. Opened) THEN
      PRINT *,'Opening file... '

      call getLun_s(unit)
      OPEN(unit, FILE=TRIM(toFile), IOSTAT=ios, ERR=99, STATUS=Status, &
            POSITION=Position, ACCESS='sequential', FORM='formatted', &
            BLANK='zero')
   ENDIF

   DO i = 1, NumWrite
      WRITE(unit, 52, IOSTAT=ios, ERR=99) TRIM(ADJUSTL(TLE_Data(i)%SatID))
      ! Attention: nDO20, nDDO60 and BSTAR are zero.
      WRITE(unit, 54, IOSTAT=ios, ERR=99) i, TLE_Data(i)%IntDesig, TLE_Data(i)%EPOCH
      WRITE(unit, 56, IOSTAT=ios, ERR=99) i, TLE_Data(i)%i0, &
           TLE_Data(i)%Node0, TLE_Data(i)%e0, TLE_Data(i)%omega0, &
           TLE_Data(i)%M0, TLE_Data(i)%n0
      WRITE (*,*) 'Processed.'
   ENDDO

   CLOSE(unit, IOSTAT=ios)
   call freeLun_s(unit)

   GOTO 99

52 FORMAT (A)
54 FORMAT ('1',1x,I5,2x,A8,1x,F14.8,2x,'.00000000',2(2x,'00000-0'),&
           1x,'0',2x,'0010')
56 FORMAT ('2',1x,I5,2(1X,F8.4),F8.7,0P,2(1x,F8.4),1X,F11.8,&
           '     ','0',T26,' ')

99 IF (PRESENT(IOSTAT)) IOSTAT = ios
END SUBROUTINE Write_TLE_Rec
!-------------------------------------------------------------------------
!
! TLEtoKepler converts two-line-element sets into Kepler-element-records (for
! use within KeplerOrbit or CircularOrbit).
!
! Input: TLE
! Output: Kepler Element set.
!
SUBROUTINE TLEtoKepler(TLE, KeplerElem)
   TYPE (TLE_Set), INTENT(INOUT) :: TLE
   TYPE (Kepler_Elements), INTENT(OUT) :: KeplerElem

   REAL( DP ) :: dummy

   ! calculate semimajor axis (out of TLE-set) if necessary:
   CALL ModifyTLESet(TLE)

   KeplerElem%New_Elements = .TRUE.
   KeplerElem%a = TLE%a
   KeplerElem%e = TLE%e0
   KeplerElem%i = TLE%i0 * RAD_TO_DEG
   KeplerElem%Node0 = TLE%Node0 * RAD_TO_DEG
   KeplerElem%omega0 = TLE%omega0 * RAD_TO_DEG
   KeplerElem%EPOCH = TLE%EPOCH
   CALL ModifyKeplerSet(KeplerElem)
   ! calculate epoch of perigee passage (mean anomaly = 0)
   dummy = TLE%M0/KeplerElem%n_pert
   KeplerElem%EPOCH = KeplerElem%EPOCH - dummy
   KeplerElem%Node0_rad = MOD(KeplerElem%Node0_rad - dummy * &
        KeplerElem%Node_dot, TWO_PI)
   IF (KeplerElem%Node0_rad < 0.0_DP) &
        KeplerElem%Node0_rad = KeplerElem%Node0_rad + TWO_PI
   KeplerElem%omega0_rad = MOD(KeplerElem%omega0_rad - dummy * &
        KeplerElem%omega_dot, TWO_PI)
   IF (KeplerElem%omega0_rad < 0.0_DP) &
        KeplerElem%omega0_rad = KeplerElem%omega0_rad + TWO_PI
   KeplerElem%Node0 = KeplerElem%Node0_rad * RAD_TO_DEG
   KeplerElem%omega0 = KeplerElem%omega0_rad * RAD_TO_DEG
END SUBROUTINE TLEtoKepler



!-------------------------------------------------------------------------
!
! KeplerToTLE converts a set of Kepler elements to a TLE-set.
!
! Input: Kepler = Kepler elements to be converted
! Output: TLE = two-line-elements corresponding to Kepler elements.
!
SUBROUTINE KeplertoTLE(Kepler, TLE)
   TYPE (Kepler_Elements), INTENT(INOUT) :: Kepler
   TYPE (TLE_Set), INTENT(OUT) :: TLE

   CALL ModifyKeplerSet(Kepler)
   TLE%M0 = 0.0_DP               ! "mean" mean anomaly at epoch
   TLE%Node0 = Kepler%Node0   ! "mean" longitude of ascending node at epoch
   TLE%omega0 = Kepler%omega0 ! "mean" argument of perigee at epoch
   TLE%e0 = Kepler%e          ! "mean" eccentricity at epoch
   TLE%i0 = Kepler%i          ! "mean" inclination at epoch
   TLE%nDO20 = 0.0_DP         ! 1/2 of the First Time Derivative of Mean
                              ! Motion at epoch
   TLE%nDDO60 = 0.0_DP        ! 1/6 of the Second Time Derivative of Mean
                              ! Motion at epoch
   TLE%BSTAR = 0.0_DP         ! the SGP4 type drag coefficient is zero
   TLE%EPOCH = Kepler%Epoch   ! the epoch (modified julian date)
   TLE%n0 = ModifiedMeanMotion(Kepler%n, TLE%i0, TLE%e0)/TWO_PI   ! SGP type
                              ! "mean" mean motion at epoch
   TLE%New_Elements = .TRUE.
END SUBROUTINE KeplertoTLE



!-------------------------------------------------------------------------
!
! ModifiedMeanMotion calculates the modified mean motion (saved in the
! TLE files and later used by routine SGPOrbit) out of the original
! mean motion.
! Input: n = original mean motion [radians/day]
!   i = inclination [degrees]
!   e = eccentricity
! Output: modified mean motion [radians/day]
!
REAL( DP ) FUNCTION ModifiedMeanMotion(n, i, e)
   REAL( DP ), INTENT(IN) :: n, i, e

   REAL( DP ) :: COSIO, DFAC, A2, D2, n1, A1, D1, AO, DO, n0

   COSIO = COS(i*DE2RA)
   DFAC = (3.0_DP*COSIO*COSIO - 1.0_DP)/(1.0_DP - e*e)**1.5_DP
   A2 = (XKE/(n/XMNPDA))**TOTHIRD
   D2 = CK2*1.5_DP/A2/A2*DFAC
   n1 = n * (1.0_DP + D2)               ! 1st guess for modified n
   A1 = (XKE/(n1/XMNPDA))**TOTHIRD
   D1 = CK2*1.5_DP/A1/A1*DFAC
   AO = A1*(1.0_DP - 1.0_DP/3.0_DP*D1 - D1*D1 - 134.0_DP/81.0_DP*D1*D1*D1)
   DO = CK2*1.5_DP/AO/AO*DFAC
   n0 = n * (1.0_DP + DO)               ! cf. [1] (KELSO report), p.21
   ModifiedMeanMotion = n0
END FUNCTION ModifiedMeanMotion

*/


/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite azimuth/elevation angle
% [func]   : calculate satellite azimuth/elevation angle
% [argin]  : spos = satellite postion [posx;posy;posz] (m) (ecef)
%            rpos = station position [posx;posy;posz] (m) (ecef)
% [argout] : azel = azimuth/elevation angle [az,el] (rad)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%            08/11/30   0.2  suppress warning
%-----------------------------------------------------------------------------*/
extern void satazel(const double *spos, const double *rpos, double *azel)
{
	double gpos[3],rs[3],rss[3],E[9];
	ecef2geod(rpos,gpos,E);
	rs[0]=spos[0]-rpos[0];
	rs[1]=spos[1]-rpos[1];
	rs[2]=spos[2]-rpos[2];
	Mv(E,rs,rss);
	azel[0]=atan2(rss[0],rss[1]);
	azel[1]=asin(rss[2]/NORM(rs));
}

#define DE2RA       0.174532925E-1
#define E6A         1.E-6
#define PIO2        1.57079633
#define QO          120.0
#define SO          78.0
#define TOTHRD      0.66666667
#define TWOPI       6.2831853
#define X3PIO2      4.71238898
#define XJ2         1.082616E-3
#define XJ3         -0.253881E-5
#define XJ4         -1.65597E-6
#define XKE         0.743669161E-1
#define XKMPER      6378.135
#define XMNPDA      1440.0
#define AE          1.0
#define CK2         5.413080E-4         /* = 0.5*XJ2*AE*AE */
#define CK4         0.62098875E-6       /* = -0.375*XJ4*AE*AE*AE*AE */
#define QOMS2T      1.88027916E-9       /* = pow((QO-SO)*AE/XKMPER,4.0) */
#define S           1.01222928          /* = AE*(1.0+SO/XKMPER) */

static void SGP4_STR3(double tsince, const /*tled_t *data*/double *data, double *rs)
{
    double xnodeo,omegao,xmo,eo,xincl,xno,xndt2o,xndd6o,bstar;
    double a1,cosio,theta2,x3thm1,eosq,betao2,betao,del1,ao,delo,xnodp,aodp,s4;
    double qoms24,perige,pinvsq,tsi,eta,etasq,eeta,psisq,coef,coef1,c1,c2,c3,c4;
    double c5,sinio,a3ovk2,x1mth2,theta4,xmdot,x1m5th,omgdot,xhdot1,xnodot;
    double omgcof,xmcof,xnodcf,t2cof,xlcof,aycof,delmo,sinmo,x7thm1,c1sq,d2,d3;
    double d4,t3cof,t4cof,t5cof,xmdf,omgadf,xnoddf,omega,xmp,tsq,xnode,delomg;
    double delm,tcube,tfour,a,e,xl,beta,xn,axn,xll,aynl,xlt,ayn,capu,sinepw;
    double cosepw,epw,ecose,esine,elsq,pl,r,rdot,rfdot,betal,cosu,sinu,u,sin2u;
    double cos2u,rk,uk,xnodek,xinck,rdotk,rfdotk,sinuk,cosuk,sinik,cosik,sinnok;
    double cosnok,xmx,xmy,ux,uy,uz,vx,vy,vz,x,y,z,xdot,ydot,zdot;
    double temp,temp1,temp2,temp3,temp4,temp5,temp6,tempa,tempe,templ;
    int i,isimp;
    
    xnodeo=data[0];//data->OMG*DE2RA;
    omegao=data[1];//data->omg*DE2RA;
    xmo=data[2];//data->M*DE2RA;
    xincl=data[3];//data->inc*DE2RA;
    temp=TWOPI/XMNPDA/XMNPDA;
    xno=data[4];//data->n*temp*XMNPDA;
    xndt2o=data[5];//data->ndot*temp;
    xndd6o=data[6];//data->nddot*temp/XMNPDA;
    bstar=data[7];//data->bstar/AE;
    eo=data[8];//data->ecc;
    /*
    * recover original mean motion (xnodp) and semimajor axis (aodp)
    * from input elements
    */
    a1=pow(XKE/xno,TOTHRD);
    cosio=cos(xincl);
    theta2=cosio*cosio;
    x3thm1=3.0*theta2-1.0;
    eosq=eo*eo;
    betao2=1.0-eosq;
    betao=sqrt(betao2);
    del1=1.5*CK2*x3thm1/(a1*a1*betao*betao2);
    ao=a1*(1.0-del1*(0.5*TOTHRD+del1*(1.0+134.0/81.0*del1)));
    delo=1.5*CK2*x3thm1/(ao*ao*betao*betao2);
    xnodp=xno/(1.0+delo);
    aodp=ao/(1.0-delo);
    /*
    * initialization
    * for perigee less than 220 kilometers, the isimp flag is set and
    * the equations are truncated to linear variation in sqrt a and
    * quadratic variation in mean anomaly. also, the c3 term, the
    * delta omega term, and the delta m term are dropped.
    */
    isimp=0;
    if ((aodp*(1.0-eo)/AE)<(220.0/XKMPER+AE)) isimp=1;
    
    /* for perigee below 156 km, the values of s and qoms2t are altered */
    s4=S;
    qoms24=QOMS2T;
    perige=(aodp*(1.0-eo)-AE)*XKMPER;
    if (perige<156.0) {
        s4=perige-78.0;
        if (perige<=98.0) s4=20.0;
        qoms24=pow((120.0-s4)*AE/XKMPER,4.0);
        s4=s4/XKMPER+AE;
    }
    pinvsq=1.0/(aodp*aodp*betao2*betao2);
    tsi=1.0/(aodp-s4);
    eta=aodp*eo*tsi;
    etasq=eta*eta;
    eeta=eo*eta;
    psisq=fabs(1.0-etasq);
    coef=qoms24*pow(tsi,4.0);
    coef1=coef/pow(psisq,3.5);
    c2=coef1*xnodp*(aodp*(1.0+1.5*etasq+eeta*(4.0+etasq))+0.75*
            CK2*tsi/psisq*x3thm1*(8.0+3.0*etasq*(8.0+etasq)));
    c1=bstar*c2;
    sinio=sin(xincl);
    a3ovk2=-XJ3/CK2*pow(AE,3.0);
    c3=coef*tsi*a3ovk2*xnodp*AE*sinio/eo;
    x1mth2=1.0-theta2;
    c4=2.0*xnodp*coef1*aodp*betao2*(eta*
            (2.0+0.5*etasq)+eo*(0.5+2.0*etasq)-2.0*CK2*tsi/
            (aodp*psisq)*(-3.0*x3thm1*(1.0-2.0*eeta+etasq*
            (1.5-0.5*eeta))+0.75*x1mth2*(2.0*etasq-eeta*
            (1.0+etasq))*cos(2.0*omegao)));
    c5=2.0*coef1*aodp*betao2*(1.0+2.75*(etasq+eeta)+eeta*etasq);
    theta4=theta2*theta2;
    temp1=3.0*CK2*pinvsq*xnodp;
    temp2=temp1*CK2*pinvsq;
    temp3=1.25*CK4*pinvsq*pinvsq*xnodp;
    xmdot=xnodp+0.5*temp1*betao*x3thm1+0.0625*temp2*betao*
            (13.0-78.0*theta2+137.0*theta4);
    x1m5th=1.0-5.0*theta2;
    omgdot=-0.5*temp1*x1m5th+0.0625*temp2*(7.0-114.0*theta2+
            395.0*theta4)+temp3*(3.0-36.0*theta2+49.0*theta4);
    xhdot1=-temp1*cosio;
    xnodot=xhdot1+(0.5*temp2*(4.0-19.0*theta2)+2.0*temp3*(3.0-
            7.0*theta2))*cosio;
    omgcof=bstar*c3*cos(omegao);
    xmcof=-TOTHRD*coef*bstar*AE/eeta;
    xnodcf=3.5*betao2*xhdot1*c1;
    t2cof=1.5*c1;
    xlcof=0.125*a3ovk2*sinio*(3.0+5.0*cosio)/(1.0+cosio);
    aycof=0.25*a3ovk2*sinio;
    delmo=pow(1.0+eta*cos(xmo),3.0);
    sinmo=sin(xmo);
    x7thm1=7.0*theta2-1.0;
    
    if (isimp!=1) {
        c1sq=c1*c1;
        d2=4.0*aodp*tsi*c1sq;
        temp=d2*tsi*c1/3.0;
        d3=(17.0*aodp+s4)*temp;
        d4=0.5*temp*aodp*tsi*(221.0*aodp+31.0*s4)*c1;
        t3cof=d2+2.0*c1sq;
        t4cof=0.25*(3.0*d3+c1*(12.0*d2+10.0*c1sq));
        t5cof=0.2*(3.0*d4+12.0*c1*d3+6.0*d2*d2+15.0*c1sq*(2.0*d2+c1sq));
    }
    else {
        d2=d3=d4=t3cof=t4cof=t5cof=0.0;
    }
    /* update for secular gravity and atmospheric drag */
    xmdf=xmo+xmdot*tsince;
    omgadf=omegao+omgdot*tsince;
    xnoddf=xnodeo+xnodot*tsince;
    omega=omgadf;
    xmp=xmdf;
    tsq=tsince*tsince;
    xnode=xnoddf+xnodcf*tsq;
    tempa=1.0-c1*tsince;
    tempe=bstar*c4*tsince;
    templ=t2cof*tsq;
    if (isimp==1) {
        delomg=omgcof*tsince;
        delm=xmcof*(pow(1.0+eta*cos(xmdf),3.0)-delmo);
        temp=delomg+delm;
        xmp=xmdf+temp;
        omega=omgadf-temp;
        tcube=tsq*tsince;
        tfour=tsince*tcube;
        tempa=tempa-d2*tsq-d3*tcube-d4*tfour;
        tempe=tempe+bstar*c5*(sin(xmp)-sinmo);
        templ=templ+t3cof*tcube+tfour*(t4cof+tsince*t5cof);
    }
    a=aodp*pow(tempa,2.0);
    e=eo-tempe;
    xl=xmp+omega+xnode+xnodp*templ;
    beta=sqrt(1.0-e*e);
    xn=XKE/pow(a,1.5);
    
    /* long period periodics */
    axn=e*cos(omega);
    temp=1.0/(a*beta*beta);
    xll=temp*xlcof*axn;
    aynl=temp*aycof;
    xlt=xl+xll;
    ayn=e*sin(omega)+aynl;
    
    /* solve keplers equation */
    capu=fmod(xlt-xnode,TWOPI);
    temp2=capu;
    for (i=0;i<10;i++) {
        sinepw=sin(temp2);
        cosepw=cos(temp2);
        temp3=axn*sinepw;
        temp4=ayn*cosepw;
        temp5=axn*cosepw;
        temp6=ayn*sinepw;
        epw=(capu-temp4+temp3-temp2)/(1.0-temp5-temp6)+temp2;
        if (fabs(epw-temp2)<=E6A) break;
        temp2=epw;
    }
    /* short period preliminary quantities */
    ecose=temp5+temp6;
    esine=temp3-temp4;
    elsq=axn*axn+ayn*ayn;
    temp=1.0-elsq;
    pl=a*temp;
    r=a*(1.0-ecose);
    temp1=1.0/r;
    rdot=XKE*sqrt(a)*esine*temp1;
    rfdot=XKE*sqrt(pl)*temp1;
    temp2=a*temp1;
    betal=sqrt(temp);
    temp3=1.0/(1.0+betal);
    cosu=temp2*(cosepw-axn+ayn*esine*temp3);
    sinu=temp2*(sinepw-ayn-axn*esine*temp3);
    u=atan2(sinu,cosu);
    sin2u=2.0*sinu*cosu;
    cos2u=2.0*cosu*cosu-1.0;
    temp=1.0/pl;
    temp1=CK2*temp;
    temp2=temp1*temp;
    
    /* update for short periodics */
    rk=r*(1.0-1.5*temp2*betal*x3thm1)+0.5*temp1*x1mth2*cos2u;
    uk=u-0.25*temp2*x7thm1*sin2u;
    xnodek=xnode+1.5*temp2*cosio*sin2u;
    xinck=xincl+1.5*temp2*cosio*sinio*cos2u;
    rdotk=rdot-xn*temp1*x1mth2*sin2u;
    rfdotk=rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);
    
    /* orientation vectors */
    sinuk=sin(uk);
    cosuk=cos(uk);
    sinik=sin(xinck);
    cosik=cos(xinck);
    sinnok=sin(xnodek);
    cosnok=cos(xnodek);
    xmx=-sinnok*cosik;
    xmy=cosnok*cosik;
    ux=xmx*sinuk+cosnok*cosuk;
    uy=xmy*sinuk+sinnok*cosuk;
    uz=sinik*sinuk;
    vx=xmx*cosuk-cosnok*sinuk;
    vy=xmy*cosuk-sinnok*sinuk;
    vz=sinik*cosuk;
    
    /* position and velocity */
    x=rk*ux;
    y=rk*uy;
    z=rk*uz;
    xdot=rdotk*ux+rfdotk*vx;
    ydot=rdotk*uy+rfdotk*vy;
    zdot=rdotk*uz+rfdotk*vz;
    
    rs[0]=x*XKMPER/AE*1E3; /* (m) */
    rs[1]=y*XKMPER/AE*1E3;
    rs[2]=z*XKMPER/AE*1E3;
    rs[3]=xdot*XKMPER/AE*XMNPDA/86400.0*1E3; /* (m/s) */
    rs[4]=ydot*XKMPER/AE*XMNPDA/86400.0*1E3;
    rs[5]=zdot*XKMPER/AE*XMNPDA/86400.0*1E3;
}
/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : spheric harmonic functions
% [func]   : spheric harmonic functions by azimuth/elevation angle
% [argin]  : az,el = azimuth/elevation angle (rad)
%            nmax  = max degrees of spheric harmonic functions
% [argout] : fc,fs = spheric harmonic functions
%                fc(n+1,m+1) = Pnm(-cos(2*el))*cos(m*az)
%                fs(n+1,m+1) = Pnm(-cos(2*el))*sin(m*az)
%                (Pnm = normalized Legendre polynomial)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/18  0.1  new
%-----------------------------------------------------------------------------*/

/* factorial -----------------------------------------------------------------*/
static double factorial(int n)
{
	double f; for (f=1.0;n>1;n--) f*=(double)n; return f;
}
/* normalized legendre functions ---------------------------------------------*/
static void LegendreFunc(double x, int nmax, double *p)
{
	int n,m;
	
	p[NM(0,0)]=1.0;
	p[NM(1,0)]=x;
	p[NM(1,1)]=sqrt(1.0-x*x);
	
	for (n=2;n<=nmax;n++) {
	    p[NM(n,n)]=(2.0*n-1.0)*sqrt(1.0-x*x)*p[NM(n-1,n-1)];
	    for (m=0;m<n;m++) {
	        p[NM(n,m)]=(x*(2.0*n-1.0)*p[NM(n-1,m)]-(n+m-1.0)*p[NM(n-2,m)])/(n-m);
	    }
	}
	for (n=1;n<=nmax;n++) {
	    p[NM(n,0)]*=sqrt(2.0*n+1.0);
	    for (m=1;m<=n;m++)
	    	p[NM(n,m)]*=sqrt(factorial(n-m)*(4.0*n+2.0)/factorial(n+m));
	}
}
/* spheric harmonic functions ------------------------------------------------*/
extern void sphfunc_azel(double az, double el, int nmax, double *fc, double *fs)
{
	double *p,*cosm,*sinm;
	int n,m;
	
	p=ZMAT(nmax+1,nmax+1); cosm=VEC(nmax+1); sinm=VEC(nmax+1);
	if (p==NULL||cosm==NULL||sinm==NULL) return;
	
	for (m=0;m<=nmax;m++) {
		cosm[m]=cos((double)m*az);
		sinm[m]=sin((double)m*az);
	}
	LegendreFunc(-cos(2.0*el),nmax,p);
	
	fc[0]=1.0;
	for (n=1;n<=nmax;n++)
	for (m=0;m<=n;m++) {
		fc[NM(n,m)]=p[NM(n,m)]*cosm[m];
		fs[NM(n,m)]=p[NM(n,m)]*sinm[m];
	}
	FreeMat(p); FreeMat(cosm); FreeMat(sinm);
}
extern void sphfunc_latlon(double lat, double lon, int nmax, double *fc, double *fs)
{
	double *p,*cosm,*sinm;
	int n,m;
	
	p=ZMAT(nmax+1,nmax+1); cosm=VEC(nmax+1); sinm=VEC(nmax+1);
	if (p==NULL||cosm==NULL||sinm==NULL) return;
	
	for (m=0;m<=nmax;m++) {
		cosm[m]=cos((double)m*lon);
		sinm[m]=sin((double)m*lon);
	}
	LegendreFunc(-sin(lat),nmax,p);
	
	fc[0]=1.0;
	for (n=1;n<=nmax;n++)
	for (m=0;m<=n;m++) {
		fc[NM(n,m)]=p[NM(n,m)]*cosm[m];
		fs[NM(n,m)]=p[NM(n,m)]*sinm[m];
	}
	FreeMat(p); FreeMat(cosm); FreeMat(sinm);
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