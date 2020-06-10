#ifndef __TIMES_H__
#define __TIMES_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define  sPresice 0.00001
static int dayofmonth[2][13]={{0,31,59,90,120,151,181,212,243,273,304,334,365},{0,31,60,91,121,152,182,213,244,274,305,335,366}};
int iauJd2cal(double dj1, double dj2,
              int *iy, int *im, int *id, double *fd);
int iauCal2jd(int iy, int im, int id, double *djm0, double *djm);
void ymd2doy(int,int,int,int*);
void ymd2mjd(int,int,int,int,int,double,double*);
void doy2ymd(int,int,int*,int*);
void doy2mjd(int,int,int,int,double,double*);
void mjd2ymd(double,int*,int*,int*,int*,int*,double*);
void mjd2doy(double,int*,int*,int*,int*,double*);
void gpswd2mjd(int,int,double,double*);
void gpsws2mjd(int,double,double*);
void mjd2gpswd(double,int*,int*,double*);
void mjd2gpsws(double ,int*,double*);
void ymd2gpswd(int ,int ,int ,int ,int ,double ,int*,int*,double*);//14
void ymd2gpsws(int ,int ,int ,int ,int ,double ,int*,double*);//15
void doy2gpswd(int ,int ,int ,int ,double ,int*,int*,double*);//24
void doy2gpsws(int ,int ,int ,int ,double ,int*,double*);//25
void gpswd2ymd(int ,int ,double ,int* ,int* ,int* ,int* ,int*,double*);//41
void gpsws2ymd(int ,double ,int* ,int* ,int* ,int* ,int* ,double*);//51
void gpswd2doy(int ,int ,double ,int* ,int* ,int* ,int* ,double*);//42
void gpsws2doy(int ,double ,int* ,int* ,int* ,int* ,double*);//52
void gpswd2gpsws(int ,double ,double*);//45
void gpsws2gpswd(double ,int*,double*);//54
#endif