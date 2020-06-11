#include "CTIMES.h"
void ymd2doy(int Yr,int Mh,int Dy,int*Do){
   if(Yr%4!=0||(Yr%100==0&&Yr%400!=0))
      *Do=dayofmonth[0][Mh-1]+Dy;
   else
      *Do=dayofmonth[1][Mh-1]+Dy;
}
void ymd2mjd(int Yr,int Mh,int Dy,int h,int m,double s,double*Mjd){
  double Mjdf=(h+(m+s/60.)/60.)/24.,Mjd1,Mjd2;
  iauCal2jd(Yr, Mh, Dy,&Mjd1,&Mjd2);
  *Mjd=Mjd1+Mjd2+Mjdf-2400000.5;
}
void doy2ymd(int Yr,int Do,int* Mh,int* Dy){
   if(Yr%4!=0||(Yr%100==0&&Yr%400!=0))
   for(int i=1;i<13;i++)
   {if(dayofmonth[0][i]>=Do){*Mh=i;*Dy=Do-dayofmonth[0][i-1];break;}}
   else
   for(int i=1;i<13;i++)
   {if(dayofmonth[1][i]>=Do){*Mh=i;*Dy=Do-dayofmonth[1][i-1];break;}}
}
void doy2mjd(int Yr,int Do,int h,int m,double s,double*Mjd){
   int Mh,Dy;
   doy2ymd(Yr,Do,&Mh,&Dy);
   ymd2mjd(Yr,Mh,Dy,h,m,s,Mjd);
}
void mjd2ymd(double Mjd,int* Yr,int* Mh,int* Dy,int* h,int* m,double*s){
   double Mjdf;
   iauJd2cal(2400000.5,Mjd+sPresice/86400,Yr,Mh,Dy,&Mjdf);
   *h=(int)(fmod(Mjdf*24+sPresice/3600,24));
   *m=(int)(fmod(Mjdf*24*60+sPresice/60,60.0));
   *s=fmod(Mjdf*24*3600+sPresice,60.0);
}
void mjd2doy(double Mjd,int* Yr,int* Do,int* h,int* m,double*s){
   int Mh,Dy;
   mjd2ymd(Mjd,Yr,&Mh,&Dy,h,m,s);
   ymd2doy(*Yr,Mh,Dy,Do);
}

void gpswd2mjd(int Gw,int Gd,double Sd,double*Mjd){
    *Mjd=44244+Gw*7+Gd+Sd/86400.0;
}
void gpsws2mjd(int Gw,double Gs,double*Mjd){
     *Mjd=44244+Gw*7+Gs/86400.0;
}
void mjd2gpswd(double Mjd,int*Gw,int*Gd,double*Sd)
{
    double DurDay=Mjd-44244;
    *Gw=floor(DurDay/7.0+sPresice/86400/7.0);
    *Gd=floor(DurDay-*Gw*7.0+sPresice/86400);
    *Sd=fmod(DurDay+sPresice/86400,1.0)*86400.0;
}
void mjd2gpsws(double Mjd,int*Gw,double*Gs)
{
    double DurDay=Mjd-44244;
    *Gw=floor(DurDay/7.0+sPresice/86400/7.0);
    *Gs=fmod(DurDay,7.0)*86400.0;
}

void ymd2gpswd(int Yr,int Mh,int Dy,int h,int m,double s,int*Gw,int*Gd,double*Sd){
    double Mjd;
    ymd2mjd(Yr,Mh,Dy,h,m,s,&Mjd);
    mjd2gpswd(Mjd,Gw,Gd,Sd);
}
void ymd2gpsws(int Yr,int Mh,int Dy,int h,int m,double s,int*Gw,double*Gs){
    double Mjd;
    ymd2mjd(Yr,Mh,Dy,h,m,s,&Mjd);
    mjd2gpsws(Mjd,Gw,Gs);
}
void doy2gpswd(int Yr,int Do,int h,int m,double s,int*Gw,int*Gd,double*Sd){
    double Mjd;
    doy2mjd(Yr,Do,h,m,s,&Mjd);
    mjd2gpswd(Mjd,Gw,Gd,Sd);
}
void doy2gpsws(int Yr,int Do,int h,int m,double s,int*Gw,double*Gs){
    double Mjd;
    doy2mjd(Yr,Do,h,m,s,&Mjd);
    mjd2gpsws(Mjd,Gw,Gs);
}
void gpswd2ymd(int Gw,int Gd,double Sd,int* Yr,int* Mh,int* Dy,int* h,int* m,double*s){
    double Mjd;
    gpswd2mjd(Gw,Gd,Sd,&Mjd);
    mjd2ymd(Mjd,Yr,Mh,Dy,h,m,s);
}
void gpsws2ymd(int Gw,double Gs,int* Yr,int* Mh,int* Dy,int* h,int* m,double*s){
    double Mjd;
    gpsws2mjd(Gw,Gs,&Mjd);
    mjd2ymd(Mjd,Yr,Mh,Dy,h,m,s);
}
void gpswd2doy(int Gw,int Gd,double Sd,int* Yr,int* Do,int* h,int* m,double*s){
    double Mjd;
    gpswd2mjd(Gw,Gd,Sd,&Mjd);
    mjd2doy(Mjd,Yr,Do,h,m,s);
}
void gpsws2doy(int Gw,double Gs,int* Yr,int* Do,int* h,int* m,double*s){
    double Mjd;
    gpsws2mjd(Gw,Gs,&Mjd);
    mjd2doy(Mjd,Yr,Do,h,m,s);
}
void gpswd2gpsws(int Gd,double Sd,double*Gs){
   *Gs=Gd*86400.0+Sd;    
}
void gpsws2gpswd(double Gs,int*Gd,double*Sd){
   *Gd=floor(Gs/86400.0+sPresice/86400.0);
   *Sd=Gs-*Gd*86400.0;
}