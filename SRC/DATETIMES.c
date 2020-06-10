
//C source of Mapoet Niphy 
/*
Author    :  Mapoet Niphy
Date      :  2020
Institude :  SHAO
Name      :  DATETIMES
Project   :  GEOTOOL

  Created by Mapoet Niphy on 2020/06/10.
  Copyright © 2020年 Mapoet Niphy. All rights reserved.

*/
#include "CTIMES.h"

int main(int argc,char**argv){
// YMD->1
// DOY->2
// MJD->3 
// GWD->4
// GWS->5
// big uint:0
// smail unit:1
// add->0
// sub->9
// dif->8
   int Yr,Mh,Dy,Do,h,m,id,Gw,Gd,Yr2,Mh2,Dy2,Do2,h2,m2,Gw2,Gd2;
   double    s,Mjd,Sd,Gs,s2,Mjd2,Sd2,Gs2,d;
   id=atoi(argv[1]);
   while(!feof(stdin))switch(id)
   {
      case 120:if(fscanf(stdin,"%d%d%d\n",&Yr,&Mh,&Dy)!=3)continue;
             ymd2doy(Yr,Mh,Dy,&Do);fprintf(stdout,"%4.4d%4.3d\n",Yr,Do);break;
      case 121:if(fscanf(stdin,"%d%d%d%d%d%lf\n",&Yr,&Mh,&Dy,&h,&m,&s)!=6)continue;
             ymd2doy(Yr,Mh,Dy,&Do);fprintf(stdout,"%4.4d%4.3d%3.2d%3.2d%6.2lf\n",Yr,Do,h,m,s);break;
      case 130:if(fscanf(stdin,"%d%d%d\n",&Yr,&Mh,&Dy)!=3)continue;
             ymd2mjd(Yr,Mh,Dy,0,0,0,&Mjd);fprintf(stdout,"%20.12lf\n",Mjd);break;
      case 131:if(fscanf(stdin,"%d%d%d%d%d%lf\n",&Yr,&Mh,&Dy,&h,&m,&s)!=6)continue;
             ymd2mjd(Yr,Mh,Dy,h,m,s,&Mjd);fprintf(stdout,"%20.12lf\n",Mjd);break;
      case 140:if(fscanf(stdin,"%d%d%d\n",&Yr,&Mh,&Dy)!=3)continue;
             ymd2gpswd(Yr,Mh,Dy,0,0,0,&Gw,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d\n",Gw,Gd);break;
      case 141:if(fscanf(stdin,"%d%d%d%d%d%lf\n",&Yr,&Mh,&Dy,&h,&m,&s)!=6)continue;
             ymd2gpswd(Yr,Mh,Dy,h,m,s,&Gw,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d%12.4lf\n",Gw,Gd,Sd);break;
      case 150:if(fscanf(stdin,"%d%d%d\n",&Yr,&Mh,&Dy)!=3)continue;
             ymd2gpsws(Yr,Mh,Dy,0,0,0,&Gw,&Gs);fprintf(stdout,"%5.4d\n",Gw);break;
      case 151:if(fscanf(stdin,"%d%d%d%d%d%lf\n",&Yr,&Mh,&Dy,&h,&m,&s)!=6)continue;
             ymd2gpsws(Yr,Mh,Dy,h,m,s,&Gw,&Gs);fprintf(stdout,"%5.4d%12.4lf\n",Gw,Gs);break;
      case 210:if(fscanf(stdin,"%d%d\n",&Yr,&Do)!=2)continue;
             doy2ymd(Yr,Do,&Mh,&Dy);fprintf(stdout,"%4.4d%3.2d%3.2d\n",Yr,Mh,Dy);break;
      case 211:if(fscanf(stdin,"%d%d%d%d%lf\n",&Yr,&Do,&h,&m,&s)!=5)continue;
             doy2ymd(Yr,Do,&Mh,&Dy);fprintf(stdout,"%4.4d%3.2d%3.2d%3.2d%3.2d%6.2lf\n",Yr,Mh,Dy,h,m,s);break;
      case 230:if(fscanf(stdin,"%d%d\n",&Yr,&Do)!=2)continue;
             doy2mjd(Yr,Do,0,0,0,&Mjd);fprintf(stdout,"%20.12lf\n",Mjd);break;
      case 231:if(fscanf(stdin,"%d%d%d%d%lf\n",&Yr,&Do,&h,&m,&s)!=5)continue;
             doy2mjd(Yr,Do,h,m,s,&Mjd);fprintf(stdout,"%20.12lf\n",Mjd);break;
      case 240:if(fscanf(stdin,"%d%d\n",&Yr,&Do)!=2)continue;
             doy2gpswd(Yr,Do,0,0,0,&Gw,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d\n",Gw,Gd);break;
      case 241:if(fscanf(stdin,"%d%d%d%d%lf\n",&Yr,&Do,&h,&m,&s)!=5)continue;
             doy2gpswd(Yr,Do,h,m,s,&Gw,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d%12.4lf\n",Gw,Gd,Sd);break;
      case 250:if(fscanf(stdin,"%d%d\n",&Yr,&Do)!=2)continue;
             doy2gpsws(Yr,Do,0,0,0,&Gw,&Gs);fprintf(stdout,"%5.4d\n",Gw);break;
      case 251:if(fscanf(stdin,"%d%d%d%d%lf\n",&Yr,&Do,&h,&m,&s)!=5)continue;
             doy2gpsws(Yr,Do,h,m,s,&Gw,&Gs);fprintf(stdout,"%5.4d%12.4lf\n",Gw,Gs);break;
      case 310:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2ymd(Mjd,&Yr,&Mh,&Dy,&h,&m,&s);fprintf(stdout,"%4.4d%3.2d%3.2d\n",Yr,Mh,Dy);break;
      case 311:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2ymd(Mjd,&Yr,&Mh,&Dy,&h,&m,&s);fprintf(stdout,"%4.4d%3.2d%3.2d%3.2d%3.2d%6.2lf\n",Yr,Mh,Dy,h,m,s);break;
      case 320:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2doy(Mjd,&Yr,&Do,&h,&m,&s);fprintf(stdout,"%4.4d%4.3d\n",Yr,Do);break;
      case 321:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2doy(Mjd,&Yr,&Do,&h,&m,&s);fprintf(stdout,"%4.4d%4.3d%3.2d%3.2d%6.2lf\n",Yr,Do,h,m,s);break;
      case 340:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2gpswd(Mjd,&Gw,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d\n",Gw,Gd);break;
      case 341:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2gpswd(Mjd,&Gw,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d%12.4lf\n",Gw,Gd,Sd);break;
      case 350:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2gpsws(Mjd,&Gw,&Gs);fprintf(stdout,"%5.4d\n",Gw);break;
      case 351:if(fscanf(stdin,"%lf\n",&Mjd)!=1)continue;
             mjd2gpsws(Mjd,&Gw,&Gs);fprintf(stdout,"%5.4d%12.4lf\n",Gw,Gs);break;
      case 410:if(fscanf(stdin,"%d%d",&Gw,&Gd)!=2)continue;
             gpswd2ymd(Gw,Gd,0,&Yr,&Mh,&Dy,&h,&m,&s);fprintf(stdout,"%4.4d%3.2d%3.2d\n",Yr,Mh,Dy);break;
      case 411:if(fscanf(stdin,"%d%d%lf",&Gw,&Gd,&Sd)!=3)continue;
             gpswd2ymd(Gw,Gd,Sd,&Yr,&Mh,&Dy,&h,&m,&s);fprintf(stdout,"%4.4d%3.2d%3.2d%3.2d%3.2d%6.2lf\n",Yr,Mh,Dy,h,m,s);break;
      case 420:if(fscanf(stdin,"%d%d",&Gw,&Gd)!=2)continue;
             gpswd2doy(Gw,Gd,0,&Yr,&Do,&h,&m,&s);fprintf(stdout,"%4.4d%4.3d\n",Yr,Do);break;
      case 421:if(fscanf(stdin,"%d%d%lf",&Gw,&Gd,&Sd)!=3)continue;
             gpswd2doy(Gw,Gd,Sd,&Yr,&Do,&h,&m,&s);fprintf(stdout,"%4.4d%4.3d%3.2d%3.2d%6.2lf\n",Yr,Do,h,m,s);break;
      case 430:if(fscanf(stdin,"%d%d",&Gw,&Gd)!=2)continue;
             gpswd2mjd(Gw,Gd,0,&Mjd);fprintf(stdout,"%20.11lf\n",Mjd);break;
      case 431:if(fscanf(stdin,"%d%d%lf",&Gw,&Gd,&Sd)!=3)continue;
             gpswd2mjd(Gw,Gd,Sd,&Mjd);fprintf(stdout,"%20.11lf\n",Mjd);break;
      case 450:if(fscanf(stdin,"%d%d",&Gw,&Gd)!=2)continue;
             gpswd2gpsws(Gd,0,&Gs);fprintf(stdout,"%5.4d%12.4lf\n",Gw,Gs);break;
      case 451:if(fscanf(stdin,"%d%d%lf",&Gw,&Gd,&Sd)!=3)continue;
             gpswd2gpsws(Gd,Sd,&Gs);fprintf(stdout,"%5.4d%12.4lf\n",Gw,Gs);break;
      case 510:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2ymd(Gw,Gs,&Yr,&Mh,&Dy,&h,&m,&s);fprintf(stdout,"%4.4d%3.2d%3.2d\n",Yr,Mh,Dy);break;
      case 511:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2ymd(Gw,Gs,&Yr,&Mh,&Dy,&h,&m,&s);fprintf(stdout,"%4.4d%3.2d%3.2d%3.2d%3.2d%6.2lf\n",Yr,Mh,Dy,h,m,s);break;
      case 520:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2doy(Gw,Gs,&Yr,&Do,&h,&m,&s);fprintf(stdout,"%4.4d%4.3d\n",Yr,Do);break;
      case 521:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2doy(Gw,Gs,&Yr,&Do,&h,&m,&s);fprintf(stdout,"%4.4d%4.3d%3.2d%3.2d%6.2lf\n",Yr,Do,h,m,s);break;
      case 530:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2mjd(Gw,Gs,&Mjd);fprintf(stdout,"%20.11lf\n",Mjd);break;
      case 531:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2mjd(Gw,Gs,&Mjd);fprintf(stdout,"%20.11lf\n",Mjd);break;
      case 540:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2gpswd(Gw,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d\n",Gw,Gd);break;
      case 541:if(fscanf(stdin,"%d%lf",&Gw,&Gs)!=2)continue;
             gpsws2gpswd(Gs,&Gd,&Sd);fprintf(stdout,"%5.4d%2.1d%12.4lf\n",Gw,Gd,Sd);break;
      case 100:if(fscanf(stdin,"%d%d%d%lf\n",&Yr,&Mh,&Dy,&d)!=4)continue;
             ymd2mjd(Yr,Mh,Dy,0,0,0,&Mjd);mjd2ymd(Mjd+d,&Yr2,&Mh2,&Dy2,&h2,&m2,&s2);
             fprintf(stdout,"%4.4d%3.2d%3.2d\n",Yr2,Mh2,Dy2);break;
      case 190:if(fscanf(stdin,"%d%d%d%lf\n",&Yr,&Mh,&Dy,&d)!=4)continue;
             ymd2mjd(Yr,Mh,Dy,0,0,0,&Mjd);mjd2ymd(Mjd-d,&Yr2,&Mh2,&Dy2,&h2,&m2,&s2);
             fprintf(stdout,"%4.4d%3.2d%3.2d\n",Yr2,Mh2,Dy2);break;
      case 180:if(fscanf(stdin,"%d%d%d%d%d%d\n",&Yr,&Mh,&Dy,&Yr2,&Mh2,&Dy2)!=6)continue;
             ymd2mjd(Yr,Mh,Dy,0,0,0,&Mjd);ymd2mjd(Yr2,Mh2,Dy2,0,0,0,&Mjd2);
             fprintf(stdout,"%20.12lf\n",Mjd-Mjd2);break;
      case 101:if(fscanf(stdin,"%d%d%d%d%d%lf%lf\n",&Yr,&Mh,&Dy,&h,&m,&s,&d)!=7)continue;
             ymd2mjd(Yr,Mh,Dy,h,m,s,&Mjd);mjd2ymd(Mjd+d/86400.0,&Yr2,&Mh2,&Dy2,&h2,&m2,&s2);
             fprintf(stdout,"%4.4d%3.2d%3.2d%3.2d%3.2d%6.2lf\n",Yr2,Mh2,Dy2,h2,m2,s2);break;
      case 191:if(fscanf(stdin,"%d%d%d%d%d%lf%lf\n",&Yr,&Mh,&Dy,&h,&m,&s,&d)!=7)continue;
             ymd2mjd(Yr,Mh,Dy,h,m,s,&Mjd);mjd2ymd(Mjd-d/86400.0,&Yr2,&Mh2,&Dy2,&h2,&m2,&s2);
             fprintf(stdout,"%4.4d%3.2d%3.2d%3.2d%3.2d%6.2lf\n",Yr2,Mh2,Dy2,h2,m2,s2);break;
      case 181:if(fscanf(stdin,"%d%d%d%d%d%lf%d%d%d%d%d%lf\n",&Yr,&Mh,&Dy,&h,&m,&s,&Yr2,&Mh2,&Dy2,&h2,&m2,&s2)!=12)continue;
             ymd2mjd(Yr,Mh,Dy,h,m,s,&Mjd);ymd2mjd(Yr2,Mh2,Dy2,h2,m2,s2,&Mjd2);
             fprintf(stdout,"%20.12lf\n",(Mjd-Mjd2)*86400.0);break;
      case 200:if(fscanf(stdin,"%d%d%lf\n",&Yr,&Do,&d)!=3)continue;
             doy2mjd(Yr,Do,0,0,0,&Mjd);mjd2doy(Mjd+d,&Yr2,&Do2,&h2,&m2,&s2);
            fprintf(stdout,"%4.4d%4.3d\n",Yr2,Do2);break;
      case 290:if(fscanf(stdin,"%d%d%lf\n",&Yr,&Do,&d)!=3)continue;
             doy2mjd(Yr,Do,0,0,0,&Mjd);mjd2doy(Mjd-d,&Yr2,&Do2,&h2,&m2,&s2);
            fprintf(stdout,"%4.4d%4.3d\n",Yr2,Do2);break;
      case 280:if(fscanf(stdin,"%d%d%d%d\n",&Yr,&Do,&Yr2,&Do2)!=4)continue;
             doy2mjd(Yr,Do,0,0,0,&Mjd);doy2mjd(Yr2,Do2,0,0,0,&Mjd2);
             fprintf(stdout,"%20.12lf\n",Mjd-Mjd2);break;
      case 201:if(fscanf(stdin,"%d%d%d%d%lf%lf\n",&Yr,&Do,&h,&m,&s,&d)!=6)continue;
             doy2mjd(Yr,Do,h,m,s,&Mjd);mjd2doy(Mjd+d/86400.0,&Yr2,&Do2,&h2,&m2,&s2);
             fprintf(stdout,"%4.4d%4.3d%3.2d%3.2d%6.2lf\n",Yr2,Do2,h2,m2,s2);break;
      case 291:if(fscanf(stdin,"%d%d%d%d%lf%lf\n",&Yr,&Do,&h,&m,&s,&d)!=6)continue;
             doy2mjd(Yr,Do,h,m,s,&Mjd);mjd2doy(Mjd-d/86400.0,&Yr2,&Do2,&h2,&m2,&s2);
             fprintf(stdout,"%4.4d%4.3d%3.2d%3.2d%6.2lf\n",Yr2,Do2,h2,m2,s2);break;
      case 281:if(fscanf(stdin,"%d%d%d%d%lf%d%d%d%d%lf\n",&Yr,&Do,&h,&m,&s,&Yr2,&Do2,&h2,&m2,&s2)!=10)continue;
             doy2mjd(Yr,Do,h,m,s,&Mjd);doy2mjd(Yr2,Do2,h2,m2,s2,&Mjd2);
             fprintf(stdout,"%20.12lf\n",(Mjd-Mjd2)*86400.0);break;
      case 300:if(fscanf(stdin,"%lf%lf\n",&Mjd,&d)!=2)continue;
             fprintf(stdout,"%20.12lf\n",Mjd+d);break;
      case 390:if(fscanf(stdin,"%lf%lf\n",&Mjd,&d)!=2)continue;
             fprintf(stdout,"%20.12lf\n",Mjd-d);break;
      case 380:if(fscanf(stdin,"%lf%lf\n",&Mjd,&Mjd2)!=2)continue;
             fprintf(stdout,"%20.12lf\n",Mjd-Mjd2);break;
      case 301:if(fscanf(stdin,"%lf%lf\n",&Mjd,&d)!=2)continue;
             fprintf(stdout,"%20.12lf\n",Mjd+d/86400.0);break;
      case 391:if(fscanf(stdin,"%lf%lf\n",&Mjd,&d)!=2)continue;
             fprintf(stdout,"%20.12lf\n",Mjd-d/86400.0);break;
      case 381:if(fscanf(stdin,"%lf%lf\n",&Mjd,&Mjd2)!=2)continue;
             fprintf(stdout,"%20.12lf\n",(Mjd-Mjd2)*86400.0);break;
      case 400:if(fscanf(stdin,"%d%d%lf",&Gw,&Gd,&d)!=3)continue;
             gpswd2mjd(Gw,Gd,0,&Mjd);mjd2gpswd(Mjd+d,&Gw2,&Gd2,&Sd2);
             fprintf(stdout,"%5.4d%2.1d\n",Gw2,Gd2);break;
      case 490:if(fscanf(stdin,"%d%d%lf",&Gw,&Gd,&d)!=3)continue;
             gpswd2mjd(Gw,Gd,0,&Mjd);mjd2gpswd(Mjd-d,&Gw2,&Gd2,&Sd2);
             fprintf(stdout,"%5.4d%2.1d\n",Gw2,Gd2);break;
      case 480:if(fscanf(stdin,"%d%d%d%d",&Gw,&Gd,&Gw2,&Gd2)!=4)continue;
             gpswd2mjd(Gw,Gd,0,&Mjd);gpswd2mjd(Gw2,Gd2,0,&Mjd2);
             fprintf(stdout,"%20.12lf\n",Mjd-Mjd2);break;
      case 401:if(fscanf(stdin,"%d%d%lf%lf",&Gw,&Gd,&Sd,&d)!=4)continue;
             gpswd2mjd(Gw,Gd,Sd,&Mjd);mjd2gpswd(Mjd+d/86400.0,&Gw2,&Gd2,&Sd2);
             fprintf(stdout,"%5.4d%2.1d%12.4lf\n",Gw2,Gd2,Sd2);break;
      case 491:if(fscanf(stdin,"%d%d%lf%lf",&Gw,&Gd,&Sd,&d)!=4)continue;
             gpswd2mjd(Gw,Gd,Sd,&Mjd);mjd2gpswd(Mjd+d/86400.0,&Gw2,&Gd2,&Sd2);
             fprintf(stdout,"%5.4d%2.1d%12.4lf\n",Gw2,Gd2,Sd2);break;
      case 481:if(fscanf(stdin,"%d%d%d%d%d%d",&Gw,&Gd,&Sd,&Gw2,&Gd2,&Sd2)!=6)continue;
             gpswd2mjd(Gw,Gd,Sd,&Mjd);gpswd2mjd(Gw2,Gd2,Sd2,&Mjd2);
             fprintf(stdout,"%20.12lf\n",(Mjd-Mjd2)*86400.0);break;
      case 500:if(fscanf(stdin,"%d%lf\n",&Gw,&d)!=2)continue;
             gpsws2mjd(Gw,0,&Mjd);mjd2gpsws(Mjd-d,&Gw2,&Gs2);
             fprintf(stdout,"%5.4d\n",Gw2);break;
      case 590:if(fscanf(stdin,"%d%lf\n",&Gw,&d)!=2)continue;
             gpsws2mjd(Gw,0,&Mjd);mjd2gpsws(Mjd-d,&Gw2,&Gs2);
             fprintf(stdout,"%5.4d\n",Gw2);break;
      case 580:if(fscanf(stdin,"%d%d\n",&Gw,&Gw2)!=2)continue;
             gpsws2mjd(Gw,0,&Mjd);gpsws2mjd(Gw2,0,&Mjd2);
             fprintf(stdout,"%20.12lf\n",Mjd-Mjd2);break;
      case 501:if(fscanf(stdin,"%d%lf%lf\n",&Gw,&Gs,&d)!=3)continue;
             gpsws2mjd(Gw,Gs,&Mjd);mjd2gpsws(Mjd+d/86400.0,&Gw2,&Gs2);
             fprintf(stdout,"%5.4d%12.4lf\n",Gw2,Gs2);break;
      case 591:if(fscanf(stdin,"%d%lf%lf\n",&Gw,&Gs,&d)!=3)continue;
             gpsws2mjd(Gw,Gs,&Mjd);mjd2gpsws(Mjd+d/86400.0,&Gw2,&Gs2);
             fprintf(stdout,"%5.4d%12.4lf\n",Gw2,Gs2);break;
      case 581:if(fscanf(stdin,"%d%lf%d%lf\n",&Gw,&Gs,&Gw2,&Gs2)!=4)continue;
             gpsws2mjd(Gw,Gs,&Mjd);gpsws2mjd(Gw2,Gs2,&Mjd2);
             fprintf(stdout,"%20.12lf\n",(Mjd-Mjd2)*86400.0);break;
      default:return 1;
   }
   return 0;
}