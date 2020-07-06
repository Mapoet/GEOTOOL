
//C source of Mapoet Niphy
/*
Author    :  Mapoet Niphy
Date      :  2020
Institude :  SHAO
Name      :  POSCONVERT
Project   :  GEOTOOL

  Created by Mapoet Niphy on 2020/06/10.
  Copyright © 2020年 Mapoet Niphy. All rights reserved.

*/
#include "CTIMES.h"
#include "XFORM.h"
/*
argv: coorType2type, coorSys2sys:Time  
coorType : e(lipise),p(olar), c(art), l(athour)
coorSys  : gei geo gse gsm mag sm
Time     : mjd;
           wk,sec;
           wk,day,sec;
           y,d,hour,min,sec;
           y,m,d,hour,min,sec;
*/

int main(int argc, char **argv)
{
   double rin[3 + 6], temp[2][3], rout[3], ep[6], sec;
   int coorType, times[2], p, dates[4];
   char *ptime, psys[30] = {0};
   if (argc < 2)
      return -1;
   coorType = atoi(argv[1]);
   if (argc > 2 && (ptime = strstr(argv[2], ":")) != NULL)
   {
      p = 0;
      strncpy(psys, argv[2], ptime - argv[2]);
      ptime++;
      p=sscanf(ptime,"%lf%lf%lf%lf%lf%lf",rin+3,rin+4,rin+5,rin+6,rin+7,rin+8);
      switch (p)
      {
      case 1:
         mjd2doy(rin[3], dates, dates + 1, dates + 2, dates + 3, &sec);
         times[0] = dates[0] * 1000 + dates[1];
         times[1] = ((dates[2] * 60 + dates[3]) * 60 + sec) * 1e3;
         break;
      case 2:
         gpsws2doy((int)rin[3], rin[4], dates, dates + 1, dates + 2, dates + 3, &sec);
         times[0] = dates[0] * 1000 + dates[1];
         times[1] = ((dates[2] * 60 + dates[3]) * 60 + sec) * 1e3;
         break;
      case 3:
         gpswd2doy((int)rin[3], (int)rin[4], rin[5], dates, dates + 1, dates + 2, dates + 3, &sec);
         times[0] = dates[0] * 1000 + dates[1];
         times[1] = ((dates[2] * 60 + dates[3]) * 60 + sec) * 1e3;
         break;
      case 5:
         times[0] = (int)rin[3] * 1000 + (int)rin[4];
         times[1] = ((rin[5] * 60 + rin[6]) * 60 + rin[7]) * 1e3;
         break;
      case 6:
         ymd2doy((int)rin[3], (int)rin[4], rin[5], dates);
         times[0] = rin[3] * 1000 + dates[0];
         times[1] = ((rin[6] * 60 + rin[7]) * 60 + rin[8]) * 1e3;
         break;
      default:
         break;
      }
   }
   while (!feof(stdin))
      switch (coorType)
      {
      case 12:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         geod2ecef(rin, temp[0], NULL);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         cart2pol(temp[1], &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]/DEG2RAD, rout[2]);
         break;
      case 13:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         geod2ecef(rin, temp[0], NULL);
         if (argc < 3 || p <= 0)
         {
            memcpy(rout, temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], rout);
         }
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 14:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         geod2ecef(rin, temp[0], NULL);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         cart2pol(temp[1], &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]*24/360/DEG2RAD, rout[2]);
         break;
      case 21:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         rin[0] *= DEG2RAD;
         rin[1] *= DEG2RAD;
         pol2cart(rin[0], rin[1], rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         ecef2geod(temp[1], rout, NULL);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 23:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         rin[0] *= DEG2RAD;
         rin[1] *= DEG2RAD;
         pol2cart(rin[0], rin[1], rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(rout, temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], rout);
         }
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 24:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         pol2cart(rin[0]*DEG2RAD,rin[1]*DEG2RAD,rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         cart2pol(temp[1], &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]*24/360/DEG2RAD, rout[2]);
         break;
      case 31:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[0], rin, sizeof(double)*3);
         }
         else
         {
            xform(psys, times, rin, temp[0]);
         }
         ecef2geod(temp[0], rout, NULL);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 32:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[0], rin, sizeof(double)*3);
         }
         else
         {
            xform(psys, times, rin, temp[0]);
         }
         cart2pol(rin, &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]/DEG2RAD, rout[2]);
         break;
      case 34:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[0], rin, sizeof(double)*3);
         }
         else
         {
            xform(psys, times, rin, temp[0]);
         }
         cart2pol(rin, &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]*24/360/DEG2RAD, rout[2]);
         break;
      case 41:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         rin[0] *= DEG2RAD;
         rin[1] /= 24/360/DEG2RAD;
         pol2cart(rin[0], rin[1], rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         ecef2geod(temp[1], rout, NULL);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 42:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         pol2cart(rin[0]*DEG2RAD,rin[1]/24*360*DEG2RAD,rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         cart2pol(temp[1], &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]/DEG2RAD, rout[2]);
         break;
      case 43:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         rin[0] *= DEG2RAD;
         rin[1] /= 24/360/DEG2RAD;
         pol2cart(rin[0], rin[1], rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(rout, temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], rout);
         }
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 1:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         geod2ecef(rin, temp[0], NULL);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         ecef2geod(temp[1], rout,NULL);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 2:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         pol2cart(rin[0]*DEG2RAD,rin[1]*DEG2RAD,rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         cart2pol(temp[1], &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]/DEG2RAD, rout[2]);
         break;
      case 3:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         if (argc < 3 || p <= 0)
         {
            memcpy(rout, rin, sizeof(double)*3);
         }
         else
         {
            xform(psys, times, rin, rout);
         }
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0], rout[1], rout[2]);
         break;
      case 4:
         if (fscanf(stdin, "%lf%lf%lf", &rin[0], &rin[1], &rin[2]) != 3)
            continue;
         pol2cart(rin[0]*DEG2RAD,rin[1]/24*360*DEG2RAD,rin[2], temp[0]);
         if (argc < 3 || p <= 0)
         {
            memcpy(temp[1], temp[0], sizeof(double)*3);
         }
         else
         {
            xform(psys, times, temp[0], temp[1]);
         }
         cart2pol(temp[1], &rout[0], &rout[1], &rout[2]);
         fprintf(stdout, "%15.3lf%15.3lf%15.4lf\n", rout[0]/DEG2RAD, rout[1]*24/360/DEG2RAD, rout[2]);
         break;
      default:
         return 1;
      }
   return 0;
}
