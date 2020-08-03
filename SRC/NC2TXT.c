
//C source of Mapoet Niphy
/*
Author    :  Mapoet Niphy
Date      :  2020
Institude :  SHAO
Name      :  NC2TXT
Project   :  GEOTOOL

  Created by Mapoet Niphy on 2020/06/10.
  Copyright © 2020年 Mapoet Niphy. All rights reserved.

*/
#include <string.h>
#include <stdio.h>
#include <netcdf.h>
#include "CTIMES.h"
#define R2D (180/3.1415926)

int  main(int argc, char**argv){
	int itype = -1, ntype[] = { 18, 6, 25, 12 };
	const char names[][25][15] = { { "time", "caL1Snr", "pL1Snr", "pL2Snr", "xLeo", "yLeo", "zLeo", "xdLeo", "ydLeo", "zdLeo", "xGps", "yGps", "zGps", "xdGps", "ydGps", "zdGps", "exL1", "exL2" }, { "MSL_alt", "GEO_lat", "GEO_lon", "OCC_azi", "TEC_cal", "ELEC_dens" }, { "time", "caL1Snr", "pL1Snr", "pL2Snr", "xLeo", "yLeo", "zLeo", "xdLeo", "ydLeo", "zdLeo", "xGps", "yGps", "zGps", "xdGps", "ydGps", "zdGps", "exL1", "exL2", "exLC", "xmdl", "xmdldd", "xrng", "xmdl2", "xmdldd2", "xrng2" }, { "Lat", "Lon", "MSL_alt", "Ref", "Azim", "Pres", "Temp", "Bend_ang", "Opt_bend_ang", "Impact_parm", "Bend_ang_stdv", "Ref_stdv"/*, "OL_par", "OL_ipar", "OL_vec1", "OL_vec2", "OL_vec3", "OL_vec4"*/ } };
	int ncid = 0, varid = 0, itemp = 0;
	int i,j,inum = 0,wk;
	double ep[6] = { 0 },mjd, **vals = NULL;
	int rt = nc_open(argv[1], 0, &ncid), nw = 0;
	if (argc != 2){
		fprintf(stderr, "help:\n %s ncfile\n", argv[0]);
		return 1;
	}
	if (strstr(argv[1], "ionPhs") != NULL)itype = 0;
	if (strstr(argv[1], "ionPrf") != NULL)itype = 1;
	if (strstr(argv[1], "atmPhs") != NULL)itype = 2;
	if (strstr(argv[1], "atmPrf") != NULL)itype = 3;
	if(itype==-1)return -2;
	rt += nc_get_att_int(ncid, NC_GLOBAL, "year", &itemp);
	ep[0] = itemp;
	rt += nc_get_att_int(ncid, NC_GLOBAL, "month", &itemp);
	ep[1] = itemp;
	rt += nc_get_att_int(ncid, NC_GLOBAL, "day", &itemp);
	ep[2] = itemp;
	rt += nc_get_att_int(ncid, NC_GLOBAL, "hour", &itemp);
	ep[3] = itemp;
	rt += nc_get_att_int(ncid, NC_GLOBAL, "minute", &itemp);
	ep[4] = itemp;
	ymd2mjd(ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], &mjd);
	rt += nc_get_att_double(ncid, NC_GLOBAL, "second", ep + 5);
	//fprintf(stdout, "time:%4.4d %2.2d %2.2d %2.2d %2.2d %4.1lf\n", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], (int)ep[5]);
	vals = (double**)malloc(sizeof(double*)*ntype[itype]);
	rt += nc_inq_varid(ncid, names[itype][0], &varid);
	rt += nc_inq_dimlen(ncid, varid, &inum);
	for (i = 0; i < ntype[itype]; i++){
		//fprintf(stdout, "%20s",names[itype][i]);
		vals[i] = (double*)malloc(sizeof(double)*inum);
		rt += nc_inq_varid(ncid, names[itype][i], &varid);
		rt += nc_get_var_double(ncid, varid, vals[i]);
	}
	rt+=nc_close(ncid);
	if (rt != 0)return -1;
	//fprintf(stdout, "\n");
	for (i = 0; i < inum; i++)
	{
		for (j = 0; j < ntype[itype]; j++)
		{
			if (j == 0){
				mjd2gpsws(mjd + vals[j][i] / 86400, &wk, &ep[5]);
				fprintf(stdout, "%4.4d %15.3lf", wk,ep[5]);
			}
			else
				fprintf(stdout, "%20.5lf", vals[j][i]);
		}
		fprintf(stdout, "\n");
	}
	for (i = 0; i < ntype[itype]; i++)free(vals[i]);
	free(vals);
	return 0;
}
