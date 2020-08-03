/*
Author    :  Mapoet Niphy
Date      :  2020
Institude :  SHAO
Name      :  OBS2STEC
Project   :  GEOTOOL

  Created by Mapoet Niphy on 2020/06/10.
  Copyright © 2020年 Mapoet Niphy. All rights reserved.

*/
#include "rtklib.h"
#define  MAXOBSEPOCH   (24 * 60 * 60/24)
#define  MAXWINDOWSIZE      20
#define  NTYPE              3// P,L,D
typedef enum {
	FIT_WS=1,
	FIT_PWS,
	FIT_OC
}FitType;
typedef enum{
	CACU_CSP,
	CACU_P,
	CACU_OC
}CacuType;
typedef struct{
	int w;
	double lambda;
	double pho;
}forsmooth;
typedef struct{
	int w;
	int k;
	double lambda;
	double pho;
}forpoly;
typedef struct{
	double sigmw;// noise of WM(m)
	double sigif;// noise of ion - free phase(m)
	double dgmax;// cycle - slip threshold of LG(m)
	double dmp1;// cycle - slip threshold of MP1(m)
	double dmp2;// cycle - slip threshold of MP2(m)
	double dwl;// cycle - slip threshold of WL(cycle)
	double outl;// outlier threshold(sigma)
	int elwe;// threshold elevation weighting(1:on, 0 : off)
	int	   wind[2];// slip detection window(points)
	// wind(1) : moving avarage window width
	//wind(2) : avaraging window width
	double reps;// repair cycle - slip flag(1:on, 0 : off)
	double npnt;// no of obs points for fitting
	double nmax;//degree of poly. for fitting
}forcombine;
union fits
{
	forsmooth	ws;
	forpoly		pws;
	forcombine	oc;
};
typedef struct{
	char				sp3fn[50];
	char				obsfn[50];
	FILE*				ftec;
	double				frequence[NFREQ];
	double				hiono;
	double              targetalt[3];//max and min,if 3,must upto targetalt[2];(for gnss-leo)
	double				cutele[3];//max and min,if 3 ,must upto cutele[2];
	double				septime;
	FitType				prft;
	FitType				tecft;
	union fits			pseudofit, tecfit;
	CacuType			cacutype;
	double tecrange;
	double tecmin;
	double tecmax;
}forparaments;
int ParserArgs		  (int argc, char**argv, forparaments*p);
void GetSatTesSeries  (int sat, sta_t sta, nav_t*navsp3, gtime_t*t, double**in, double**out, short*hlt, int num, forparaments*p);
void WritetoFile	  (int sat, sta_t sta, gtime_t*t, double**out, short*hlt, int num, FILE*fid);
void SpanProceeding   (int sat, sta_t sta, nav_t*navsp3, gtime_t*t, double**in, double**out, short*hlt, int s, int e, forparaments*p);
void TecProceeding    (gtime_t*t, double **in, double**out, short*hlt, int s, int e, forparaments*p);
void CheckandRepairPLD(gtime_t*t, double**in, double**out, short*hlt, int s, int e, forparaments*p);
void CheckandRepairTec(gtime_t*t, double **in, double**out, short*hlt, int s, int e, forparaments*p);
void Fit_PWS          (gtime_t*t, double **in, double**out, short*hlt, int s, int e, int k, union fits p);
void Fit_WS           (gtime_t*t, double **in, double**out, short*hlt, int s, int e, int k, union fits p);
void Fit_OC           (gtime_t*t, double **in, double**out, short*hlt, int s, int e, union fits p, double*frequence);
void Cacu_byPseudoRangeWithCSP(gtime_t*t, double **in, double**out, short*hlt, int s, int e, double*frequence);
void Cacu_byPhase     (gtime_t*t, double **in, double**out, short*hlt, int s, int e, double*frequence);
void Cacu_byOC        (gtime_t*t, double **in, double**out, short*hlt, int s, int e, double*frequence);
/*
OBS2TEC:-sp3:igs1923?.sp3 -obs:20na3180.16o -tec:20na3180_cut20.stec  -fpt:WS   "20 2.5 0.2 1"
OBS2TEC:-sp3:peph\gfz19230.sp3 -obs:day\ARUB3180.16o -tec:stec\ARUB3180.16o -prft:PWS   "15 3 0.4 3" -tft:WS   "15 0.4 3"
CUTELEV:1lsu3350.14n 1lsu3350.14o 1lsu3350.14s 15 350 rnx.14o
*/
int main(int argc, char**argv)
{
	int num = 0;
	int sys = 0;
	obs_t obs = { 0 };
	nav_t navsp3 = { 0 };
	sta_t sta = { 0 };
	gtime_t*t;
	double**in;
	double**out;
	short*hlt;
	forparaments p;
	if(ParserArgs(argc, argv, &p)!=0)
		return -1;
	readsp3(p.sp3fn, &navsp3, 1);
	if(navsp3.ne==0)
		return 1;
	if(readrnx(p.obsfn, 1, "", &obs, NULL, &sta)!=1)
		return 2;
	for (int i = 0; i < MAXSAT; i++)
	{
		sys=satsys(i, NULL);
		if (sys != SYS_GPS&&sys != SYS_GLO)
			continue;
		num = 0;
		for (int j = 0; j < obs.n; j++)
		if (obs.data[j].sat == i)
			num++;
		if (!num)
			continue;
		if (!(t = (gtime_t*)malloc(sizeof(gtime_t)* num)))
		{
			fprintf(stderr, "malloc t[%5.5d] of gtime_t err!\n", num);
			system("pause");
			return -2;
		}
		if (!(hlt = (short*)malloc(sizeof(short)* num)))
		{
			fprintf(stderr, "malloc hlt[%5.5d] of short err!\n", num);
			system("pause");
			return -2;
		}
		if (!(in = (double**)malloc(sizeof(double*)* num)))
		{
			fprintf(stderr, "malloc in[%5.5d] of double* err!\n", num);
			system("pause");
			return -2;
		}
		if (!(out = (double**)malloc(sizeof(double*)* num)))
		{
			fprintf(stderr, "malloc out[%5.5d] of double* err!\n", num);
			system("pause");
			return -2;
		}
		for (int j = 0, n = 0; j < obs.n; j++)
			if (obs.data[j].sat == i)
			{
				if (!(in[n] = (double*)malloc(sizeof(double)* NTYPE*NFREQ)))
				{
					fprintf(stderr, "malloc in[%5.5d][%2.2d] err!\n", n, NTYPE * NFREQ);
					system("pause");
					return -2;
				}
				if (!(out[n] = (double*)malloc(sizeof(double)* 12)))
				{
					fprintf(stderr, "malloc out[%5.5d][10] err!\n", n);
					system("pause");
					return -2;
				}
				t[n] = obs.data[j].time;
				for (int ifr = 0; ifr < NFREQ; ifr++)
				{
					in[n][ifr * NTYPE + 0] = obs.data[j].P[ifr];
					in[n][ifr * NTYPE + 1] = obs.data[j].L[ifr];
					in[n][ifr * NTYPE + 2] = obs.data[j].D[ifr];
				}
				n++;
			}
		GetSatTesSeries(i, sta, &navsp3, t, in, out, hlt, num,&p);
		WritetoFile(i, sta, t, out, hlt, num, p.ftec);
		for (int i = 0; i < num; i++){
			free(in[i]);
			free(out[i]);
		}
		free(in);
		free(out);
	}
	if(p.ftec!=stdout)
		fclose(p.ftec);
	return 0;
}
int ParserArgs(int argc, char**argv, forparaments*p){
	p->cutele[0] = 90 * D2R;
	p->cutele[1] = -10 * D2R;
	p->cutele[2] = 0;
	p->frequence[0] = FREQ1;
	p->frequence[1] = FREQ2;
	p->frequence[2] = FREQ5;
	p->frequence[3] = FREQ7;
	p->ftec = stdout;
	p->hiono = 350e3;
	strcpy(p->obsfn, "");
	p->septime = 30;
	p->prft = FIT_OC;
	p->pseudofit.oc.wind[0] = 32;
	p->pseudofit.oc.wind[0] = 16;
	p->tecft = 0;
	strcpy(p->sp3fn, "*.sp3");
	p->targetalt[0] = 80e3;
	p->targetalt[1] = 800e3;
	p->targetalt[2] = 8000e3;
	p->tecmax = 500;
	p->tecmin = -500;
	p->tecrange = 1000;
	p->cacutype = CACU_P;
	p->cacutype = CACU_CSP;
	p->cacutype = CACU_OC;
	for (int i = 1; i < argc; i++)
	{
		if (strncmp(argv[i], "-sp3", 4) == 0 && (argv[i][4] == ':' || argv[i][4] == '='))
			strcpy(p->sp3fn, argv[i] + 5);
		if (strncmp(argv[i], "--Sp3File", 9) == 0 && (argv[i][9] == ':' || argv[i][9] == '='))
			strcpy(p->sp3fn, argv[i] + 10);
		if (strncmp(argv[i], "-obs", 4) == 0 && (argv[i][4] == ':' || argv[i][4] == '='))
			strcpy(p->obsfn, argv[i] + 5);
		if (strncmp(argv[i], "--ObsFile", 9) == 0 && (argv[i][9] == ':' || argv[i][9] == '='))
			strcpy(p->obsfn, argv[i] + 10);
		if (strncmp(argv[i], "-tec", 4) == 0 && (argv[i][4] == ':' || argv[i][4] == '='))
			p->ftec = fopen(argv[i]+5, "w");
		if (strncmp(argv[i], "--TecFile", 9) == 0 && (argv[i][9] == ':' || argv[i][9] == '='))
			p->ftec = fopen(argv[i] + 10, "w");
		if (strncmp(argv[i], "-nt", 3) == 0 && (argv[i][3] == ':' || argv[i][3] == '='))
			p->septime = atoi(argv[i] + 4);
		if (strncmp(argv[i], "--NTIME", 7) == 0 && (argv[i][7] == ':' || argv[i][7] == '='))
			p->septime = atoi(argv[i] + 8);
		if (strncmp(argv[i], "-prft", 5) == 0 && (argv[i][5] == ':' || argv[i][5] == '='))
		{
			if (strcmp(argv[i] + 6, "no") == 0 || strcmp(argv[i] + 6, "") == 0 || strcmp(argv[i] + 6, "N") == 0 || strcmp(argv[i] + 6, "No") == 0)
			{
				p->prft = 0;
			}
			if (strcmp(argv[i] + 6, "WindowSmooth") == 0 || strcmp(argv[i] + 6, "WS") == 0)
			{
				p->prft = FIT_WS;
				sscanf(argv[++i], "%d%lf%lf", &p->pseudofit.ws.w, &p->pseudofit.ws.lambda, &p->pseudofit.ws.pho);
			}
			if (strcmp(argv[i] + 6, "PolyWindowSmooth") == 0 || strcmp(argv[i] + 6, "PWS") == 0)
			{
				p->prft = FIT_PWS;
				sscanf(argv[++i], "%d%d%lf%lf", &p->pseudofit.pws.w, &p->pseudofit.pws.k, &p->pseudofit.pws.lambda, &p->pseudofit.pws.pho);
			}
			if (strcmp(argv[i] + 6, "ObsCombine") == 0 || strcmp(argv[i] + 6, "OC") == 0)
			{
				p->prft = FIT_OC;
				//sscanf(argv[++i], "%lf", &p->pseudofit.oc.pho);
			}
		}
		if (strncmp(argv[i], "--PseudoRangeFitType", 20) == 0 && (argv[i][20] == ':' || argv[i][20] == '='))
		{
			if (strcmp(argv[i] + 21, "no") == 0 || strcmp(argv[i] + 21, "") == 0 || strcmp(argv[i] + 21, "N") == 0 || strcmp(argv[i] + 21, "No") == 0)
			{
				p->prft = 0;
			}
			if (strcmp(argv[i] + 21, "WindowSmooth") == 0 || strcmp(argv[i] + 21, "WS") == 0)
			{
				p->prft = FIT_WS;
				sscanf(argv[++i], "%d%lf%lf", &p->pseudofit.ws.w, &p->pseudofit.ws.lambda, &p->pseudofit.ws.pho);
			}
			if (strcmp(argv[i] + 21, "PolyWindowSmooth") == 0 || strcmp(argv[i] + 21, "PWS") == 0)
			{
				p->prft = FIT_PWS;
				sscanf(argv[++i], "%d%d%lf%lf", &p->pseudofit.pws.w, &p->pseudofit.pws.k, &p->pseudofit.pws.lambda, &p->pseudofit.pws.pho);
			}
			if (strcmp(argv[i] + 21, "ObsCombine") == 0 || strcmp(argv[i] + 21, "OC") == 0)
			{
				p->prft = FIT_OC;
				//sscanf(argv[++i], "%lf",&p->pseudofit.oc.pho);
			}
		}
		if (strncmp(argv[i], "-tft", 4) == 0 && (argv[i][4] == ':' || argv[i][4] == '='))
		{
			if (strcmp(argv[i] + 5, "no") == 0 || strcmp(argv[i] + 5, "") == 0 || strcmp(argv[i] + 5, "N") == 0 || strcmp(argv[i] + 5, "No") == 0)
			{
				p->tecft = 0;
			}
			if (strcmp(argv[i] + 5, "WindowSmooth") == 0 || strcmp(argv[i] + 5, "WS") == 0)
			{
				p->tecft = FIT_WS;
				sscanf(argv[++i], "%d%lf%lf", &p->tecfit.ws.w, &p->tecfit.ws.lambda, &p->tecfit.ws.pho);
			}
			if (strcmp(argv[i] + 5, "PolyWindowSmooth") == 0 || strcmp(argv[i] + 5, "PWS") == 0)
			{
				p->tecft = FIT_PWS;
				sscanf(argv[++i], "%d%d%lf%lf", &p->tecfit.pws.w, &p->tecfit.pws.k, &p->tecfit.pws.lambda, &p->tecfit.pws.pho);
			}
			if (strcmp(argv[i] + 5, "ObsCombine") == 0 || strcmp(argv[i] + 5, "OC") == 0)
			{
				p->tecft = FIT_OC;
				//sscanf(argv[++i], "%lf", &p->tecfit.oc.pho);
			}
		}
		if (strncmp(argv[i], "--TecFitType", 12) == 0 && (argv[i][12] == ':' || argv[i][12] == '='))
		{
			if (strcmp(argv[i] + 13, "no") == 0 || strcmp(argv[i] + 13, "") == 0 || strcmp(argv[i] + 13, "N") == 0 || strcmp(argv[i] + 13, "No") == 0)
			{
				p->tecft = 0;
			}
			if (strcmp(argv[i] + 13, "WindowSmooth") == 0 || strcmp(argv[i] + 13, "WS") == 0)
			{
				p->tecft = FIT_WS;
				sscanf(argv[++i], "%d%lf%lf", &p->tecfit.ws.w, &p->tecfit.ws.lambda, &p->tecfit.ws.pho);
			}
			if (strcmp(argv[i] + 13, "PolyWindowSmooth") == 0 || strcmp(argv[i] + 13, "PWS") == 0)
			{
				p->tecft = FIT_PWS;
				sscanf(argv[++i], "%d%d%lf%lf", &p->tecfit.pws.w, &p->tecfit.pws.k, &p->tecfit.pws.lambda, &p->tecfit.pws.pho);
			}
			if (strcmp(argv[i] + 13, "ObsCombine") == 0 || strcmp(argv[i] + 13, "OC") == 0)
			{
				p->tecft = FIT_OC;
				//sscanf(argv[++i], "%lf", &p->tecfit.oc);
			}
		}
		if (strcmp(argv[i], "-tecrg") == 0)
			sscanf(argv[++i], "%lf%lf%lf", &p->tecmin, &p->tecmax, &p->tecrange);
		if (strcmp(argv[i], "--TecRange") == 0)
			sscanf(argv[++i], "%lf%lf%lf", &p->tecmin, &p->tecmax, &p->tecrange);
	}
	return 0;
}
void GetSatTesSeries(int sat, sta_t sta, nav_t*navsp3, gtime_t*t, double**in, double**out, short*hlt, int num, forparaments*p){
	char satn[4];
	satno2id(sat, satn);
	switch (satn[0])
	{
	case 'C':p->frequence[0] = FREQ1_CMP; p->frequence[1] = FREQ2_CMP; p->frequence[2] = FREQ3_CMP; break;
	case 'R':p->frequence[0] = FREQ1_GLO; p->frequence[1] = FREQ2_GLO; p->frequence[2] = FREQ3_GLO; break;
	default:
		p->frequence[0] = FREQ1; p->frequence[1] = FREQ2, p->frequence[2] = FREQ5; break;
	}
	double *dt = (double*)malloc(sizeof(double)*(num - 1)), tmin = 3000, sep = 0;
	int *span = (int*)malloc(sizeof(int)*(num - 1) / 5), nspan = 1;
	span[0] = 0;
	for (int i = 0; i < num - 1; i++)
	{
		dt[i] = (t[i + 1].time - t[i].time) - (t[i + 1].sec - t[i].sec);
		tmin = (tmin<fabs(dt[i]) ? tmin : fabs(dt[i]));
	}
	if (fabs(tmin) < 1e-3)
		tmin = dt[0];
	if (p->septime>0)
		sep = p->septime;
	else
		sep = p->septime*tmin;
	for (int i = 0; i< num - 1; i++)
	if (fabs(dt[i])>sep)
		span[nspan++] = i + 1;
	span[nspan++] = num;
	free(dt);
	for (int i = 0; i < nspan-1;i++)
		SpanProceeding(sat, sta, navsp3, t, in, out, hlt, span[i], span[i + 1], p);
	free(span);
}
void SpanProceeding(int sat, sta_t sta, nav_t*navsp3, gtime_t*t, double**in, double**out, short*hlt, int s, int e, forparaments*p){
	double tx[6] = { 0 }, rt[2] = { 0 }, var = 0,
		dr[3] = { 0 }, dl = 0, dphi[NFREQ] = { 0 }, pos[3] = { 0 },
		azel[2] = { 0 }, pppos[3] = { 0 };
	for(int i = s; i < e; i++)
	{
		if ((!peph2pos(t[i], sat, navsp3, 1, tx, rt, &var)) || fabs(in[i][0])<1e-3 || fabs(in[i][1])<1e-3
			|| fabs(in[i][3])<1e-3 || fabs(in[i][4])<1e-3)
		{
			memset(out[i], 0, sizeof(double)* 10);
			hlt[i] = 0;
			continue;
		}
		dl=geodist(tx, sta.pos, dr);
		ecef2pos(sta.pos, pos);
		satazel(pos, dr, azel);
		if (azel[1] < p->cutele[1]){
			memset(out[i], 0, sizeof(double)* 10);
			hlt[i] = 0;
			continue;
		}
		
		for (int ifr = 0; ifr < NFREQ; ifr++)
		{
			in[i][ifr * 3 + 0] += rt[0] * CLIGHT;
			in[i][ifr * 3 + 1] += rt[0] * p->frequence[0];
		}
		out[i][0] = tx[0];
		out[i][1] = tx[1];
		out[i][2] = tx[2];
		out[i][3] = sta.pos[0];
		out[i][4] = sta.pos[1];
		out[i][5] = sta.pos[2];
		out[i][10] = ionppp(pos, azel, RE_WGS84, p->hiono, pppos);
		out[i][6] = pppos[0];
		out[i][7] = pppos[1];
		out[i][7] = out[i][7] < 0 ? out[i][7] + 2 * PI : out[i][7];
		out[i][8] = azel[0];
		out[i][8] = out[i][8] < 0 ? out[i][8] + 2 * PI : out[i][8];
		out[i][9] = azel[1];
		hlt[i] = 1;
	}
	CheckandRepairPLD(t, in, out, hlt, s, e, p);
	TecProceeding(t, in, out, hlt, s, e, p);
	CheckandRepairTec(t, in, out, hlt, s, e, p);
}
void TecProceeding(gtime_t*t, double **in, double**out, short*hlt, int s, int e, forparaments*p){
	switch (p->cacutype)
	{
	case CACU_CSP:Cacu_byPseudoRangeWithCSP(t, in, out, hlt, s, e, p->frequence); break;
	case CACU_P:Cacu_byPhase(t, in, out, hlt, s, e, p->frequence); break;
	case CACU_OC:Cacu_byOC(t, in, out, hlt, s, e, p->frequence); break;
	default:
		Cacu_byPseudoRangeWithCSP(t, in, out, hlt, s, e, p->frequence); break;
		break;
	}
}
void WritetoFile(int sat, sta_t sta, gtime_t*t, double**out, short*hlt, int num, FILE*fid){
	char satn[4], stan[5];
	satno2id(sat, satn);
	strncpy(stan, sta.name, 4);
	stan[4] = '\0';
	double ep[6], ep0[6];
	time2epoch(t[0], ep0);
	for (int i = 0; i < num; i++)
	if (hlt[i]){
		time2epoch(t[i], ep);
		//fprintf(fid, "%4d%3d%3d%3d%3d%7.3lf%5s%5s%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf\n",
		//	(int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5], satn, sta.name, out[i][0], out[i][1], out[i][2], out[i][3], out[i][4], out[i][5]);
		//fprintf(fid, "%15.3lf%5s%5s%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf\n",
		//	ep[5] + 60 * (ep[4] + 60 * (ep[3] + 24 * (ep[2] - ep0[2]))), satn, stan, out[i][0], out[i][1], out[i][2], out[i][3], out[i][4], out[i][5]);
		for (int j = 0; j < 12; j++)
		{
			if (fabs(out[i][j])>1e10)
				out[i][j] = INFINITY;
			if (isinf(out[i][j])){
				if (i == 0){
					if (!isinf(out[i + 1][j]))
						out[i][j] = out[i + 1][j];
					else
						out[i][j] = 0;
				}
				else if (i == num - 1){
					if (!isinf(out[i - 1][j]))
						out[i][j] = out[i - 1][j];
					else
						out[i][j] = 0;
				}
				else if (!isinf(out[i - 1][j]) && !isinf(out[i + 1][j]))
					out[i][j] = (out[i - 1][j] + out[i + 1][j]) / 2;
				else if (!isinf(out[i - 1][j]))
					out[i][j] = out[i - 1][j];
				else if (!isinf(out[i + 1][j]))
					out[i][j] = out[i + 1][j];
				else out[i][j] = 0;
			}
		}
		fprintf(fid, "%4d%3d%3d%3d%3d%7.3lf%5s%5s%15.3lf%15.3lf%15.3lf%15.3lf%15.3lf%15.3lf%12.6lf%12.6lf%12.6lf%12.6lf%15.6lf%15.6lf\n",
			(int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5], satn, stan, out[i][0], out[i][1], out[i][2], out[i][3], out[i][4], out[i][5], out[i][6] * R2D, out[i][7] * R2D, out[i][8] * R2D, out[i][9] * R2D, out[i][10], out[i][11]);
	}
}
void CheckandRepairPLD(gtime_t*t, double**in, double** out, short*hlt, int s, int e, forparaments*p)
{
	if(!p->prft)
		return;
	if (p->prft == FIT_WS)
	for (int k = 0; k < NTYPE*NFREQ; k++)
		Fit_WS(t, in, out, hlt, s, e, k, p->pseudofit);
	if (p->prft == FIT_PWS)
	for (int k = 0; k < NTYPE*NFREQ; k++)
		Fit_PWS(t, in, out, hlt, s, e, k, p->pseudofit);
	if (p->prft == FIT_OC)
	//for (int k = 0; k < NTYPE*NFREQ; k++)
		Fit_OC(t, in, out, hlt, s, e, p->pseudofit,p->frequence);
}
void CheckandRepairTec(gtime_t*t, double **in, double**out, short*hlt, int s, int e, forparaments*p){
	if (e - s < 3)
		return;
	double*polyA, *polyL, *polyX, *polyArg, *polyN, *polyATL, tm, dt, dtec = 0;
	static double pre=0;
	int n = 4,m=e-s,c=0,ic;
	/*
	polyArg = zeros(1 + 1 + n, 1 + 1 + n);
	polyN = zeros(1 + 1 + n, 1 + 1 + n);
	polyATL = zeros(1 + 1 + n, 1);
	for (int it = s; it < e; it += m)
	{
		tm = 0;
		ic = 0;
		for (int i = it; i < it + m; i++)
		if(hlt[i]){
			tm +=t[i].time + t[i].sec;
			ic++;
			if (fabs(pre) < 1e-6)pre = out[i][11];
		}
		c = ic;
		tm /= c;
		polyA = zeros(1 + 1 + n, c);
		polyL = zeros(c, 1);
		polyX = zeros(c, 1);
		ic = 0;
		for (int i = it; i < it + m; i++)
		if(hlt[i]){
			dt = (t[i].time + t[i].sec - tm);// / (t[it].time + t[it].sec - t[it + m - 1].time - t[it + m - 1].sec) * 2;
			polyA[ic] = out[i][10];
			for (int j = 1; j <= n; j++)polyA[j*c + ic] = polyA[(j-1)*c + ic] * dt;
			polyA[(n+1) * c + ic] = -1;
			if (i>0 && fabs(out[i][11] - pre)>5.){
				if (hlt[i - 1])pre = out[i - 1][11];
				dtec+= pre - out[i][11];
				pre = out[i][11];
			}
				polyL[ic] = out[i][11]+dtec;
			ic++;
		}
		matmul("TN", n + 2, n + 2, c, 1., polyA, polyA, 0., polyN);
		matmul("TN", n + 2, 1, c, 1., polyA, polyL, 0., polyATL);
		//FILE*polyAf = fopen("polyAf.dat", "w");
		//matfprint(polyA, c, n + 2, 30, 6, polyAf);
		//fclose(polyAf);
		//FILE*polyLf = fopen("polyLf.dat", "w");
		//matfprint(polyL, c, 1, 30, 6, polyLf);
		//fclose(polyLf);
		//FILE*polyNf = fopen("polyNf.dat", "w");
		//matfprint(polyN, n + 2, n + 2, 30, 6, polyNf);
		//fclose(polyNf);
		//FILE*polyATLf = fopen("polyATLf.dat", "w");
		//matfprint(polyATL, n + 2, 1, 30, 6, polyATLf);
		//fclose(polyATLf);
		if (lusolve(polyN, polyATL, n+2 , polyArg) == -1){
			memset(polyArg, 0, sizeof(double)*(n+2));
		}
		for (int i = it; i < it + m; i++)
		if(hlt[i])out[i][11] += polyArg[n+1];
		fprintf(stdout, "DCB:%15.6lf\n", polyArg[n + 1]);
		fprintf(stdout, "VTEC:%15.6lf+%15.6lfdt+%15.6lfdt^2+%15.6lfdt^3+%15.6lfdt^4\n", polyArg[0], polyArg[1], polyArg[2], polyArg[3], polyArg[4]);
		//matmul("NN", c, 1, n+2, 1.0, polyA, polyArg, 0.0, polyX);
		//for (int ii = 0; ii < c; ii++)
		//	printf("%lf-%lf=%lf\n", polyL[ii], polyX[ii], polyL[ii] - polyX[ii]);
		if(c>0)pre = polyX[c-1];
		//FILE*polyXf = fopen("polyXf.dat", "w");
		//matfprint(polyX, c, 1, 30, 6, polyXf);
		//fclose(polyXf);
		free(polyA);
		free(polyL);
		//free(polyX);
	}
	free(polyArg);
	free(polyATL);
	free(polyN);
	*/
	//if (p->tecft == FIT_WS)
	//	Fit_WS (t, in, out,hlt, s, e, -11, p->tecfit);
	//if (p->tecft == FIT_PWS)
	//	Fit_PWS(t, in, out, hlt, s, e, -11, p->tecfit);
	
}
void Fit_PWS(gtime_t*t, double **in, double**out, short*hlt, int s, int e, int k, union fits p){
	double **data=in,tm, pm, pt, coff, val, *polyA, *polyL, *polyN, *polyATL, *polyArg, *polyLambda;
	int mark[MAXWINDOWSIZE] = { 0 }, pn, pn0;
	int cn, cm, ck, itype = 0,flag = 0;
	double ep[6];
	cn = (int)(p.pws.w*(1 - p.pws.lambda));
	ck = p.pws.k + 1;
	cm = p.pws.w - cn;
	polyA = zeros(cn, ck);
	polyL = zeros(cn, 1);
	polyN = zeros(ck, ck);
	polyATL = zeros(ck, 1);
	polyArg = zeros(ck, 1);
	polyLambda = zeros(ck, 1);
	if (k < 0)
	{
		k = -k;
		data = out;
		itype = 1;
	}
	if (e - s < ck + 3)
		return;
	for (int i = s + p.pws.w / 2; i < e - p.pws.w / 2; i += cm)
	{
		memset(mark, 0, cn);
		pn0 = -1;
		flag = 0;
		while (1){
			tm = 0; pn = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i + p.pws.w / 2; c++)
			if (abs(c - i)>cm / 2.&&hlt[c]&& (!mark[u++]))
			{
				tm += (t[c].time + t[c].sec);
				pn++;
			}
			if (pn < ck){
				flag = 1;
				break;
			}
			tm /= pn;
			if (pn0 == pn)
				break;
			pn0 = pn; pn = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i + p.pws.w / 2; c++)
			if (fabs(c - i)>cm / 2. &&hlt[c] && (!mark[u++]))
			{
				polyA[pn] = 1;
				for (int pk = 1; pk <= p.pws.k; pk++)
					polyA[pk*pn0 + pn] = polyA[(pk - 1)*pn0 + pn] * (t[c].time + t[c].sec - tm);
				polyL[pn] = data[c][k];
				pn++;
			}
			for (int pk = 0; pk <= p.pws.k; pk++){
				if (fabs(polyLambda[pk] = sqrt(dot(polyA + pk*pn0, polyA + pk*pn0, pn0))) < 1e-16)
					polyLambda[pk] = 1;
				for (int c = 0; c < pn0; c++)
					polyA[pk*pn0 + c] /= polyLambda[pk];
			}
			matmul("TN", ck, ck, pn, 1., polyA, polyA, 0., polyN);
			matmul("TN", ck,  1, pn, 1., polyA, polyL, 0., polyATL);
			if (lusolve(polyN, polyATL, ck, polyArg) == -1){
				memset(polyArg, 0, sizeof(double)*ck);
				memcpy(polyArg, polyATL, sizeof(double)* 2);
			}
			matmul("NN", pn, 1, ck, 1.0, polyA, polyArg, 0.0, polyL);
			pn = 0; pm = 0; pt = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i + p.pws.w / 2; c++)
			if (abs(c - i)>cm / 2.&&hlt[c] && (!mark[u++]))
			{
				pm += data[c][k] - polyL[pn];
				pt += (data[c][k] - polyL[pn])*(data[c][k] - polyL[pn]);
				pn++;
			}
			pm /= pn;
			pt = sqrt((pt - pn*pm*pm) / (pn - 1));
			pn = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i + p.pws.w / 2; c++)
			if (abs(c - i)>cm / 2.&&hlt[c] && (!mark[u++]))
			{
				if (fabs(data[c][k] - polyL[pn] - pm) > p.pws.pho * pt)
					mark[pn] = 1;
				pn++;
			}
		}
		time2epoch(t[i], ep);
		if (itype){
			fprintf(stderr, "TEC:%5d%3.2d%3.2d%3.2d%3.2d%6.2lf", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
			for (int pk = 0; pk < ck; pk++)
				fprintf(stderr, " %19.6lf %19.6lf", polyArg[pk], polyLambda[pk]);
			fprintf(stderr, "\n");
		}
		else {
			fprintf(stderr, "PR:%5d%3.2d%3.2d%3.2d%3.2d%6.2lf", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
			for (int pk = 0; pk < ck; pk++)
				fprintf(stderr, " %19.6lf %19.6lf", polyArg[pk], polyLambda[pk]);
			fprintf(stderr, "\n");
		}
		for (int c = (int)(i - cm / 2.); c < i + cm / 2.; c++)
		{
			coff = 1;
			val = coff / polyLambda[0] * polyArg[0];
			for (int pk = 1; pk <= p.pws.k; pk++)
			{
				coff = coff*(t[c].time + t[c].sec - tm);
				val += coff / polyLambda[pk] * polyArg[pk];
			}
			if (fabs(data[c][k] - val - pm)>p.pws.pho*pt)
			{
				//fprintf(stderr, "%15.4lf%15.4lf%15.4lf%15.4lf%15.4lf\n", data[c][k], val, pm, pt, fabs(data[c][k] - val - pm) / pt);
				//system("pause");
				data[c][k] = val;
			}
		}
	}
	free(polyA);
	free(polyArg);
	free(polyATL);
	free(polyL);
	free(polyLambda);
	free(polyN);
}
void Fit_P2WS(gtime_t*t, double **in, double**out, short*hlt, int s, int e, int k, union fits p){
	double **data = in, tm1,tm2, pm, pt, coff, val1,val2, *polyA, *polyL, *polyN, *polyATL, *polyArg , *polyLambda;
	int mark[MAXWINDOWSIZE] = { 0 }, pn, pn0,flag=0;
	int cn, cm, ck, itype = 0;
	double ep[6];
	cn = (int)(p.pws.w*(1 - p.pws.lambda)/2);
	ck = p.pws.k + 1;
	cm = p.pws.w - cn;
	polyA = zeros(cn, ck);
	polyL = zeros(cn, 1);
	polyN = zeros(ck, ck);
	polyATL = zeros(ck, 1);
	polyArg = zeros(2*ck, 1);
	polyLambda = zeros(2*ck, 1);
	if (k < 0)
	{
		k = -k;
		data = out;
		itype = 1;
	}
	if (e - s < 2*ck + cm)
		return;
	for (int i = s + p.pws.w / 2; i < e - p.pws.w / 2; i += cm)
	{
		memset(mark, 0, cn);
		pn0 = -1;
		flag = 0;
		while (1){
			tm1 = 0; pn = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i - cm / 2; c++)
			if (!mark[u++])
			{
				tm1 += (t[c].time + t[c].sec);
				pn++;
			}
			if (pn < ck){
				flag = 1;
				break;
			}
			tm1 /= pn;
			if (pn0 == pn)
				break;
			pn0 = pn; pn = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i - cm / 2; c++)
			if (hlt[c] && (!mark[u++]))
			{
				polyA[pn] = 1;
				for (int pk = 1; pk <= p.pws.k; pk++)
					polyA[pk*pn0 + pn] = polyA[(pk - 1)*pn0 + pn] * (t[c].time + t[c].sec - tm1);
				polyL[pn] = data[c][k];
				pn++;
			}
			for (int pk = 0; pk <= p.pws.k; pk++){
				if (fabs(polyLambda[pk] = sqrt(dot(polyA + pk*pn0, polyA + pk*pn0, pn0))) < 1e-16)
					polyLambda[pk] = 1;
				for (int c = 0; c < pn0; c++)
					polyA[pk*pn0 + c] /= polyLambda[pk];
			}
			matmul("TN", ck, ck, pn, 1., polyA, polyA, 0., polyN);
			matmul("TN", ck,  1, pn, 1., polyA, polyL, 0., polyATL);
			if (lusolve(polyN, polyATL, ck, polyArg) == -1){
				memset(polyArg, 0, sizeof(double)*ck);
				memcpy(polyArg, polyATL, sizeof(double)* 2);
			}
			matmul("NN", pn, 1, ck, 1.0, polyA, polyArg, 0.0, polyL);
			pn = 0; pm = 0; pt = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i - cm / 2; c++)
			if (hlt[c] && (!mark[u++]))
			{
				pm += data[c][k] - polyL[pn];
				pt += (data[c][k] - polyL[pn])*(data[c][k] - polyL[pn]);
				pn++;
			}
			pm /= pn;
			pt = sqrt((pt - pn*pm*pm) / (pn - 1));
			pn = 0;
			for (int c = i - p.pws.w / 2, u = 0; c < i - cm / 2; c++)
			if (hlt[c] && (!mark[u++]))
			{
				if (fabs(data[c][k] - polyL[pn] - pm) > p.pws.pho * pt)
					mark[pn] = 1;
				pn++;
			}
		}
		memset(mark, 0, cn);
		pn0 = -1;
		while (1){
			tm2 = 0; pn = 0;
			for (int c = i + p.pws.w / 2, u = 0; c > i + cm / 2; c--)
			if (!mark[u++])
			{
				tm2 += (t[c].time + t[c].sec);
				pn++;
			}
			if (pn < ck){
				flag = 1;
				break;
			}
			tm2 /= pn;
			if (pn0 == pn)
				break;
			pn0 = pn; pn = 0;
			for (int c = i + p.pws.w / 2, u = 0; c > i + cm / 2; c--)
			if (hlt[c] && (!mark[u++]))
			{
				polyA[pn] = 1;
				for (int pk = 1; pk <= p.pws.k; pk++)
					polyA[pk*pn0 + pn] = polyA[(pk - 1)*pn0 + pn] * (t[c].time + t[c].sec - tm2);
				polyL[pn] = data[c][k];
				pn++;
			}
			for (int pk = 0; pk <= p.pws.k; pk++){
				if (fabs(polyLambda[ck+pk] = sqrt(dot(polyA + pk*pn0, polyA + pk*pn0, pn0))) < 1e-16)
					polyLambda[ck+pk] = 1;
				for (int c = 0; c < pn0; c++)
					polyA[pk*pn0 + c] /= polyLambda[ck+pk];
			}
			matmul("TN", ck, ck, pn, 1., polyA, polyA, 0., polyN);
			matmul("TN", ck, 1, pn, 1., polyA, polyL, 0., polyATL);
			if (lusolve(polyN, polyATL, ck, polyArg+ck) == -1){
				memset(polyArg, 0, sizeof(double)*ck);
				memcpy(polyArg, polyATL, sizeof(double)* 2);
			}
			matmul("NN", pn, 1, ck, 1.0, polyA, polyArg+ck, 0.0, polyL);
			pn = 0; pm = 0; pt = 0;
			for (int c = i + p.pws.w / 2, u = 0; c > i + cm / 2; c--)
			if (hlt[c] && (!mark[u++]))
			{
				pm += data[c][k] - polyL[pn];
				pt += (data[c][k] - polyL[pn])*(data[c][k] - polyL[pn]);
				pn++;
			}
			pm /= pn;
			pt = sqrt((pt - pn*pm*pm) / (pn - 1));
			pn = 0;
			for (int c = i + p.pws.w / 2, u = 0; c > i + cm / 2; c--)
			if (hlt[c] && (!mark[u++]))
			{
				if (fabs(data[c][k] - polyL[pn] - pm) > p.pws.pho * pt)
					mark[pn] = 1;
				pn++;
			}
		}
		if (flag)
			continue;
		time2epoch(t[i], ep);
		if (itype){
			fprintf(stderr, "TEC1:%5d%3.2d%3.2d%3.2d%3.2d%6.2lf", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
			for (int pk = 0; pk < ck; pk++)
				fprintf(stderr, " %19.6lf %19.6lf", polyArg[pk], polyLambda[pk]);
			fprintf(stderr, "\n");
			fprintf(stderr, "TEC2:%5d%3.2d%3.2d%3.2d%3.2d%6.2lf", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
			for (int pk = 0; pk < ck; pk++)
				fprintf(stderr, " %19.6lf %19.6lf", polyArg[ck+pk], polyLambda[ck+pk]);
			fprintf(stderr, "\n");
		}
		else {
			fprintf(stderr, "PR1:%5d%3.2d%3.2d%3.2d%3.2d%6.2lf", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
			for (int pk = 0; pk < ck; pk++)
				fprintf(stderr, " %19.6lf %19.6lf", polyArg[pk], polyLambda[pk]);
			fprintf(stderr, "\n");
			fprintf(stderr, "PR2:%5d%3.2d%3.2d%3.2d%3.2d%6.2lf", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
			for (int pk = 0; pk < ck; pk++)
				fprintf(stderr, " %19.6lf %19.6lf", polyArg[ck + pk], polyLambda[ck + pk]);
			fprintf(stderr, "\n");
		}
		for (int c =(int)( i - cm / 2.), u = 0; c < i + cm / 2.; c++)
		{
			coff = 1;
			val1 = coff / polyLambda[0] * polyArg[0];
			for (int pk = 1; pk <= p.pws.k; pk++)
			{
				coff = coff*(t[c].time + t[c].sec - tm1);
				val1 += coff / polyLambda[pk] * polyArg[pk];
			}
			coff = 1;
			val2 = coff / polyLambda[ck] * polyArg[ck];
			for (int pk = 1; pk <= p.pws.k; pk++)
			{
				coff = coff*(t[c].time + t[c].sec - tm2);
				val2 += coff / polyLambda[ck+pk] * polyArg[ck+pk];
			}
			if (fabs(data[c][k] - val1 - pm)>p.pws.pho*pt)
			{
				//fprintf(stderr, "%15.4lf%15.4lf%15.4lf%15.4lf%15.4lf\n", data[c][k], val, pm, pt, fabs(data[c][k] - val - pm) / pt);
				//system("pause");
				data[c][k] = val1;
			}
		}
	}
	free(polyA);
	free(polyArg);
	free(polyATL);
	free(polyL);
	free(polyLambda);
	free(polyN);
}
void Fit_WS (gtime_t*t, double **in, double**out, short*hlt, int s, int e, int k, union fits p){
	int mark[MAXWINDOWSIZE] = { 0 }, pn, pn0,cn, cm;
	double **data=in,pm, pt;
	cn = (int)(p.ws.w*(1 - p.ws.lambda));
	cm = p.ws.w - cn;
	if (k < 0)
	{
		k = -k;
		data = out;
	}
	for (int i = s + p.ws.w / 2; i < e - p.ws.w / 2; i += cm)
	{
		memset(mark, 0, cn);
		pn0 = -1;
		while (1){
			pm = 0; pt = 0; pn = 0;
			for (int c = i - p.ws.w / 2, u = 0; c < i + p.ws.w / 2; c++)
			if (abs(c - i)>cm / 2.&&hlt[c] && (!mark[u++]))
			{
				pm += data[c][k];
				pt += data[c][k] * data[c][k];
				pn++;
			}
			if (pn){
				pm /= pn; pt = sqrt((pt - pn*pm*pm) / (pn - 1));
			}
			if ((pn0 == pn) || pn == 0)
				break;
			pn0 = pn;			
			for (int c = i - p.ws.w / 2, u = 0; c < i + p.ws.w / 2; c++)
			if (abs(c - i)>cm / 2.&&hlt[c] && (!mark[u++]))
			if (fabs(data[c][k] - pm)>p.ws.pho * pt)
				mark[u] = 1;
		}
		for (int c = i - p.ws.w / 2; c < i + p.ws.w / 2; c++)
		if (fabs(c - i) < cm / 2. && fabs(data[i][k] - pm)>p.ws.pho * pt)
			data[c][k] = pm;
	}
}
int pow10(int k){
	int re=1;
	for (int i = 0; i < k; i++)
		re *= 10;
	return re;
}
void Fit_OC(gtime_t*t, double **in, double**out, short*hlt, int s, int e, union fits p, double*frequence){
	double **data = in;
	if (e - s < 3)
		return;
	int nhlt = 0, n = 0,repeat=0;
	double lamb[NFREQ + 2],c[2];
	double **CM/*WM,IF,LG,MP1,MP2,WL*/, thres[5] = { 0 },nthres=0, aveDelta[5] = { 0 }, stdDelta[5] = { 0 }, ave, std;
	unsigned short *mark,pmark[5];
	char label[15];
	double ionocoff = 40.28 *(1 / (frequence[1] * frequence[1]) - 1 / (frequence[0] * frequence[0]));
	for (int i = 0; i < NFREQ; i++)
		lamb[i] = CLIGHT / frequence[i];
	lamb[3] = 1 / (1 / lamb[0] - 1 / lamb[1]);
	lamb[4] = 1 / (1 / lamb[0] + 1 / lamb[1]);
	c[0] = lamb[1] * lamb[1] / (lamb[1] * lamb[1] - lamb[0] * lamb[0]);
	c[1] = lamb[0] * lamb[0] / (lamb[1] * lamb[1] - lamb[0] * lamb[0]);
	for (int i = s, u = 0; i < e; i++, u++)
	if (hlt[i])nhlt++;
	CM = (double**)malloc(sizeof(double*)* 5);
	for (int i = 0; i < 5;i++)
		CM[i] = (double*)malloc(sizeof(double)*nhlt);
	mark = (unsigned short*)malloc(sizeof(unsigned short)*nhlt);
	memset(mark, 0, sizeof(unsigned short)*nhlt);
	for (int i = s, u = 0; i < e; i++)if (hlt[i]){
		CM[0][u] = lamb[3] * (in[i][1] - in[i][NTYPE + 1]) - lamb[4] * (in[i][0] / lamb[0] + in[i][NTYPE] / lamb[1]);
		CM[1][u] = lamb[0] * in[i][1] - lamb[1] * in[i][NTYPE + 1];
		CM[2][u] = in[i][0] - lamb[0] * in[i][1] - 2 * c[0] * (lamb[0] * in[i][1] - lamb[1] * in[i][NTYPE + 1]);
		CM[3][u] = in[i][NTYPE] - lamb[1] * in[i][NTYPE + 1] - 2 * c[1] * (lamb[0] * in[i][1] - lamb[1] * in[i][NTYPE + 1]);
		CM[4][u] = c[0] * lamb[0] * in[i][1] - c[1] * lamb[1] * in[i][NTYPE + 1] - c[0] * in[i][0] + c[1] * in[i][NTYPE];
		u++;
	}
	for (int k = 0; k < 5; k++)
	{
		thres[k] = (k==0?-2.5:-2.5);// p.oc.dgmax;
		n = 0; ave = 0; std = 0;
		for (int i = s, u = 0; i < e - 1; i++)
		if (hlt[i]){
			if (mark[u + 1]<pow10(k))
			{
				ave += (CM[k][u + 1] - CM[k][u]) / (t[i+1].time + t[i+1].sec - t[i].time - t[i].sec);
				std += (CM[k][u + 1] - CM[k][u])*(CM[k][u + 1] - CM[k][u]) / (t[i + 1].time + t[i + 1].sec - t[i].time - t[i].sec) / (t[i + 1].time + t[i + 1].sec - t[i].time - t[i].sec);
				n++;
			}
			u++;
		}
		repeat = 0;
		stdDelta[k] = 1e5;
		while (repeat<50){
			if (n == 0)
				break;
			ave /= n;
			std = sqrt((std - n*ave * ave) / (n - 1));
			if (fabs(stdDelta[k] - std) < 1e-6)
				break;
			aveDelta[k]=ave;
			stdDelta[k] = std;
			nthres = thres[k]>0 ? thres[k] : -thres[k] * stdDelta[k];
			n = 0; ave = 0; std = 0;
			for (int i = s, u = 0; i < e-1; i++)if (hlt[i])
			{
				if (mark[u + 1]<pow10(k))
					if(fabs(CM[k][u + 1] - CM[k][u] - aveDelta[k]) < nthres)
					{
						ave += CM[k][u + 1] - CM[k][u];
						std += (CM[k][u + 1] - CM[k][u])*(CM[k][u + 1] - CM[k][u]);
						n++;
					}
					else 
						mark[u + 1] += pow10(k);
				u++;
			}
			repeat++;
		}
	}
	for (int i = s, u = 0; i < e - 1; i++)if (hlt[i]){
		if (mark[u + 1] != 0){
			sprintf(label, "%5.5d", mark[u + 1]);
			printf("%15s%5d\n", label, i);
			for (int k = 0; k < 5; k++)
				pmark[k] = label[k] - '0';
		}
		u++;
	}
	for (int i = 0; i < 5; i++)
		free(CM[i]);
	free(CM);
}
void Cacu_byPseudoRangeWithCSP(gtime_t*t, double **in, double**out, short*hlt, int s, int e, double*frequence){
	double lambda = 40.28 *(1 / (frequence[1] * frequence[1]) - 1 / (frequence[0] * frequence[0]));
	for (int i = s + 1; i < e; i++)
	for (int k = 0; k < NFREQ; k++)
		in[i][NFREQ*k] = in[i][NFREQ*k] / i + (i - 1.) / i*(in[i - 1][NFREQ*k] - CLIGHT*(in[i][NFREQ*k + 1] / frequence[k] - in[i - 1][NFREQ*k+ 1] / frequence[k]));
	for (int i = s; i < e; i++)
	if (hlt[i]){
		out[i][11] = (in[i][0] - in[i][NTYPE]) / lambda *1e-16;;
	}
}
void Cacu_byPhase(gtime_t*t, double **in, double**out, short*hlt, int s, int e, double*frequence){
	double lambda = 40.28 *(1 / (frequence[1] * frequence[1]) - 1 / (frequence[0] * frequence[0]));
	for (int i = s; i < e; i++)
	if (hlt[i]){
		out[i][11] = ((in[i][1]) / frequence[0] - (in[i][NTYPE + 1]) / frequence[1])*CLIGHT / lambda *1e-16;;
	}
}
void Cacu_byOC(gtime_t*t, double **in, double**out, short*hlt, int s, int e, double*frequence){
	if (e - s < 3)
		return;
	int nhlt = 0, n = 0;
	double *PFG, *LFG, mFG = 0, tFG = 0, mean = 0, theta = 1e5;
	double lambda = 40.28 *(1 / (frequence[1] * frequence[1]) - 1 / (frequence[0] * frequence[0]));
	for (int i = s, u = 0; i < e; i++, u++)
	if (hlt[i])nhlt++;
	PFG = (double*)malloc(sizeof(double)*nhlt);
	LFG = (double*)malloc(sizeof(double)*nhlt);
	for (int i = s,u=0; i < e; i++)if(hlt[i]){
		PFG[u] = in[i][0] - in[i][3];
		LFG[u] = (in[i][1] / frequence[0] - in[i][1+NTYPE] / frequence[1])*CLIGHT;
		mFG += PFG[u] + LFG[u];
		tFG += (PFG[u] + LFG[u])*(PFG[u] + LFG[u]);
		u++;
	}
	n = nhlt;
	while (1){
		if (n == 0)
			break;
		mFG /= n;
		tFG = sqrt((tFG - n*mFG*mFG) / (n - 1));
		if (fabs(tFG - theta) < 1e-6)
			break;
		mean = mFG;
		theta = tFG; 
		n = 0; mFG = 0; tFG = 0;
		for (int i = s, u = 0; i < e; i++)if (hlt[i]){
			if (fabs(PFG[u] + LFG[u] - mean) < 2.*theta)
			{
				mFG += PFG[u] + LFG[u];
				tFG += (PFG[u] + LFG[u])*(PFG[u] + LFG[u]);
				n++;
			}
			u++;
		}
	}
	for (int i = s+1, u = 0; i < e; i++)if (hlt[i])
	{
		out[i][11] = (LFG[u] - mFG) / lambda *1e-16;
		u++;
	}
	free(PFG);
	free(LFG);
}
