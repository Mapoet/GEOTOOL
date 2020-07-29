#include "srcrtklib\rtklib.h"
#define  MAXOBSEPOCH   (24 * 60 * 60/24)
#define  MAXWINDOWSIZE      20
enum pfit{
	FIT_CSP = 1,
	FIT_Width,
};
void gettec(int sat, sta_t sta, nav_t*navsp3, gtime_t t[MAXOBSEPOCH], double in[MAXOBSEPOCH][4], double out[MAXOBSEPOCH][6], short hlt[MAXOBSEPOCH], int num, int w, double hion, int l2p);
void fitp(gtime_t t[MAXOBSEPOCH], double in[MAXOBSEPOCH][4], int s, int e, int w, double f1, double f2, double f3, int l2p);
void cacuspan(int sat, sta_t sta, nav_t*navsp3, gtime_t t[MAXOBSEPOCH], double in[MAXOBSEPOCH][4], double out[MAXOBSEPOCH][6], short hlt[MAXOBSEPOCH], int s, int e, int w, double f1, double f2, double f3, double hion, int l2p);
void writetofile(FILE*fid, int sat, sta_t sta, gtime_t t[MAXOBSEPOCH], double out[MAXOBSEPOCH][6], short hlt[MAXOBSEPOCH],int num);
int main(int argc, char**argv)
{
	char* frnx = "";
	char* fsp3 = "";
	FILE* ftec = stdout;
	int l2p = 0;
	int w = 10;
	double hion = 350e3;
	if (argc < 2||argc >7){
		fprintf(stderr, "need at least one parament as sp3 file.\n");
		return -1;
	}
	if (argc == 2){
		fsp3 = argv[1];
	}
	if (argc == 3){
		fsp3 = argv[1];
		frnx = argv[2];
	}
	if (argc == 4){
		fsp3 = argv[1];
		frnx = argv[2];
		if (!(ftec = fopen(argv[3], "w")))
			return -2;
	}
	if (argc == 5){
		fsp3 = argv[1];
		frnx = argv[2];
		if (!(ftec = fopen(argv[3], "w")))
			return -2;
		l2p = atoi(argv[4]);
	}
	if (argc == 6){
		fsp3 = argv[1];
		frnx = argv[2];
		if (!(ftec = fopen(argv[3], "w")))
			return -2;
		l2p = atoi(argv[4]);
		w = atoi(argv[5]);
	}
	if (argc == 7){
		fsp3 = argv[1];
		frnx = argv[2];
		if (!(ftec = fopen(argv[3], "w")))
			return -2;
		l2p = atoi(argv[4]);
		w = atoi(argv[5]);
		hion = atof(argv[6]);
	}
	obs_t obs = { 0 };
	nav_t navsp3 = { 0 };
	nav_t navobs = { 0 };
	sta_t sta = { 0 };
	gtime_t t[MAXOBSEPOCH];
	double in[MAXOBSEPOCH][2 + 2];
	double out[MAXOBSEPOCH][2 + 2 + 2];
	short hlt[MAXOBSEPOCH];
	int num = 0;
	int sys = 0;
	readsp3(fsp3, &navsp3, 1);
	if(navsp3.ne==0)
		return 1;
	if(readrnx(frnx, 1, "", &obs, &navobs, &sta)!=1)
		return 2;
	for (int i = 0; i < MAXSAT; i++)
	{
		num = 0;
		sys=satsys(i, NULL);
		if (sys != SYS_GPS&&sys != SYS_GLO)
			continue;
		for (int j = 0; j < obs.n; j++)
		{
			if (obs.data[j].sat == i)
			{
				t[num] = obs.data[j].time;
				in[num][0] = obs.data[j].P[0];
				in[num][1] = obs.data[j].P[1];
				in[num][2] = obs.data[j].L[0];
				in[num++][3] = obs.data[j].L[1];
			}
		}
		if (num){
			gettec(i, sta, &navsp3, t, in, out, hlt, num,w, hion, l2p);
			writetofile(ftec, i, sta, t, out, hlt, num);
		}
	}
	if(ftec!=stdout)
		fclose(ftec);
	return 0;
}
void gettec(int sat, sta_t sta, nav_t*navsp3, gtime_t t[MAXOBSEPOCH], double in[MAXOBSEPOCH][4], double out[MAXOBSEPOCH][6], short hlt[MAXOBSEPOCH], int num, int w, double hion, int l2p){
	char satn[4];
	satno2id(sat, satn);
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;
	switch (satn[0])
	{
	case 'C':f1 = FREQ1_CMP; f2 = FREQ2_CMP; f3 = FREQ3_CMP; break;
	case 'R':f1 = FREQ1_GLO; f2 = FREQ2_GLO; f3 = FREQ3_GLO; break;
	default:
		f1 = FREQ1; f2 = FREQ2, f3 = FREQ5; f4 = FREQ6; break;
	}
	double dt[MAXOBSEPOCH] = { 0 },mt=0,st=0;
	short span[MAXOBSEPOCH / 5] = { 0 }, nspan = 1;
	for (int i = 0; i < num - 1; i++)
	{
		dt[i] = (t[i + 1].time - t[i].time) - (t[i + 1].sec - t[i].sec);
		mt += dt[i] / (num - 1);
		st += dt[i] * dt[i];
	}
	st = sqrt((st-(num-1)*mt*mt)/(num-2));
	for (int i = 0; i< num - 1; i++)
	if (fabs(dt[i] - mt)>3 * mt&&fabs(dt[i] - mt)>3 * st)
		span[nspan++] = i + 1;
	span[nspan++] = num;
	for (int i = 0; i < nspan-1;i++)
		cacuspan(sat, sta, navsp3, t, in, out,hlt, span[i], span[i + 1],w, f1,f2,f3,hion,l2p);
}
void cacuspan(int sat, sta_t sta, nav_t*navsp3, gtime_t t[MAXOBSEPOCH], double in[MAXOBSEPOCH][4], double out[MAXOBSEPOCH][6], short hlt[MAXOBSEPOCH], int s, int e, int w, double f1, double f2, double f3, double hion, int l2p){
	double rx[6] = { 0 }, rt[2] = { 0 }, var = 0,
		dr[3] = { 0 }, dl = 0, pos[3] = { 0 },
		azel[2] = { 0 }, pppos[3] = { 0 }, maxmin[2] = { -HUGE_VALD, HUGE_VALD }, mtec = 0, stec = 0, ntec = 0, mtec0 = 0, stec0 = 0;
	int flag;
	fitp(t, in, s, e, w, f1, f2, f3, l2p);
	for(int i = s; i < e; i++)
	{
		if ((!peph2pos(t[i], sat, navsp3, 1, rx, rt, &var)) || fabs(in[i][0])<1e-3 || fabs(in[i][1])<1e-3
			|| fabs(in[i][2])<1e-3 || fabs(in[i][3])<1e-3)
		{
			memset(out[i], 0, sizeof(double)* 6);
			hlt[i] = 0;
			continue;
		}
		dl = 0;
		for (int j = 0; j < 3; j++){
			dr[j] = rx[j] - sta.pos[j];
			dl += dr[j] * dr[j];
		}
		dl = sqrt(dl);
		for (int j = 0; j < 3; j++)
			dr[j] /= dl;
		ecef2pos(sta.pos, pos);
		satazel(pos, dr, azel);
		ionppp(pos, azel, RE_WGS84, hion, pppos);
		out[i][0] = pppos[0] * R2D;
		out[i][1] = pppos[1] * R2D;
		out[i][2] = azel[0] * R2D;
		out[i][3] = azel[1] * R2D;
		out[i][4] = ionmapf(pos,azel);
		out[i][5] = 1e-16 / 40.28*(f1*f1*f2*f2) / (f1*f1 - f2*f2)*(in[i][0] - in[i][1]);//-9.5196*(in[i][2]-1.28333*in[i][3]);
		hlt[i] = 1;
		maxmin[0] = out[i][5] > maxmin[0] ? out[i][5] : maxmin[0];
		maxmin[1] = out[i][5] < maxmin[1] ? out[i][5] : maxmin[1];
		mtec += out[i][5];
		stec += out[i][5] * out[i][5];
		ntec++;
	}
	if (fabs(maxmin[0]-maxmin[1]) < 3e2)
	{
		if (maxmin[1]>0)
			return;
		for (int i = s; i < e; i++)
		if (hlt[i])
		{
			out[i][5] -= maxmin[1];
		}
		return;
	}
	while (1)
	{
		mtec0 = mtec/ntec;
		stec0 = sqrt((stec - ntec*mtec*mtec) / (ntec - 1));
		mtec = 0;stec = 0;ntec = 0;flag = 0;
		for (int i = s; i < e; i++)
		if (hlt[i])
		{
			if (fabs(out[i][5] - mtec0)>3 * stec0)
			{
				if (i == s&&fabs(out[i + 1][5] - mtec0) < 3 * stec0)
					out[i][5] = out[i + 1][5];
				else if (i == e&&fabs(out[i - 1][5] - mtec0) < 3 * stec0)
					out[i][5] = out[i - 1][5];
				else if (fabs(out[i - 1][5] - mtec0) < 3 * stec0&&fabs(out[i + 1][5] - mtec0) < 3 * stec0)
					out[i][5] = (out[i - 1][5] + out[i - 1][5]) / 2;
				else
				{
					hlt[i] = 0;
					flag = 1;
					continue;
				}
			}
			mtec += out[i][5];
			stec += out[i][5] * out[i][5];
			ntec++;
		}
		if (!flag)
			break;
	}
}
void writetofile(FILE*fid, int sat, sta_t sta, gtime_t t[MAXOBSEPOCH], double out[MAXOBSEPOCH][6], short hlt[MAXOBSEPOCH], int num){
	char satn[4],stan[5];
	satno2id(sat, satn);
	memcpy(stan, sta.name, sizeof(char)* 4);
	stan[4] = '\0';
	double ep[6],ep0[6];
	time2epoch(t[0], ep0);
	for (int i = 0; i < num; i++)
	if(hlt[i]){
		time2epoch(t[i], ep);
		//fprintf(fid, "%4d%3d%3d%3d%3d%7.3lf%5s%5s%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf\n",
		//	(int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5], satn, sta.name, out[i][0], out[i][1], out[i][2], out[i][3], out[i][4], out[i][5]);
		fprintf(fid, "%15.3lf%5s%5s%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf%12.3lf\n",
			ep[5] + 60 * (ep[4] + 60 * (ep[3] + 24 * (ep[2] - ep0[2]))), satn, stan, out[i][0], out[i][1], out[i][2], out[i][3], out[i][4], out[i][5]);
	}
}
void fitp(gtime_t t[MAXOBSEPOCH], double in[MAXOBSEPOCH][4], int s, int e, int w, double f1, double f2, double f3, int l2p)
{
	if (!l2p)
		return;
	double p[MAXOBSEPOCH], pm, pn, pt;
	int mark[MAXWINDOWSIZE] = { 0 }, flag = 0;
	if (l2p == 1){
		for (int i = s + 1; i < e; i++)
		for (int k = 0; k<2; k++)
			in[i][k] = in[i][k] / i + (i - 1) / i*(in[i - 1][k] + CLIGHT*(in[i][k + 2] / f1 - in[i - 1][k + 2] / f2));
		return;
	}
	if (l2p == 2)
	{
		for (int k = 0; k < 4; k++)
		{
			for (int i = s; i < e; i++)
				p[i] = in[i][k];
			for (int i = s + w / 2; i < e - w / 2.; i++)
			{
				memset(mark, 0, sizeof(int)* MAXWINDOWSIZE);
				while (1){
					pm = 0;
					pn = 0;
					pt = 0;
					flag = 0;
					for (int c = i - w / 2, u = 0; c < i + w / 2.; c++,u++)
					if (mark[u] == 0){ 
						pm += p[c]; 
						pt += p[c] * p[c];
						pn++;
					}
					pm /= pn;
					pt = sqrt((pt - pn*pm*pm) / (pn - 1));
					for (int c = i - w / 2, u = 0; c < i + w / 2; c++, u++)
					if (mark[u] == 0 && fabs(p[i] - pm)>3 * pt)
					{
						mark[u] = 1;
						flag = 1;
					}
					if (!flag)
						break;
				}
				in[i][k] = pm;
			}
		}
		return;
	}
}
