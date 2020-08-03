
/*
Author    :  Mapoet Niphy
Date      :  2020
Institude :  SHAO
Name      :  TXRX2IPP
Project   :  GEOTOOL

  Created by Mapoet Niphy on 2020/06/10.
  Copyright © 2020年 Mapoet Niphy. All rights reserved.

*/
#include "rtklib.h"
int main(int argc,char** argv){
	char line[800];
	double Tx[3],Px[3], Rx[3],Dx[3],Nx[3],azel[2],ipp[3],hion=HION;
	if (argc >= 2){
		sscanf(argv[1],"%lf", &hion);// hion *= 1e3;
	}
	while (~feof(stdin))
	{
		memset(line, 0, 800);
		if (fgets(line, 800, stdin) == EOF)
			break;
		if (strlen(line) == 0)
			break;
		//if (sscanf(line, "%lf%lf%lf%lf%lf%lf", Tx, Tx + 1, Tx + 2, Rx, Rx + 1, Rx + 2) != 6)
		//	continue;
		if (sscanf(line, "%lf%lf%lf%lf%lf%lf", Rx, Rx + 1, Rx + 2, Tx, Tx + 1, Tx + 2) != 6)
			continue;
		//for (int i = 0; i < 3; i++)Tx[i] *= 1e3, Rx[i] *= 1e3;
		for (int i = 0; i < 3; i++)Tx[i] *= 1e3;
		Dx[0] = Tx[0] - Rx[0];
		Dx[1] = Tx[1] - Rx[1];
		Dx[2] = Tx[2] - Rx[2];
		normv3(Dx,Nx);
		ecef2pos(Rx, Px);
		satazel(Px, Nx, azel);
		ionppp(Px, azel, RE_WGS84/1e3, hion, ipp);
		ipp[2] = azel[1] * R2D;// Px[2];
		fprintf(stdout, "%20.6lf %20.6lf %20.6lf\n", ipp[0] * R2D, ipp[1] * R2D, ipp[2]);
	}
	return 0;
}
