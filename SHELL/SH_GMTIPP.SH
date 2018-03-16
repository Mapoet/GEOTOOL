#!/bin/bash
SLAT=5;
ELAT=55;
SLON=75;
ELON=145;
NYEAR=2016;
NDOY=141;
B=10g10/5g5;
cd ../SIMCWORK/;
pscoast -R${SLON}/${ELON}/${SLAT}/${ELAT} \
-JN110/15c -B${B} -A250 -Glightgray -Wthinnest -K -X8c -Y3c > base.ps;
for IH in `seq 0 23`;do
   TSTR=`printf "%4.4d-%3.3d-%2.2d" ${NYEAR} ${NDOY} ${IH}`;
   PSFN=`printf "Distr%3.3d%2.2d.ps" ${NDOY} ${IH}`;
   cp base.ps ${PSFN};
   psxy -R${SLON}/${ELON}/${SLAT}/${ELAT} \
-JN110/15c -B${B}:."Distr. of IPP and Fixed NE at ${TSTR}": -T -K -O >>${PSFN};
   OBSFN=`printf "SIM%3.3d0_%2.2d00.SIMI" ${NDOY} ${IH}`;
   COSFN=`printf "COS%3.3d0_%2.2d00.SIMI" ${NDOY} ${IH}`;
   FNEFN=`printf "COS%3.3d0_%2.2d00.FNES" ${NDOY} ${IH}`;
   cat ${OBSFN} | awk '{print$(14),$(13);}'|psxy -R -J -B -K -O \
   -St0.1c -Gyellow >>${PSFN};
   cat ${COSFN} | awk '{print$(14),$(13);}'|psxy -R -J -B -K -O \
   -Sa0.15c -Ggreen >>${PSFN};
   cat ${FNEFN} | awk '{print$(2),$(1);}'|psxy -R -J -B -K -O \
   -Ss0.2c -Gred >>${PSFN};
   echo -12 0 | psxy -R-15/15/-5/5 -JX15c/1c -Bwnse -K -O -St0.1c -Gyellow -Y-2c >>${PSFN};
   echo  -2 0 | psxy -R -J -B -K -O -Sa0.15c -Ggreen >>${PSFN};
   echo   8 0 | psxy -R -J -B -K -O -Ss0.2c -Gred >>${PSFN};
   echo -11 0 12 0 0 LM GNSS   | pstext -R -J -B -K -O -W  >>${PSFN};
   echo  -1 0 12 0 0 LM COSMIC | pstext -R -J -B -K -O -W  >>${PSFN};
   echo   9 0 12 0 0 LM IONPRF | pstext -R -J -B -K -O -W >>${PSFN};
   psxy -R -J -B -T -O >>${PSFN};
   ps2raster -Tt -E180 ${PSFN};
   rm ${PSFN};
done;
rm base.ps;
