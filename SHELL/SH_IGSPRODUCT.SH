#!/bin/bash
#Author :Mapoet
#Date   :2016.12.28
#Loc    :SHAO

BINPATH="../BIN";
WORKPATH="../TEMP";
IGSPATH="../DATA/REALTIME/igs";
if [  $# -ne 3 ];then
	echo "usage: ./SH_IGSPRODUCT.SH time dday type;";
	echo "          for type is ionMap,zpdRnx,sp3Eph,dcbRnx,erpRnx,clkRnx,resRnx,sscRnx,snxRnx;";
	echo "       ./SH_IGSPRODUCT.SH  mjd  dday  type;";
	echo "       ./SH_IGSPRODUCT.SH \"year doy\"  dday  type;";
	echo "       ./SH_IGSPRODUCT.SH \"year month day\" dday  type;";
	echo "       ./SH_IGSPRODUCT.SH \"year doy h m s\" dday  type;";
	echo "       ./SH_IGSPRODUCT.SH \"year month day h m s\" dday  type;";
	exit;
fi
WORKPATH="${WORKPATH}${3}${RANDOM}";
mkdir ${WORKPATH};
cd ${WORKPATH};
ln -s ${BINPATH}/DATETIMES DATETIMES;

DATEOPT=`echo $1 | awk '{if(NF==1)print "3";\
                         if(NF==2)print "230";\
                         if(NF==3)print "130";\
                         if(NF==5)print "231";\
                         if(NF==6)print "131";}'`;
if [ ! "${DATEOPT}" == "3" ];then
   SMJD=`echo $1 | ./DATETIMES "${DATEOPT}"`;
else
   SMJD=$1;
fi
EMJD=`echo ${SMJD} $2 |./DATETIMES 300 `;
FILETYPE="$3";
if [ ! -e ${IGSPATH} ];then
    mkdir ${IGSPATH};
fi
if [ ! -e ${IGSPATH}/${FILETYPE} ];then
    mkdir ${IGSPATH}/${FILETYPE};
fi

for NMJD in `seq $SMJD 1 $EMJD`;do
	YEAR=` echo ${NMJD} |./DATETIMES 320 |awk '{print $1;}'`;
	YR=`   echo ${YEAR} | cut -c3-4`;
	DOY=`  echo ${NMJD} |./DATETIMES 320 |awk '{print $2;}'`;
	GW=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%4.4d",$1);}'`;
	GD=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%1.1d",$2);}'`;
	if [ ${FILETYPE} == "ionMap" ];then
	   #url="ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/2007/005/igsg0050.07i.Z"
	   FN="codg${DOY}0.${YR}i";
	   URL="ftp://198.118.242.40/gnss/products/ionex/${YEAR}/${DOY}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="igsg${DOY}0.${YR}i";
	   URL="ftp://198.118.242.40/gnss/products/ionex/${YEAR}/${DOY}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
	if [ ${FILETYPE} == "zpdRnx" ];then
	   #url="ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/2007/005/igsg0050.07i.Z"
         mkdir "${IGSPATH}/${FILETYPE}/${YEAR}.${DOY}";
	   FN="????${DOY}0.${YR}zpd";
	   URL="ftp://198.118.242.40/gnss/products/troposphere/zpd/${YEAR}/${DOY}/${FN}.gz";
	   wget -q ${URL};
	   for FILE in `ls ${FN}.gz`;do
         gunzip "${FILE}";
         mv "${FILE%.*}" "${IGSPATH}/${FILETYPE}/${YEAR}.${DOY}/${FILE%.*}";
         done
      fi
	if [ ${FILETYPE} == "clkRnx" ];then
	   FN="cod${GW}${GD}.clk";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="igs${GW}${GD}.clk";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
	if [ ${FILETYPE} == "resRnx" ];then
	   FN="cod${GW}${GD}.res";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="igs${YR}P${GW}${GD}.res";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
	if [ ${FILETYPE} == "sscRnx" ];then
	   FN="cod${GW}${GD}.ssc";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="igs${YR}P${GW}${GD}.ssc";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
	if [ ${FILETYPE} == "snxRnx" ];then
	   FN="cod${GW}${GD}.snx";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="igs${YR}P${GW}${GD}.snx";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
	if [ ${FILETYPE} == "sp3Eph" ];then
	   #url="ftp://cddis.gsfc.nasa.gov/gnss/products/1914/*";
	   FN="cod${GW}${GD}.eph";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget  -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="igs${GW}${GD}.sp3";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget  -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
done
for NMJD in `seq $SMJD 7 $EMJD`;do
	YEAR=` echo ${NMJD} |./DATETIMES 320 |awk '{print $1;}'`;
	GW=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%4.4d",$1);}'`;
	GD=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%1.1d",$2);}'`;
	if [ ${FILETYPE} == "erpRnx" ];then
	   FN="igs"${GW}"7.erp";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="igs"${GW}"7.sum";
	   URL="ftp://198.118.242.40/gnss/products/${GW}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
done
for NMJD in `seq $SMJD 30 $EMJD`;do
	YEAR=` echo ${NMJD} |./DATETIMES 320 |awk '{print $1;}'`;
	YR=`   echo ${YEAR} | cut -c3-4`;
	MT=`  echo ${NMJD} |./DATETIMES 310 |awk '{print $2;}'`;
	if [ ${FILETYPE} == "dcbRnx" ];then
	   FN="P1C1${YR}${MT}.DCB";
	   URL="ftp://ftp.aiub.unibe.ch/CODE/${YEAR}/${FN}.Z";
	   wget -q ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="P1C1"${YR}${MT}"_RINEX.DCB";
	   URL="ftp://ftp.aiub.unibe.ch/CODE/${YEAR}/${FN}.Z";
	   wget  ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="P1P2${YR}${MT}.DCB";
	   URL="ftp://ftp.aiub.unibe.ch/CODE/${YEAR}/${FN}.Z";
	   wget  ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="P1P2"${YR}${MT}"_ALL.DCB";
	   URL="ftp://ftp.aiub.unibe.ch/CODE/${YEAR}/${FN}.Z";
	   wget  ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
	   FN="P2C2"${YR}${MT}"_RINEX.DCB";
	   URL="ftp://ftp.aiub.unibe.ch/CODE/${YEAR}/${FN}.Z";
	   wget  ${URL};
	   gunzip "${FN}.Z";
         mv "${FN}" "${IGSPATH}/${FILETYPE}/${FN}";
      fi
done
rm -f  DATETIMES;
cd ../SHELL;
rm -rf ${WORKPATH};





