#!/bin/bash
#Author :Mapoet
#Date   :2016.12.28
#Loc    :SHAO

BINPATH="../BIN";
WORKPATH="../TEMP";
COSMICPATH="../DATA/REALTIME/cosmic";
if [  $# -ne 3 ];then
	echo "usage: ./sh_COSMICDATA.sh time dday type;";
	echo "          for type is ionPrf,ionPhs,podTec;";
	echo "       ./SH_COSMICDATA.SH  mjd  dday  type;";
	echo "       ./SH_COSMICDATA.SH \"year doy\"  dday  type;";
	echo "       ./SH_COSMICDATA.SH \"year month day\" dday  type;";
	echo "       ./SH_COSMICDATA.SH \"year doy h m s\" dday  type;";
	echo "       ./SH_COSMICDATA.SH \"year month day h m s\" dday  type;";
	exit;
fi
WORKPATH="${WORKPATH}${3}${RANDOM}";
mkdir ${WORKPATH};
cd ${WORKPATH};
ln -sf ${BINPATH}/DATETIMES DATETIMES;

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
if [ ! -e ${COSMICPATH} ];then
    mkdir ${COSMICPATH};
fi
SRE=`echo 2006 112 |./DATETIMES 230`; 
ERE=`echo 2014 120 |./DATETIMES 230`; 
SPO=`echo 2014 121 |./DATETIMES 230`; 
EPO=`echo 2020 069 |./DATETIMES 230`;
SRT=`echo 2019 001 |./DATETIMES 230`; 
ERT=`echo 2020 112 |./DATETIMES 230`;
WGETAUTHORITY=" --http-user=guopeng --http-password=760807";
CURLAUTHORITY=" -k -s -u guopeng:760807";
# HOSTURL2="https://data.cosmic.ucar.edu/gnss-ro/cosmic2/provisional/spaceWeather/level1b/${year}";
# ${CDAAC}/${day}/podTc2_prov1_${daydir}.tar.gz;
HOSTURL="https://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data";
for NMJD in `seq ${SMJD} 1 ${EMJD}`;do
   FILE=`echo ${NMJD} |./DATETIMES 320| awk '{printf("%4.4d.%3.3d",$1,$2);}' `;
   PROJECTTYPE=cosmic2013;
   if [ `echo ${NMJD} ${SPO}|awk '{if($1>$2)print "1";else print "0"; }' ` == "1" ] ;then
      PROJECTTYPE=cosmic;
   fi
   if [ `echo ${NMJD} ${SRT}|awk '{if($1>$2)print "1";else print "0"; }' ` == "1" ] ;then
      PROJECTTYPE=cosmicrt;
   fi
   URL="${HOSTURL}/${PROJECTTYPE}/${FILETYPE}/${FILE}";
   #wget -q  ${AUTHORITY} ${URL};
   curl  ${CURLAUTHORITY} ${URL} -o ${FILE}.tar;
   tar x -f ${FILE}.tar;
   cp -rf ${PROJECTTYPE}/* ${COSMICPATH};
   rm -rf ${PROJECTTYPE};
   rm -rf ${FILE} ${FILE}.tar;
done

rm -f  DATETIMES;
cd ../SHELL;
rm -rf ${WORKPATH};
