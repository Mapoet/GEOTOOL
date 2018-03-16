#!/bin/bash
#Author :Mapoet
#Date   :2016.12.28
#Loc    :SHAO

BINPATH="../BIN";
WORKPATH="../WORK";
COSMICPATH="../DATA/REALTIME";
if [  $# -ne 3 ];then
	echo "usage: ./sh_COSMICDATA.sh time dday type;";
	echo "          for type is ionPrf,ionPhs,podTec;";
	echo "       ./sh_COSMICDATA.sh  mjd  dday  type;";
	echo "       ./sh_COSMICDATA.sh \"year doy\"  dday  type;";
	echo "       ./sh_COSMICDATA.sh \"year month day\" dday  type;";
	echo "       ./sh_COSMICDATA.sh \"year doy h m s\" dday  type;";
	echo "       ./sh_COSMICDATA.sh \"year month day h m s\" dday  type;";
	exit;
fi
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
if [ ! -e ${COSMICPATH} ];then
    mkdir ${COSMICPATH};
fi
SRE=`echo 2006 112 |./DATETIMES 230`; 
ERE=`echo 2014 120 |./DATETIMES 230`; 
SPO=`echo 2014 121 |./DATETIMES 230`; 
EPO=`echo 2017 180 |./DATETIMES 230`;
SRT=`echo 2017 001 |./DATETIMES 230`; 
ERT=`echo 2017 250 |./DATETIMES 230`;
AYTHORITY=" --http-user=guopeng --http-password=760807";
HOSTURL="http://cdaac-www.cosmic.ucar.edu/cdaac/rest/tarservice/data";
for NMJD in `seq -f%30.15f ${SMJD} 1 ${EMJD}`;do
   FILE=`echo ${NMJD} |./DATETIMES 320| awk '{printf("%4.4d.%3.3d",$1,$2);}' `;
   PROJECTTYPE=cosmic2013;
   if [ `echo "${NMJD} > ${SPO}" | bc ` == "1" ] ;then
      PROJECTTYPE=cosmic;
   fi
   if [ `echo "${NMJD} > ${SRT}" | bc ` == "1" ] ;then
      PROJECTTYPE=cosmicrt;
   fi
   URL="${HOSTURL}/${PROJECTTYPE}/${FILETYPE}/${FILE}";
   wget  -q ${AYTHORITY} ${URL};
   tar x -f ${FILE};
   cp -rf ${PROJECTTYPE}/* ${COSMICPATH};
   rm -rf ${PROJECTTYPE};
   rm -rf ${FILE};
done

rm -f  DATETIMES;
