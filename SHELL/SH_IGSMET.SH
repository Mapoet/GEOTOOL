#!/bin/bash
#Author :Mapoet
#Date   :2016.12.28
#Loc    :SHAO

BINPATH="../BIN";
WORKPATH="../TEMP";
IGSPATH="../DATA/REALTIME/igs";
if [  $# -ne 3 ];then
	echo "usage: ./SH_IGSMET.SH time dday type;";
	echo "          for type is daily,hourly,highrate;";
	echo "       ./SH_IGSMET.SH  mjd  dday  type;";
	echo "       ./SH_IGSMET.SH \"year doy\"  dday  type;";
	echo "       ./SH_IGSMET.SH \"year month day\" dday  type;";
	echo "       ./SH_IGSMET.SH \"year doy h m s\" dday  type;";
	echo "       ./SH_IGSMET.SH \"year month day h m s\" dday  type;";
	exit;
fi
WORKPATH="${WORKPATH}MET${RANDOM}";
mkdir ${WORKPATH};
cd ${WORKPATH};
ln -s ${BINPATH}/DATETIMES DATETIMES;
ln -s ${BINPATH}/CRX2RNX CRX2RNX;
ln -s ${BINPATH}/GFZRNX GFZRNX;

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
if [ ! -e ${IGSPATH}/metRnx ];then
    mkdir ${IGSPATH}/metRnx;
fi

for NMJD in `seq $SMJD 1 $EMJD`;do
	YEAR=` echo ${NMJD} |./DATETIMES 320 |awk '{print $1;}'`;
	YR=`   echo ${YEAR} | cut -c3-4`;
	DOY=`  echo ${NMJD} |./DATETIMES 320 |awk '{print $2;}'`;
	GW=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%4.4d",$1);}'`;
	GD=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%1.1d",$2);}'`;
	if [ ${FILETYPE} == "highrate" ];then
	   #url="ftp://198.118.242.40/gnss/data/highrate/2016/345/16m/06/*";
	   METPATH="${IGSPATH}/metRnx/${YEAR}.${DOY}.highrate";
	   if [ ! -e ${METPATH} ];then
	      mkdir ${METPATH};
	   fi
	   for H in `seq 0 1 23`;do
	      HH=`printf "%2.2d" ${H}`;
	      URL=ftp://198.118.242.40/gnss/data/highrate/${YEAR}/${DOY}/${YR}m/${HH}/*.Z;
	      wget -q ${URL};
	      #gunzip *.${YR}d.Z;
	      for FILE in `ls ????${DOY}?.${YR}d.Z`;do
	         gunzip ${FILE};
	      done
	      for FILE in `ls *.${YR}d`;do
	         cat ${FILE} | ./CRX2RNX > NEWFILE;
	         rm -f ${FILE};
	         mv -f NEWFILE ${METPATH}/${FILE%.*}.${YR}o ;
	      done
	      URL=ftp://198.118.242.40/gnss/data/highrate/${YEAR}/${DOY}/${YR}m/${HH}/*.crx.gz;
	      wget -q ${URL};
	      for FILE in `ls *.crx.gz`;do
	         gunzip ${FILE};
	      done
	      for FILE in `ls *.crx`;do
	         cat ${FILE} | ./CRX2RNX > NEWFILE;
	         rm -f ${FILE};
	         NEWNAME=`echo ${FILE} ${H}| awk \
	         '{asic="abcdefghijklmnopqrstuvwxyz";\
	           printf("%s%c%s\n",substr($1,1,7),substr(asic,$2,$2+1),substr($1,22,2));}'`;
	         ./GFZRNX -finp NEWFILE -fout NEWFILE2 -vo 2;
	         rm -f NEWFILE;
	         mv -f NEWFILE2 ${METPATH}/${NEWNAME}.${YR}o ;
	      done
	   done
      fi
	if [ ${FILETYPE} == "hourly" ];then
	   #url="ftp://198.118.242.40/gnss/data/hourly/2016/039/21/*";
	   METPATH="${IGSPATH}/metRnx/${YEAR}.${DOY}.hourly";
	   if [ ! -e ${METPATH} ];then
	      mkdir ${METPATH};
	   fi
	   for H in `seq 0 1 23`;do
	      HH=`printf "%2.2d" ${H}`;
	      URL=ftp://198.118.242.40/gnss/data/hourly/${YEAR}/${DOY}/${HH}/*.${YR}m.Z;
	      wget -q ${URL};
	      #gunzip *.${YR}d.Z;
	      for FILE in `ls ????${DOY}?.${YR}m.Z`;do
	         gunzip ${FILE};
	      done
	      for FILE in `ls *.${YR}m`;do
	         mv -f ${FILE} ${METPATH}/${FILE%.*}.${YR}m;
	      done
	      URL=ftp://198.118.242.40/gnss/data/hourly/${YEAR}/${DOY}/${HH}/*.rnx.gz;
	      wget -q ${URL};
	      for FILE in `ls *.rnx.gz`;do
	         gunzip ${FILE};
	      done
	      for FILE in `ls *.rnx`;do
	         NEWNAME=`echo ${FILE} ${H}| awk \
	         '{asic="abcdefghijklmnopqrstuvwxyz";\
	           printf("%s%c%s\n",substr($1,1,7),substr(asic,$2,$2+1),substr($1,22,2));}'`;
	         mv -f ${FILE} ${METPATH}/${NEWNAME}.${YR}m;
	      done
	   done
      fi
	if [ ${FILETYPE} == "daily" ];then
	   #url="ftp://198.118.242.40/gnss/data/daily/2016/075/16m/*";
	   METPATH="${IGSPATH}/metRnx/${YEAR}.${DOY}.daily";
	   if [ ! -e ${METPATH} ];then
	      mkdir ${METPATH};
	   fi
	   URL=ftp://198.118.242.40/gnss/data/daily/${YEAR}/${DOY}/${YR}m/*.Z;
	   wget -q ${URL};
	   #gunzip *.${YR}d.Z;
	   for FILE in `ls ????${DOY}?.${YR}m.Z`;do
	      gunzip ${FILE};
	   done
	   for FILE in `ls *.${YR}m`;do
	      mv -f ${FILE} ${METPATH}/${FILE%.*}.${YR}m ;
	   done
	   URL=ftp://198.118.242.40/gnss/data/daily/${YEAR}/${DOY}/${YR}m/*.rnx.gz;
	   wget -q ${URL};
	   for FILE in `ls *.rnx.gz`;do
	      gunzip ${FILE};
	   done
	   for FILE in `ls *.rnx`;do
	      NEWNAME=`echo ${FILE}| awk \
	         '{printf("%s0\n",substr($1,1,7));}'`;
	      mv -f ${FILE} ${METPATH}/${NEWNAME}.${YR}m;
	   done
      fi
done

rm -f  DATETIMES;
rm -f  CRX2RNX;
rm -f  GFZRNX;
cd ../SHELL;
rm -rf ${WORKPATH};





