#!/bin/bash
#Author :Mapoet
#Date   :2016.12.28
#Loc    :SHAO

BINPATH="../BIN";
WORKPATH="../WORK";
IGSPATH="../DATA/REALTIME";
if [  $# -ne 3 ];then
	echo "usage: ./sh_IGSDATA.sh time dday type;";
	echo "          for type is daily,hourly,highrate;";
	echo "       ./sh_IGSDATA.sh  mjd  dday  type;";
	echo "       ./sh_IGSDATA.sh \"year doy\"  dday  type;";
	echo "       ./sh_IGSDATA.sh \"year month day\" dday  type;";
	echo "       ./sh_IGSDATA.sh \"year doy h m s\" dday  type;";
	echo "       ./sh_IGSDATA.sh \"year month day h m s\" dday  type;";
	exit;
fi
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
if [ ! -e ${IGSPATH}/obsRnx ];then
    mkdir ${IGSPATH}/obsRnx;
fi

for NMJD in `seq -f%30.15f $SMJD 1 $EMJD`;do
	YEAR=` echo ${NMJD} |./DATETIMES 320 |awk '{print $1;}'`;
	YR=`   echo ${YEAR} | cut -c3-4`;
	DOY=`  echo ${NMJD} |./DATETIMES 320 |awk '{print $2;}'`;
	GW=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%4.4d",$1);}'`;
	GD=`   echo ${NMJD} |./DATETIMES 340 |awk '{printf("%1.1d",$2);}'`;
	if [ ${FILETYPE} == "highrate" ];then
	   #url="ftp://198.118.242.40/gnss/data/highrate/2016/345/16d/06/*";
	   RNXPATH="${IGSPATH}/obsRnx/${YEAR}.${DOY}.highrate";
	   if [ ! -e ${RNXPATH} ];then
	      mkdir ${RNXPATH};
	   fi
	   for H in `seq 0 1 23`;do
	      HH=`printf "%2.2d" ${H}`;
	      URL=ftp://198.118.242.40/gnss/data/highrate/${YEAR}/${DOY}/${YR}d/${HH}/*.Z;
	      wget -q ${URL};
	      #gunzip *.${YR}d.Z;
	      for FILE in `ls ????${DOY}?.${YR}d.Z`;do
	         gunzip ${FILE};
	      done
	      for FILE in `ls *.${YR}d`;do
	         cat ${FILE} | ./CRX2RNX > NEWFILE;
	         rm -f ${FILE};
	         mv -f NEWFILE ${RNXPATH}/${FILE%.*}.${YR}o ;
	      done
	      URL=ftp://198.118.242.40/gnss/data/highrate/${YEAR}/${DOY}/${YR}d/${HH}/*.crx.gz;
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
	         mv -f NEWFILE2 ${RNXPATH}/${NEWNAME}.${YR}o ;
	      done
	   done
      fi
	if [ ${FILETYPE} == "hourly" ];then
	   #url="ftp://198.118.242.40/gnss/data/hourly/2016/039/21/*";
	   RNXPATH="${IGSPATH}/obsRnx/${YEAR}.${DOY}.hourly";
	   if [ ! -e ${RNXPATH} ];then
	      mkdir ${RNXPATH};
	   fi
	   for H in `seq 0 1 23`;do
	      HH=`printf "%2.2d" ${H}`;
	      URL=ftp://198.118.242.40/gnss/data/hourly/${YEAR}/${DOY}/${HH}/*.${YR}d.Z;
	      wget -q ${URL};
	      #gunzip *.${YR}d.Z;
	      for FILE in `ls ????${DOY}?.${YR}d.Z`;do
	         gunzip ${FILE};
	      done
	      for FILE in `ls *.${YR}d`;do
	         cat ${FILE} | ./CRX2RNX > NEWFILE;
	         rm -f ${FILE};
	         mv -f NEWFILE ${RNXPATH}/${FILE%.*}.${YR}o ;
	      done
	      URL=ftp://198.118.242.40/gnss/data/hourly/${YEAR}/${DOY}/${HH}/*.crx.gz;
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
	         mv -f NEWFILE2 ${RNXPATH}/${NEWNAME}.${YR}o ;
	      done
	   done
      fi
	if [ ${FILETYPE} == "daily" ];then
	   #url="ftp://198.118.242.40/gnss/data/daily/2016/075/16d/*";
	   RNXPATH="${IGSPATH}/obsRnx/${YEAR}.${DOY}.daily";
	   if [ ! -e ${RNXPATH} ];then
	      mkdir ${RNXPATH};
	   fi
	   URL=ftp://198.118.242.40/gnss/data/daily/${YEAR}/${DOY}/${YR}d/*.Z;
	   wget -q ${URL};
	   #gunzip *.${YR}d.Z;
	   for FILE in `ls ????${DOY}?.${YR}d.Z`;do
	      gunzip ${FILE};
	   done
	   for FILE in `ls *.${YR}d`;do
	      cat ${FILE} | ./CRX2RNX > NEWFILE;
	      rm -f ${FILE};
	      mv -f NEWFILE ${RNXPATH}/${FILE%.*}.${YR}o ;
	   done
	   URL=ftp://198.118.242.40/gnss/data/daily/${YEAR}/${DOY}/${YR}d/*.crx.gz;
	   wget -q ${URL};
	   #gunzip *.crx.gz;
	   for FILE in `ls *.crx.gz`;do
	      gunzip ${FILE};
	   done
	   for FILE in `ls *.crx`;do
	      cat ${FILE} | ./CRX2RNX > NEWFILE;
	      rm -f ${FILE};
	         NEWNAME=`echo ${FILE}| awk \
	         '{printf("%s0\n",substr($1,1,7));}'`;
	      ./GFZRNX -finp NEWFILE -fout NEWFILE2 -vo 2;
	      rm -f NEWFILE;
	      mv -f NEWFILE2 ${RNXPATH}/${NEWNAME}.${YR}o ;
	   done
      fi
done

rm -f  DATETIMES;
rm -f  CRX2RNX;
rm -f  GFZRNX;





