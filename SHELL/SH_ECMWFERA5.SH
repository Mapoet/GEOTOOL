#!/bin/bash
#Author :Mapoet
#Date   :2020.05.28
#Loc    :SHAO

BINPATH="../BIN";
WORKPATH="../TEMP";
IGSPATH="../DATA/REALTIME";
if [  $# -ne 2 ];then
	echo "usage: ./SH_ECMWFERA5.SH time dday;";
	echo "       ./SH_ECMWFERA5.SH  mjd  dday;";
	echo "       ./SH_ECMWFERA5.SH \"year doy\"  dday ;";
	echo "       ./SH_ECMWFERA5.SH \"year month day\" dday ;";
	echo "       ./SH_ECMWFERA5.SH \"year doy h m s\" dday ;";
	echo "       ./SH_ECMWFERA5.SH \"year month day h m s\" dday ;";
	exit;
fi
WORKPATH="${WORKPATH}ERA5${RANDOM}";
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
FILETYPE="era5";
if [ ! -e ${IGSPATH} ];then
    mkdir ${IGSPATH};
fi
if [ ! -e ${IGSPATH}/${FILETYPE} ];then
    mkdir ${IGSPATH}/${FILETYPE};
fi

for NMJD in `seq $SMJD 1  $EMJD`;do
	Y=`echo $NMJD|./DATETIMES 311|awk '{print $1;}'`;
	M=`echo $NMJD|./DATETIMES 311|awk '{print $2;}'`;
	D=`echo $NMJD|./DATETIMES 311|awk '{print $3;}'`;
	ERA5NAME=`echo $NMJD|./DATETIMES 311|awk '{printf("era5-%4.4d-%2.2d-%2.2d.nc",$1,$2,$3)}'`;
      if [ ! -e ${IGSPATH}/${FILETYPE}/${ERA5NAME} ]||[ `ls -l ${IGSPATH}/${FILETYPE}/${ERA5NAME}|awk '{print $5}'` -lt 2769714732 ];then
	   python  << EOF &
#!/usr/bin/env python
import cdsapi
import sys
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'geopotential', 'relative_humidity', 'specific_humidity', 'specific_cloud_ice_water_content',
            'specific_cloud_liquid_water_content',
            'temperature',
        ],
        'grid':'0.5/0.5',
        'pressure_level': [
            '1', '2', '3',
            '5', '7', '10',
            '20', '30', '50',
            '70', '100', '125',
            '150', '175', '200',
            '225', '250', '300',
            '350', '400', '450',
            '500', '550', '600',
            '650', '700', '750',
            '775', '800', '825',
            '850', '875', '900',
            '925', '950', '975',
            '1000',
        ],
        'year': '%s' % (${Y}),
        'month': '%s' % (${M}),
        'day': '%s' % (${D}),
        'time': [
            '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
        ],
    },
    '%s' % ("${IGSPATH}/${FILETYPE}/${ERA5NAME}"))
EOF
      fi
done
rm -f  DATETIMES;
cd ../SHELL;
rm -rf ${WORKPATH};

