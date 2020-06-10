#!/bin/bash
#Author :Mapoet
#Date   :2020.06.10
#Loc    :SHAO
date=`date`;
time=`echo $date|awk '{print substr($0,1,4),substr($0,7,2),substr($0,11,2),substr($0,19,2),substr($0,22,2),substr($0,25,2)}'`;
ymd=`echo $time|awk '{print $1,$2,$3}'`;
tyago=`echo $time|awk '{print $1-10,$2,$3}'`;
tyagos=`echo $time|awk '{print $1-10,$2,$3,$4,$5,$6}'`;
echo time: ${time};
echo ydoy:`echo ${time} |../BIN/DATETIMES 121`;
mjd=`echo ${time} |../BIN/DATETIMES 131`;
echo mjd:$mjd;
echo wds:`echo ${time} |../BIN/DATETIMES 141`;
echo ws:`echo ${time} |../BIN/DATETIMES 151`;
echo next day:`echo ${ymd} 1|../BIN/DATETIMES 100`;
echo next sec:`echo ${time} 1|../BIN/DATETIMES 101`;
echo before day:`echo ${ymd} 1|../BIN/DATETIMES 190`;
echo before sec:`echo ${time} 1|../BIN/DATETIMES 191`;
echo days between ${tyago}  and  ${ymd} :`echo ${ymd}  ${tyago} |../BIN/DATETIMES 180`;
echo secs between ${tyagos} and ${time} :`echo ${time}  ${tyagos} |../BIN/DATETIMES 181`;
echo mjd2ydmd $mjd:`echo $mjd |../BIN/DATETIMES 310`;

echo "test pos(34 112 50) in geo..."
echo "to geod in cart:" `echo 34 112 50|../BIN/POSCONVERT 13`;
echo "to geod in polar:" `echo 34 112 50|../BIN/POSCONVERT 12`;
echo "invert geo in cart to in elipise:" `echo -1982905.001    4907862.100   3546474.5234|../BIN/POSCONVERT 31`;
echo "to gei in cart:" `echo 34 112 50 |../BIN/POSCONVERT 13 "geo2gei:2019 334 00 00 00"`;
echo "invert gei in cart to geo in elipise:" `echo  -5291319.984    -144760.021   3546474.5234|../BIN/POSCONVERT 31 "gei2geo:2019 334 00 00 00"`;