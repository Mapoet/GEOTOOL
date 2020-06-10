#!/bin/bash
#Author :Mapoet
#Date   :2020.06.10
#Loc    :SHAO
echo time: 2020 6 10 22 39 21;
echo ydoy:`echo 2020 6 10 22 39 21 |../BIN/DATETIMES 121`;
mjd=`echo 2020 6 10 22 39 21 |../BIN/DATETIMES 131`;
echo mjd:$mjd;
echo wds:`echo 2020 6 10 22 39 21 |../BIN/DATETIMES 141`;
echo ws:`echo 2020 6 10 22 39 21 |../BIN/DATETIMES 151`;
echo next day:`echo 2020 6 10 1|../BIN/DATETIMES 100`;
echo next sec:`echo 2020 6 10 22 39 21 1|../BIN/DATETIMES 101`;
echo before day:`echo 2020 6 10 1|../BIN/DATETIMES 190`;
echo before sec:`echo 2020 6 10 22 39 21 1|../BIN/DATETIMES 191`;
echo days between 2000 6 10 22 39 21 and 2020 6 10 22 39 21 :`echo 2020 6 10 2000 6 10 |../BIN/DATETIMES 180`;
echo secs between 2000 6 10  and 2020 6 10 :`echo 2020 6 10 22 39 21 2000 6 10 22 39 21 |../BIN/DATETIMES 181`;
echo mjd2ydmd $mjd:`echo $mjd |../BIN/DATETIMES 310`;
echo "test pos..."
echo -1982905.001    4907862.100   3546474.5234|../BIN/POSCONVERT 31;
echo 34 112 50|../BIN/POSCONVERT 13;
echo 34 112 50 |../BIN/POSCONVERT 13 "geo2gei:2019 334 00 00 00";
echo  -5291319.984    -144760.021   3546474.5234|../BIN/POSCONVERT 31 "gei2geo:2019 334 00 00 00";