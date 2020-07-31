#!/usr/bin/make
CC=gcc
FC=gfortran
CPP=g++

INSTALLPATH=..
RTKOBJ=${RTKSRC:%.c=%.o}
GFLAG= -O3
RTKFLAG= -DDLL -DENAGLO -DENAGAL -DENACMP -DNFREQ=3
INCS=
LIBS= -lblas -llapack -lpthread -lnetcdf -lm

AR=ar rv


RTKSRC=binex.c\
convkml.c\
convrnx.c\
crescent.c\
datum.c\
download.c\
ephemeris.c\
geoid.c\
gw10.c\
ionex.c\
javad.c\
lambda.c\
novatel.c\
nvs.c\
options.c\
pntpos.c\
postpos.c\
ppp_ar.c\
ppp.c\
preceph.c\
qzslex.c\
rcvlex.c\
rcvraw.c\
rinex.c\
rt17.c\
rtcm2.c\
rtcm3.c\
rtcm3e.c\
rtcm.c\
rtkcmn.c\
rtkpos.c\
rtksvr.c\
sbas.c\
septentrio.c\
skytraq.c\
solution.c\
ss2.c\
stream.c\
streamsvr.c\
tle.c\
ublox.c


all: CRX2RNX DATETIMES POSCONVERT NC2TXT librtk.a simobs simatmion obs2tec txrx2ipp


prepare:
	@echo "preparing...";
ifeq ($(CC),pgcc)
	@echo "source ~/.bash_profile";
	@source ~/.bash_profile;
endif

.c.o:$(RTKSRC)
	$(CC) -c $(GFLAG)  $(RTKFLAG) $(INCS) $<

SOFA.o: SOFA.h SOFA.c
	$(CC) -c SOFA.c;

CTIMES.o:CTIMES.c CTIMES.h SOFA.h
	$(CC) -c CTIMES.c;

XFORM.o:XFORM.c XFORM.h
	$(CC) -c XFORM.c;

CRX2RNX:CRX2RNX.c
	$(CC) CRX2RNX.c -o CRX2RNX  $(GFLAG) $(INCS) $(LIBS);

DATETIMES:DATETIMES.c CTIMES.o  SOFA.o
	$(CC) DATETIMES.c -o DATETIMES CTIMES.o SOFA.o $(GFLAG) $(INCS) $(LIBS);

POSCONVERT:POSCONVERT.c CTIMES.o XFORM.o  SOFA.o
	$(CC) POSCONVERT.c -o POSCONVERT CTIMES.o XFORM.o SOFA.o $(GFLAG) $(INCS) $(LIBS);

NC2TXT:NC2TXT.c  CTIMES.o	
	$(CC) NC2TXT.c -o NC2TXT CTIMES.o XFORM.o SOFA.o \
	$(GFLAG) $(INCS) $(LIBS);

librtk.a:$(RTKOBJ)
	$(AR) librtk.a $(RTKOBJ)

simobs:librtk.a simobs.o
	$(CC) simobs.o librtk.a $(GFLAG)  $(RTKFLAG) $(INCS) $(LIBS) -o simobs

simatmion:librtk.a simatmion.o
	$(CC) simatmion.o librtk.a $(GFLAG)  $(RTKFLAG) $(INCS) $(LIBS) -o simatmion

obs2tec:librtk.a obs2tec.o
	$(CC) obs2tec.o librtk.a $(GFLAG)  $(RTKFLAG) $(INCS) $(LIBS) -o obs2tec

txrx2ipp:librtk.a txrx2ipp.o
	$(CC) txrx2ipp.o librtk.a $(GFLAG)  $(RTKFLAG) $(INCS) $(LIBS) -o txrx2ipp



clean:
	$(if *.o,rm -f *.o)
	$(if CRX2RNX,rm -f CRX2RNX)
	$(if DATETIMES,rm -f DATETIMES)
	$(if POSCONVERT,rm -f POSCONVERT)
	$(if NC2TXT,rm -f NC2TXT)
	$(if librtk.a,rm -f librtk.a)
	$(if simobs,rm -f simobs)
	$(if simatmion,rm -f simatmion)
	$(if obs2tec,rm -f obs2tec)
	$(if txrx2ipp,rm -f txrx2ipp)

install:all
	$(if $(INSTALLPATH)/BIN,rm -rf $(INSTALLPATH)/BIN)
	$(if $(INSTALLPATH)/LIB,rm -rf $(INSTALLPATH)/LIB)
	$(if $(INSTALLPATH)/INC,rm -rf $(INSTALLPATH)/INC)
	mkdir "$(INSTALLPATH)/BIN";
	cp -f GFZRNX     $(INSTALLPATH)/BIN/GFZRNX;
	cp -f TEQC     $(INSTALLPATH)/BIN/TEQC;
	cp -f CRX2RNX     $(INSTALLPATH)/BIN/CRX2RNX;
	cp -f DATETIMES     $(INSTALLPATH)/BIN/DATETIMES;
	cp -f POSCONVERT     $(INSTALLPATH)/BIN/POSCONVERT;
	cp -f NC2TXT     $(INSTALLPATH)/BIN/NC2TXT;
	mkdir "$(INSTALLPATH)/INC";
	cp -f rtklib.h     $(INSTALLPATH)/INC/rtklib.h;
	cp -f SOFA.h     $(INSTALLPATH)/INC/SOFA.h;
	cp -f XFORM.h     $(INSTALLPATH)/INC/XFORM.h;
	mkdir "$(INSTALLPATH)/LIB";
	cp -f librtk.a     $(INSTALLPATH)/LIB/librtk.a;
	cp -f simobs     $(INSTALLPATH)/BIN/simobs;
	cp -f simatmion     $(INSTALLPATH)/BIN/simatmion;
	cp -f obs2tec     $(INSTALLPATH)/BIN/obs2tec;
	cp -f txrx2ipp     $(INSTALLPATH)/BIN/txrx2ipp;
	
.PHONY:prepare

