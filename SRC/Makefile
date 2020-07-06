#!/bin/make

#FC=pgf90
#CPP=pgc++
#CC=pgcc
FC=gfortran
CPP=g++
CC=gcc
INSTALLPATH=../BIN
CLIB= -lm
NCINC= -I/usr/include 
NCLIB= -lnetcdf -L/usr/include 
#PTHEADLIB= -lpthread 
#RTKLIBFLAG= -DDLL -DENAGLO -DENAGAL -DENACMP -DNFREQ=5 
BLASLAPACK= -lblas -llapack 
GFLAG= -O4 -g
OPMFLAG=  
AR=ar rv

all:CRX2RNX DATETIMES POSCONVERT NC2TXT


prepare:
	@echo "preparing...";
ifeq ($(CC),pgcc)
	@echo "source ~/.bash_profile";
	@source ~/.bash_profile;
endif

SOFA.o: SOFA.h SOFA.c
	$(CC) -c SOFA.c;

CTIMES.o:CTIMES.c CTIMES.h SOFA.h
	$(CC) -c CTIMES.c;

XFORM.o:XFORM.c XFORM.h
	$(CC) -c XFORM.c;

CRX2RNX:CRX2RNX.c
	$(CC) CRX2RNX.c -o CRX2RNX  $(CLIB) $(GFLAG) $(OPMFLAG) $(CLIB);

DATETIMES:DATETIMES.c CTIMES.o  SOFA.o
	$(CC) DATETIMES.c -o DATETIMES CTIMES.o SOFA.o $(GFLAG) $(OPMFLAG) $(CLIB);

POSCONVERT:POSCONVERT.c CTIMES.o XFORM.o  SOFA.o
	$(CC) POSCONVERT.c -o POSCONVERT CTIMES.o XFORM.o SOFA.o $(GFLAG) $(OPMFLAG) $(CLIB);

NC2TXT:NC2TXT.c  CTIMES.o	
	$(CC) NC2TXT.c -o NC2TXT CTIMES.o XFORM.o SOFA.o \
	$(GFLAG) $(OPMFLAG) $(CLIB) $(NCINC) $(NCLIB);

clean:
	$(if *.o,rm -f *.o)
	$(if CRX2RNX,rm -f CRX2RNX)
	$(if DATETIMES,rm -f DATETIMES)
	$(if POSCONVERT,rm -f POSCONVERT)
	$(if NC2TXT,rm -f NC2TXT)

install:
	$(if $(INSTALLPATH),rm -rf $(INSTALLPATH))
	mkdir "$(INSTALLPATH)";
	cp -f GFZRNX     $(INSTALLPATH)/GFZRNX;
	cp -f TEQC     $(INSTALLPATH)/TEQC;
	cp -f CRX2RNX     $(INSTALLPATH)/CRX2RNX;
	cp -f DATETIMES     $(INSTALLPATH)/DATETIMES;
	cp -f POSCONVERT     $(INSTALLPATH)/POSCONVERT;
	cp -f NC2TXT     $(INSTALLPATH)/NC2TXT;

.PHONY:prepare
