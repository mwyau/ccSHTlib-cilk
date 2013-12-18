#  Define the C compiler for serial code:
CC = icc

#  Define the C compiler for parallel code:
MPCC = mpicc

#  Define the directory which contains the 
#  FFTW header files:
FFTWINCLUDE = ../fftw/include

#  Define the directory which contains the 
#  FFTW library files:
FFTWLIBRARY = ../fftw/lib

#  If FORTRAN adds an underscore to the object names
#  or uses upper case object names one or both of 
#  the following should be added to the CFLAGS and 
#  MPCFLAGS definitions: 
#
#  -DF77_APPEND_UNDERSCORE
#  -DF77_UPPERCASE_NAMES 

#  Define the flags to be passed to the 
#  serial C complier.  
CFLAGS = -O3

#  Define the flags to be passed to the 
#  parallel C compiler.  
MPCFLAGS = -O3

##################################################
#    DO NOT CHANGE MAKEFILE BELOW THIS POINT     #
##################################################

.PHONY : all serial library executable clean tidy

all : libccSHT.a ccSHTx libccSHTmpi.a ccSHTmpix libccSHTcilk.a ccSHTcilkx tidy

serial : libccSHT.a ccSHTx tidy

library : libccSHT.a libccSHTmpi.a tidy

executable : ccSHTx ccSHTmpix tidy

clean : 
	rm -f lib/libccSHT.a lib/libccSHTmpi.a bin/ccSHTx bin/ccSHTmpix *.o makeSizeType include/sizedType.h

tidy :
	rm -f *.o makeSizedType


FFTW = -I$(FFTWINCLUDE) -L$(FFTWLIBRARY)
CFLAGS2 = $(CFLAGS)
CILKFLAGS2 = $(CFLAGS) -D USE_CILK -std=c99
MPCFLAGS2 = $(MPCFLAGS) -D USE_MPI
INCLUDE = -I./include
LINK = -L./lib
VPATH = src:include:lib

sizedType.h : makeSizedType.c
	$(CC) $(CFLAGS2) src/makeSizedType.c -o makeSizedType
	./makeSizedType
	rm -f makeSizedType

generalTools_ser.o : generalTools.c generalTools.h sizedType.h
	$(CC) $(CFLAGS2) $(INCLUDE) -c src/generalTools.c 
	mv generalTools.o generalTools_ser.o

generalTools_par.o : generalTools.c generalTools.h sizedType.h
	$(MPCC) $(MPCFLAGS2) $(INCLUDE) -c src/generalTools.c
	mv generalTools.o generalTools_par.o

fftw_complex_helper.o : fftw_complex_helper.c fftw_complex_helper.h generalTools.h sizedType.h
	$(CC) $(CFLAGS2) $(INCLUDE) $(FFTW) -c src/fftw_complex_helper.c 

ccSHT.o : ccSHT.c ccSHT.h generalTools.h fftw_complex_helper.h sizedType.h
	$(CC) $(CFLAGS2) $(INCLUDE) $(FFTW) -c src/ccSHT.c

ccSHTcilk.o : ccSHTcilk.c ccSHT.h generalTools.h fftw_complex_helper.h sizedType.h
	$(CC) $(CILKFLAGS2) $(INCLUDE) $(FFTW) -c src/ccSHTcilk.c

ccSHTmpi.o : ccSHTmpi.c ccSHTmpi.h ccSHT.h generalTools.h fftw_complex_helper.h sizedType.h
	$(MPCC) $(MPCFLAGS2) $(INCLUDE) $(FFTW) -c src/ccSHTmpi.c

ccSHTfortran.o : ccSHTfortran.c ccSHT.h generalTools.h fftw_complex_helper.h sizedType.h ccSHTfortran.h
	$(CC) $(CFLAGS2) $(INCLUDE) $(FFTW) -c src/ccSHTfortran.c

ccSHTmpiFortran.o : ccSHTmpiFortran.c ccSHTmpi.h ccSHT.h generalTools.h fftw_complex_helper.h sizedType.h ccSHTfortran.h
	$(MPCC) $(MPCFLAGS2) $(INCLUDE) $(FFTW) -c src/ccSHTmpiFortran.c

libccSHT.a : generalTools_ser.o fftw_complex_helper.o ccSHT.o ccSHTfortran.o
	ar rc lib/libccSHT.a generalTools_ser.o fftw_complex_helper.o ccSHT.o ccSHTfortran.o
	ranlib lib/libccSHT.a

libccSHTcilk.a : generalTools_ser.o fftw_complex_helper.o ccSHTcilk.o ccSHTfortran.o
	ar rc lib/libccSHTcilk.a generalTools_ser.o fftw_complex_helper.o ccSHTcilk.o ccSHTfortran.o
	ranlib lib/libccSHTcilk.a

libccSHTmpi.a : generalTools_par.o fftw_complex_helper.o ccSHT.o ccSHTmpi.o ccSHTfortran.o ccSHTmpiFortran.o
	ar rc lib/libccSHTmpi.a generalTools_par.o fftw_complex_helper.o ccSHT.o ccSHTmpi.o ccSHTfortran.o ccSHTmpiFortran.o
	ranlib lib/libccSHTmpi.a

ccSHTx : ccSHTx.c libccSHT.a 
	$(CC) $(CFLAGS2) $(FFTW) $(INCLUDE) $(LINK) src/ccSHTx.c -o bin/ccSHTx -lm -lfftw -lccSHT -lfftw -lm

ccSHTmpix : ccSHTx.c libccSHTmpi.a
	$(MPCC) $(MPCFLAGS2) $(FFTW) $(INCLUDE) $(LINK) src/ccSHTx.c -o bin/ccSHTmpix -lm -lfftw -lccSHTmpi -lfftw -lm 

ccSHTcilkx : ccSHTx.c libccSHTcilk.a
	$(CC) $(CFLAGS2) $(FFTW) $(INCLUDE) $(LINK) src/ccSHTx.c -o bin/ccSHTcilkx -lm -lfftw -lccSHTcilk -lfftw -lm -lpthread 

