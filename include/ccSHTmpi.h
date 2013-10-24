/*******************************************************************************
*  ccSHTmpi:                                                                   *
*    Some functions which do spherical harmonic transforms (SHT's) in parallel *
*    via the MPI library.  The two principal functions are the following:      *
*                                                                              *
*    forwardSHTmpi()                                                           *
*      parallel forward SHT                                                    *
*    backwardSHTmpi()                                                          *
*      parallel backward SHT                                                   *
*                                                                              *
*  C.M. Cantalupo 11/20/01 last updated 7/14/03                                *
*******************************************************************************/


/*******************************************************************************
*   Version 1.03 July 2003                                                     *
*                                                                              *
*   Copyright (C) 2003  C.M. Cantalupo                                         *
*                                                                              *
*   This program is free software; you can redistribute it and/or modify       *
*   it under the terms of the GNU General Public License as published by       *
*   the Free Software Foundation; either version 2 of the License, or          *
*   (at your option) any later version.                                        *
*                                                                              *
*   This program is distributed in the hope that it will be useful,            *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
*   GNU General Public License for more details.                               *
*                                                                              *
*   You should have received a copy of the GNU General Public License          *
*   along with this program; if not, write to the Free Software                *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  *
*                                                                              *
*******************************************************************************/

#ifndef _H_ccSHTmpi
#define _H_ccSHTmpi

#include "mpi.h"
#include "ccSHT.h"

/*********************************************************** 
*  almStruct is a data structure used to store spherical   *
*  harmonic coefficients that have been distributed over   *
*  a group of processors.                                  *
***********************************************************/
typedef struct
{
  int lmax;
  int nm;
  int *numData;
  int *mList;
  fftw_complex **A;
  void *allocPtrs[4];
} almStruct;


/***********************************************************
*  rowFFTstruct is a data structure used to store Fourier  *
*  transforms of rows of pixels which have been            *
*  distributed over a group of processors.                 *
***********************************************************/
typedef struct 
{
  int numRows;
  int *numData;
  int *rowList;
  int *rowLength;
  fftw_complex **data;
} rowFFTstruct;

/***********************************************************
*                   LIBRARY FUNCTIONS                      *
***********************************************************/
void forwardSHTmpi( void *map, int inputIsComplex, coordStruct coords, int lmax, 
                    fftw_complex *outputAlms, int outputRoot, MPI_Comm theComm);
void backwardSHTmpi( fftw_complex *alm, coordStruct coords, int lmax, void *map, 
                     int outputIsComplex, int outputRoot, MPI_Comm theComm);


/***********************************************************
*                   INTERNAL FUNCTIONS                     *
***********************************************************/

/* Functions for breaking up gathering pixel domain objects */
coordStruct breakUpCoords( coordStruct coords, int numProc, int myRank);
void breakUpVector( void *map, coordStruct coords, coordStruct myCoords, 
                    int numProc, int myRank, void *myMap, int isComplex);
void assembleVector( void *myMap, coordStruct coords, coordStruct myCoords, 
                     void *map, int isComplex, int outputRoot, MPI_Comm theComm);
void gatherSmallQlm(int *lStartLocal, int *lStartAll, int **lStart, 
                    flt8 *QstartLocal, flt8 *QstartAll, flt8 **Qstart, 
                    int nThetaValsAll, int lmax, MPI_Comm theComm );

/* Functions for breaking up and gathering spherical harmonic domain objects */
almStruct createAlmStruct( int lmax, int numProc, int myRank);
void destroyAlmStruct(almStruct *A );
void gatherAlm(almStruct *myAlm, fftw_complex *allAlm, int outputRoot, MPI_Comm theComm);

/* Functions for breaking up, gathering and computing linear fourier domain objects */
rowFFTstruct createRowFFTstruct(int numRow, flt8 *deltaPhi, int numProc, int myRank);
void destroyRowFFTstruct(rowFFTstruct *theStruct);
void calculateRowFFTstruct(coordStruct splitCoords, void *splitMap, 
                           int mapIsComplex, rowFFTstruct *splitFFT);
void gatherRowFFTstruct(rowFFTstruct *myFFT, rowFFTstruct *allFFT, MPI_Comm theComm);



#endif
