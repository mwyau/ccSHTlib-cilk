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

/*******************************************************************************
These are a set of functions which can be easily linked to a fortran code. 
These are simply wrappers for the functions found in ccSHT.c.  In fortran 
the subroutines would look as follows:  

       forwardsht_real_fortran( map, coords, lmax, alm )
         REAL*8, DIMENSION(:), INTENT(IN) :: map
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         COMPLEX*8, DIMENSION(:), INTENT(OUT) :: alm

       forwardsht_complex_fortran( map, coords, lmax, alm )
         COMPLEX*8, DIMENSION(:), INTENT(IN) :: map
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         COMPLEX*8, DIMENSION(:), INTENT(OUT) :: alm

       backwardsht_real_fortran( alm, coords, lmax, map)
         COMPLEX*8, DIMENSION(:), INTENT(IN) :: alm
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         REAL*8, DIMENSION(:), INTENT(OUT) :: map

       backwardsht_complex_fortran( alm, coords, lmax, map)
         COMPLEX*8, DIMENSION(:), INTENT(IN) :: alm
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         COMPLEX*8, DIMENSION(:), INTENT(OUT) :: map

       convertradec2coords_fortran(  pixSize,  nPix,  raDecArray,
       ndra, dra, coords)
         REAL*8, INTENT(IN) :: pixSize
         INTEGER*4, INTENT(IN) :: nPix
         REAL*8, DIMENSION(:), INTENT(IN) :: raDecArray
         INTEGER*4, INTENT(IN) :: ndra
         REAL*8, DIMENSION(:), INTENT(IN) :: dra
         INTEGER*8, INTENT(OUT) :: coords

       destroycoords_fortran( coords)
         INTEGER*8, INTENT(INOUT) :: coords

       calculateqlm_fortran( lmax, x, calc1, calc2, calc3, calc4, Q )
         INTEGER*4, INTENT(IN) :: lmax 
         REAL*8, INTENT(IN) :: x
         REAL*8, INTENT(IN) :: calc1 
         REAL*8, INTENT(IN) :: calc2 
         REAL*8, INTENT(IN) :: calc3 
         REAL*8, INTENT(IN) :: calc4 
         REAL*8, INTENT(OUT) :: Q

       lmcalculations_fortran( lmax, calc1, calc2, calc3, calc4 )
         INTEGER*4, INTENT(IN) :: lmax
         REAL*8, INTENT(OUT) :: calc1 
         REAL*8, INTENT(OUT) :: calc2 
         REAL*8, INTENT(OUT) :: calc3 
         REAL*8, INTENT(OUT) :: calc4

       calcuatesmallqlm_fortran( lmax, x, calc1, calc2, calc3, calc4, Qstart, lStart )
         INTEGER*4, INTENT(IN) :: lmax 
         REAL*8, INTENT(IN) :: x 
         REAL*8, INTENT(IN) :: calc1
         REAL*8, INTENT(IN) :: calc2 
         REAL*8, INTENT(IN) :: calc3 
         REAL*8, INTENT(IN) :: calc4 
         REAL*8, INTENT(OUT) :: Qstart 
         INTEGER*4, INTENT(OUT) :: lStart


*******************************************************************************/

#include "ccSHT.h"
#include "ccSHTfortran.h"


void forwardsht_real_fortran( FORTRAN_REAL_8 *map, 
                              FORTRAN_INTEGER_8 *coords, 
                              FORTRAN_INTEGER_4 *lmax, 
                              FORTRAN_COMPLEX_8 *alm )
{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  forwardSHT( map, 0, *coordsPtr, *lmax, alm );
}


void forwardsht_complex_fortran( FORTRAN_COMPLEX_8 *map, 
                                 FORTRAN_INTEGER_8 *coords, 
                                 FORTRAN_INTEGER_4 *lmax, 
                                 FORTRAN_COMPLEX_8 *alm )
{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  forwardSHT( map, 1, *coordsPtr, *lmax, alm );
}


void backwardsht_real_fortran( FORTRAN_COMPLEX_8 *alm, 
                               FORTRAN_INTEGER_8 *coords,
                               FORTRAN_INTEGER_4 *lmax, 
                               FORTRAN_REAL_8 *map)
{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  backwardSHT( alm, *coordsPtr, *lmax, map, 0);
}


void backwardsht_complex_fortran( FORTRAN_COMPLEX_8 *alm, 
                                  FORTRAN_INTEGER_8 *coords,
                                  FORTRAN_INTEGER_4 *lmax, 
                                  FORTRAN_COMPLEX_8 *map)
{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  backwardSHT( alm, *coordsPtr, *lmax, map, 1);
}


void convertradec2coords_fortran( FORTRAN_REAL_8 *pixSize, 
                                  FORTRAN_INTEGER_4 *nPix, 
                                  FORTRAN_REAL_8 *raDecArray,
                                  FORTRAN_INTEGER_4 *ndra, 
                                  FORTRAN_REAL_8 *dra,
                                  FORTRAN_INTEGER_8 *output )
{
  coordStruct *coords;
  
  coords = (coordStruct*)malloc(sizeof(coordStruct));
  
  *coords = convertRaDec2Coords( *pixSize, *nPix, raDecArray, *ndra, dra);

  *((coordStruct**)output) = coords;
}


void destroycoords_fortran( FORTRAN_INTEGER_8 *coords)
{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  destroyCoords(coordsPtr);
  
  free(coordsPtr);

  coordsPtr = NULL;
}

void calculateqlm_fortran( FORTRAN_INTEGER_4 *lmax, 
                           FORTRAN_REAL_8 *x, 
                           FORTRAN_REAL_8 *calc1, 
                           FORTRAN_REAL_8 *calc2, 
                           FORTRAN_REAL_8 *calc3, 
                           FORTRAN_REAL_8 *calc4, 
                           FORTRAN_REAL_8 *Q)
{
  calculateQlm( *lmax, *x, calc1, calc2, calc3, calc4, Q);
}


void lmcalculations_fortran( FORTRAN_INTEGER_4 *lmax, 
                             FORTRAN_REAL_8 *calc1, 
                             FORTRAN_REAL_8 *calc2, 
                             FORTRAN_REAL_8 *calc3, 
                             FORTRAN_REAL_8 *calc4)
{
  lmCalculations( *lmax, calc1, calc2, calc3, calc4 );
}

void calcuatesmallqlm_fortran( FORTRAN_INTEGER_4 *lmax, 
                               FORTRAN_REAL_8 *x, 
                               FORTRAN_REAL_8 *calc1, 
                               FORTRAN_REAL_8 *calc2, 
                               FORTRAN_REAL_8 *calc3, 
                               FORTRAN_REAL_8 *calc4, 
                               FORTRAN_REAL_8 *Qstart, 
                               FORTRAN_INTEGER_4 *lStart )
{
  calculateSmallQlm( *lmax, *x, calc1, calc2, calc3, calc4, Qstart, lStart );
}
