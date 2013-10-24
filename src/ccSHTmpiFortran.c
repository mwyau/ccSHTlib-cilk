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
These are simply wrappers for the functions found in ccSHTmpi.c.  In 
fortran the subroutines would look as follows: 

       forwardshtmpi_real_fortran( map, coords, lmax, alm, outpu-
       tRoot, theComm )
         REAL*8, DIMENSION(:), INTENT(IN) :: map
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         COMPLEX*8, DIMENSION(:), INTENT(OUT) :: alm
         INTEGER*4, INTENT(IN) :: outputRoot
         INTEGER, INTENT(IN) :: theComm

       forwardshtmpi_complex_fortran(  map,  coords,  lmax,  alm,
       outputRoot, theComm )
         COMPLEX*8, DIMENSION(:), INTENT(IN) :: map
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         COMPLEX*8, DIMENSION(:), INTENT(OUT) :: alm
         INTEGER*4, INTENT(IN) :: outputRoot
         INTEGER, INTENT(IN) :: theComm

       backwardshtmpi_real_fortran( alm, coords, lmax, map,  out-
       putRoot, theComm )
         COMPLEX*8, DIMENSION(:), INTENT(IN) :: alm
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         REAL*8, DIMENSION(:), INTENT(OUT) :: map
         INTEGER*4, INTENT(IN) :: outputRoot
         INTEGER, INTENT(IN) :: theComm

       backwardshtmpi_complex_fortran(  alm,  coords,  lmax, map,
       outputRoot, theComm )
         COMPLEX*8, DIMENSION(:), INTENT(IN) :: alm
         INTEGER*8, INTENT(IN) :: coords
         INTEGER*4, INTENT(IN) :: lmax
         COMPLEX*8, DIMENSION(:), INTENT(OUT) :: map
         INTEGER*4, INTENT(IN) :: outputRoot
         INTEGER, INTENT(IN) :: theComm

*******************************************************************************/

#include "ccSHTmpi.h"
#include "ccSHTfortran.h"



void forwardshtmpi_real_fortran( FORTRAN_REAL_8 *map, 
                                 FORTRAN_INTEGER_8 *coords, 
                                 FORTRAN_INTEGER_4 *lmax, 
                                 FORTRAN_COMPLEX_8 *alm,
                                 FORTRAN_INTEGER_4 *outputRoot,
                                 MPI_Comm *theComm )

{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  forwardSHTmpi( map, 0, *coordsPtr, *lmax, alm, *outputRoot, *theComm );
}


void forwardshtmpi_complex_fortran( FORTRAN_COMPLEX_8 *map, 
                                    FORTRAN_INTEGER_8 *coords, 
                                    FORTRAN_INTEGER_4 *lmax, 
                                    FORTRAN_COMPLEX_8 *alm,
                                    FORTRAN_INTEGER_4 *outputRoot, 
                                    MPI_Comm *theComm )
{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  forwardSHTmpi( map, 1, *coordsPtr, *lmax, alm, *outputRoot, *theComm );
}


void backwardshtmpi_real_fortran( FORTRAN_COMPLEX_8 *alm, 
                                  FORTRAN_INTEGER_8 *coords,
                                  FORTRAN_INTEGER_4 *lmax, 
                                  FORTRAN_REAL_8 *map,
                                  FORTRAN_INTEGER_4 *outputRoot, 
                                  MPI_Comm *theComm )
{
  coordStruct *coordsPtr;
  
  coordsPtr = *((coordStruct**)coords);

  backwardSHTmpi( alm, *coordsPtr, *lmax, map, 0, *outputRoot, *theComm);
}


void backwardshtmpi_complex_fortran( FORTRAN_COMPLEX_8 *alm, 
                                     FORTRAN_INTEGER_8 *coords,
                                     FORTRAN_INTEGER_4 *lmax, 
                                     FORTRAN_COMPLEX_8 *map,
                                     FORTRAN_INTEGER_4 *outputRoot, 
                                     MPI_Comm *theComm )
{
  coordStruct *coordsPtr;

  coordsPtr = *((coordStruct**)coords);

  backwardSHTmpi( alm, *coordsPtr, *lmax, map, 1, *outputRoot, *theComm);
}
