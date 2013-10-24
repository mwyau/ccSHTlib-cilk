/*******************************************************************************
*  fftw_complex_helper                                                         *
*    A set of functions to do length and angle computations as well as         *
*    complex arithmetic.  The complex data type is defined by fftw.h           *
*    version 2.  The definition of the complex number in fftw changed with     *
*    version 3 (it is not a struct, but a length two vector).                  *
*    fftw_complex_helper is not compatible with fftw version 3.                *
*                                                                              *
*    comPlusCom()                                                              *
*      complex plus a complex.                                                 *
*    scalTimesCom()                                                            *
*      scalar times a complex.                                                 *
*    comTimesCom()                                                             *
*      complex times complex.                                                  *
*    comDivideCom()                                                            *
*      complex divided by complex.                                             *
*    ccSHT_conj()                                                              *
*      complex conjugate.                                                      *
*    length()                                                                  *
*      complex modulus.                                                        *
*    angle()                                                                   *
*      complex arg, angle in radians.                                          *
*                                                                              *
*  C.M. Cantalupo  10/29/01 last updated 7/11/03                               *
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

#ifndef _H_fftw_complex_helper
#define _H_fftw_complex_helper

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw.h>
#include "generalTools.h"


fftw_complex comPlusCom(fftw_complex A, fftw_complex B);
fftw_complex scalTimesCom(flt8 A, fftw_complex B);
fftw_complex comTimesCom(fftw_complex A, fftw_complex B);
fftw_complex comDivideCom(fftw_complex A, fftw_complex B);
fftw_complex ccSHT_conj( fftw_complex A );
double length( fftw_complex A );
double angle( fftw_complex A );

#endif
