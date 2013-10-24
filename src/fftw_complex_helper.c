/*******************************************************************************
*   Version 1.03 July 2003                                                     *
*                                                                              *
*   Copyright (C) 2003  C.M. Cantalupo                                         *
*                                                                              *
*   fftw_complex_helper is free software; you can redistribute it and/or modify*
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

#include "fftw_complex_helper.h"

fftw_complex comPlusCom(fftw_complex A, fftw_complex B)
{
  A.re = A.re+B.re;
  A.im = A.im+B.im;
  return(A);
}


fftw_complex scalTimesCom(flt8 A, fftw_complex B)
{
  B.re *= A;
  B.im *= A;
  return(B);
}


fftw_complex comTimesCom(fftw_complex A, fftw_complex B)
{
  fftw_complex temp;
  temp.re = A.re*B.re - A.im*B.im;
  temp.im = A.re*B.im + A.im*B.re;
  return(temp);
}


fftw_complex comDivideCom(fftw_complex A, fftw_complex B)
{
  double length2;
  length2 = length(B);
  length2 *= length2;

  return(comTimesCom(A,scalTimesCom(1/length2,ccSHT_conj(B))));
}


fftw_complex ccSHT_conj( fftw_complex A )
{
  A.im = -A.im;
  return A;
}


double length( fftw_complex A )
{
  return(sqrt(A.re*A.re+A.im*A.im));
}


double angle( fftw_complex A )
{
  double theAngle = atan2(A.im, A.re);
  if(theAngle >= 0)
    return(theAngle);
  else
    return(theAngle + 2*PI);
}
