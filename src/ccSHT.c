/*******************************************************************************
*   Version 1.031 October 2003                                                 *
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

#include "ccSHT.h"

void forwardSHT( void *mapv, int inputIsComplex, coordStruct coords, int lmax, fftw_complex *alm )
{
/*******************************************************************************

       These functions are the members of the ccSHT library which
       perform  the forward discrete spherical harmonic transform
       (SHT).  forwardSHT() performs the  computation  in  serial
       and  forwardSHTmpi() performs the computation in parallel.
       For  some  background,  and  more  information  about  the
       library see the man page for ccSHT.

       forwardSHT() and forwardSHTmpi perform a band limited dis-
       crete SHT, where the band limit is set by the input param-
       eter  called  lmax.   More exactly, the transform computes
       the spherical harmonic coefficients (a = a_{l,m}) for  all
       l  in  the  set {0, 1, 2, ..., lmax} and m in the set {-l,
       -l+1, -l+2, ..., l}.  The spherical harmonic  coefficients
       are defined to be

              /
             |      _
        a =  |  f * Y dA
             |
            /

       where

        a = a_{l,m}

        f = f(theta,phi)
        _   _
        Y = Y_{l,m}(theta, phi)

       and  the  integral is over the surface of the unit sphere.
       In the above formulation, Y bar are the complex conjugates
       of  the spherical harmonics.  More explicitly, the complex
       conjugates of the spherical harmonics are defined

                ________________
        _      / (2*l+1)*(l-m)!
        Y =   / ----------------*P(cos(theta))*exp(-i*m*phi)
            \/     4*pi*(l+m)!

       where P(cos(theta)) = P_{l,m}(cos(theta)) are the  associ-
       ated  Legendre polynomials, and i is the imaginary number.
       The above formulation is only  valid  for  non-negative  m
       since  P  is  only defined for non-negative m.  Y_{l,m} is
       related to Y_{l,-m} as follows:

                       m _
        Y_{l,-m} = (-1) *Y_{l,m}

       so we can infer the values of Y for negative m.

       The input to the transform is a vector of evaluations of f
       at  a set of positions on the sphere.  The location of the
       positions at which the function f is evaluated are  stored
       in  the  structure  called  coords (see man page ccSHT for
       more details about this structure).  Let  us  refer  to  a
       vector  of values of the spherical harmonics for a fixed l
       and m over all of the pixels on the sphere as a  spherical
       harmonic vector.  The integral is then calculated as a dot
       product in pixel space of  the  function  vector  and  the
       spherical  harmonic  vector  with  a  constant  quadrature
       weight (coords.pixSize).  If a more complicated quadrature
       scheme  is preferred, simply set coords.pixSize to one and
       pass f multiplied by the preferred quadrature weighting to
       forwardSHT() in place of f.

       For  more information about the calculation of the associ-
       ated Legendre polynomials see the man page for  calculate-
       Qlm().   The  transform  does not actually call calculate-
       Qlm() for optimization purposes.  Instead most of the cod-
       ing  required to do the calculation of the associated Leg-
       endre polynomials is done in the functions calculateSmall-
       Qlm()  and lmCalculations() which the transform does call.

       Below is a table describing the parameters passed to  for-
       wardSHT() and forwardSHTmpi().

       map    A  vector  of length number of pixels (coords.nPix)
              of function values to be transformed.  map  can  be
              either  real  or  complex  (flt8  or  fftw_complex)
              depending on the value of inputIsComplex.

       inputIsComplex
              A  logical  indicator  of  the  type  of  map.   If
              inputIsComplex  is  not zero then map is assumed to
              be  a   vector   of   length   number   of   pixels
              (coords.nPix)  eight  byte  floating  point complex
              pairs (fftw_complex).  If  inputIsComplex  is  zero
              then map is assumed to be a vector of length number
              of pixels (coords.nPix) eight byte  floating  point
              numbers (flt8).

       coords A  structure  which  contains the information about
              the pixelization of the sphere (see  man  page  for
              ccSHT  for  a  detailed  description of this struc-
              ture).

       lmax   Band limit which defines the maximum  l  value  for
              which  the spherical harmonic coefficients are com-
              puted.

       alm    The spherical harmonic  coefficients  corresponding
              to the input vector.  This is the primary output of
              forwardSHT().  Note that these are complex, and the
              FFTW  structure  for complex numbers has been used.
              This should be  a  vector  with  enough  space  for
              (lmax+1)^2   fftw_complex   values.   The  indexing
              scheme is discussed in the man page for ccSHT.

*******************************************************************************/

  int l, m, i, j, k, maxJ, nPhi;
  fftw_complex CC, tempc1, tempc2;
  flt8 x, minDeltaPhi, temp, Q0, Q1, Q2, nQ;
  int *lStart;
  flt8 *Qstart, *mapd, *norm, *calc1, *calc2, *calc3, *calc4, *calc3ptr, *calc4ptr;
  fftw_complex *mapc, *mapRow, *mapRowTrans, *almPosPtr, *almNegPtr;
  fftw_plan **plans;

  if(inputIsComplex)
  {
    mapc = (fftw_complex*)mapv;
    mapd = NULL;
  }
  else
  {
    mapd = (flt8*)mapv;
    mapc = NULL;
  }

  /* figure out the smallest deltaPhi */

  minDeltaPhi = findSmallestEl(coords.deltaPhi,coords.nThetaVals);

  /* Allocate arrays */

  maxJ = ccSHT_round((2*PI)/fabs(minDeltaPhi));
  mapRow = (fftw_complex *)calloc(maxJ,sizeof(fftw_complex));
  mapRowTrans = (fftw_complex *)calloc(maxJ,sizeof(fftw_complex));
  norm = (flt8 *)malloc(sizeof(flt8)*(lmax + 1));
  calc1 = (flt8 *)malloc(sizeof(flt8)*(lmax+1));
  calc2 = (flt8 *)malloc(sizeof(flt8)*(lmax+1));
  calc3 = (flt8 *)malloc(sizeof(flt8)*((lmax*lmax + 3*lmax)/2+1));
  calc4 = (flt8 *)malloc(sizeof(flt8)*((lmax*lmax + 3*lmax)/2+1));
  Qstart = (flt8 *)malloc(sizeof(flt8)*2*(lmax+1));
  lStart = (int *)malloc(sizeof(int)*(lmax+1));

  errorCheck(-1, "forwardSHT", (mapRow && mapRowTrans && norm && calc1 && calc2 && calc3 && calc4 && Qstart && lStart), 1);
  
  /* Create the plans for fftw */
  plans = createPlans(coords);

  /* calculate sqrt((2*l + 1)/(4*PI)) */
  temp = coords.pixSize*sqrt(1/(4*PI));
  for( l = 0; l <= lmax; l++)
    norm[l] = temp*sqrt(2*l+1);

  /* initalize alm's to zero */
  memset( alm, 0, sizeof(fftw_complex)*(lmax+1)*(lmax+1));

  /* calculate the functions of l and m */
  lmCalculations( lmax, calc1, calc2, calc3, calc4);

  /* Now that all the dirty work is done, do the transform */

  k = 0;
  for( i = 0; i < coords.nThetaVals; i++)
  {
    /* Figure out the number of pixels in the padded row */
    maxJ = ccSHT_round((2*PI)/fabs(coords.deltaPhi[i])); 
 
    x = cos(coords.thetaVals[i]);

    /*  Calculate the first (in terms of l) two non-zero Q_lm(x) values for all m  */
    calculateSmallQlm(lmax, x, calc1, calc2, calc3, calc4, Qstart, lStart);

    /* Figure out number of pixels sampled in this row */
    if(i < coords.nThetaVals -1)
      nPhi = coords.thetaBreaks[i+1] - coords.thetaBreaks[i];
    else
      nPhi = coords.nPix - coords.thetaBreaks[i];

    /* Make a zero padded fftw_complex array corresponding to this row */
    j = 0;
    l = 0;
    while( j < nPhi)
    {
      if( k < coords.nGaps && coords.thetaBreaks[i] + j == coords.gaps[2*k] )
      {
        memset(mapRow + l, 0, sizeof(fftw_complex)*coords.gaps[2*k+1]);
        l += coords.gaps[2*k+1];
        k++;
      }
      if(inputIsComplex)
        mapRow[l] = mapc[coords.thetaBreaks[i]+j];
      else
      {
        mapRow[l].re = mapd[coords.thetaBreaks[i]+j];
        mapRow[l].im = 0;
      }
      j++;
      l++;      
    }
    if( maxJ-l )
      memset(mapRow + l, 0, sizeof(fftw_complex)*(maxJ-l));
    
    /* Do the Fourier transform */

    fftw_one(*(plans[i]),mapRow,mapRowTrans);
   
    /* Do the Spherical Harmonic transform */

    /* Calculate for m == 0 */

    almPosPtr = alm + almIndex(0,0,lmax);
    calc3ptr = calc3 + Qindex(0,0,lmax);
    calc4ptr = calc4 + Qindex(0,0,lmax);

    Q2 = Qstart[0];
    Q1 = Qstart[1];

    l = lStart[0];  /* note that this should be zero */

    /*************************************************************
    *  Solve for real and imaginary component of a_{l,0} if      *
    *  input vector is complex                                   *
    *************************************************************/
    if( inputIsComplex )
    {
      if( l <= lmax )
      {
        nQ = norm[l]*Q2;
        almPosPtr[l].re += nQ*mapRowTrans->re; 
        almPosPtr[l].im += nQ*mapRowTrans->im;      
        l++;
        if( l <= lmax )
        {
          nQ = norm[l]*Q1;
          almPosPtr[l].re += nQ*mapRowTrans->re;
          almPosPtr[l].im += nQ*mapRowTrans->im;
      
          l++;
          for( ; l <= lmax; l++)
          {
            Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
            nQ = norm[l]*Q0;
            almPosPtr[l].re += nQ*mapRowTrans->re;
            almPosPtr[l].im += nQ*mapRowTrans->im;
            Q2 = Q1;
            Q1 = Q0;
          }
        }
      }
    }
    /**********************************************************************
    *  If the input vector is real valued then we know that the zero      *
    *  mode of the fft will be real valued as well, so only calculate     *
    *  for the real component of a_{l,0}.                                 *
    **********************************************************************/
    else 
    {
      if( l <= lmax )
      {
        nQ = norm[l]*Q2;
        almPosPtr[l].re += nQ*mapRowTrans->re; 
        l++;
        if( l <= lmax )
        {
          nQ = norm[l]*Q1;
          almPosPtr[l].re += nQ*mapRowTrans->re;
      
          l++;
          for( ; l <= lmax; l++)
          {
            Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
            nQ = norm[l]*Q0;
            almPosPtr[l].re += nQ*mapRowTrans->re;
            Q2 = Q1;
            Q1 = Q0;
          }
        }
      }
    }

    /* Calculate a_{l,m} for m != 0*/
    for( m = 1; m <= lmax; m++ )
    {
      /* Set the phi0 rotation variable for positive m */
      CC.re = cos(m*coords.phi0[i]);
      CC.im = -sin(m*coords.phi0[i]);

      /*****************************************************************
      *  Pick out the Fourier component of intrest for positive m and  *
      *  rotate it.                                                    *
      *****************************************************************/
      j = myMod((m),maxJ);
      tempc1 = comTimesCom(CC,mapRowTrans[j]);

      if( inputIsComplex )
      {
        /* Set the phi0 rotation variable for negative m */
        CC.im *= -1;
        /*****************************************************************
        *  Pick out the Fourier component of intrest for negetive m and  *
        *  rotate it.                                                    *
        *****************************************************************/
        j = myMod((-m),maxJ);
        tempc2 =  scalTimesCom(pown1(m),comTimesCom(CC,mapRowTrans[j]));
      }

      /* Set a pointer which will index a_{l,m} for this m, and -m */
      almPosPtr = alm + almIndex(m,m,lmax) - m;
      almNegPtr = alm + almIndex(m,-m,lmax) - m; 

      /* Set pointers which will index the Legendre recursion coefficients */
      calc3ptr = calc3 + Qindex(m,m,lmax) - m;
      calc4ptr = calc4 + Qindex(m,m,lmax) - m;
      
      /* Set the lowest l for which the Q's are non-zero */
      l = lStart[m];
      /* Set the starting values of the Legendre recursion */
      Q2 = Qstart[2*m];
      Q1 = Qstart[2*m+1];

      /************************************************************
      *  If the input vector is complex then we must calculate    *
      *  a_{l,m} for both positive and negetive m.                *
      ************************************************************/ 
      if( inputIsComplex )
      {
        if( l <= lmax )
        {
          nQ = norm[l]*Q2;
          almPosPtr[l].re += nQ*tempc1.re;
          almPosPtr[l].im += nQ*tempc1.im;
          almNegPtr[l].re += nQ*tempc2.re;
          almNegPtr[l].im += nQ*tempc2.im;

          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;
            almPosPtr[l].re += nQ*tempc1.re;
            almPosPtr[l].im += nQ*tempc1.im;
            almNegPtr[l].re += nQ*tempc2.re;
            almNegPtr[l].im += nQ*tempc2.im;

            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2; 
              nQ = norm[l]*Q0;
              almPosPtr[l].re += nQ*tempc1.re;
              almPosPtr[l].im += nQ*tempc1.im;
              almNegPtr[l].re += nQ*tempc2.re;
              almNegPtr[l].im += nQ*tempc2.im;
              Q2 = Q1;
              Q1 = Q0;
            }
          }
        }
      }
      /*************************************************************
      *  If the input vector is real valued then we can            *
      *  calculate a_{l,m} for m > 0, and then infer the values    *
      *  for m < 0 later.                                          *
      *************************************************************/
      else
      {
        if( l <= lmax )
        {
          nQ = norm[l]*Q2;
          almPosPtr[l].re += nQ*tempc1.re;
          almPosPtr[l].im += nQ*tempc1.im;

          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;
            almPosPtr[l].re += nQ*tempc1.re;
            almPosPtr[l].im += nQ*tempc1.im;

            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2; 
              nQ = norm[l]*Q0;
              almPosPtr[l].re += nQ*tempc1.re;
              almPosPtr[l].im += nQ*tempc1.im;
              Q2 = Q1;
              Q1 = Q0;
            }
          }
        }
      }  
    }
  }

  /* if input is real, then fill in the negetive m component of the alm's */
  if( !inputIsComplex )
  {
    for( m = 1; m <= lmax; m++ )
    {
      almPosPtr = alm + almIndex(m,m,lmax) - m;
      almNegPtr = alm + almIndex(m,-m,lmax) - m;      
      
      if( m%2 )
      {
        for( l = m; l <= lmax; l++ )
	{
          almNegPtr[l].re = -almPosPtr[l].re;
          almNegPtr[l].im = almPosPtr[l].im;
        }
      }
      else
      {
        for( l = m; l <= lmax; l++ )
	{
          almNegPtr[l].re = almPosPtr[l].re;
          almNegPtr[l].im = -almPosPtr[l].im;
        }
      }
    }
  }
  /* free up memory */

  destroyPlans(coords,  plans);
  free(lStart);
  free(Qstart);
  free(calc4);
  free(calc3);
  free(calc2);
  free(calc1);
  free(norm);
  free(mapRowTrans);
  free(mapRow);
}

void backwardSHT( fftw_complex *alm, coordStruct coords, int lmax, void *mapv, int outputIsComplex)
{
/*******************************************************************************

       These functions are the members of the ccSHT library which
       perform the backward discrete spherical harmonic transform
       (SHT).  backwardSHT() performs the computation  in  serial
       and backwardSHTmpi() performs the computation in parallel.
       For  some  background,  and  more  information  about  the
       library see the man page for ccSHT.

       backwardSHT()  and backwardSHTmpi() perform a band limited
       discrete SHT, where the band limit is  set  by  the  input
       parameter  called lmax.  That is to say that the transform
       assumes a_{l,m} is zero for all l greater than lmax.   For
       this reason the user supplies the spherical harmonic coef-
       ficients only for values of l less than or equal  to  lmax
       (l  in  the set {0, 1, 2, ..., lmax} and m in the set {-l,
       -l+1, -l+2, ..., l}).

       The backward discrete spherical harmonic transform creates
       a  pixel domain vector (g) from a vector of spherical har-
       monic coefficients (a).  This is  done  by  the  following
       computation:

             lmax     l
             -----  -----
              \      \
        g =    )      )    a*Y
              /      /
             -----  -----
             l = 0  m = -l

       where

        g = g(theta,phi)

        a = a_{l,m}

        Y = Y_{l,m}(theta,phi)

       In  the  above formulation, Y are the spherical harmonics.
       More explicitly, the spherical harmonics are defined

                ________________
               / (2*l+1)*(l-m)!
        Y =   / ----------------*P(cos(theta))*exp(i*m*phi)
            \/     4*pi*(l+m)!

       where P(cos(theta)) = P_{l,m}(cos(theta)) are the  associ-
       ated  Legendre polynomials, and i is the imaginary number.
       The above formulation is only  valid  for  non-negative  m
       since  P  is  only defined for non-negative m.  Y_{l,m} is
       related to Y_{l,-m} as follows:
                          _
        Y_{l,-m} = (-1)^m*Y_{l,m}

       so we can infer the values of Y for negative m.

       The input to the transform is a vector of  spherical  har-
       monic  coefficients  (alm).  The indexing of the spherical
       harmonic coefficient vector is described in the  man  page
       for  ccSHT(1).   The  output of the transform is the pixel
       domain vector g (map).  The locations of the positions  on
       the sphere where the function g is evaluated are stored in
       the structure called coords (see man page ccSHT  for  more
       details about this structure).

       For  more information about the calculation of the associ-
       ated Legendre polynomials see the man page for  calculate-
       Qlm().  The transform doesn't actually call calculateQlm()
       for optimization purposes.  Instead  most  of  the  coding
       required  to do the calculation of the associated Legendre
       polynomials is done in the  functions  calculateSmallQlm()
       and lmCalculations() which the transform does call.

       Below is a table describing the parameters passed to back-
       wardSHT() and backwardSHTmpi().

       alm    The coefficients of the band limited spherical har-
              monic  expansion  of  the  desired output function.
              Note that these are complex, and we have  used  the
              FFTW  structure for complex numbers.  alm should be
              a length (lmax+1)^2 fftw_complex vector. The index-
              ing scheme is defined in in the man page for ccSHT.

       coords A structure which contains  the  information  about
              the  pixelization  of  the sphere (see man page for
              ccSHT for a detailed  description  of  this  struc-
              ture).

       lmax   Band  limit  which  defines the maximum l value for
              which the spherical harmonic coefficients  are  not
              assumed to be zero.

       map    A  pointer  to  the  output (g).  The output can be
              either real or complex depending on  the  value  of
              outputIsComplex.   If  outputIsComplex is zero then
              the imaginary part of g is ignored.  In this  case,
              if the imaginary part of g is non-trivial a warning
              message  is   printed   to   stderr.    There   are
              coords.nPix  elements (real or complex) in the out-
              put vector g.

       outputIsComplex
              A logical  indicator  as  to  the  type  of  output
              desired  if  outputIsComplex  is  non-zero then the
              output is complex, otherwise the output is just the
              real  part  of  g.  See description of map for more
              information.


*******************************************************************************/


  int i,j,k,l,m, maxJ, maxM, nBigI;
  flt8 x, minDeltaPhi, sq4pin1, nQ, Q0, Q1, Q2;
  fftw_complex CC, cTemp1, cTemp2;
  int *lStart;
  flt8 *mapd, *norm, *calc1, *calc2, *calc3, *calc4, *calc3ptr, *calc4ptr, *Qstart;
  fftw_complex *mapc, *mapRow, *mapRowTrans, *almPosPtr, *almNegPtr;
  fftw_plan **plans;

  if( outputIsComplex )
  {
    mapc = (fftw_complex*)mapv;
    mapd = NULL;
  }
  else
  {
    mapd = (flt8*)mapv;
    mapc = NULL;
  }

  /* figure out the smallest deltaPhi */
  minDeltaPhi = findSmallestEl(coords.deltaPhi,coords.nThetaVals);

  /* Allocate arrays */
  maxJ = ccSHT_round((2*PI)/fabs(minDeltaPhi));
  mapRow = (fftw_complex *)calloc(maxJ,sizeof(fftw_complex));
  mapRowTrans = (fftw_complex *)calloc(maxJ,sizeof(fftw_complex));
  norm = (flt8*)malloc(sizeof(flt8)*(lmax+1));
  calc1 = (flt8 *)malloc(sizeof(flt8)*(lmax+1));
  calc2 = (flt8 *)malloc(sizeof(flt8)*(lmax+1));
  calc3 = (flt8 *)malloc(sizeof(flt8)*((lmax*lmax + 3*lmax)/2+1));
  calc4 = (flt8 *)malloc(sizeof(flt8)*((lmax*lmax + 3*lmax)/2+1));
  Qstart = (flt8 *)malloc(sizeof(flt8)*2*(lmax+1));
  lStart = (int *)malloc(sizeof(int)*(lmax+1));

  
  errorCheck(-1, "backwardSHT", ( mapRow && mapRowTrans && norm && calc1 && calc2 && calc3 && calc4 && Qstart && lStart), 1);

  /* Create the plans for fftw */
  plans = createPlans(coords);

  /* calculate sqrt((2*l + 1)/(4*PI)) */
  sq4pin1 = 1/sqrt(4*PI); 
  for( l = 0; l <= lmax; l++)
    norm[l] = sq4pin1*sqrt(2*l+1);

  /* calculate the functions of l and m */
  lmCalculations( lmax, calc1, calc2, calc3, calc4);

  k = 0;
  nBigI = 0;
  for( i = 0; i < coords.nThetaVals; i++)
  {
    maxJ = ccSHT_round((2*PI)/fabs(coords.deltaPhi[i]));

    x = cos(coords.thetaVals[i]);

    /*  Calculate the first (in terms of l) two non-zero Q_lm(x) values for all m  */
    calculateSmallQlm(lmax, x, calc1, calc2, calc3, calc4, Qstart, lStart);

    /* make vectors for the Fourier transforms */

    memset(mapRow, 0, sizeof(fftw_complex)*maxJ);

    /* calculate contribution to mapRow[0] by a_{l,0} */
    almPosPtr = alm + almIndex(0,0,lmax);
    calc3ptr = calc3 + Qindex(0,0,lmax);
    calc4ptr = calc4 + Qindex(0,0,lmax);

    Q2 = Qstart[0];
    Q1 = Qstart[1];

    if( outputIsComplex )
    {
      l = lStart[0]; /* note that this should be zero */
      if( l <= lmax )
      {
        nQ = norm[l]*Q2;
        mapRow->re += nQ*almPosPtr[l].re;
        mapRow->im += nQ*almPosPtr[l].im;
      
        l++;
        if( l <= lmax )
        {
          nQ = norm[l]*Q1;
          mapRow->re += nQ*almPosPtr[l].re;
          mapRow->im += nQ*almPosPtr[l].im;
          l++;
          for( ; l <= lmax; l++ )
          {
            Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
            nQ = norm[l]*Q0;
            mapRow->re += nQ*almPosPtr[l].re;
            mapRow->im += nQ*almPosPtr[l].im;
            Q2 = Q1;
            Q1 = Q0;
          }
        }
      }
    }
    else
    {
      l = lStart[0]; /* note that this should be zero */
      if( l <= lmax )
      {
        nQ = norm[l]*Q2;
        mapRow->re += nQ*almPosPtr[l].re;
      
        l++;
        if( l <= lmax )
        {
          nQ = norm[l]*Q1;
          mapRow->re += nQ*almPosPtr[l].re;
          l++;
          for( ; l <= lmax; l++ )
          {
            Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
            nQ = norm[l]*Q0;
            mapRow->re += nQ*almPosPtr[l].re;
            Q2 = Q1;
            Q1 = Q0;
          }
        }
      }
    }


    /* calculate the rest of mapRow */
    for( m = 1; m <= lmax; m++)
    {
      cTemp1.re = cTemp1.im = cTemp2.re = cTemp2.im = 0;

      almPosPtr = alm + almIndex(m,m,lmax) - m;
      almNegPtr = alm + almIndex(m,-m,lmax) - m;
      calc3ptr = calc3 + Qindex(m,m,lmax) - m;
      calc4ptr = calc4 + Qindex(m,m,lmax) - m;

      Q2 = Qstart[2*m];
      Q1 = Qstart[2*m+1];

      if( outputIsComplex )
      {
        l = lStart[m];
        if( l <= lmax )
        {
          nQ = norm[l]*Q2;
          cTemp1.re += nQ*almPosPtr[l].re;
          cTemp1.im += nQ*almPosPtr[l].im;
          cTemp2.re += nQ*almNegPtr[l].re;
          cTemp2.im += nQ*almNegPtr[l].im;
        
          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;
            cTemp1.re += nQ*almPosPtr[l].re;
            cTemp1.im += nQ*almPosPtr[l].im;
            cTemp2.re += nQ*almNegPtr[l].re;
            cTemp2.im += nQ*almNegPtr[l].im;

            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
              nQ = norm[l]*Q0;
              cTemp1.re += nQ*almPosPtr[l].re;
              cTemp1.im += nQ*almPosPtr[l].im;
              cTemp2.re += nQ*almNegPtr[l].re;
              cTemp2.im += nQ*almNegPtr[l].im;
              Q2 = Q1;
              Q1 = Q0;
            }
          }
        }
      }
      else
      {
        l = lStart[m];
        if( l <= lmax )
        {
          nQ = norm[l]*Q2;
          cTemp1.re += nQ*almPosPtr[l].re;
          cTemp1.im += nQ*almPosPtr[l].im;
        
          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;
            cTemp1.re += nQ*almPosPtr[l].re;
            cTemp1.im += nQ*almPosPtr[l].im;

            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
              nQ = norm[l]*Q0;
              cTemp1.re += nQ*almPosPtr[l].re;
              cTemp1.im += nQ*almPosPtr[l].im;
              Q2 = Q1;
              Q1 = Q0;
            }
          }
        }
      }

      CC.re = cos(m*coords.phi0[i]);
      CC.im = sin(m*coords.phi0[i]);

      j = myMod(-m,maxJ);
      mapRow[j] = comPlusCom(mapRow[j], comTimesCom(CC,cTemp1));

      if( !outputIsComplex )
      {
        if( m%2 )
        {
          cTemp2.re = -cTemp1.re;
          cTemp2.im = cTemp1.im;
        }
        else
        {
          cTemp2.re = cTemp1.re;
          cTemp2.im = -cTemp1.im;
        }
      }

      CC.im *= -1;    
      j = myMod(m,maxJ);
      mapRow[j] = comPlusCom(mapRow[j],scalTimesCom(pown1(m),comTimesCom(CC,cTemp2)));
    }
    /* Do Fourier transform */

    fftw_one(*(plans[i]),mapRow,mapRowTrans);

    if( i < coords.nThetaVals - 1 )
      maxM = coords.thetaBreaks[i+1] - coords.thetaBreaks[i];
    else
      maxM = coords.nPix - coords.thetaBreaks[i];

    /* remove padded parts */
    m = 0; /* m is the index for the unpadded pixel "column" */
    j = 0; /* j is the index for the padded pixel "column" */
    while( m < maxM )
    {
      if( k < coords.nGaps && coords.thetaBreaks[i] + m == coords.gaps[2*k] )
      {
        j += coords.gaps[2*k+1];
        k++;
      }

      if(outputIsComplex)
        mapc[coords.thetaBreaks[i] + m] = mapRowTrans[j];
      else
      {
        if(fabs(mapRowTrans[j].im) > fabs(mapRowTrans[j].re)*MY_EPSILON )
          nBigI++;
        mapd[coords.thetaBreaks[i] + m] = mapRowTrans[j].re;
      }

      j++;
      m++;
    }
  }

  if(nBigI)
  {
    fprintf(stderr, "backwardSHT:\n");
    fprintf(stderr, "%.2f percent of output pixels have non-trivial imaginary component.\n",(100.0*nBigI)/coords.nPix);
    fprintf(stderr, "WARNING: Will only return real part\n");
  }


  destroyPlans(coords, plans);
  free(lStart);
  free(Qstart);
  free(calc4);
  free(calc3);
  free(calc2);
  free(calc1);
  free(norm);
  free(mapRowTrans);
  free(mapRow); 

}

void calculateQlm(int lmax, flt8 x, flt8 *calc1, flt8 *calc2, flt8 *calc3, flt8 *calc4, flt8 *Q)
{
/*******************************************************************************

       This function is the member of  the  ccSHT  library  which
       calculates  the  scaled  associated  Legendre polynomials.
       For  some  background,  and  more  information  about  the
       library see the man page for ccSHT.

       The  associated  Legendre  polynomials  are functions of a
       real variable, x, and two integer indexes, l and m.   cal-
       culateQlm()  solves  for  the  scaled  associated Legendre
       polynomials for a fixed x, all l in the set  {0,  1,  ...,
       lmax},  and  all m in the set {0, 1, ..., l}.  The scaling
       of the polynomials is given here:

                      ________
                     / (l-m)!
        Q_{l,m} =   / -------- P_{l,m}
                  \/   (l+m)!


       where P_{l,m} are the associated Legendre polynomials  and
       Q_{l,m}  are  the  scaled associated Legendre polynomials.
       The calculation of the polynomials is done with  the  fol-
       lowing three term recursion on l:

                                  2  m/2
        Q_{m,m}(x) = c1(m)*( 1 - x  )

        Q_{m+1,m}(x) = c2(m)*x*Q_{m,m}(x)

        Q_{l,m}(x) =  c3(l,m)*x*Q_{l-1,m}(x)
                          - c4(l,m)*Q_{l-2,m}(x)

       where the functions of l and m are given below:

                           ____________________
                          /    m               |
                     m   /  --------'
        c1(m)  = (-1) * /  '  |  |          1
                       /      |  |    1 - -----
                      /       |  |         2*i
                    \/       i = 1

                    _____
                   /     |
        c2(m) =   / 2*m+1
                \/

                    ( 2*l-1 )
        c3(l,m) = --------------
                      2   2  1/2
                   ( l - m  )

                       _____________
                      /      2   2  |
                     /  (l-1) - m
        c4(l,m) =   / -------------
                   /       2   2
                 \/       l - m

       It  can  be  seen in the above equations that Q_{m,m} gets
       exponentially small with large m.  For x close to -1 or  1
       and  large  m,  Q_{m,m}  is  not representable as a double
       because  the  absolute  value  of  Q_{m,m}  is  less  than
       2^(-1074).  This algorithm rescales Q when it would under-
       flow and keeps track of  the  scaling.   This  is  because
       although Q_{m,m}(x) may be virtually zero, Q_{m+500,m} may
       be substantial but unless the first  two  numbers  in  the
       recursion   over   l   are  non-zero  the  calculation  of
       Q_{l,m}(x) would remain 0 for all l.

       Before calling calculateQlm, some preliminary calculations
       must  be  done  with the function called lmCalculations().
       See the man page of lmCalculations() for more information.
       The preliminary calculations are passed to calculateQlm by
       way of the inputs calc1, calc2, calc3, and calc4.

       Below is a table describing the parameters passed to  cal-
       culateQlm()

       lmax   Maximum l value for the computation

       x      The  fixed  x  value (between negative one and one)
              which determines the position for which the polyno-
              mials are computed.

       calc1, calc2, calc3, calc4
              Calculations  which  are functions only of l and m.
              These are to be calculated by the  function  called
              lmCalculations().   See  man  page  for  lmCalcula-
              tions() for more details.

       Q      This is the vector where the computed scaled  asso-
              ciated  Legendre polynomials will be stored.  There
              should be space for (lmax*lmax + 3*lmax)/2+1  eight
              byte  floating  point  numbers  in  the array.  The
              indexing scheme for this vector is explained in the
              ccSHT man page.

*******************************************************************************/

  int l,m;
  int *lStart;
  flt8 *Qstart, *calc3ptr, *calc4ptr, *Qptr;

  if(lmax < -1)
    return;

  /* Allocate some memory */
  lStart = (int*)malloc(sizeof(int)*(lmax+1));
  Qstart = (flt8*)malloc(sizeof(flt8)*2*(lmax+1));
  errorCheck(-1, "calculateQlm", (lStart && Qstart), 1);

  /* Calculate first represenably non-zero Q values for each m */
  calculateSmallQlm( lmax, x, calc1, calc2, calc3, calc4, Qstart, lStart);

  for( m = 0; m <= lmax; m++ )
  {
    /* Set some pointers for easy array indexing */
    Qptr = Q + Qindex(m,m,lmax) - m;
    calc3ptr = calc3 + Qindex(m,m,lmax) - m;
    calc4ptr = calc4 + Qindex(m,m,lmax) - m;

    if( lStart[m] == lmax + 1)
    {
      /* Q's are all representably zero for this m */ 
      memset( Q + Qindex(m,m,lmax), 0, sizeof(flt8)*(lmax - m + 1));
    }
    else    
    {
      /* Set representably zero Q's to zero */ 
      memset( Q + Qindex(m,m,lmax), 0, sizeof(flt8)*(lStart[m] - m ));
        
      /* Set inital conditions for the recurrence */
      Qptr[lStart[m]] = Qstart[2*m];
      if( lStart[m] < lmax )
      {
        Qptr[lStart[m]+1] = Qstart[2*m+1];

        /* Calculate the recurrence */
        for( l = lStart[m]+2; l <= lmax; l++)
          Qptr[l] = x*calc3ptr[l]*Qptr[l-1] + calc4ptr[l]*Qptr[l-2];
      }
    }
  }
  free(Qstart);
  free(lStart);
}

void lmCalculations(int lmax, flt8 *calc1, flt8 *calc2, flt8 *calc3, flt8 *calc4)
{
/******************************************************************************* 

       This function is a member of the ccSHT library which  does
       preliminary calculations used by calculateQlm() to compute
       the associated Legendre polynomials.   These  calculations
       are  only  functions  of  l and m.  The user must allocate
       memory for  each  of  these  calculations  before  calling
       lmCalculations().    The  amount  of  memory  required  is
       explained below in the description of  the  input  parame-
       ters.

       lmax   This  is  the maximum l value for which the associ-
              ated Legendre polynomials will be  calculated.   It
              is  important  that  this  is  the same as the lmax
              which is passed to calculateQlm().

       calc1  A vector of length  (lmax+1)  eight  byte  floating
              point numbers.

       calc2  A  vector  of  length  (lmax+1) eight byte floating
              point numbers.

       calc3  A vector of length ((lmax*lmax + 3*lmax)/2+1) eight
              byte floating point numbers.

       calc4  A vector of length ((lmax*lmax + 3*lmax)/2+1) eight
              byte floating point numbers.

*******************************************************************************/
  int l,m;
  flt8 dlm0, dlm1, temp;
  flt8 *calc3ptr, *calc4ptr;

  calc1[0] = 1.0;
  calc2[0] = 1.0;
  for(l = 1; l <= lmax; l++)
  {
    calc3[Qindex(l,0,lmax)] = 2 - 1.0/l;
    calc4[Qindex(l,0,lmax)] = 1.0/l - 1;
  }

  for( m = 1; m <= lmax; m++)
  {
    calc1[m] = -sqrt(1 - 1.0/(2*m));
    calc2[m] = sqrt(2*m+1);
      
    calc3ptr = calc3 + Qindex(m,m,lmax) - m;
    calc4ptr = calc4 + Qindex(m,m,lmax) - m;

    dlm0 = 0;
    for( l = m+1; l <= lmax; l++)
    {
      dlm1 = dlm0;
      dlm0 = sqrt(l*l - m*m);
      temp = 1/dlm0;
      calc3ptr[l] = temp*(2*l-1);
      calc4ptr[l] = -temp*dlm1;
    }
  }
}


void calculateSmallQlm( int lmax, flt8 x, flt8 *calc1, flt8 *calc2, flt8 *calc3, flt8 *calc4, flt8 *Qstart, int *lStart)
{
  int checkExp, thisgExp, thisQexp, Q1exp, Q2exp, xExp, rr, dr;
  int i, l, m, lmaxFile;
  flt8 sint, thisgFrac, thisQfrac, Qmm, Q0, Q1, Q2, smallestDouble;
  flt8 *calc3ptr, *calc4ptr;
  char preCalcFileName[STRLEN];
  FILE *preCalcFile;

  static int printWarningMessage = 0;
  
/*******************************************************************************
*  Set the name of the file where we will look for solution.  These            *
*  files are assumed to be in the directory FILE_DIR (defined in ccSHT.h).     *
*  Note that if the user has created a directory called FILE_DIR               *
*  then the the solution files will be created if they don't already exist.    *
*******************************************************************************/
    
  sprintf(preCalcFileName, "%s/QlmStart_%.16e", FILE_DIR, x);
  preCalcFile = fopen(preCalcFileName, "r");

  /* Check to see if the file exists. */
  if(preCalcFile)
  {
    /* If it does and the lmax of the file is large enough, 
       read it in and return.  */
    fread(&lmaxFile, sizeof(int), 1, preCalcFile);
    if( lmaxFile >= lmax )
    {
      fread( lStart, sizeof(int), lmax+1, preCalcFile);
      if( lmaxFile > lmax )
        fseek( preCalcFile, sizeof(int)*(lmaxFile-lmax), SEEK_CUR);
      fread( Qstart, sizeof(flt8), 2*(lmax+1), preCalcFile);

      fclose(preCalcFile);
      return;
    }
    fclose(preCalcFile);
  }


  /*************************************************************
  *  Calculate Qlm = sqrt( (l-m)!/(l+m)! ) * Plm(x)            *
  *  Underflows at the beginning of the recurrence are avoided *
  *  by rescaling the exponent of the floating point number    *
  *************************************************************/

  /* Set sin(theta), where x := cos(theta) */
  sint = sqrt(1-x*x);
  /* Find out the exponent in the floating point number x */
  frexp(x, &xExp);

  /* Calculate Qstart and lStart for m == 0 */
  Qstart[0] = Q2 = Qmm = 1;
  Qstart[1] = Q1 = calc2[0]*x;
  lStart[0] = 0;

  /* Calculate for the m where no scaling needs to be done */

  smallestDouble = ldexp(64, MIN_EXP );
  for( m = 1; (m <= lmax) && (fabs(Qstart[2*m-2]) > smallestDouble); m++ )
  {
    Qstart[2*m] = sint*calc1[m]*Qstart[2*m-2];
    Qstart[2*m+1] = x*calc2[m]*Qstart[2*m];
    lStart[m] = m;
  }

  /* Now calculate for the m where scaling must be done */

  Qmm = Qstart[2*m-2];
  rr = 0;

  for( ; m <= lmax; m++)
  {
    /* thisQ is Q_{m-1,m-1}, thisg is calc1[m]*sint  */
    thisQfrac = frexp(Qmm,&thisQexp);
    thisgFrac = frexp(calc1[m]*sint,&thisgExp);          

    checkExp = thisgExp + thisQexp + xExp - 2;

    /* check to see if we need to rescale */
    if( checkExp >= MIN_EXP )
    {
      Qmm = calc1[m]*sint*Qmm;
    }
    else
    {
      /* increase scale factor ( 2^rr ) */
      dr = MIN_EXP - checkExp;
      rr += dr;
      Qmm = ldexp(thisgFrac*thisQfrac,thisgExp + thisQexp + dr);
    }

    Q2 = Qmm;
    Q1 = calc2[m]*x*Qmm;

    /* If no rescaling has been done */
    if( rr == 0 )
    {
      Qstart[2*m] = Q2;
      Qstart[2*m+1] = Q1;
      lStart[m] = m;
    }
    else 
    {
      /* Calculate Q for increasing l until they become representable 
         as non-zero doubles. */

      calc3ptr = calc3 + Qindex(m,m,lmax) - m;
      calc4ptr = calc4 + Qindex(m,m,lmax) - m;

      frexp(Q1, &Q1exp);
      frexp(Q2, &Q2exp);

      /* If we need to rescale more than once for this m use 
         a loop that checks for overflow. */
      if( rr > MAX_EXP - MIN_EXP - 32)
      {
        for( l = m+2; (l <= lmax) && (Q2exp < MIN_EXP + rr); l++ )
        {
          if( Q1exp >= MAX_EXP )
          {
            /***************************************************
            *  rescale again, but our scale factor must still  *
            *  be greater than 1.                              *
            ***************************************************/
            dr = MIN_EXP - Q2exp;
            Q2 = ldexp(Q2,dr);
            Q1 = ldexp(Q1,dr);
            rr += dr;
          }
          /* calculate Q_{l,m}  */
          Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
          Q2 = Q1;
          Q1 = Q0;

          Q2exp = Q1exp;
          frexp(Q1, &Q1exp);
        }
      }
      else
      {
        /* Use a loop that doesn't check for overflow */
        for( l = m+2; (l <= lmax) && (Q2exp < MIN_EXP + rr); l++ )
        {
          /* calculate Q_{l,m}  */
          Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
          Q2 = Q1;
          Q1 = Q0;

          frexp(Q2, &Q2exp);
        }        
      }
      /* If we found non-zero Q's record them */
      if( l <= lmax )
      {
        Qstart[2*m] = ldexp(Q2,-rr);
        Qstart[2*m+1] = ldexp(Q1,-rr);
        lStart[m] = l-2;
      }
      /* Otherwise we are done */
      else
      {
        memset(Qstart+2*m, 0, sizeof(flt8)*2*(lmax+1-m));
        for( ; m <= lmax; m++)
          lStart[m] = lmax+1;
      }
    }
  }

  preCalcFile = fopen(preCalcFileName, "w");

  if(preCalcFile)
  {
    fwrite( &lmax, sizeof(int), 1, preCalcFile);
    fwrite( lStart, sizeof(int), lmax+1, preCalcFile);
    fwrite( Qstart, sizeof(flt8), 2*(lmax+1), preCalcFile);

    fclose(preCalcFile);
  }
  else if( printWarningMessage )
  {
    fprintf( stderr, "\n\n");
    fprintf( stderr, "WARNING: calculateSmallQlm() did not write speed up information to disk.\n");
    fprintf( stderr, "         You can allow calculateSmallQlm() to write these files to disk\n");
    fprintf( stderr, "         by creating a directory called \"ccSHT_files\".  You can free up all\n");
    fprintf( stderr, "         the disk space used by this function by deleting \"ccSHT_files/QlmStart_*\"\n");
    fprintf( stderr, "         The speed up files take up (16+sizeof(int))*(lmax+1) bytes for\n");
    fprintf( stderr, "         each x value passed to calculateSmallQlm() (i.e. for all theta).\n");
    fprintf( stderr, "         Therefore the disk use for an all sky map without a truncated\n");
    fprintf( stderr, "         spectrum is proportional to the number of pixels in the map.\n\n\n");
    printWarningMessage = 0;
  } 
}

fftw_plan **createPlans(coordStruct coords)
{
/*******************************************************************************
* createPlans:                                                                 *
*                                                                              *
*   This function creates an array of pointers to fftw_plan's to be used in    *
*   an SHT.  There is a pointer for every row of pixels in coords, but only    *
*   one plan for every particular length of row.  That is to say that the      *
*   pointers corresponding to rows which are the same length will point to     *
*   the same plan.  To free up the plans call destroyPlans().                  *
*                                                                              *
*******************************************************************************/


  int i, j, nPhi;
  int *n, *s;
  fftw_plan **plans;

  n = (int*)malloc(coords.nThetaVals*sizeof(int));
  s = (int*)malloc(coords.nThetaVals*sizeof(int));
  plans = (fftw_plan**)malloc(coords.nThetaVals*sizeof(fftw_plan*));  

  errorCheck(-1, "createPlans", (n && plans), 1);

  for( i = 0; i < coords.nThetaVals; i++ )
  {
    n[i] = ccSHT_round((2*PI)/fabs(coords.deltaPhi[i]));
    s[i] = coords.deltaPhi[i]/fabs(coords.deltaPhi[i]);
    
    for( j = 0; j < i; j++)
      if(n[j] == n[i] && s[i] == s[j])
      {
        plans[i] = plans[j];
        break;
      }
    if(j == i)
    {
      plans[i] = (fftw_plan*)malloc(sizeof(fftw_plan));
      if( coords.deltaPhi[i] > 0 )
        *(plans[i]) = fftw_create_plan(n[i],FFTW_FORWARD,FFTW_ESTIMATE);
      else
        *(plans[i]) = fftw_create_plan(n[i],FFTW_BACKWARD,FFTW_ESTIMATE);
    }
  }
  free(s);
  free(n);
  return(plans);
}


void destroyPlans(coordStruct coords, fftw_plan **plans)
{
/*******************************************************************************
* destroyPlans:                                                                *
*                                                                              *
*   This function frees up memory allocated by createPlans.                    *
*                                                                              *
*******************************************************************************/


  int i, j, nPhi;
  int *n;
  
  n = (int*)malloc(coords.nThetaVals*sizeof(int));
  for( i = 0; i < coords.nThetaVals; i++ )
  {
    n[i] = ccSHT_round((2*PI)/fabs(coords.deltaPhi[i]));
    
    for( j = 0; j < i & n[j] != n[i]; j++);
    if(j == i)
    {
      fftw_destroy_plan(*(plans[i]));
      free(plans[i]);
    }  
  }
  free(plans);
  free(n);
}


coordStruct convertRaDec2Coords(flt8 pixSize, int nPix, flt8 *raDecArray, int ndra, flt8 *dra)
{
/*******************************************************************************

       This function is a member of the ccSHT library which  cre-
       ates  the  structure  that  the transform functions use to
       describe the pixelization of the sphere.  In the ccSHT man
       page there is a subsection on the pixelization which gives
       some general information about the discretization  of  the
       sphere.   The  documentation for convertRaDec2Coors() will
       be clearer if the reader first looks at  the  pixelization
       subsection in the man page for ccSHT.

       This  function  makes  use  of the astronomer's coordinate
       system of right  ascension  (ra)  and  declination  (dec).
       Right  ascension  is  a  measure  of longitudinal angle in
       degrees from 0 to 360 increasing counter clockwise (to the
       east)  around the north pole.  Declination is a measure of
       latitudinal angle in degrees from 90 at the north pole  to
       -90  at the opposite (south) pole.  This is converted into
       standard mathematical spherical coordinates theta and  phi
       as follows:

        theta = (90 - dec)*pi/180
        phi = ra*pi/180

       convertRaDec2Coords() takes an array of pixel positions in
       right ascension and declination coordinates  and  converts
       the  positions into a coordStruct.  This structure is used
       when calling both the forward and backward transform.  The
       positions  of  the pixels within the array raDecArray must
       conform to the assumptions described in  the  pixelization
       subsection  of the man page for ccSHT.  The positions must
       also be ordered in  a  specific  way.   This  ordering  is
       important,  as  all  pixel  domain vectors used with ccSHT
       must also have this ordering.  It is discussed briefly  in
       the man page for ccSHT, but we will discuss it again here.

       convertRaDec2Coords() assumes a  sequential  ordering  for
       the  pixels in the ra dec list.  The ordering must conform
       to the following specification: two pixels which have  the
       same declination and are contiguous on the sphere (or only
       separated by a gap in the pixelization) must be contiguous
       in  the  pixel ordering unless they are the first and last
       pixels listed for a particular value of declination.  This
       specification  does  not  completely determine an ordering
       for a pixelization, none the less, any ordering  which  is
       used must conform to the above specification.

       The  ordering  of  the ra dec array is the ordering of all
       pixel domain vectors used with the coordStruct produced by
       calling convertRaDec2Coords().

       This  function  allocates  memory to create the structure.
       To set all of the fields in the structure to zero and free
       up  the  memory  associated  with a structure created with
       converRaDec2Coords() pass a pointer to  the  structure  to
       the function called destroyCoords().  For more information
       about this function see the man page for  destroyCoords().

       Below  is a table describing the parameters passed to con-
       vertRaDec2Coords()

       pixSize
              Quadrature weight for use in the forward  transform
              (i.e.  area of pixels in radians^2).  If quadrature
              weight is not a constant, set pixSize  to  one  (in
              which  case the integral in the forward SHT is just
              a sum) and point wise multiply your input vector by
              the desired quadrature weighting scheme before for-
              ward transforming.

       nPix   Number of elements in a pixel vector not  including
              any inferred pixels or gaps.

       raDecArray
              An  array  of right ascension and declination pairs
              for the pixelization.   The  vector  should  be  of
              length  nPix*2  and every even indexed entry is the
              ra value corresponding to the next entry  which  is
              the  dec  value  of  a position on the sphere.  The
              ordering of this array is important,  see  descrip-
              tion above.

       ndra   Because  the  constraints on the pixelization allow
              for missing pixels there is  no  way  to  determine
              directly  from a list of ra and dec values the lon-
              gitudinal pixel separation in each row  of  pixels.
              For  this reason this function also takes a list of
              pixel widths.  ndra is the number of widths in this
              list.   If  no  pixel widths are specified (i.e. if
              ndra == 0 or dra == NULL) then the pixel widths are
              inferred  from the raDecArray assuming that no gaps
              exist between the first and second  pixel  in  each
              row of the list.

       dra    An  array of declination and delta ra pairs for the
              pixelization.  If this pointer  is  NULL  then  the
              list  will be inferred from the values in the raDe-
              cArray assuming that  no  gaps  exist  between  the
              first and second pixel in each row of the list.  If
              the vector is given it should be of  length  ndra*2
              and  every even indexed entry is the dec value cor-
              responding to the next entry which is the delta  ra
              value  for  the  row.   Note that there may be more
              rows in this list than there are in the actual pix-
              elization.   The  list  however  must have the same
              ordering of dec values as in raDecArray.  Note that
              the  sign  of dra should correspond to the ordering
              of ra  values  in  raDecArray  (i.e.  if  they  are
              increasing for a given row then dra should be posi-
              tive).

*******************************************************************************/

  flt8 lastDec, thisDec, deg2rad, thisDeltaPhi, temp;
  int i, j, k, jmax, countedGaps, parsedGaps;
  coordStruct coords;

  deg2rad = PI/180.0;
  coords.pixSize = pixSize;
  coords.nPix = nPix;
  coords.nThetaVals = 0;
  coords.nGaps = 0;
  
  /* Figure out how many rows are in data set. */
  lastDec = -100;
  for( i = 1; i < nPix*2; i += 2 )
  {
    thisDec = raDecArray[i];
    if(fabs(thisDec - lastDec) > MY_EPSILON)
      coords.nThetaVals++;
    lastDec = thisDec;
  }
  

  /* Allocate memory for coords. */
  coords.thetaVals = (flt8*)malloc(coords.nThetaVals*sizeof(flt8));
  coords.thetaBreaks = (int*)malloc(coords.nThetaVals*sizeof(int));
  coords.deltaPhi = (flt8*)malloc(coords.nThetaVals*sizeof(flt8));
  coords.phi0 = (flt8*)malloc(coords.nThetaVals*sizeof(flt8));

  errorCheck( -1, "convertRaDec2Coords", 
               (coords.thetaVals && coords.thetaBreaks && coords.deltaPhi && coords.phi0), 1);

  /* Figure out where rows begin, the dec of each row. */
  lastDec = -100;
  j = 0;
  for( i = 0; i < nPix*2; i += 2)
  {
    thisDec = raDecArray[i+1];
  
    if(fabs(thisDec - lastDec) > MY_EPSILON)
    {
      coords.thetaVals[j] = (90 - thisDec)*deg2rad;
      coords.thetaBreaks[j] = i/2;
      coords.phi0[j] = raDecArray[i]*deg2rad;
      j++;
    }

    lastDec = thisDec;
  }

  /* Figure out deltaPhi for each row.  */

  if( ndra && dra )
  {
    j = 0;
    for( i = 0; i < coords.nThetaVals; i++)
    {
      while( j < ndra && fabs((90 - dra[2*j])*deg2rad - coords.thetaVals[i]) > MY_EPSILON )
        j++;

      errorCheck(-1, "dra vector is inconsistent with pixel vector, exiting\n\n", j < ndra, 0);     
      coords.deltaPhi[i] = dra[2*j+1]*deg2rad;
    }
  }
  else
  {
    for( i = 0; i < coords.nThetaVals; i++ )
      coords.deltaPhi[i] = (raDecArray[2*(coords.thetaBreaks[i]+1)]-raDecArray[2*coords.thetaBreaks[i]])*deg2rad;
  }

  /* Figure out the gaps in the pixelization. */

  /****************************************************************
  *  The loop below is executed once or twice.  The first time    *
  *  through the number of gaps is counted.  If there are         *
  *  any gaps in the pixelization then then the loop is executed  *
  *  again.  On the second pass through memory is allocated       *
  *  and the information about the gaps is stored in this memory. *
  ****************************************************************/

  countedGaps = 0;
  parsedGaps = 0;
  while( !countedGaps || (countedGaps && coords.nGaps != 0 && !parsedGaps ) )
  {
    if( countedGaps )
    {    
      /* If we have already counted the gaps then allocate memory.  */
      coords.gaps = (int *)malloc(sizeof(int)*coords.nGaps*2);
      errorCheck(-1, "convertRaDec2Coords", (int)coords.gaps, 1);
    }
  
    k = 0;
    for( i = 0; i < coords.nThetaVals; i++)
    {
      if( i != coords.nThetaVals - 1 )
        jmax = coords.thetaBreaks[i+1];
      else
        jmax = nPix;
      for( j = coords.thetaBreaks[i]; j < jmax - 1; j++)
      {
        /* Check the change in phi for every pixel. */
        thisDeltaPhi = (raDecArray[2*(j+1)] - raDecArray[2*j])*deg2rad;

        /* Account for going accross phi == 0. */
        if( fabs(thisDeltaPhi - 2*PI) < fabs(thisDeltaPhi) )
          thisDeltaPhi -= 2*PI;
        if( fabs(thisDeltaPhi + 2*PI) < fabs(thisDeltaPhi) )
          thisDeltaPhi += 2*PI;

        /* If this change does not equal the expected change we have found a gap. */
        if( fabs(thisDeltaPhi - coords.deltaPhi[i]) > MY_EPSILON )
	{
          /* If this is the first time through count the gap. */
          if( !countedGaps )
            coords.nGaps++;
          else
          {
            /* If this is the second time through record the gap location */ 
            temp = thisDeltaPhi/coords.deltaPhi[i];
           
            errorCheck(-1, "Pixels in a row are not evenly spaced",  fabs(ccSHT_round(temp) - temp) < 0.001, 0 );

            coords.gaps[2*k] = j+1;
            coords.gaps[2*k+1] = ccSHT_round(temp) - 1;
            k++;
          }
        }
      }
    }

    if(countedGaps)
      parsedGaps = 1;
    countedGaps = 1;
  }

  if( coords.nGaps == 0 )
    coords.gaps = NULL;
  
  return(coords);
}


void destroyCoords(coordStruct *coords)
{
/*******************************************************************************

      This is the destructor function for the coordStruct struc-
      ture.  It frees up the vectors, and sets all of the fields
      to  zero of the structure pointed to by the input argument
      coords.    

*******************************************************************************/
 
  coords->pixSize = 0;
  coords->nPix = 0;
  coords->nThetaVals = 0;
  coords->nGaps = 0;
  free(coords->thetaVals);
  free(coords->thetaBreaks);
  free(coords->phi0);
  free(coords->deltaPhi);
  free(coords->gaps);
  coords->thetaVals = NULL;
  coords->thetaBreaks = NULL;
  coords->phi0 = NULL;
  coords->deltaPhi = NULL;
  coords->gaps = NULL;
}

