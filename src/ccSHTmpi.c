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

#include "ccSHTmpi.h"

void forwardSHTmpi( void *map, int inputIsComplex, coordStruct coords, int lmax, 
                    fftw_complex *outputAlms, int outputRoot, MPI_Comm theComm)
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

       outputRoot
              This is only required  for  forwardSHTmpi().   This
              input  indicates the rank of the processor on which
              the output will be stored.  If  outputRoot  ==  -1,
              then the output will be stored on all processors.

       theComm
              This is only required for forwardSHTmpi().  This is
              the MPI communicator on which the computation is to
              be done.

*******************************************************************************/


  int numProc, myRank, i, j, k, l, m;
  int *lStartBuffer, *lStartPtr;
  int **lStart;
  flt8 temp, x, Q0, Q1, Q2, nQ;
  fftw_complex tempc1, tempc2, CC;
  flt8 *norm, *calc1, *calc2, *calc3, *calc4, *calc3ptr, *calc4ptr;
  flt8 *QstartBuffer, *QstartPtr;
  fftw_complex *summedAlm;
  void *myRowMap, *myLocalMap;
  flt8 **Qstart;
  almStruct myAlm;
  rowFFTstruct myLocalFFT, myRowFFT;
  coordStruct myLocalCoords, myRowCoords;
  int wrapAround[2];
  int procMatrixSize[2];
  int procMatrixPos[2];
  int freeCoords[2];
  int myCartGridRank;
  MPI_Comm cartGridComm, myProcRowComm, myProcColComm;
  MPI_Status theStatus;

  MPI_Comm_size(theComm, &numProc);
  MPI_Comm_rank(theComm, &myRank);

  /* Create processor grid and communicators */
  procMatrixSize[0] = (int)(sqrt(numProc));
  procMatrixSize[1] = (int)(numProc/procMatrixSize[0]);
  while( procMatrixSize[0]*procMatrixSize[1] != numProc )
  {
    procMatrixSize[0]--;
    procMatrixSize[1] = (int)(numProc/procMatrixSize[0]);
  }

  wrapAround[0] = wrapAround[1] = 0;

  MPI_Cart_create( theComm, 2, procMatrixSize, wrapAround, 1, &cartGridComm );
  MPI_Comm_rank(cartGridComm, &myCartGridRank);
  MPI_Cart_coords(cartGridComm, myCartGridRank, 2, procMatrixPos );

  freeCoords[0] = 0;
  freeCoords[1] = 1;
  MPI_Cart_sub(cartGridComm, freeCoords, &myProcRowComm);

  freeCoords[0] = 1;
  freeCoords[1] = 0;
  MPI_Cart_sub(cartGridComm, freeCoords, &myProcColComm);


  /* Create structure for local alm's distributing over processor columns */
  myAlm = createAlmStruct(lmax, procMatrixSize[1], procMatrixPos[1]);

  /* Break up coords so that the theta values are cyclicly distributed over processor rows */
  myRowCoords = breakUpCoords( coords, procMatrixSize[0], procMatrixPos[0]);

  /* Break up myRowCoords so that this rows theta values are cyclicly distributed over processor collumns */
  myLocalCoords = breakUpCoords( myRowCoords, procMatrixSize[1], procMatrixPos[1]);

  /* Allocate some memory */
  if( inputIsComplex )
  {
    myRowMap = malloc(sizeof(fftw_complex)*myRowCoords.nPix);
    myLocalMap = malloc(sizeof(fftw_complex)*myLocalCoords.nPix);
  }
  else
  {
    myRowMap = malloc(sizeof(flt8)*myRowCoords.nPix);
    myLocalMap = malloc(sizeof(flt8)*myLocalCoords.nPix);
  }
  errorCheck(myRank, "forwardSHTmpi", (!(myRowCoords.nPix && !myRowMap) && !(myLocalCoords.nPix && !myLocalMap)), 1);

  /* Break up input map so that the theta values are cyclicly distributed over the processor rows */
  breakUpVector( map, coords, myRowCoords, procMatrixSize[0], procMatrixPos[0], myRowMap, inputIsComplex);

  /* Break up this processor row's map for each column.  */
  breakUpVector( myRowMap, myRowCoords, myLocalCoords, procMatrixSize[1], procMatrixPos[1], myLocalMap, inputIsComplex);
  
  /* Create fft structures */
  myLocalFFT = createRowFFTstruct( myRowCoords.nThetaVals, myRowCoords.deltaPhi, procMatrixSize[1], procMatrixPos[1]);
  
  myRowFFT = createRowFFTstruct( myRowCoords.nThetaVals, myRowCoords.deltaPhi, 1, 0);

  /* Do the FFT's for this processor */
  calculateRowFFTstruct( myLocalCoords, myLocalMap, inputIsComplex, &myLocalFFT);

  free(myLocalMap);

  /* Gather the distributed FFT's */
  gatherRowFFTstruct( &myLocalFFT, &myRowFFT, myProcRowComm);

  /* Free up memory for local FFT's */
  destroyRowFFTstruct(&myLocalFFT);

  /* Calculate small Qlm's in parallel and distibute */

  /* Allocate memory for (l,m) calculations */
  norm = (flt8 *)malloc(sizeof(flt8)*(lmax+1));
  calc1 = (flt8 *)malloc(sizeof(flt8)*(lmax+1));
  calc2 = (flt8 *)malloc(sizeof(flt8)*(lmax+1));
  calc3 = (flt8 *)malloc(sizeof(flt8)*((lmax*lmax + 3*lmax)/2 + 1));
  calc4 = (flt8 *)malloc(sizeof(flt8)*((lmax*lmax + 3*lmax)/2 + 1));


  /* Allocate memory for small Qlm calculations */
  lStartBuffer = (int *)malloc(sizeof(int)*(lmax+1)*(myLocalCoords.nThetaVals+myRowCoords.nThetaVals));
  lStart = (int **)malloc(sizeof(int*)*myRowCoords.nThetaVals);
  QstartBuffer = (flt8 *)malloc(sizeof(flt8)*2*(lmax+1)*(myLocalCoords.nThetaVals + myRowCoords.nThetaVals));
  Qstart = (flt8 **)malloc(sizeof(flt8*)*2*myRowCoords.nThetaVals);

  errorCheck(myRank, "forwardSHTmpi", (norm && calc1 && calc2 && calc3 && calc4 && 
                                         lStartBuffer && lStart && QstartBuffer && Qstart), 1);

  /* Calculate normalization term (pixSize*sqrt((2*l + 1)/(4*PI))) */
  temp = coords.pixSize*sqrt(1/(4*PI));
  for( l = 0; l <= lmax; l++)
    norm[l] = temp*sqrt(2*l+1);

  /* Do other (l,m) calculations */
  lmCalculations( lmax, calc1, calc2, calc3, calc4);


  /* Calculate small Qlm for this processor's theta values */
  for( i = 0; i < myLocalCoords.nThetaVals; i++)
  {
    calculateSmallQlm( lmax, cos(myLocalCoords.thetaVals[i]), calc1, calc2, calc3, calc4, 
                       QstartBuffer + 2*(lmax+1)*i, lStartBuffer + (lmax+1)*i);
  }

  /* Gather small Qlms */
  gatherSmallQlm(lStartBuffer, lStartBuffer + (lmax+1)*myLocalCoords.nThetaVals, lStart, 
                 QstartBuffer, QstartBuffer + 2*(lmax+1)*myLocalCoords.nThetaVals, Qstart, 
                 myRowCoords.nThetaVals, lmax, myProcRowComm );

  /* Set the alm's to zero */
  memset(myAlm.allocPtrs[3], 0, sizeof(fftw_complex)*myAlm.numData[procMatrixPos[1]]);


  /* Do forward SHT for this processor column's m's */
  for( i = 0; i < myRowCoords.nThetaVals; i++)
  {
    x = cos(myRowCoords.thetaVals[i]);
    if( procMatrixPos[1] == 0 )
    {
      /* Calculate for m == 0 */
      calc3ptr = calc3 + Qindex(0,0,lmax);
      calc4ptr = calc4 + Qindex(0,0,lmax);

      Q2 = Qstart[i][0];
      Q1 = Qstart[i][1];

      l = lStart[i][0];  /* note that this should be zero */

      /*************************************************************
      *  Solve for real and imaginary component of a_{l,0} if      *
      *  input vector is complex                                   *
      *************************************************************/
      if( inputIsComplex )
      {
        if( l <= lmax )
        {
          nQ = norm[l]*Q2;
          myAlm.A[0][l].re = nQ*myRowFFT.data[i][0].re;
          myAlm.A[0][l].im = nQ*myRowFFT.data[i][0].im;

          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;
            myAlm.A[0][l].re += nQ*myRowFFT.data[i][0].re;
            myAlm.A[0][l].im += nQ*myRowFFT.data[i][0].im;
        
            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
              nQ = norm[l]*Q0;
              myAlm.A[0][l].re += nQ*myRowFFT.data[i][0].re;
              myAlm.A[0][l].im += nQ*myRowFFT.data[i][0].im;
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
          myAlm.A[0][l].re += nQ*myRowFFT.data[i][0].re;

          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;
            myAlm.A[0][l].re += nQ*myRowFFT.data[i][0].re;
        
            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2;
              nQ = norm[l]*Q0;
              myAlm.A[0][l].re += nQ*myRowFFT.data[i][0].re;
              Q2 = Q1;
              Q1 = Q0;
            }
          }
        }
      }
    }

    /* Calculate a_{l,m} for m != 0*/
    for( k = 0; k < myAlm.nm; k++ )
    {
      /* Skip over any negative m values (they will be computed with the positive ones). */
      while( myAlm.mList[k] <= 0 )
        k++;
      if( k >= myAlm.nm)
        break;

      /* Set m to the first positive m value in the list */
      m = myAlm.mList[k];

      /* Set the phi0 rotation variable for positive m */ 
      CC.re = cos(m*myRowCoords.phi0[i]);
      CC.im = -sin(m*myRowCoords.phi0[i]);

      /*****************************************************************
      *  Pick out the Fourier component of intrest for positive m and  *
      *  rotate it.                                                    *
      *****************************************************************/
      j = myMod(m, myRowFFT.rowLength[i]);
      tempc1 = comTimesCom(CC, myRowFFT.data[i][j]);

      if( inputIsComplex )
      {
        /* Set the phi0 rotation variable for negative m */
        CC.im *= -1;
        /*****************************************************************
        *  Pick out the Fourier component of intrest for negetive m and  *
        *  rotate it.                                                    *
        *****************************************************************/
        j = myMod(-m,myRowFFT.rowLength[i]);
        tempc2 =  scalTimesCom(pown1(m),comTimesCom(CC,myRowFFT.data[i][j]));
      }

      /* Set pointers which will index the Legendre recursion coefficients */
      calc3ptr = calc3 + Qindex(m,m,lmax) - m;
      calc4ptr = calc4 + Qindex(m,m,lmax) - m;
      
      /* Set the lowest l for which the Q's are non-zero */
      l = lStart[i][m];
      /* Set the starting values of the Legendre recursion */
      Q2 = Qstart[i][2*m];
      Q1 = Qstart[i][2*m+1];

      /************************************************************
      *  If the input vector is complex then we must calculate    *
      *  a_{l,m} for both positive and negetive m.                *
      ************************************************************/ 
      if( inputIsComplex )
      {
        if( l <= lmax )
        {
          nQ = norm[l]*Q2;

          myAlm.A[m][l].re += nQ*tempc1.re;
          myAlm.A[m][l].im += nQ*tempc1.im;
          myAlm.A[-m][l].re += nQ*tempc2.re;
          myAlm.A[-m][l].im += nQ*tempc2.im;

          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;

            myAlm.A[m][l].re += nQ*tempc1.re;
            myAlm.A[m][l].im += nQ*tempc1.im;
            myAlm.A[-m][l].re += nQ*tempc2.re;
            myAlm.A[-m][l].im += nQ*tempc2.im;

            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2; 
              nQ = norm[l]*Q0;

              myAlm.A[m][l].re += nQ*tempc1.re;
              myAlm.A[m][l].im += nQ*tempc1.im;
              myAlm.A[-m][l].re += nQ*tempc2.re;
              myAlm.A[-m][l].im += nQ*tempc2.im;
              
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

          myAlm.A[m][l].re += nQ*tempc1.re;
          myAlm.A[m][l].im += nQ*tempc1.im;

          l++;
          if( l <= lmax )
          {
            nQ = norm[l]*Q1;

            myAlm.A[m][l].re += nQ*tempc1.re;
            myAlm.A[m][l].im += nQ*tempc1.im;

            l++;
            for( ; l <= lmax; l++)
            {
              Q0 = x*calc3ptr[l]*Q1 + calc4ptr[l]*Q2; 
              nQ = norm[l]*Q0;

              myAlm.A[m][l].re += nQ*tempc1.re;
              myAlm.A[m][l].im += nQ*tempc1.im;

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
    for( k = 0; k < myAlm.nm; k++ )
    {
      /* Skip over any negative m values in the list. */
      while( myAlm.mList[k] <= 0 )
        k++;
      if( k >= myAlm.nm)
        break;

      /* Set m to the first positive m value in the list */
      m = myAlm.mList[k];

      /* Set the m < 0 alm values */
      if( m%2 )
      {
        for( l = m; l <= lmax; l++ )
	{
          myAlm.A[-m][l].re = -myAlm.A[m][l].re;
          myAlm.A[-m][l].im = myAlm.A[m][l].im;
        }
      }
      else
      {
        for( l = m; l <= lmax; l++ )
	{
          myAlm.A[-m][l].re = myAlm.A[m][l].re;
          myAlm.A[-m][l].im = -myAlm.A[m][l].im;
        }
      }
    }
  }

  /* do communication to sum and gather alm's */
  
  if( outputRoot == -1 )
  {
    /* Allocate some memory for sum which must be done out of place */
    summedAlm = (fftw_complex *)malloc(sizeof(fftw_complex)*myAlm.numData[procMatrixPos[1]]);
    errorCheck(myRank, "forwardSHTmpi(): summedAlm ", (int)summedAlm, 1);

    MPI_Allreduce( myAlm.allocPtrs[3], summedAlm, myAlm.numData[procMatrixPos[1]]*2, MPI_DOUBLE, MPI_SUM, myProcColComm);

    memcpy(myAlm.allocPtrs[3], summedAlm, sizeof(fftw_complex)*myAlm.numData[procMatrixPos[1]]);
    free(summedAlm);    

    gatherAlm(&myAlm, outputAlms, -1, myProcRowComm);

  }
  else
  {
    if(procMatrixPos[0] == 0 )
    {
      summedAlm = (fftw_complex *)malloc(sizeof(fftw_complex)*myAlm.numData[procMatrixPos[1]]);
      errorCheck(myRank, "forwardSHTmpi(): summedAlm ", (int)summedAlm, 1);
    }
    else
      summedAlm = NULL;
 
    /* Sum up the alm's from processors in the same column and store the result on processors in row zero */
    MPI_Reduce( myAlm.allocPtrs[3], summedAlm, myAlm.numData[procMatrixPos[1]]*2, MPI_DOUBLE, MPI_SUM, 0, myProcColComm);

    if( procMatrixPos[0] == 0 )
    {
      memcpy(myAlm.allocPtrs[3], summedAlm, sizeof(fftw_complex)*myAlm.numData[procMatrixPos[1]]);
      free(summedAlm);
    }

    /* In case the upper left hand corner of our processor matrix is not the output root 
       we must allocate some memory */
    if( procMatrixPos[0] == 0 && procMatrixPos[1] == 0 && myRank != outputRoot )
    {
      outputAlms = (fftw_complex*)malloc(sizeof(fftw_complex)*(lmax+1)*(lmax+1));
      errorCheck(myRank, "forwardSHTmpi(): outputAlms", (int)outputAlms, 1);
    }

    /* Gather up the alm's from all of the processors in row zero and store it on procesor (0,0)*/
    if( procMatrixPos[0] == 0 )
      gatherAlm(&myAlm, outputAlms, 0, myProcRowComm);

    /* In case the upper left hand corner of out processor matrix is not the root of theComm
       we must send the root the solution */
    if( procMatrixPos[0] == 0 && procMatrixPos[1] == 0 && myRank != outputRoot )
    {
      MPI_Send(outputAlms, 2*(lmax+1)*(lmax+1), MPI_DOUBLE, outputRoot, 1111, theComm);
      free(outputAlms);
    }
    else if( ( procMatrixPos[0] != 0 || procMatrixPos[1] != 0 ) && myRank == outputRoot )
      MPI_Recv(outputAlms, 2*(lmax+1)*(lmax+1), MPI_DOUBLE, MPI_ANY_SOURCE, 1111, theComm, &theStatus);
    
  }

  /* Clean up */
  free(Qstart);
  free(QstartBuffer);
  free(lStart);
  free(lStartBuffer);
  free(calc4);
  free(calc3);
  free(calc2);
  free(calc1);
  free(norm); 
  free(myLocalMap);
  free(myRowMap);


  destroyRowFFTstruct(&myRowFFT);
  destroyRowFFTstruct(&myLocalFFT);
  destroyCoords(&myRowCoords);
  destroyCoords(&myLocalCoords);
  destroyAlmStruct(&myAlm);

  MPI_Comm_free(&myProcColComm);
  MPI_Comm_free(&myProcRowComm);
  MPI_Comm_free(&cartGridComm);
}


void backwardSHTmpi( fftw_complex *alm, coordStruct coords, int lmax, void *map, 
                     int outputIsComplex, int outputRoot, MPI_Comm theComm) 
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

       outputRoot
              This is only required for  backwardSHTmpi().   This
              input  indicates the rank of the processor on which
              the output will be stored.  If  outputRoot  ==  -1,
              then the output will be stored on all processors.

       theComm
              This  parameter  is  only  required  for  backward-
              SHTmpi().  This is the MPI  communicator  on  which
              the computation is to be done.

*******************************************************************************/

  int numProc, myRank;
  void *myOut;
  coordStruct myCoords;
  fftw_complex *myAlm;

  MPI_Comm_size(theComm, &numProc);
  MPI_Comm_rank(theComm, &myRank);

  if(numProc > coords.nThetaVals)
  {
    if(!myRank)
    {
      fprintf(stderr, "ERROR:  There are more processors than there are rows of pixels,\n");
      fprintf(stderr, "returning without doing computation.\n\n");
    }
    return;
  }

  myCoords = breakUpCoords( coords, numProc, myRank);

  if( outputIsComplex)
    myOut = calloc(myCoords.nPix, sizeof(fftw_complex));
  else
    myOut = calloc(myCoords.nPix, sizeof(flt8));

  errorCheck(myRank, "backwardSHTmpi", (int)myOut, 1);
    
  backwardSHT( alm, myCoords, lmax, myOut, outputIsComplex);
  
  assembleVector( myOut, coords, myCoords, map, outputIsComplex, outputRoot, theComm);

  free(myOut);
  destroyCoords(&myCoords);
}


coordStruct breakUpCoords( coordStruct coords, int numProc, int myRank)
{
/*******************************************************************************
* breakUpCoords:                                                               *
*                                                                              *
*   This function breaks down the pixelization so that each processor has a    *
*   smaller SHT to do.  Essentially the rows of constant latitude are          *
*   distributed cyclically to the processors.  If there were four processors   *
*   and six rows then processor 0 would get rows 0 and 4, processor 1 would    *
*   get rows 1 and 5, processor 2 would get row 2, and processor 3 would get   *
*   row 3.  This function takes as input:                                      *
*                                                                              *
*   coords: The coordStruct to be broken up.                                   *
*   numProc:  The number of processors to which the SHT will be distributed.   *
*   myRank:  The rank of the processor who's coordStruct is to be              *
*   constructed.                                                               *
*                                                                              *
*   This function returns a new coordStruct which should be destroyed with     *
*   the function destroyCoords().                                              *
*                                                                              *
*******************************************************************************/

  int i, j, k, l, nPixThisRow, firstI, lastI;
  coordStruct myCoords; 

 
  myCoords.pixSize = coords.pixSize;

  /* figure out how many rows this processor will get */
  myCoords.nThetaVals = coords.nThetaVals/numProc;
  if( myRank < coords.nThetaVals%numProc )
    myCoords.nThetaVals++;

  /* figure out how many gaps will be in this processor's rows */
  myCoords.nGaps = 0;
  myCoords.gaps = NULL;

  if( coords.nGaps )
  {
    k = 0;
    for( j = myRank; j < coords.nThetaVals; j += numProc )
    {
      firstI = coords.thetaBreaks[j];
      if( j < coords.nThetaVals - 1 )
        lastI = coords.thetaBreaks[j+1];
      else
        lastI = coords.nPix;
      while( k < coords.nGaps && coords.gaps[2*k] < firstI )
        k++;
      if( k == coords.nGaps )
        break;
      while( k < coords.nGaps && coords.gaps[2*k] < lastI )
      {
        myCoords.nGaps++;
        k++;
      }
      if( k == coords.nGaps)
        break;
    }
    if( myCoords.nGaps )
    {
      myCoords.gaps = (int*)malloc(sizeof(int)*myCoords.nGaps*2);
      errorCheck( myRank, "breakUpCoords",(int)myCoords.gaps, 1);
    }
  }

  /* allocate memory for myCoords. */
  myCoords.thetaVals = (flt8*)malloc(sizeof(flt8)*myCoords.nThetaVals);
  myCoords.thetaBreaks = (int*)malloc(sizeof(int)*myCoords.nThetaVals);
  myCoords.phi0 = (flt8*)malloc(sizeof(flt8)*myCoords.nThetaVals);  
  myCoords.deltaPhi = (flt8*)malloc(sizeof(flt8)*myCoords.nThetaVals);

  errorCheck( myRank, "breakUpCoords",(myCoords.thetaVals && myCoords.thetaBreaks && myCoords.phi0 && myCoords.deltaPhi), 1);

  /**********************************************************************
  *   Break up coords.  Note that there are a lot of indexes.           *
  *     i is the index for the pixel row for this processor's coords    *
  *     j is the index for the pixel row for the original coords        *
  *     k is the index for the gaps for the original coords             *
  *     l is the index for the gaps for this processor's coords         *
  **********************************************************************/

  myCoords.thetaBreaks[0] = 0;
  myCoords.nPix = 0;
  i = 0;
  k = 0;
  l = 0;

  for( j = myRank; j < coords.nThetaVals; j += numProc )
  {
    myCoords.thetaVals[i] = coords.thetaVals[j];
    myCoords.phi0[i] = coords.phi0[j];
    myCoords.deltaPhi[i] = coords.deltaPhi[j]; 

    firstI = coords.thetaBreaks[j];

    if( j < coords.nThetaVals -1 )
       lastI = coords.thetaBreaks[j+1];
     else
      lastI = coords.nPix;

    nPixThisRow = lastI - firstI;
    
    myCoords.nPix += nPixThisRow;

    if( i < myCoords.nThetaVals - 1) 
      myCoords.thetaBreaks[i+1] = myCoords.thetaBreaks[i] + nPixThisRow;

    while( k < coords.nGaps && coords.gaps[2*k] < firstI )
      k++;
    while( k < coords.nGaps && coords.gaps[2*k] < lastI )
    {
      myCoords.gaps[2*l] = myCoords.thetaBreaks[i] + (coords.gaps[2*k] - firstI);
      myCoords.gaps[2*l+1] = coords.gaps[2*k+1];
      k++;
      l++;
    }
    i++;
  }

  return(myCoords);
}

void breakUpVector( void *mapv, coordStruct coords, coordStruct myCoords, int numProc, int myRank, void *mymapv, int isComplex)
{
/*******************************************************************************
* breakUpVector:                                                               *
*                                                                              *
*   This is analogous to breakUpCoords, except that this function breaks up    *
*   the data vector associated with the pixelization which was broken down by  *
*   breakUpCoords().  This function works for both complex and real data       *
*   vectors and is used by forwardSHT().  This function takes as input:        *
*                                                                              *
*   mapv:  Pointer to the data vector to be broken up.  This data vector can   *
*        be either real or complex.  The type of data is determined by the     *
*        value of isComplex.                                                   *
*   coords:  The original coordStruct corresponding to mapv.                   *
*   myCoords:  The coordStruct which is a subset of the rows in coords.  It    *
*              is assumed that this structure was created by calling           *
*              breakUpCoords(coords,numProc,myRank).                           *
*   numProc:  The number of processors to which the data is to be distributed. *
*   myRank:  The rank of the processor for which myCoords will be created.     *
*   mymapv:  This is the output which is the subset of mapv which is to be     *
*          processed by the processor specified by myRank.                     *
*   isComplex:  A logical indicator of the type of data vector (real or        *
*               complex).                                                      *
*                                                                              *
*   The inverse of this function is assembleVector().                          *
*                                                                              *
*******************************************************************************/


  int i,j, nPixThisRow;
  flt8 *mapd, *mymapd;
  fftw_complex *mapc, *mymapc;

  if(isComplex)
  {
    mapc = mapv;
    mymapc = mymapv;
    mapd = NULL;
    mymapd = NULL;
  }
  else
  {
    mapd = mapv;
    mymapd = mymapv;
    mapc = NULL;
    mymapc = NULL;
  }

  i = 0;
  for( j = myRank; j < coords.nThetaVals; j += numProc )
  {
    if( j < coords.nThetaVals -1 )
      nPixThisRow = coords.thetaBreaks[j+1] - coords.thetaBreaks[j];
    else
      nPixThisRow = coords.nPix - coords.thetaBreaks[j];

    if( isComplex )
      memcpy(mymapc+myCoords.thetaBreaks[i], mapc + coords.thetaBreaks[j], sizeof(fftw_complex)*nPixThisRow);
    else
      memcpy(mymapd+myCoords.thetaBreaks[i], mapd + coords.thetaBreaks[j], sizeof(flt8)*nPixThisRow);
    i++;
  }
}

void assembleVector( void *mymap, coordStruct coords, coordStruct myCoords, void *mapv, int isComplex, int outputRoot, MPI_Comm theComm)
{
/*******************************************************************************
* assembleVector:                                                              *
*                                                                              *
*   This is associated with breakUpCoords, and is the inverse function to      *
*   breakUpVector().  This function collects pieces of a data vector           *
*   distributed over all of the processors and constructs the complete vector  *
*   on each processor.  The data is assumed to be distributed as described in  *
*   breakUpCoords(), and can be either complex or real depending on the value  *
*   of isComplex.  This function is called by backwardSHT.  The function       *
*   takes as inputs:                                                           *
*                                                                              *
*   mymap:  The input vector that is stored on the processor running the code. *
*   coords:  The coordStruct for the complete pixelization.                    *
*   myCoords:  The coordStruct for the data in mymap.                          *
*   theComm:  The communicator on which the data is distributed.               *
*   mapv:  The output which is a complete data vector corresponding to coords. *
*   isComplex:  A logical indicator of the type (real or complex) of data      *
*               stored in mymap and mapv.                                      *
*                                                                              *
*******************************************************************************/

  int numProc, myRank, i, j, chunkSize, istart, temp, place;
  int *nRecv, *recvOff;
  flt8 *mapd, *outputS;
  fftw_complex *mapc;

  if( isComplex )
  {
    mapc = mapv;
    mapd = NULL;
  }
  else
  {
    mapd = mapv;
    mapc = NULL;
  }

  /* find out about processors */
  MPI_Comm_size(theComm, &numProc);
  MPI_Comm_rank(theComm, &myRank);

  /* allocate memory */
  nRecv = (int*)calloc(numProc, sizeof(int));
  recvOff = (int*)calloc(numProc, sizeof(int));
  errorCheck( myRank, "assembleVector", (nRecv && recvOff), 1);

  if( outputRoot == -1 || outputRoot == myRank )
  {
    if( isComplex )
      outputS = (flt8*)malloc(sizeof(fftw_complex)*coords.nPix);
    else
      outputS = (flt8*)malloc(sizeof(flt8)*coords.nPix);
    errorCheck( myRank, "assembleVector", (int)outputS, 1);
  }
  else
    outputS = NULL;

  /* generate arrays for MPI_Allgatherv */
  for( i = 0; i < coords.nThetaVals - 1; i++ )
  {
    j = i % numProc;
    nRecv[j] += coords.thetaBreaks[i + 1] - coords.thetaBreaks[i];
  }
  j = i % numProc;
  nRecv[j] += coords.nPix - coords.thetaBreaks[i];

  if(isComplex)
    for( i = 0; i < numProc; i++)
      nRecv[i] *= 2;

  recvOff[0] = 0;
  for( i = 1; i < numProc; i++ )
    recvOff[i] = recvOff[i-1] + nRecv[i-1];

  if( outputRoot == -1 )
    MPI_Allgatherv( mymap, (1+(isComplex != 0))*myCoords.nPix, flt8mpiType, outputS, nRecv, recvOff, flt8mpiType, theComm );
  else
    MPI_Gatherv( mymap, (1+(isComplex != 0))*myCoords.nPix, flt8mpiType, outputS, nRecv, recvOff, flt8mpiType, outputRoot, theComm );


  if( outputRoot == -1 || outputRoot == myRank )
  {
    /* put output in correct order */
    place = 0;
    while( place < coords.nThetaVals )
    {
      for( i = 0; i < numProc && place < coords.nThetaVals; i++ )
      {
        if( place < coords.nThetaVals - 1 )
          temp = coords.thetaBreaks[place+1] - coords.thetaBreaks[place];
        else
          temp = coords.nPix - coords.thetaBreaks[place];

        if(isComplex)
        {
          memcpy(mapc + coords.thetaBreaks[place], outputS + recvOff[i], temp*sizeof(fftw_complex));
          recvOff[i] += 2*temp;
        }
        else
        {
          memcpy(mapd + coords.thetaBreaks[place], outputS + recvOff[i], temp*sizeof(flt8));
          recvOff[i] += temp;
        }
        place++;
      }
    }
    free(outputS);
  }
  free(recvOff);
  free( nRecv );
}


almStruct createAlmStruct( int lmax, int numProc, int myRank)
{
  int i, j, k, l, m, n;
  fftw_complex *tempPtr;
  almStruct myAlm;
  
  myAlm.lmax = lmax;

  /* Figure out how many m's this processor will deal with. */
  /* Divide up lmax+1 non-negative m values. */
  myAlm.nm = (lmax+1)/numProc + ( ( myRank < (lmax+1)%numProc ) ? 1 : 0 );
  /* Now account for the negative m values. */
  if( myRank )
    myAlm.nm = 2*myAlm.nm;
  else 
    myAlm.nm = 2*myAlm.nm - 1;

  /* If this processor has no alms then return */
  if( !myAlm.nm )
  {
    myAlm.lmax = lmax;
    myAlm.nm = 0;
    myAlm.mList = NULL;
    myAlm.A = NULL;
    myAlm.allocPtrs[0] = myAlm.allocPtrs[1] = myAlm.allocPtrs[2] = myAlm.allocPtrs[3] = 0;
    return;
  }

  /* Allocate some memory */
  myAlm.mList = myAlm.allocPtrs[0] = (int *)malloc(sizeof(int)*myAlm.nm);
  myAlm.numData = myAlm.allocPtrs[1] = (int *)malloc(sizeof(int)*numProc);
  errorCheck(myRank, "createAlmStruct", (myAlm.allocPtrs[0] && myAlm.allocPtrs[1]), 1);

  
  /*******************************************************************************
  * Generate the list of m's which this processor will deal with.  Note if a     *
  * processor is computing for m, then it is also computing for -m, and the      *
  * non-negetive m's are distributed cyclicly.                                   *
  *******************************************************************************/
  if( myRank )
  {
    myAlm.mList[myAlm.nm/2] = myRank;
    myAlm.mList[myAlm.nm/2-1] = -myRank;
    i = myAlm.nm/2-2;
    j = myAlm.nm/2+1;
  }
  else
  {
    myAlm.mList[myAlm.nm/2] = 0;
    i = myAlm.nm/2-1;
    j = myAlm.nm/2+1;
  }

  for( ; j < myAlm.nm; i--, j++)
  {
    myAlm.mList[j] = myAlm.mList[j-1] + numProc;
    myAlm.mList[i] = -myAlm.mList[j];
  }
  /***********************************************************
  * Count number of alm's on every processor (this is used   *
  * in the call to MPI_Allgatherv()).                        *
  ***********************************************************/

  /* j is the number of left over non-negetive m's after dividing up evenly */
  j = (lmax+1)%numProc;

  /* Count for processor zero */
  k = (lmax+1)/numProc + (j ? 1 : 0);
  l = 2*k - 1;
  myAlm.numData[0] = l*(lmax+1) - numProc*(k*k - k);

  /* Count for processors between 1 and j-1 */
  k = (lmax+1)/numProc + 1;
  l = 2*k;
  m = l*(lmax+1) - numProc*(k*k - k);
  n = -2*k;
  for(i = 1; i < j; i++)
    myAlm.numData[i] = m + n*i;

  if( j == 0 ) 
    j = 1;
  /* Count for processors between j and numProc - 1 */
  k = (lmax+1)/numProc;
  l = 2*k;
  m = l*(lmax+1) - numProc*(k*k - k);
  n = -2*k;
  for(i = j; i < numProc; i++)
    myAlm.numData[i] = m + n*i;
   

  /* Allocate memory for alms so that you can access a_{l,m} with myAlm.A[m][l]. */
  myAlm.allocPtrs[2] = calloc((2*lmax+1), sizeof(fftw_complex *));
  errorCheck(myRank, "createAlmStruct", (int)myAlm.allocPtrs[2], 1);
  myAlm.A = ((fftw_complex **)(myAlm.allocPtrs[2])) + lmax;

  myAlm.allocPtrs[3] = calloc(myAlm.numData[myRank], sizeof(fftw_complex));
  errorCheck(myRank, "createAlmStruct", (int)myAlm.allocPtrs[3], 1);
 
  tempPtr = (fftw_complex*)myAlm.allocPtrs[3];
  for( i = 0; i < myAlm.nm; i++)
  {
    myAlm.A[myAlm.mList[i]] = tempPtr - abs(myAlm.mList[i]);
    tempPtr += lmax+1-abs(myAlm.mList[i]);
  }

  return(myAlm);
}

void destroyAlmStruct(almStruct *A )
{
  free(A->allocPtrs[3]);
  free(A->allocPtrs[2]);
  free(A->allocPtrs[1]);
  free(A->allocPtrs[0]);
  A->allocPtrs[0] = A->allocPtrs[1] = A->allocPtrs[2] = A->allocPtrs[3] = NULL;

  A->lmax = 0;
  A->nm = 0;
  A->mList = NULL;
  A->A = NULL;
}

void gatherAlm(almStruct *myAlm, fftw_complex *allAlm, int outputRoot, MPI_Comm theComm)
{
  int i, j, k, myRank, numProc, imax, lmax, nm;
  int *tempiPtr, *recvCounts, *displacements, *mListJumble;
  fftw_complex *almJumble, *tempcPtr;

  /* find out about processors */
  MPI_Comm_size(theComm, &numProc);
  MPI_Comm_rank(theComm, &myRank);

  lmax = myAlm->lmax;

  /* Allocate some memory. */
  recvCounts = (int*)malloc(sizeof(int)*numProc);
  displacements = (int*)malloc(sizeof(int)*numProc);
  errorCheck(myRank, "gatherAlm", (recvCounts && displacements), 1);

  if(outputRoot == -1 || myRank == outputRoot )
  {
    mListJumble = (int *)malloc((2*lmax+1)*sizeof(int));
    almJumble = (fftw_complex *)malloc((lmax+1)*(lmax+1)*sizeof(fftw_complex));
    errorCheck(myRank, "gatherAlm", (mListJumble && almJumble), 1);
  }
  else
  {
    mListJumble = NULL;
    almJumble = NULL;
  }

  /* calculate vectors for MPI_Gatherv() */
  recvCounts[0] = 2*myAlm->numData[0]; /* note that the factor of two is there because alm is complex */
  displacements[0] = 0;
  for( i = 1; i < numProc; i++ )
  {
    recvCounts[i] = 2*myAlm->numData[i];
    displacements[i] = displacements[i-1] + recvCounts[i-1];
  }

  if( outputRoot == -1 )
  {
    MPI_Allgatherv(myAlm->allocPtrs[3], recvCounts[myRank], MPI_DOUBLE, 
                   almJumble, recvCounts, displacements, MPI_DOUBLE, theComm);
  }
  else
  {
    MPI_Gatherv(myAlm->allocPtrs[3], recvCounts[myRank], MPI_DOUBLE, 
                almJumble, recvCounts, displacements, MPI_DOUBLE, outputRoot, theComm);
  }

  if( outputRoot == -1 || outputRoot == myRank )
  {
    tempiPtr = mListJumble;
    for( k = 0; k < numProc; k++ )
    {
      /* Figure out how many m's processor k will deal with. */
      /* Divide up lmax+1 non-negative m values. */
      nm = (lmax+1)/numProc + ( ( k < (lmax+1)%numProc ) ? 1 : 0 );
      /* Now account for the negative m values. */
      if( k )
        nm = 2*nm;
      else 
        nm = 2*nm - 1;

      if( k )
      {
        tempiPtr[nm/2] = k;
        tempiPtr[nm/2-1] = -k;
        i = nm/2-2;
        j = nm/2+1;
      }
      else
      {
        tempiPtr[nm/2] = 0;
        i = nm/2-1;
        j = nm/2+1;
      }

      for( ; j < nm; i--, j++)
      {
        tempiPtr[j] = tempiPtr[j-1] + numProc;
        tempiPtr[i] = -tempiPtr[j];
      }
      tempiPtr += nm; 
    }

    tempcPtr = almJumble;
    imax = 2*lmax + 1;
    for( i = 0; i < imax; i++ )
    {
      memcpy(allAlm+almIndex(abs(mListJumble[i]), mListJumble[i], lmax), tempcPtr, 
             sizeof(fftw_complex)*(lmax + 1 - abs(mListJumble[i])));
      tempcPtr += lmax + 1 - abs(mListJumble[i]);
    }
    free(mListJumble);
    free(almJumble);
  }
  
  free(displacements);
  free(recvCounts);
}


rowFFTstruct createRowFFTstruct(int numRows, flt8 *deltaPhi, int numProc, int myRank)
{
  /* numRows and deltaPhi are for the full pixelization */

  int i, j, k, thisRowLength;
  rowFFTstruct splitFFT;
  
  splitFFT.numRows = numRows/numProc + (myRank < numRows%numProc ? 1 : 0);

  splitFFT.numData = (int *)calloc(numProc, sizeof(int));
  splitFFT.rowList = (int *)malloc(sizeof(int)*splitFFT.numRows);
  splitFFT.rowLength = (int *)malloc(sizeof(int)*splitFFT.numRows);
  errorCheck(myRank, "createRowFFTstruct", (splitFFT.numData && splitFFT.rowList && splitFFT.rowLength), 1);

  for( i = 0; i < numRows; i++ )
  {
    j = i%numProc;
    k = i/numProc;
    thisRowLength = ccSHT_round(2*PI/fabs(deltaPhi[i]));
    if( j == myRank)
    {
      splitFFT.rowLength[k] = thisRowLength;
      splitFFT.rowList[k] = i;
    }
    splitFFT.numData[j] += thisRowLength;
  }

  splitFFT.data = (fftw_complex **)malloc(sizeof(fftw_complex*)*splitFFT.numRows);
  errorCheck(myRank, "createRowFFTstruct", (int)splitFFT.data, 1);
  splitFFT.data[0] = (fftw_complex *)malloc(sizeof(fftw_complex)*splitFFT.numData[myRank]);
  errorCheck(myRank, "createRowFFTstruct", (int)splitFFT.data[0], 1);

  for( i = 1; i < splitFFT.numRows; i++ )
    splitFFT.data[i] = splitFFT.data[i-1] + splitFFT.rowLength[i-1];

  return(splitFFT);
}


void destroyRowFFTstruct(rowFFTstruct *theStruct)
{
  theStruct->numRows = 0;
  free(theStruct->data[0]);
  free(theStruct->data);
  free(theStruct->rowLength);
  free(theStruct->rowList);
  free(theStruct->numData);
  theStruct->data = NULL;
  theStruct->rowLength = NULL;
  theStruct->rowList = NULL;
  theStruct->numData = NULL;
}


void calculateRowFFTstruct(coordStruct splitCoords, void *splitMap, int mapIsComplex, rowFFTstruct *splitFFT)
{
  int i, j, k, l, maxJ, nPhi;
  flt8 minDeltaPhi;
  flt8 *splitMapd;
  fftw_complex *splitMapc;
  fftw_complex *mapRow;
  fftw_plan **plans;

  if( mapIsComplex )
  {
    splitMapc = (fftw_complex *)splitMap;
    splitMapd = NULL;
  }
  else
  {
    splitMapd = (flt8 *)splitMap;
    splitMapc = NULL;
  }  
  
  /* figure out the smallest deltaPhi */

  minDeltaPhi = findSmallestEl(splitCoords.deltaPhi,splitCoords.nThetaVals);

  /* Allocate array */

  maxJ = ccSHT_round(2*PI/fabs(minDeltaPhi));
  mapRow = (fftw_complex *)calloc(maxJ,sizeof(fftw_complex));
  errorCheck(-1, "calculateRowFFTstruct", (int)mapRow, 1);

  plans = createPlans(splitCoords);
  k = 0;
  for( i = 0; i < splitCoords.nThetaVals; i++ )
  {
    maxJ = ccSHT_round(2*PI/fabs(splitCoords.deltaPhi[i]));
    
    /* Figure out number of pixels sampled in this row */
    if(i < splitCoords.nThetaVals - 1)
      nPhi = splitCoords.thetaBreaks[i+1] - splitCoords.thetaBreaks[i];
    else
      nPhi = splitCoords.nPix - splitCoords.thetaBreaks[i];

    /* Make a zero padded fftw_complex array corresponding to this row */
    j = 0;
    l = 0;
    while( j < nPhi)
    {
      if( k < splitCoords.nGaps && splitCoords.thetaBreaks[i] + j == splitCoords.gaps[2*k] )
      {
        memset(mapRow + l, 0, sizeof(fftw_complex)*splitCoords.gaps[2*k+1]);
        l += splitCoords.gaps[2*k+1];
        k++;
      }
      if(mapIsComplex)
        mapRow[l] = splitMapc[splitCoords.thetaBreaks[i]+j];
      else
      {
        mapRow[l].re = splitMapd[splitCoords.thetaBreaks[i]+j];
        mapRow[l].im = 0;
      }
      j++;
      l++;      
    }
    memset(mapRow + l, 0, sizeof(fftw_complex)*(maxJ-l));
    
    /* Do the Fourier transform */

    fftw_one(*(plans[i]), mapRow, splitFFT->data[i]);
  }
  destroyPlans(splitCoords, plans);
  free(mapRow);

}


void gatherRowFFTstruct(rowFFTstruct *myFFT, rowFFTstruct *allFFT, MPI_Comm theComm)
{
  int myRank, numProc, i, j;
  int *recvCounts, *displacements;
  fftw_complex *tempFFTbuffer, *tempPtr;

  MPI_Comm_size(theComm, &numProc);
  MPI_Comm_rank(theComm, &myRank);

  /* Allocate some memory */
  tempFFTbuffer = (fftw_complex *)malloc(sizeof(fftw_complex)*(allFFT->numData[0]));
  recvCounts = (int *)malloc(sizeof(int)*numProc);
  displacements = (int *)malloc(sizeof(int)*numProc);
  errorCheck(myRank, "assembleRowFFTstruct", (tempFFTbuffer && recvCounts && displacements), 1);

  /* Create vectors for MPI_Allgatherv */  
  recvCounts[0] = myFFT->numData[0]*2;
  displacements[0] = 0;
  for( i = 1; i < numProc; i++)
  {
    recvCounts[i] = myFFT->numData[i]*2;
    displacements[i] = displacements[i-1] + recvCounts[i-1];
  }

  /* Gather the FFT's from each processor */
  MPI_Allgatherv( myFFT->data[0], recvCounts[myRank], MPI_DOUBLE, 
                  tempFFTbuffer, recvCounts, displacements, MPI_DOUBLE, theComm);

  /* Put them in the correct order */
  tempPtr = tempFFTbuffer;
  for( i = 0; i < numProc; i++ )
  {
    for( j = i; j < allFFT->numRows; j += numProc)
    {
      memcpy(allFFT->data[j], tempPtr, sizeof(fftw_complex)*allFFT->rowLength[j]);
      tempPtr += allFFT->rowLength[j];
    }
  }

  /* Free up memory */
  free(recvCounts);
  free(displacements);
  free(tempFFTbuffer);
}

void gatherSmallQlm(int *lStartLocal, int *lStartAll, int **lStart, flt8 *QstartLocal, 
                    flt8 *QstartAll, flt8 **Qstart, int nThetaValsAll, int lmax, MPI_Comm theComm )
{
  /* 
     Note that this function could be improved if the lStart data, and
     the Qstart data were packed up and communicated together, instead of 
     in two seperate calls to MPI_Allgatherv.  
  */

  int myRank, numProc, i, j;
  int *lStartPtr, *recvCounts, *displacements;
  flt8 *QstartPtr;

  MPI_Comm_size(theComm, &numProc);
  MPI_Comm_rank(theComm, &myRank);

  /* Set up Qstart and lStart data structures */
  memset(Qstart, 0, sizeof(flt8*)*nThetaValsAll);
  memset(lStart, 0, sizeof(int*)*nThetaValsAll);

  QstartPtr = QstartAll;
  lStartPtr = lStartAll;
  for( i = 0; i < numProc; i++ )
  {
    for(j = i; j < nThetaValsAll; j += numProc )
    {
      Qstart[j] = QstartPtr;
      QstartPtr += 2*(lmax+1);
      lStart[j] = lStartPtr;
      lStartPtr += lmax+1;
    }
  }

  /* Allocate memory for MPI_Allgatherv call */
  recvCounts = (int*)malloc(sizeof(int)*numProc);
  displacements = (int*)malloc(sizeof(int)*numProc);

  errorCheck(myRank, "gatherSmallQlm", (recvCounts && displacements), 1);
 
  /* Set up vectors for MPI_Allgatherv() */
  j = (nThetaValsAll/numProc + 1)*(lmax+1);
  for( i = 0; i < nThetaValsAll%numProc; i++)
    recvCounts[i] = j;

  j = (nThetaValsAll/numProc)*(lmax+1);
  for( i = nThetaValsAll%numProc; i < numProc; i++ )
    recvCounts[i] = j;

  displacements[0] = 0;
  for( i = 1; i < numProc; i++ )
    displacements[i] = displacements[i-1] + recvCounts[i-1];

  /* Gather the lStart vectors from all of the processors */
  MPI_Allgatherv(lStartLocal, recvCounts[myRank], MPI_INT, lStartAll, recvCounts, displacements, MPI_INT, theComm);

  /* Fix MPI_Allreducev vectors for the Qstart vectors which are twice as long */
  recvCounts[0] *= 2;
  for( i = 1; i < numProc; i++ )
  {
    recvCounts[i] *=2;
    displacements[i] = displacements[i-1] + recvCounts[i-1];
  }

  /* Gather Qstart vectors from all of the processors */
  MPI_Allgatherv(QstartLocal, recvCounts[myRank], flt8mpiType, QstartAll, recvCounts, displacements, flt8mpiType, theComm);

  /* Free up some memory */
  free(displacements);
  free(recvCounts);
}
