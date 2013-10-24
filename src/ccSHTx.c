/*******************************************************************************

     ccSHTx(1)                                               ccSHTx(1)
     
     
     
     NAME
            ccSHTx, ccSHTmpix - Spherical Harmonic Transform
     
     SYNOPSIS
            ccSHTx -n numPix -l lmax -f fFile -a almFile -p pixFile -d
            directionFlag [-w pixArea] [-q quadWeightFile] [-c]
     
            ccSHTmpix -n numPix -l lmax -f fFile -a almFile -p pixFile
            -d directionFlag [-w pixArea] [-q quadWeightFile] [-c]
     
     
     DESCRIPTION
            These  are  executable  applications which make use of the
            ccSHT spherical harmonic transform library.  ccSHTx is the
            executable  binary  which  performs  a  spherical harmonic
            transform in  serial,  and  ccSHTmpix  is  the  executable
            binary  which does the transform in parallel.  The calling
            sequence for these two executables is  exactly  the  same.
            There  are six required parameters which must be specified
            and two optional parameters.  Please see the man page  for
            ccSHT(1) for more information about the ccSHT package.
     
            -n numPix
                   Specifies  the  number of pixels in the discretiza-
                   tion of the sphere.
     
            -l lmax
                   Sets the band cut off for the  SHT.   This  is  the
                   largest value of l for which the spherical harmonic
                   coefficients are not assumed to be zero.
     
            -f mapFile
                   Specifies the name of the file which is to  contain
                   the  pixel  domain  vector.   When  doing a forward
                   transform this file is read in, and  when  doing  a
                   backward  transform  this  file  is created.  It is
                   formatted as a raw binary file of eight byte float-
                   ing point reals or complex pairs.  The pixel domain
                   vector is assumed to be real valued unless  the  -c
                   option  is  specified  in  which case the vector is
                   interpreted as an  array  of  complex  pairs.   The
                   length of the vector is specified by the -n parame-
                   ter.
     
            -a almFile
                   Specifies the name of the file which is to  contain
                   the  spherical harmonic coefficients.  When doing a
                   forward transform this file is  created,  and  when
                   doing  a  backward  transform this file is read in.
                   It is formatted as a raw binary file of eight  byte
                   floating  point  complex pairs.  The number of com-
                   plex pairs is (lmax+1)*(lmax+1) where lmax is spec-
                   ified by the -l parameter.
     
            -p pixFile
                   Specifies  the  name of the file which contains the
                   right ascension (RA) and declination (DEC)  coordi-
                   nates  for each of the pixels in the discretization
                   of the sphere.  This file is  formatted  as  a  raw
                   binary file of four byte floating point pairs of RA
                   and DEC coordinates.  There are as many  coordinate
                   pairs   in  this  file  as  the  number  of  pixels
                   specified by the -n parameter.  If there  are  gaps
                   in  the pixelization then an additional file may be
                   necessary.  This file must have the  same  name  as
                   the  pixel  file,  only with a .dra appended to it.
                   For instance if the pixel file was  called  pixels,
                   then  the  additional  file  would  be  called pix-
                   els.dra.  This file specifies the latitudinal spac-
                   ing  between  pixels  on every row.  This file is a
                   raw binary file.  The first number in the file is a
                   four  byte  integer  which  specifies the number of
                   rows which are described in the file.  This integer
                   is  followed by as many pairs of four byte floating
                   point numbers.  The first of  each  pair  specifies
                   the  declination  of the row, and the second speci-
                   fies the change in right ascension  between  neigh-
                   boring  pixels  in the given row.  These pairs must
                   be ordered in the same way as the rows in the pixel
                   location  file,  but  there may be more rows in the
                   .dra file than in the pixel location file.
     
            -d directionFlag
                   Specifies the direction of the transform to be per-
                   formed.   If  (f,  F, or 1) are given then the code
                   performs a forward  transform,  deriving  spherical
                   harmonic coefficients from a pixel domain function.
                   If (b, B, or -1) are given then the code performs a
                   backward  transform,  deriving a pixel domain func-
                   tion from spherical harmonic coefficients.
     
            -w pixArea
                   This optional parameter specifies the area  of  one
                   pixel.  This is the constant quadrature weight used
                   to compute the integral in the  forward  transform.
                   If  this number is not specified, than the integral
                   is computed as a sum,  (i.e.  the  weight  used  is
                   one).   If  a  weight file is specified with the -q
                   option then the -w option should not be used.
     
            -q quadWeightFile
                   This optional parameter specifies a file  which  is
                   to  contain  the  quadrature weights to compute the
                   integral in the forward transform.  This file is  a
                   raw  binary  file of eight byte floating point num-
                   bers that specify a weight for every pixel.
     
            -c     If this option is given then the pixel domain func-
                   tion  is  assumed  to be complex valued.  The pixel
                   domain function is either the input to the  forward
                   transform  or  the  output from the backward trans-
                   form.
     
     COPYRIGHT
            Version 1.03 July 2003
     
            Copyright (C) 2003  Christopher M. Cantalupo
     
            ccSHTx and ccSHTmpix are free  software;  you  can  redis-
            tribute them and/or modify them under the terms of the GNU
            General Public License as published by the  Free  Software
            Foundation;  either  version 2 of the License, or (at your
            option) any later version.
     
            This program is distributed in the hope that  it  will  be
            useful, but WITHOUT ANY WARRANTY; without even the implied
            warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR
            PURPOSE.   See  the  GNU  General  Public License for more
            details.
     
            Christopher Cantalupo <cmc@nersc.gov>
     
            Send bug reports or comments to the above address.
     
     
     
                                July 15 2003                 ccSHTx(1)
     
*******************************************************************************/
        

/*******************************************************************************
*   Version 1.03 July 2003                                                     *
*                                                                              *
*   Copyright (C) 2003  C.M. Cantalupo                                         *
*                                                                              *
*   ccSHTx and ccSHTmpix are free software; you can redistribute it and/or     *
*   modify it under the terms of the GNU General Public License as published by*
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

#include <sys/time.h>

#ifdef USE_MPI
#include "ccSHTmpi.h"
#else
#include "ccSHT.h"
#endif

int main(int argc, char **argv)
{
  FILE *fFile, *pixFile, *draFile, *almFile, *weightFile;
  char fFileName[STRLEN];
  char pixFileName[STRLEN];
  char draFileName[STRLEN];
  char almFileName[STRLEN];
  char weightFileName[STRLEN];
  char *strPtr;
  int dirFlag, numPix, lmax, fIsComplex, i;
  int4 ndra;
  flt8 *fr, *raDecArray, *dra, *weights;
  flt8 pixSize, timeLaps, flopRate;
  fftw_complex *fc, *alm;
  void *fv;
  struct timeval tic, toc;
  struct timezone tz;
  coordStruct coords;

  int myRank = 0;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

  /* Initialize the input parameters. */
  numPix = 0;
  lmax = 0;
  fFileName[0] = 0;
  almFileName[0] = 0;
  pixFileName[0] = 0;
  dirFlag = 0;
  pixSize = 1;
  weightFileName[0] = 0;
  fIsComplex = 0;  

  /* Parse input parameters from command line. */
  for( i = 0; i < argc; i++ )
  {
    if( argv[i][0] == '-')
    {
      if(argv[i][2])
        strPtr = argv[i] + 2;
      else if( i+1 < argc )
        strPtr = argv[i+1];
      else 
        strPtr = NULL;

      switch( argv[i][1] )
      {
        case 'n':
          numPix = atoi(strPtr);
	  break;
        case 'l':
          lmax = atoi(strPtr);
          break;
        case 'f':
          strcpy(fFileName, strPtr);
          break;
        case 'a':
          strcpy(almFileName, strPtr);
          break;
        case 'p':
          strcpy(pixFileName, strPtr);
          break;
        case 'd':
          if( strPtr[0] == 'f' || strPtr[0] == 'F')
            dirFlag = 1;
          else if( strPtr[0] == 'b' || strPtr[0] == 'B')
            dirFlag = -1;
          else
            dirFlag = atoi(strPtr);
          break;
        case 'w':
          pixSize = atof(strPtr);
          break;
        case 'q':
          strcpy(weightFileName, strPtr);
          break;
        case 'c':
          fIsComplex = 1;
          break;
      }
    }
  }

  if( !myRank )
  {
    fprintf( stdout, "  numPix = %i\n", numPix);
    fprintf( stdout, "  lmax = %i\n", lmax);
    fprintf( stdout, "  constant quadrature weight = %e\n", pixSize);

    fprintf( stdout, "  Pixel space function file name is %s.\n", fFileName);
    fprintf( stdout, "  Spherical harmonic coefficient file name is %s.\n", almFileName);
    fprintf( stdout, "  Pixel location file name is %s.\n", pixFileName);
    if( weightFileName[0])
      fprintf( stdout, "  Variable quadrature weight file name is %s.\n", weightFileName);
    else
      fprintf( stdout, "  No variable quadrature weight file specified.\n");

    if( fIsComplex )
      fprintf( stdout, "  Pixel space function is complex valued.\n");
    else
      fprintf( stdout, "  Pixel space function is real valued.\n");

    if( dirFlag == -1 )
      fprintf( stdout, "  Computing backward transform.\n");
    else if( dirFlag == 1 )
      fprintf( stdout, "  Computing forward transform.\n");
    else
      fprintf( stdout, "  Direction flag does not make sense.\n");

  }

  /* Make sure that all of the required parameters were specified. */
  errorCheck(myRank, "Please specify the number of pixels with the -n arguement", (numPix > 0), 0);
  errorCheck(myRank, "Please specify the maximum l value (lmax) with the -l arguement", (lmax > 0), 0);
  errorCheck(myRank, "Please specify the pixel space function file name with the -f arguement", fFileName[0], 0);
  errorCheck(myRank, "Please specify the spherical harmonic coefficient file name with the -a arguement", almFileName[0], 0);
  errorCheck(myRank, "Please specify the pixel location file name with the -p arguement", pixFileName[0], 0);
  errorCheck(myRank, "Please specify the direction of the transform to be computed with the -d arguement (backward: B, b, -1; forward: F, f, 1)", (dirFlag == -1 || dirFlag == 1), 0);


  /* Allocate memory for raDecArray. */
  raDecArray = (flt8*)malloc(numPix*2*sizeof(flt8));
  errorCheck(myRank,"main",(int)raDecArray, 1);


  /* Read in raDecArray from file. */
  if( !myRank )
  {
    pixFile = fopen(pixFileName,"r");
    errorCheck(myRank,pixFileName,(int)pixFile,2);
  
    input_ftod(pixFile,raDecArray,numPix*2);
    fclose(pixFile);
  }

#ifdef USE_MPI
  errorCheck(myRank, " ", 1, 0);
  MPI_Bcast( raDecArray, numPix*2, flt8mpiType, 0, MPI_COMM_WORLD );
#endif

  /* Deal with dra file */
  ndra = 0;
  dra = NULL;
  if( !myRank )
  {
    /* Look for dra file. */
    sprintf(draFileName, "%s.dra", pixFileName);
    draFile = fopen(draFileName,"r");
    if( draFile )
    {
      /* If dra file exists read in the number of rows of pixels. */
      fread(&ndra, sizeof(int4), 1, draFile);
    }
  }

#ifdef USE_MPI
  MPI_Bcast( &ndra, 1,  int4mpiType, 0, MPI_COMM_WORLD );
#endif

  if( ndra )
  {
    dra = malloc(sizeof(flt8)*ndra*2);
    errorCheck(myRank, "main", (int)dra, 1);
    if( !myRank )  
    {  
      input_ftod(draFile, dra, ndra*2);
      fclose(draFile);
    }

#ifdef USE_MPI
    errorCheck(myRank, " ", 1, 0);
    MPI_Bcast( dra, ndra*2, flt8mpiType, 0, MPI_COMM_WORLD );
#endif

  }
   
  /* create coordStruct for the pixelization */
  coords = convertRaDec2Coords(pixSize, numPix, raDecArray, ndra, dra);

  /* Free memory. */
  free(raDecArray);
  if(dra)
    free(dra);
  
  if( myRank == 0 || dirFlag == 1 )
  {
    /* Allocate memory for pixel function. */
    if( fIsComplex )
    {
      fc = (fftw_complex*)malloc(numPix*sizeof(fftw_complex));
      fr = NULL;
      errorCheck(myRank, "main", (int)fc, 1);
    }
    else
    {
      fr = (flt8 *)malloc(numPix*sizeof(flt8));
      fc = NULL;
      errorCheck(myRank, "main", (int)fr, 1);
    }
  }
  else
  {
    fr = NULL;
    fc = NULL;
  }

  if( myRank == 0 || dirFlag == -1 )
  {
    /* Allocate memory for alms. */
    alm = (fftw_complex *)malloc((lmax+1)*(lmax+1)*sizeof(fftw_complex));
    errorCheck(myRank, "main", (int)alm, 1);
  }
  else
    alm = NULL;

  if(dirFlag == 1)
  { 
    /* FORWARD TRANSFORM */

    /* Read in pixel domain function. */
    if( !myRank )
    {
      fFile = fopen(fFileName,"r");
      errorCheck(myRank, fFileName, (int)fFile, 2);

      if(fIsComplex)
        fread(fc, sizeof(fftw_complex), numPix, fFile);
      else
        fread(fr, sizeof(flt8), numPix, fFile);

      fclose(fFile);
    }

#ifdef USE_MPI
    if(fIsComplex)
      MPI_Bcast( fc, numPix*2, flt8mpiType, 0, MPI_COMM_WORLD );
    else
      MPI_Bcast( fr, numPix, flt8mpiType, 0, MPI_COMM_WORLD );
#endif

    /* See if should use quadrature weights. */
    if( weightFileName[0] )
    {
      /* Allocate memory for weights. */
      weights = (flt8 *)malloc(numPix*sizeof(flt8));
      errorCheck(myRank, "main", (int)weights, 1);

      /* Read in quadrature weights. */
      if( !myRank )
      {
        weightFile = fopen(weightFileName, "r");
        errorCheck(myRank, weightFileName, (int)weightFile, 2);
        fread(weights, sizeof(flt8), numPix, weightFile);
        fclose(weightFile);
      }
#ifdef USE_MPI
      errorCheck(myRank, " ", 1, 0);
      MPI_Bcast( weights, numPix, flt8mpiType, 0, MPI_COMM_WORLD );
#endif
     
      /* Multiply by quadrature weighting scheme. */
      if(fIsComplex)
        for( i = 0; i < numPix; i++ )
        {
          fc[i].re *= weights[i];
          fc[i].im *= weights[i];
        }
      else
        for( i = 0; i < numPix; i++ )
          fr[i] *= weights[i];

      free(weights);
    }


    /* Do transform and time the operation */
    if(fIsComplex)
    {
      gettimeofday(&tic, &tz);
#ifdef USE_MPI
      forwardSHTmpi(fc, 1, coords, lmax, alm, 0, MPI_COMM_WORLD);
#else
      forwardSHT(fc,1, coords, lmax, alm);
#endif
      gettimeofday(&toc, &tz);
    }
    else
    {
      /* Do transform and time the operation */
      gettimeofday(&tic, &tz);
#ifdef USE_MPI
      forwardSHTmpi(fr, 0, coords, lmax, alm, 0, MPI_COMM_WORLD);
#else
      forwardSHT(fr,0, coords, lmax, alm);
#endif
      gettimeofday(&toc, &tz);
    }

    /* Write alm's to output file. */
    if( !myRank )
    {
      almFile = fopen( almFileName, "w");
      errorCheck(myRank, almFileName, (int)almFile,3);
      fwrite((flt8*)alm, sizeof(fftw_complex), (lmax+1)*(lmax+1), almFile); 
      fclose(almFile);
    }
  }
  else if(dirFlag == -1)
  {
    /* BACKWARD TRANSFORM */
    /* Read in alm's. */
    if( !myRank )
    {
      almFile = fopen(almFileName, "r"); 
      errorCheck(myRank, almFileName, (int)almFile, 2);
      fread(alm, sizeof(fftw_complex), (lmax+1)*(lmax+1), almFile);
      fclose(almFile);
    }
#ifdef USE_MPI
    MPI_Bcast( alm, 2*(lmax+1)*(lmax+1), flt8mpiType, 0, MPI_COMM_WORLD );
#endif

    if( fIsComplex )
    {
      /* Do backward transform and time the computation. */
      gettimeofday(&tic, &tz);
#ifdef USE_MPI
      backwardSHTmpi(alm, coords, lmax, fc, 1, 0, MPI_COMM_WORLD);
#else
      backwardSHT( alm, coords, lmax, fc, 1);
#endif
      gettimeofday(&toc, &tz);

      if( !myRank )
      {
        /* Write the output pixel domain vector to file */
        fFile = fopen(fFileName, "w");
        errorCheck(myRank, fFileName, (int)fFile, 3);
        fwrite(fc, sizeof(fftw_complex), numPix, fFile);
        fclose(fFile);
      }
    }
    else
    {
      gettimeofday(&tic, &tz);
#ifdef USE_MPI
      backwardSHTmpi( alm, coords, lmax, fr, 0, 0, MPI_COMM_WORLD);
#else
      backwardSHT( alm, coords, lmax, fr, 0);
#endif
      gettimeofday(&toc, &tz);

      /* Write the output pixel domain vector to file */
      if( !myRank )
      {
        fFile = fopen(fFileName, "w");
        errorCheck(myRank, fFileName, (int)fFile, 3);
        fwrite(fr, sizeof(flt8), numPix, fFile);
        fclose(fFile);
      }
    }
  }
  else
    errorCheck(myRank, "direction flag must be (f, F, or 1)  for forward, and (b, B, or -1) for backward\n\n",0,0);

  /* Record and print the computation time. */
  if( !myRank )
  {
    timeLaps = (toc.tv_sec - tic.tv_sec) + (toc.tv_usec - tic.tv_usec)*1e-6;
    fprintf(stderr, "Did transform, computation took %f seconds.\n", timeLaps);
    if( fIsComplex )
      flopRate = (13.0/2.0*coords.nThetaVals*lmax*lmax)/(timeLaps)*1e-6;
    else
      flopRate = (9.0/2.0*coords.nThetaVals*lmax*lmax)/(timeLaps)*1e-6;
    fprintf(stderr, "This implies a computation rate of %f Mflop/s.\n\n", flopRate);
  }

  /* Clean up. */
  if( fIsComplex )
    free(fc);
  else
    free(fr);
  free(alm);
  destroyCoords(&coords);
#ifdef USE_MPI
  MPI_Finalize();
#endif


}
