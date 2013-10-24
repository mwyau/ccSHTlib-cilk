#include <stdlib.h>
#include <stdio.h>

/*******************************************************************************
*  makeSizedType                                                               *
*    This simple executable creates a header file which explicitly defines     *
*    integer and floating point types with specific sizes (i.e. four or eight  *
*    bytes).  The file will be called "include/sizedType.h".                   *
*******************************************************************************/

/*******************************************************************************
*   Version 1.03 July 2003                                                     *
*                                                                              *
*   Copyright (C) 2003  C.M. Cantalupo                                         *
*                                                                              *
*   makeSizedType is free software; you can redistribute it and/or modify      *
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



int main(int argc, char **argv)
{
  FILE *fid;
  char fileName[256] = "include/sizedType.h";

  fid = fopen( fileName , "w");

  /* Print a comment */
  fprintf(fid, "/* File created automatically by makeSizedType */\n\n");
  
  /* Print no repeat include header */
  fprintf(fid, "#ifndef _H_sizedType\n");
  fprintf(fid, "#define _H_sizedType\n\n\n");


  /* Find four byte integer */
  if(sizeof(int) == 4)
  {
    fprintf(fid, "#define int4 int\n");
    fprintf(fid, "#define int4mpiType MPI_INT\n");
  }
  else if(sizeof(short int) == 4)
  {
    fprintf(fid, "#define int4 short int\n");
    fprintf(fid, "#define int4mpiType MPI_SHORT\n");
  }
  else if(sizeof(long int) == 4)
  {
    fprintf(fid, "#define int4 long int\n");
    fprintf(fid, "#define int4mpiType MPI_LONG\n");
  }
  else
  {
    fclose(fid);
    remove(fileName);
    fprintf(stderr, "ERROR:  Can not find a four byte integer type.\n\n");
    return(1);
  }


  /* Find eight byte integer */
  if(sizeof(long int) == 8)
  {
    fprintf(fid, "#define int8 long int\n");
    fprintf(fid, "#define int8mpiType MPI_LONG\n");
  }
  else if(sizeof(long long int) == 8)
  {
    fprintf(fid, "#define int8 long long int\n");
    fprintf(fid, "#define int8mpiType MPI_LONG_LONG_INT\n");
  }
  else if(sizeof(int) == 8)
  {
    fprintf(fid, "#define int8 int\n");
    fprintf(fid, "#define int8mpiType MPI_INT\n");
  }
  else
  {
    fclose(fid);
    remove(fileName);
    fprintf(stderr, "ERROR:  Can not find a eight byte integer type.\n\n");
    return(1);
  }


  /* Find four byte float */
  if(sizeof(float) == 4)
  {
    fprintf(fid, "#define flt4 float\n");
    fprintf(fid, "#define flt4mpiType MPI_FLOAT\n");
  }
  else if(sizeof(double) == 4)
  {
    fprintf(fid, "#define flt4 double\n");
    fprintf(fid, "#define flt4mpiType MPI_DOUBLE\n");
  }
  else
  {
    fclose(fid);
    remove(fileName);
    fprintf(stderr, "ERROR:  Can not find four byte floating point number type.\n\n");
    return(1);
  }


  /* Find eight byte float */
  if( sizeof(double) == 8 )
  {
    fprintf(fid, "#define flt8 double\n");
    fprintf(fid, "#define flt8mpiType MPI_DOUBLE\n\n\n");
  }
  else if( sizeof(long double) == 8 )
  {
    fprintf(fid, "#define flt8 long double\n");
    fprintf(fid, "#define flt8mpiType MPI_LONG_DOUBLE\n\n\n");
  }
  else
  {
    fclose(fid);
    remove(fileName);
    fprintf(stderr, "ERROR:  Can not find eight byte floating point type.\n\n");
    return(1);
  }

  fprintf(fid, "#endif\n\n");


  fclose(fid);
}
