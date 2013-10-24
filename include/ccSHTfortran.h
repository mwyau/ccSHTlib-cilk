
#define FORTRAN_INTEGER_4 int4
#define FORTRAN_INTEGER_8 int8
#define FORTRAN_REAL_8 flt8
#define FORTRAN_COMPLEX_8 fftw_complex


#ifndef f77func
#if defined (F77_APPEND_UNDERSCORE)
#  if defined (F77_UPPERCASE_NAMES)
#    define f77func(f, F) F##_
#  else
#    define f77func(f, F) f##_
#  endif
#else
#  if defined (F77_UPPERCASE_NAMES)
#    define f77func(f, F) F
#  else
#    define f77func(f, F) f
#  endif
#endif
#endif


#define forwardsht_real_fortran          f77func(forwardsht_real_fortran,        FORWARDSHT_REAL_FORTRAN)
#define forwardsht_complex_fortran       f77func(forwardsht_complex_fortran,     FORWARDSHT_COMPLEX_FORTRAN)
#define backwardsht_real_fortran         f77func(backwardsht_real_fortran,       BACKWARDSHT_REAL_FORTRAN)
#define backwardsht_complex_fortran      f77func(backwardsht_complex_fortran,    BACKWARDSHT_COMPLEX_FORTRAN)
#define convertradec2coords_fortran      f77func(convertradec2coords_fortran,    CONVERTRADEC2COORDS_FORTRAN)
#define destroycoords_fortran            f77func(destroycoords_fortran,          DESTROYCOORDS_FORTRAN)

#define forwardshtmpi_real_fortran       f77func(forwardshtmpi_real_fortran,     FORWARDSHTMPI_REAL_FORTRAN)
#define forwardshtmpi_complex_fortran    f77func(forwardshtmpi_complex_fortran,  FORWARDSHTMPI_COMPLEX_FORTRAN)
#define backwardshtmpi_real_fortran      f77func(backwardshtmpi_real_fortran,    BACKWARDSHTMPI_REAL_FORTRAN)
#define backwardshtmpi_complex_fortran   f77func(backwardshtmpi_complex_fortran, BACKWARDSHTMPI_COMPLEX_FORTRAN)
#define calculateqlm_fortran             f77func(calculateqlm_fortran,           CALCULATEQLM_FORTRAN)
#define lmcalculations_fortran           f77func(lmcalculations_fortran,         LMCALCULATIONS_FORTRAN)
#define calculatesmallqlm_fortran        f77func(calculatesmallqlm_fortran,      CALCULATESMALLQLM_FORTRAN)

