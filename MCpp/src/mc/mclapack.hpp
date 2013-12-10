// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

#ifndef MC__MCLAPACK_H
#define MC__MCLAPACK_H

#include <iostream>
#include <iomanip>
#include <string>

#undef  MC__DEBUG_EIGEN

extern "C" void dsyev_
( const char*jobz, const char*uplo, const unsigned int*n, double*a,
  const unsigned int*lda, double*w, double*work, const int*lwork, int*info );

namespace mc
{
 
void pause()
{
  int tmp;
  std::cout << "ENTER <1> TO CONTINUE" << std::endl;
  std::cin  >> tmp;
}

template< typename U > inline void display
( const unsigned int m, const unsigned int n, U*&a, const unsigned int lda,
  const std::string&stra, std::ostream&os )
{
  os << stra << " =" << std::endl << std::scientific
     << std::setprecision(5);
  for( unsigned int im=0; a && im<m; im++ ){
    for( unsigned int in=0; in<n; in++ ){
      os << a[in*lda+im] << "  ";
    }
    os << std::endl;
  }
  os << std::endl;

  if( os == std::cout || os == std::cerr ) pause();
}

//! @brief Wrapper to LAPACK function <tt>_dsyev<tt> doing eigenvalue decomposition of symmetric <tt>n</tt>-by-<tt>n</tt> matrix <tt>A</tt>. A pointer to an array containing the eignenvalues is returned - NULL pointer if the eigenvalue decomposition was unsuccessful. If <tt>eigv</tt> is set to <tt>true</tt> eigenvector too are computed and returned in <tt>A</tt>, which is therefore alterred. Note that <tt>A</tt> is altered even if eigenvector are not computed.
inline double* dsyev_wrapper
( const unsigned int n, double*A, const bool eigv=false )
{
  int info;
  double*D = new double[n];
#ifdef MC__DEBUG_EIGEN
  display( n, n, A, n, "Matrix A", std::cout );
#endif

  // get optimal size
  char JOBZ = (eigv?'V':'N'), UPLO = 'U';
  double worktmp;
  int lwork = -1;
  dsyev_( &JOBZ, &UPLO, &n, A, &n, D, &worktmp, &lwork, &info );

  // perform eigenvalue decomposition
  lwork = (int)worktmp;
  double*work = new double[lwork];
  dsyev_( &JOBZ, &UPLO, &n, A, &n, D, work, &lwork, &info );
#ifdef MC__DEBUG_EIGEN
  if( eigv ) TModel<T>::_display( n, n, A, n, "Matrix U", std::cout );
  display( 1, n, D, 1, "Matrix D", std::cout );
#endif
  delete[] work;

#ifdef MC__DEBUG_EIGEN
  std::cout << "INFO: " << info << std::endl;
  pause;
#endif
  if( info ){ delete[] D; return 0; }
  return D;
}

} // namespace mc

#endif
