#define TEST_EXP	// <-- select test function here
const int NTE = 4;	// <-- select Taylor expansion order here
const int NX = 500;	// <-- select X discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic

////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

#include "tmodel.hpp"
typedef mc::TModel<MC> TMMC;
typedef mc::TVar<MC> TVMC;

#include "mcfadbad.hpp"
#include "fadiff.h"
typedef fadbad::F<TVMC> FTVMC;

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_POW )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return 1. - 5.*x - pow(x,2)/2. + pow(x,3)/3.;
}

#elif defined( TEST_EXP )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
const double Xref =  0.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-pow(x,2));
}

#elif defined( TEST_TRIG )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
const double Xref =  PI/4.;	// <-- linearization point
template <class T>
T myfunc
( const T&x )
{
  return tan(cos(x*atan(x)));
}

#elif defined( TEST_TRIG2 )
const double XL   = PI/6.;	// <-- X range lower bound
const double XU   = PI/3.;	// <-- X range upper bound
const double Xref = PI/4.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return sin(pow(x,-3))*cos(sqrt(x));
}
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{

#ifdef SAVE_RESULTS
  ofstream res( "TM-1D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 

    // Define Taylor model environment
    TMMC mod( 1, NTE );

    // <-- set options here -->
    MC::options.MVCOMP_USE = true;
    mod.options.BOUNDER_TYPE = TMMC::Options::BERNSTEIN;
    mod.options.BOUNDER_ORDER = 0;
    mod.options.BERNSTEIN_USE = true;

    // Define variable X, and evaluate Taylor model
    TVMC TVMCX( &mod, 0, MC(I(XL,XU),Xref).sub(1,0) );
    TVMC TVMCF = myfunc( TVMCX );
    std::cout << "\nMcCormick-Taylor model of f(x):" << TVMCF;

    // Evaluate Taylor models of first derivatives using FADBAD++
    FTVMC FTVMCX = TVMCX;
    FTVMCX.diff(0,1);
    FTVMC FTVMCF = myfunc( FTVMCX );
    std::cout << "\nMcCormick-Taylor model of df/dx: " << FTVMCF.d(0) << std::endl;

    // Repeated calculations at grid points (for display)
    for( int iX=0; iX<NX; iX++ ){ 

      double DX = XL+iX*(XU-XL)/(NX-1.);
      double DF = myfunc( DX );
      TVMCX.set( &mod, 0, MC(I(XL,XU),DX).sub(1,0) );
      TVMCF = myfunc( TVMCX );
      double PF = TVMCF.P( &DX );
      MC MCF = PF + TVMCF.R();

#ifdef SAVE_RESULTS
      res << std::setw(14) << DX << std::setw(14) << std::setw(14) << DF
          << std::setw(14) << PF << std::setw(14) << MCF.l()  << setw(14) << MCF.u()
          << std::setw(14) << MCF.cv() << setw(14) << MCF.cc()
          << std::setw(14) << TVMCF.B().l()  << std::setw(14) << TVMCF.B().u()
          << std::setw(14) << TVMCF.B().cv() << std::setw(14) << TVMCF.B().cc()
          << std::setw(14) << TVMCF.B().cvsub(0) << std::setw(14) << TVMCF.B().ccsub(0)
          << std::endl;
#endif
    }

  }
  
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( MC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in McCormick relaxation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
  catch( TMMC::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in Taylor model computation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

#ifdef SAVE_RESULTS
  res.close();
#endif
  return 0;
}
