// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

#ifndef MC__MCPROFIL_HPP
#define MC__MCPROFIL_HPP

#include "mcop.hpp"
#include "mcfunc.hpp"
#include <Interval.h>
#include <Functions.h>
#include <Constants.h>
namespace mc
{
//! @brief Specialization of the structure mc::Op to allow usage of the type INTERVAL of <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> as a template parameter in the classes mc::McCormick, mc::TModel, mc::TVar, and mc::SpecBnd
template <> struct Op< ::INTERVAL >
{
  typedef ::INTERVAL T;
  static T point( const double c ) { return T(c); }
  static T zeroone() { return T(0.,1.); }
  static void I(T& x, const T& y) { x = y; }
  static double l(const T& x) { return ::Inf(x); }
  static double u(const T& x) { return ::Sup(x); }
  static double abs (const T& x) { return ::Abs(x);  }
  static double mid (const T& x) { return ::Mid(x);  }
  static double diam(const T& x) { return ::Diam(x); }
  static T inv (const T& x) { return T(1.)/x;  }
  static T sqr (const T& x) { return ::Sqr(x);  }
  static T sqrt(const T& x) { return ::Sqrt(x); }
  static T log (const T& x) { return ::Log(x);  }
  static T xlog(const T& x) { return T( ::Pred(mc::xlog(mc::mid(::Inf(x),::Sup(x),std::exp(-1.)))), ::Succ(std::max(mc::xlog(::Inf(x)),mc::xlog(::Sup(x)))) ); }
  static T fabs(const T& x) { return T(::Pred(mc::mid(::Inf(x),::Sup(x),0.)),::Succ(::Abs(x))); }
  static T exp (const T& x) { return ::Exp(x);  }
  static T sin (const T& x) { return ::Sin(x);  }
  static T cos (const T& x) { return ::Cos(x);  }
  static T tan (const T& x) { return ::Tan(x);  }
  static T asin(const T& x) { return ::ArcSin(x); }
  static T acos(const T& x) { return ::ArcCos(x); }
  static T atan(const T& x) { return ::ArcTan(x); }
  static T hull(const T& x, const T& y) { return ::Hull(x,y); }
  static T min (const T& x, const T& y) { return T( ::Pred(std::min(::Inf(x),::Inf(y))), ::Succ(std::min(::Sup(x),::Sup(y))) ); }
  static T max (const T& x, const T& y) { return T( ::Pred(std::max(::Inf(x),::Inf(y))), ::Succ(std::max(::Sup(x),::Sup(y))) ); }
  static T arh (const T& x, const double k) { return ::Exp(-x/k); }
  template <typename X> static T pow(const X& x, const int n) { return( (n>=3&&n%2)? ::Hull(::Power(Inf(x),n),::Power(Sup(x),n)): ::Power(x,n) ); }
  template <typename X, typename Y> static T pow(const X& x, const Y& y) { return ::Power(x,y); }
  static T monomial (const unsigned int n, const T* x, const int* k) { return n? ::Power(x[0], k[0]) * monomial( n-1, x+1, k+1 ): 1.; }
  static bool inter(T& xIy, const T& x, const T& y) { return ::Intersection(xIy,x,y); }
  static bool eq(const T& x, const T& y) { return x==y; }
  static bool ne(const T& x, const T& y) { return x!=y; }
  static bool lt(const T& x, const T& y) { return x<y;  }
  static bool le(const T& x, const T& y) { return x<=y; }
  static bool gt(const T& x, const T& y) { return y<x;  }
  static bool ge(const T& x, const T& y) { return y<=x; }
};

} // namespace mc

#endif
