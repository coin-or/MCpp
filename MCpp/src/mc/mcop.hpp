// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

#ifndef MC__MCOP_HPP
#define MC__MCOP_HPP

namespace mc
{
  //! @brief C++ structure to allow usage of various template parameters in the types mc::McCormick, mc::TModel, mc::TVar, and mc::SpecBnd of MC++ _ Specialization of this structure is required for the template parameters can be found in the header files defining the types mc::McCormick, mc::TModel, mc::TVar, and mc::SpecBnd.
  template <typename T>
  struct Op
  {
    static T point( const double c ) { return T(c); }
    static T zeroone() { return T(0,1); }
    static void I(T& x, const T& y) { x = y; }
    static double l(const T& x) { return x.l(); }
    static double u(const T& x) { return x.u(); }
    static double abs (const T& x) { return abs(x);  }
    static double mid (const T& x) { return mid(x);  }
    static double diam(const T& x) { return diam(x); }
    static T inv (const T& x) { return inv(x);  }
    static T sqr (const T& x) { return sqr(x);  }
    static T sqrt(const T& x) { return sqrt(x); }
    static T log (const T& x) { return log(x);  }
    static T xlog(const T& x) { return xlog(x); }
    static T fabs(const T& x) { return fabs(x); }
    static T exp (const T& x) { return exp(x);  }
    static T sin (const T& x) { return sin(x);  }
    static T cos (const T& x) { return cos(x);  }
    static T tan (const T& x) { return tan(x);  }
    static T asin(const T& x) { return asin(x); }
    static T acos(const T& x) { return acos(x); }
    static T atan(const T& x) { return atan(x); }
    static T erf (const T& x) { return erf(x);  }
    static T erfc(const T& x) { return erfc(x); }
    static T hull(const T& x, const T& y) { return hull(x,y); }
    static T min (const T& x, const T& y) { return min(x,y);  }
    static T max (const T& x, const T& y) { return max(x,y);  }
    static T arh (const T& x, const double k) { return arh(x,k); }
    template <typename X, typename Y> static T pow(const X& x, const Y& y) { return pow(x,y); }
    static T monomial (const unsigned int n, const T* x, const int* k) { return monomial(n,x,k); }
    static bool inter(T& xIy, const T& x, const T& y) { return inter(xIy,x,y); }
    static bool eq(const T& x, const T& y) { return x==y; }
    static bool ne(const T& x, const T& y) { return x!=y; }
    static bool lt(const T& x, const T& y) { return x<y;  }
    static bool le(const T& x, const T& y) { return x<=y; }
    static bool gt(const T& x, const T& y) { return x>y;  }
    static bool ge(const T& x, const T& y) { return x>=y; }
  };
}
#endif
