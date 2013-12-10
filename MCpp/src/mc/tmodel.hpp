// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

/*!
\page page_TAYLOR Taylor Model Arithmetic for Factorable Functions
\author Beno&icirc;t Chachuat

A \f$q\f$th-order Taylor model of a multivariate function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ that is (at least) \f$(q+1)\f$-times continuously differentiable on the domain \f$D\f$, consists of the \f$q^{\rm th}\f$-order multivariate Taylor polynomial \f$\mathcal P\f$, expanded around a point \f$\hat{x}\in D\f$, plus a remainder term \f$\mathcal R\f$, so that
\f{align*}
  f({x}) \in \mathcal P({x}-\hat{x}) \oplus \mathcal R, \quad \forall {x}\in D.
\f}
The polynomial part \f$\mathcal P\f$ is propagated symbolically and accounts for functional dependencies. The remainder term \f$\mathcal R\f$, on the other hand, is traditionally computed using interval analysis [Neumaier, 2002; Makino & Berz, 2003]; see figure below. More generally, convex/concave bounds for the remainder term can be propagated using McCormick relaxations, leading to so-called McCormick-Taylor models [Sahlodin & Chachuat, 2011], and ellipsoidal enclosure can be computed for the remainder term of vector-valued functions too [Houska <I>et al.</I>, 2013]. In particular, it can be established that the remainder term has convergence order (no less than) \f$q+1\f$ with respect to the diameter of the domain set \f$D\f$ under mild conditions [Bompadre <I>et al.</I>, 2012].

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html Taylor_model.png</TD>
</TR>
</TABLE></CENTER>

The classes mc::TModel and mc::TVar provide an implementation of Taylor model arithmetic. We note that mc::McCormick is <b>not a verified implementation</b> in the sense that rounding errors are not accounted for in propagating the coefficients in the multivariate polynomial part, which are treated as floating-point numbers.

The implementation of mc::TModel and mc::TVar relies on the operator/function overloading mechanism of C++. This makes the computation of Taylor models both simple and intuitive, similar to computing function values in real arithmetics or function bounds in interval arithmetic (see \ref page_INTERVAL). Moreover, mc::TVar can be used as the template parameter of other available types in MC++; for instance, mc::McCormick can be used in order to propagate the underlying interval bounds. Likewise, mc::McCormick can be used as the template parameter of the types fadbad::F, fadbad::B and fadbad::T of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing Taylor models of either the partial derivatives or the Taylor coefficients of a factorable function (see \ref sec_TAYLOR_fadbad).

mc::TModel and mc::TVar themselves are templated in the type used to propagate bounds on the remainder term. By default, mc::TModel and mc::TVar can be used with the non-verified interval type mc::Interval of MC++. For reliability, however, it is strongly recommended to use verified interval arithmetic such as <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> (header file <tt>mcprofil.hpp</tt>) or <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> (header file <tt>mcfilib.hpp</tt>). As already noted, convex/concave bounds on the remainder term can also be propagated by using the type mc::McCormick of MC++, thereby enabling McCormick-Taylor models.

As well as propagating Taylor models for factorable functions, mc::TModel and mc::TVar provide support for computing bounds on the Taylor model range (multivariate polynomial part). We note that computing exact bounds for multivariate polynomials is a hard problem in general. Instead, a number of computationally tractable, yet typically conservative, bounding approaches are implemented in mc::TModel and mc::TVar, which include:
- Bounding every monomial term independently and adding these bounds;
- Bounding the first- and diagonal second-order terms exactly and adding bounds for the second-order off-diagonal and higher-order terms computed independently [Lin & Stadtherr, 2007];
- Bounding the terms up to order 2 based on an eigenvalue decomposition of the corresponding Hessian matrix and adding bounds for the higher-order terms computed independently;
- Rewriting the multivariate polynomial in Bernstein form, thereby providing bounds as the minimum/maximum among all Bernstein coefficients [Lin & Rokne, 1995; 1996].
.

Examples of Taylor and McCormick-Taylor models (blue lines) constructed with mc::TModel and mc::TVar are shown on the left and right plots of the figure below, respectively, for the factorable function \f$f(x)=x \exp(-x^2)\f$ (red line) for \f$x\in [-0.5,1]\f$. Also shown on these plots are the bounds, either interval or convex/concave bounds, computed from the Taylor models.

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html TM-1D.png</TD>
<TD>\image html MCTM-1D.png</TD>
</TR>
</TABLE></CENTER>


\section sec_TAYLOR_I How do I compute a Taylor model with interval remainder bound of a factorable function?

Suppose we want to compute a 4th-order Taylor model for the real-valued function \f$f(x,y)=x\exp(x+y^2)-y^2\f$ with \f$(x,y)\in [1,2]\times[0,1]\f$. For simplicity, bounds on the remainder terms are computed using the default interval type mc::Interval here:

\code
      #include "interval.hpp"
      #include "tmodel.hpp"
      typedef mc::Interval I;
      typedef mc::TModel<I> TM;
      typedef mc::TVar<I> TV;
\endcode

First, the number of independent variables in the factorable function (\f$x\f$ and \f$y\f$ here) as well as the order of the Taylor model (4th order here) are specified by defining an mc::TModel object as:

\code
      TM mod( 2, 4 );
\endcode

Next, the variables \f$x\f$ and \f$y\f$ are defined as follows:

\code
      TV X( &mod, 0, I(1.,2.) );
      TV Y( &mod, 1, I(0.,1.) );
\endcode

Essentially, the first line means that <tt>X</tt> is a variable of class mc::TVar, participating in the Taylor model <tt>mod</tt>, belonging to the interval \f$[1,2]\f$, and having index 0 (indexing in C/C++ start at 0 by convention!). The same holds for the Taylor variable <tt>Y</tt>, participating in the model <tt>mod</tt>, belonging to the interval \f$[0,1]\f$, and having index 1.

Having defined the variables, a Taylor model of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ on \f$[1,2]\times[0,1]\f$ at the mid-point \f$(\frac{3}{2},\frac{1}{2})\f$ is simply computed as:

\code
      TV F = X*exp(X+pow(Y,2))-pow(Y,2);
\endcode

This model can be displayed to the standard output as:

\code
      std::cout << "f Taylor model: " << F << std::endl;
\endcode

which produces the following output:

\verbatim
f Taylor model: 
   a0    =  8.38199e+00     0  0
   a1    =  1.90755e+00     0  1
   a2    =  3.59621e+00     1  0
   a3    =  7.47482e-01     0  2
   a4    =  9.00782e-01     1  1
   a5    =  6.30186e-01     2  0
   a6    =  1.56945e-01     0  3
   a7    =  3.35238e-01     1  2
   a8    =  1.55141e-01     2  1
   a9    =  6.67468e-02     3  0
   a10   =  3.49519e-02     0  4
   a11   =  6.58449e-02     1  3
   a12   =  6.04330e-02     2  2
   a13   =  1.80397e-02     3  1
   a14   =  5.41191e-03     4  0
   R     =  [ -2.09182e+00 :  2.22652e+00 ]
   B     =  [ -1.02564e+01 :  3.93973e+01 ]
\endverbatim

<tt>a0</tt>,...,<tt>a14</tt> refer to the coefficients of the monomial terms in the Taylor model, with the corresponding variable orders given in the subsequent columns. The remainder term as well as the Taylor model range estimator are reported next.

Other operations involve retreiving the remainder bound, centering the remainder term in a Taylor model, or computing the value of its polynomial part at a given point:

\code
      I B = F.B();
      F.C();
      double x[2] = { 0.5, 1.5 };
      double Pval = F.P( x );
\endcode

See the documentations of mc::TModel and mc::TVar for a complete list of member functions. 


\section sec_TAYLOR_MC How do I compute a Taylor model with convex/concave remainder bounds of a factorable function?

Instead of using standard interval arithmetic for propagating bounds on the remainder term of the Taylor expansion, or for bounding the range of the multivariate polynomial part, one can also make use of the McCormick relaxation technique to propagate convex/concave bounds. Such McCormick-Taylor models are enabled in MC++ simply by selecting mc::McCormick as the template parameter in mc::TModel and mc::TVar:

\code
      #include "mccormick.hpp"
      typedef mc::McCormick<I> MC;
      typedef mc::TModel<MC> TMMC;
      typedef mc::TVar<MC> TVMC;
\endcode

Then, the procedure for computing the Taylor model remains essentially the same as described in the previous section. In order to compute a 4th-order McCormick-Taylor model of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ for \f$(x,y)\in[1,2]\times[0,1]\f$ and obtain convex/concave bounds for the resulting McCormick-Taylor model at \f$(1.5,0.5)\f$, we proceed as follows:

\code
      TMMC model( 2, 4 );
      TVMC X( &model, 0, MC( I(1.,2.), 1.5 ) );
      TVMC Y( &model, 1, MC( I(0.,1.), 0.5 ) );
      TVMC F = X*exp(X+pow(Y,2))-pow(Y,2);
      MC B = F.B();
\endcode

The only difference here concerns the initialization of the McCormick variables inside the Taylor variables <tt>X</tt> and <tt>Y</tt>, which passes the bounds for the variable as well as the point at which the convex/concave bounds are computed&mdash;See \ref sec_MCCORMICK_use

More information on the Taylor model and the corresponding bounds can again be obtained as:

\code
      std::cout << "f McCormick-Taylor model: " << F << std::endl;
\endcode

In this case, the following information is displayed:

\verbatim
f McCormick-Taylor model: 
   a0    =  8.38199e+00     0  0
   a1    =  1.90755e+00     0  1
   a2    =  3.59621e+00     1  0
   a3    =  7.47482e-01     0  2
   a4    =  9.00782e-01     1  1
   a5    =  6.30186e-01     2  0
   a6    =  1.56945e-01     0  3
   a7    =  3.35238e-01     1  2
   a8    =  1.55141e-01     2  1
   a9    =  6.67468e-02     3  0
   a10   =  3.49519e-02     0  4
   a11   =  6.58449e-02     1  3
   a12   =  6.04330e-02     2  2
   a13   =  1.80397e-02     3  1
   a14   =  5.41191e-03     4  0
   R     =  [ -2.09182e+00 :  2.22652e+00 ] [ -1.70505e+00 :  1.80226e+00 ]
   B     =  [ -1.02564e+01 :  3.93973e+01 ] [ -2.63875e+00 :  2.66234e+01 ]
\endverbatim

As expected, the McCormick-Taylor model coefficients are identical to their Taylor model counterparts, and so are the interval bounds of the remainder term and Taylor model range estimator. The convex/concave bounds at \f$(1.5,0.5)\f$ come as additional information and are seen to be tighter than the interval bounds. Naturally, subgradients of the resulting convex and concave relaxations can also be propagated as explained at \ref sec_MCCORMICK_sub.

Finally, the McCormick relaxation of the Taylor model can also be obtained as:

\code
      MC B = F.B();
\endcode

Naturally, subgradients of the convex/concave remainder bounds can also be propagated, e.g., by seeding the subgradients of <tt>X</tt> and <tt>Y</tt> as follows:

\code
      TMMC model( 2, 4 );
      TVMC X( &model, 0, MC( I(1.,2.), 1.5 ).sub(2,0) );
      TVMC Y( &model, 1, MC( I(0.,1.), 0.5 ).sub(2,1) );
      TVMC F = X*exp(X+pow(Y,2))-pow(Y,2);
      MC B = F.B();
\endcode

See \ref sec_MCCORMICK_sub for more details.


\section sec_TAYLOR_fadbad How do I compute Taylor models of the partial derivatives or the Taylor coefficients of a factorable function using FADBAD++?

The combination of mc::TVar with the classes fadbad::F, and fadbad::B of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> to compute Taylor models of the partial derivatives or the Taylor coefficients of a factorable function is essentially the same as with mc::McCormick (see \ref sec_MCCORMICK_fadbad).

We present the case of fadbad::F only.  Continuing the previous example, McCormick-Taylor models of the partial derivatives of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ for \f$(x,y)\in [1,2]\times[0,1]\f$ can be computed as follows:

\code
      #include "mcfadbad.hpp" // available in MC++
      #include "fadiff.h"     // available in FADBAD++
      typedef fadbad::F<TVMC> FTVMC;
\endcode

\code
      FTVMC FX = X;           // initialize FX with McCormick-Taylor variable X
      FX.diff(0,2);           // differentiate with respect to x (index 0 of 2)

      FTVMC FY = Y;           // initialize FY with McCormick-Taylor variable Y
      FY.diff(1,2);           // differentiate with respect to y (index 1 of 2)

      FTVMC FF = FX*exp(FX+pow(FY,2))-pow(FY,2);
      std::cout << "f McCormick-Taylor model with subgradients at (1.5,0.5): " << FF.x() << std::endl;
      std::cout << "df/dx McCormick-Taylor model with subgradients at (1.5,0.5): " << FF.d(0) << std::endl;
      std::cout << "df/dy McCormick-Taylor model with subgradients at (1.5,0.5): " << FF.d(1) << std::endl;
\endcode

producing the output:

\verbatim
f McCormick-Taylor model with subgradients at (1.5,0.5): 
   a0    =  8.38199e+00     0  0
   a1    =  1.90755e+00     0  1
   a2    =  3.59621e+00     1  0
   a3    =  7.47482e-01     0  2
   a4    =  9.00782e-01     1  1
   a5    =  6.30186e-01     2  0
   a6    =  1.56945e-01     0  3
   a7    =  3.35238e-01     1  2
   a8    =  1.55141e-01     2  1
   a9    =  6.67468e-02     3  0
   a10   =  3.49519e-02     0  4
   a11   =  6.58449e-02     1  3
   a12   =  6.04330e-02     2  2
   a13   =  1.80397e-02     3  1
   a14   =  5.41191e-03     4  0
   R     =  [ -2.09182e+00 :  2.22652e+00 ] [ -1.70505e+00 :  1.80226e+00 ] [ (-7.35405e-01, 1.52576e-02) : ( 8.33268e-01, 1.52576e-02) ]
   B     =  [ -1.02564e+01 :  3.93973e+01 ] [ -2.63875e+00 :  2.66234e+01 ] [ ( 7.24415e+00, 2.32253e+00) : ( 8.81282e+00, 1.67350e+01) ]

df/dx McCormick-Taylor model with subgradients at (1.5,0.5): 
   a0    =  1.43867e+01     0  0
   a1    =  3.59591e+00     0  1
   a2    =  5.03458e+00     1  0
   a3    =  1.34997e+00     0  2
   a4    =  1.26158e+00     1  1
   a5    =  8.10583e-01     2  0
   a6    =  2.61575e-01     0  3
   a7    =  4.68731e-01     1  2
   a8    =  1.98437e-01     2  1
   a9    =  8.11786e-02     3  0
   a10   =  5.82532e-02     0  4
   a11   =  9.23031e-02     1  3
   a12   =  7.84726e-02     2  2
   a13   =  2.28503e-02     3  1
   a14   =  6.61455e-03     4  0
   R     =  [ -1.34116e+00 :  1.40851e+00 ] [ -1.32972e+00 :  1.39325e+00 ] [ ( 1.52576e-02, 1.52576e-02) : ( 1.52576e-02, 1.52576e-02) ]
   B     =  [ -1.11441e+01 :  5.89599e+01 ] [  1.45326e-01 :  3.96272e+01 ] [ ( 1.10351e+01, 7.44520e+00) : ( 1.10351e+01, 2.76304e+01) ]

df/dy McCormick-Taylor model with subgradients at (1.5,0.5): 
   a0    =  7.63199e+00     0  0
   a1    =  5.97354e+00     0  1
   a2    =  3.59621e+00     1  0
   a3    =  1.88876e+00     0  2
   a4    =  2.69889e+00     1  1
   a5    =  6.30186e-01     2  0
   a6    =  5.61936e-01     0  3
   a7    =  7.85628e-01     1  2
   a8    =  4.70235e-01     2  1
   a9    =  6.67468e-02     3  0
   a10   =  1.13425e-01     0  4
   a11   =  2.33464e-01     1  3
   a12   =  1.38004e-01     2  2
   a13   =  5.14131e-02     3  1
   a14   =  5.41191e-03     4  0
   R     =  [ -7.13855e+00 :  7.40794e+00 ] [ -5.04673e+00 :  5.18142e+00 ] [ (-2.68232e+00,-1.50133e+00) : ( 2.81702e+00, 1.63602e+00) ]
   B     =  [ -3.93313e+01 :  7.87946e+01 ] [ -2.40726e+01 :  5.36565e+01 ] [ (-9.08762e+00, 7.54498e+00) : (-3.58827e+00, 5.38645e+01) ]
\endverbatim


\section sec_TAYLOR_fct Which functions are overloaded for Taylor model arithmetic?

mc::TVar overloads the usual functions <tt>exp</tt>, <tt>log</tt>, <tt>sqr</tt>, <tt>sqrt</tt>, <tt>pow</tt>, <tt>inv</tt>, <tt>cos</tt>, <tt>sin</tt>, <tt>tan</tt>, <tt>acos</tt>, <tt>asin</tt>, <tt>atan</tt>. Unlike mc::Interval and mc::McCormick, the functions <tt>min</tt>, <tt>max</tt> and <tt>fabs</tt> are not overloaded in mc::TVar as they are nonsmooth. Moreover, mc::TVar defines the following functions:
- <tt>inter(x,y,z)</tt>, computing a Taylor model of the intersection \f$x = y\cap z\f$ of two Taylor models and returning true/false if the intersection is nonempty/empty. With Taylor models \f$\mathcal P_y\oplus\mathcal R_y\f$ and \f$\mathcal P_z\oplus\mathcal R_z\f$, this intersection is computed as follows:
\f{align*}
  \mathcal P_{x} =\ & (1-\eta) \mathcal P_y^{\rm C} + \eta \mathcal P_z^{\rm C}\\
  \mathcal R_{x} =\ & [\mathcal R_y^{\rm C}\oplus\eta\mathcal{B}(\mathcal P_y^{\rm C}-\mathcal P_z^{\rm C})] \cap [\mathcal R_z^{\rm C}\oplus (1-\eta)\mathcal{B}(\mathcal P_z^{\rm C}-\mathcal P_y^{\rm C})]\,.
\f}
with \f$\mathcal{B}(\cdot)\f$ the Taylor model range bounder, and \f$\eta\f$ a real scalar in \f$[0,1]\f$. Choosing \f$\eta=1\f$ amounts to setting the polynomial part \f$\mathcal P_{x}\f$ as \f$\mathcal P_y\f$, whereas \f$\eta=0\f$ sets \f$\mathcal P_{x}\f$ as \f$\mathcal P_z\f$. The parameter \f$\eta\f$ can be defined in mc::TModel::Options::REF_POLY.
- <tt>hull(x,y)</tt>, computing a Taylor model of the union \f$x = y\cup z\f$ of two Taylor models. With Taylor models \f$\mathcal P_y\oplus\mathcal R_y\f$ and \f$\mathcal P_z\oplus\mathcal R_z\f$, this union is computed as follows:
\f{align*}
  \mathcal P_{x} =\ & (1-\eta) \mathcal P_y^{\rm C} + \eta \mathcal P_z^{\rm C}\\
  \mathcal R_{x} =\ & {\rm hull}\{\mathcal R_y^{\rm C}\oplus\eta\mathcal{B}(\mathcal P_y^{\rm C}-\mathcal P_z^{\rm C}), \mathcal R_z^{\rm C}\oplus (1-\eta)\mathcal{B}(\mathcal P_z^{\rm C}-\mathcal P_y^{\rm C})\}\,.
\f}
with \f$\mathcal{B}(\cdot)\f$ and \f$\eta\f$ as previously.


\section sec_TAYLOR_opt How are the options set for the computation of a Taylor model?

The class mc::TModel has a public member called mc::TModel::options that can be used to set/modify the options; e.g.,

\code
      model.options.BOUNDER_TYPE = TM::Options::EIGEN;
      model.options.SCALE_VARIABLES = true;
\endcode

The available options are the following:

<TABLE border="1">
<CAPTION><EM>Options in mc::TModel::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>BOUNDER_TYPE</tt> <TD><tt>mc::TModel::Options::BOUNDER</tt> <TD>mc::TModel::Options::LSB
         <TD>Taylor model range bounder.
     <TR><TH><tt>BOUNDER_ORDER</tt> <TD><tt>unsigned int</tt> <TD>0
         <TD>Order of Bernstein polynomial for Taylor model range bounding, when mc::TModel::options::BOUNDER_TYPE = mc::TModel::options::BERNSTEIN is selected. Only values greater than the actual Taylor model order are accounted for; see [Lin & Rokne, 1996].
     <TR><TH><tt>PROPAGATE_BNDT</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to propagate bounds in arithmetic of the template parameter along with Taylor model arithmetic.
     <TR><TH><tt>INTER_WITH_BNDT</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to intersect bounds derived from the Taylor model with those from the template parameter arithmetic. Only if mc::TModel::options::PROPAGATE_BNDT is set to true.
     <TR><TH><tt>SCALE_VARIABLES</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to scale the variable ranges to [-1,1] internally.
     <TR><TH><tt>CENTER_REMAINDER</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to center the remainder term during Taylor model propagation.
     <TR><TH><tt>REF_MIDPOINT</tt> <TD><tt>bool</tt> <TD>true
         <TD>Whether to take the midpoint of the inner range as the reference in the outer composition with a univariate function (true), as opposed to taking the constant coefficient of the inner Taylor model (false).
     <TR><TH><tt>REF_POLY</tt> <TD><tt>double</tt> <TD>0.
         <TD>Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull (see \ref sec_TAYLOR_fct). A value of 0. amounts to selecting the polynomial part of the left operand, whereas a value of 1. selects the right operand.
     <TR><TH><tt>BERNSTEIN_USE</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to compute a Berstein model [Stancu, 1963] of the outer function in a univariate composition and use it instead of the Taylor model when its remainder term is smaller. This enables tighter Taylor models on wider variable ranges, while retaining the high-order convergence order of Taylor models, but incurs a computational overhead.
     <TR><TH><tt>BERNSTEIN_OPT</tt> <TD><tt>bool</tt> <TD>true
         <TD>Whether to compute exact remainder bounds for Berstein models of convex/concave univariates exp, log, inv and sqrt. This provides tighter estimators, but incurs an extra computational overhead. The exact remainder bounds are computed using the golden section search method.
     <TR><TH><tt>BERNSTEIN_TOL</tt> <TD><tt>double</tt> <TD>1e-10
         <TD>Termination tolerance for determination of the exact remainder bounds in a Berstein model of convex/concave univariates exp, log, inv and sqrt.
     <TR><TH><tt>BERNSTEIN_MAXIT</tt> <TD><tt>int</tt> <TD>100
         <TD>Maximum number of iterations for determination of the exact remainder bounds in a Berstein model of convex/concave univariates exp, log, inv and sqrt.
     <TR><TH><tt>DISPLAY_DIGITS</tt> <TD><tt>unsigned int</tt> <TD>5
         <TD>Number of digits in output stream for Taylor model coefficients.
</TABLE>


\section sec_TM_err Errors What errors can I encounter during computation of a Taylor model?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::TModel::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the computation of a Taylor model, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Possible errors encountered during the computation of a Taylor model are:

<TABLE border="1">
<CAPTION><EM>Errors during the Computation of a Taylor Model</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Division by zero
     <TR><TH><tt>2</tt> <TD>Failed to compute eigenvalue decomposition in range bounder TModel::Options::EIGEN
     <TR><TH><tt>3</tt> <TD>Failed to compute the maximum gap between a univariate term and its Bernstein model
     <TR><TH><tt>-1</tt> <TD>Number of variable in Taylor model must be nonzero
     <TR><TH><tt>-2</tt> <TD>Failed to construct Taylor variable
     <TR><TH><tt>-3</tt> <TD>Taylor model bound does not intersect with bound in template parameter arithmetic
     <TR><TH><tt>-4</tt> <TD>Operation between Taylor variables linked to different Taylor models
     <TR><TH><tt>-5</tt> <TD>Maximum size of Taylor model reached (monomials indexed as unsigned int)
     <TR><TH><tt>-33</tt> <TD>Feature not yet implemented in mc::Specbnd
</TABLE>

Moreover, exceptions may be thrown by the template parameter class itself.


\section sec_TM_refs References

- Berz, M., and G. Hoffstaetter, <A href="http://dx.doi.org/10.1023/A:1009958918582">Computation and Application of Taylor Polynomials with Interval Remainder Bounds</A>, <i>Reliable Computing</i>, <b>4</b>:83-97, 1998
- Bompadre, A., A. Mitsos, and B. Chachuat, <A href="http://dx.doi.org/10.1007/s10898-012-9998-9">Convergence analysis of Taylor models and McCormick-Taylor models</A>, <i>Journal of Global Optimization</i>, <b>in press</b>, October 2012
- Houska, B., M.E. Villanueva, B. Chachuat, <A href="http://cdc2013.units.it/index.php">A validated integration algorithm for nonlinear ODEs using Taylor models and ellipsoidal calculus</A>, <I>52nd IEEE Conference on Decision and Control (CDC)</I>, December 10-13, 2013, Florence, Italy
- Lin, Q., and J.G, Rokne, <A href="http://dx.doi.org/10.1016/0377-0427(93)E0270-V">Methods for bounding the range of a polynomial</A>, <i>Journal of Computational & Applied Mathematics</i>, <b>58</b>:193-199, 1995
- Lin, Q., and J.G, Rokne, <A href="http://dx.doi.org/10.1016/0898-1221(96)00020-X">Interval approximation of higher order to the ranges of functions</A>, <i>Computers & Mathematics with Applications</i>, <b>31</b>(7):101-109, 1996
- Lin, Y., and M.A. Stadtherr, <A href="http://dx.doi.org/10.1016/j.apnum.2006.10.006">Validated solutions of initial value problems for parametric ODEs</A>, <i>Applied Numerical Mathematics</i>, <b>57</b>(10):1145-1162, 2007
- Makino, K., and M. Berz, <A href="http://www.ijpam.eu/contents/2003-6-3/1/">Taylor models and other validated functional inclusion methods</A>, <i>International Journal of Pure & Applied Mathematics</i>, <b>6</b>(3):239-312, 2003
- Neumaier, A., <A href="http://dx.doi.org/10.1023/A:1023061927787">Taylor forms--Use and limits</A>, <i>Reliable Computing</i>, <b>9</b>(1):43-79, 2002
- Sahlodin, M.A., and B. Chachuat, <A href="http://dx.doi.org/10.1016/j.compchemeng.2011.01.031">Convex/concave relaxations of parametric ODEs using Taylor models</A>, <i>Computers & Chemical Engineering</i>, <b>35</b>(5):844-857, 2011
- Stancu, D.D., <A href="http://www.jstor.org/stable/2003844">Evaluation of the remainder term in approximation formulas by Bernstein polynomials</A>, <i>Mathematics of Computation</i>, <b>17</b>(83):270-278, 1963
.
*/

#ifndef MC__TMODEL_H
#define MC__TMODEL_H

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <string>
#include <stdarg.h>
#include <cassert>
#include <climits>
#include <stdlib.h>
#include <complex>

#include "mcfunc.hpp"
#include "mcop.hpp"
#include "mclapack.hpp"

#undef  MC__TMODEL_DEBUG
#undef  MC__TMODEL_DEBUG_SCALE
#define MC__TMODEL_CHECK
#undef  MC__TVAR_DEBUG_EXP
#undef  MC__TVAR_DEBUG_BERSTEIN
#undef  MC__TVAR_HYBRID_EIGEN

namespace mc
{

template <typename T> class TVar;
typedef unsigned long long TM_size;

//! @brief C++ class for Taylor model computation of factorable function - Taylor model environment
////////////////////////////////////////////////////////////////////////
//! mc::TModel is a C++ class for definition of Taylor model
//! environment. Propagation of Taylor models for factorable functions
//! is via the C++ class mc::TVar. The template parameter corresponds
//! to the type used to propagate the remainder bound.
////////////////////////////////////////////////////////////////////////
template <typename T>
class TModel
////////////////////////////////////////////////////////////////////////
{
  friend class TVar<T>;
  template <typename U> friend class TModel;

  template <typename U> friend TVar<U> inv
    ( const TVar<U>& );
  template <typename U> friend TVar<U> sqrt
    ( const TVar<U>& );
  template <typename U> friend TVar<U> log
    ( const TVar<U>& );
  template <typename U> friend TVar<U> exp
    ( const TVar<U>& );
  template <typename U> friend TVar<U> pow
    ( const TVar<U>&, const int );
    
public:

  /** @addtogroup TAYLOR Taylor Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Constructor of Taylor model environment for <tt>nvar</tt> variables and order <tt>nord</tt>
  TModel
    ( const unsigned int nvar, const unsigned int nord )
    { _size( nvar, nord ); }
  //! @brief Destructor of Taylor model environment

  ~TModel()
    { _cleanup(); }

  //! @brief Number of variables in Taylor model environment
  unsigned int nvar() const
    { return _nvar; };

  //! @brief Order of Taylor model environment
  unsigned int nord() const
    { return _nord; };

  //! @brief Total number of monomial terms in Taylor variable
  unsigned int nmon() const
    { return _nmon; };

  //! @brief Const pointer to array of size <tt>nmon()</tt> with bounds on each monomial term
  const T* bndmon() const
    { return _bndmon; };

  //! @brief Const pointer to array of size <tt>nmon()*nvar()</tt> with variable exponents on each monomial term. The exponent for variable <tt>ivar</tt> in monomial term <tt>imon</tt> is at position <tt>imon*nvar()+ivar</tt>.
  const unsigned int* expmon() const
    { return _expmon; };

  //! @brief Index of monomial term whose variable exponents are the same as those in array <tt>iexp</tt> (of size <tt>nvar()</tt>)
  unsigned int loc_expmon
    ( const unsigned int *iexp ) const
    { return _loc_expmon( iexp ); }

  //! @brief Exceptions of mc::TModel
  class Exceptions
  {
  public:
    //! @brief Enumeration type for TModel exception handling
    enum TYPE{
      DIV=1,	//!< Division by zero scalar
      EIGEN,	//!< Failed to compute eigenvalue decomposition in range bounder TModel::Options::EIGEN
      BERNSTEIN,//!< Failed to compute the maximum gap between a univariate term and its Bernstein model
      SIZE=-1,	//!< Number of variable in Taylor model must be nonzero
      INIT=-2,	//!< Failed to construct Taylor variable
      INCON=-3, //!< Taylor model bound does not intersect with bound in template parameter arithmetic
      TMODEL=-4,//!< Operation between Taylor variables linked to different Taylor models
      MAXSIZE=-5,//!< Maximum size of Taylor model reached (monomials indexed as unsigned int)
      UNDEF=-33 //!< Feature not yet implemented in mc::Specbnd
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DIV:
        return "mc::TModel\t Division by zero scalar";
      case EIGEN:
        return "mc::TModel\t Range bounder with eigenvalue decomposition failed";
      case BERNSTEIN:
        return "mc::TModel\t Bernstein remainder evalution failed";
      case SIZE:
        return "mc::TModel\t Inconsistent Taylor model dimension";
      case INIT:
        return "mc::TModel\t Taylor variable initialization failed";
      case INCON:
        return "mc::TModel\t Inconsistent bounds with template parameter arithmetic";
      case TMODEL:
        return "mc::TModel\t Operation between Taylor variables in different Taylor model environment not allowed";
      case MAXSIZE:
        return "mc::TModel\t Maximum size in Taylor model reached";
      case UNDEF:
        return "mc::TModel\t Feature not yet implemented in mc::Specbnd class";
      default:
        return "mc::TModel\t Undocumented error";
      }
    }

  private:
    TYPE _ierr;
  };

  //! @brief Options of mc::TModel
  struct Options
  {
    //! @brief Constructor of mc::TModel::Options
    Options():
      BOUNDER_TYPE(LSB), BOUNDER_ORDER(0), PROPAGATE_BNDT(false),
      INTER_WITH_BNDT(false), SCALE_VARIABLES(false), CENTER_REMAINDER(false),
      REF_MIDPOINT(true), REF_POLY(0.), BERNSTEIN_USE(false),
      BERNSTEIN_OPT(true), BERNSTEIN_MAXIT(100), BERNSTEIN_TOL(1e-10),
      DISPLAY_DIGITS(5)
      {}
    //! @brief Copy constructor of mc::TModel::Options
    template <typename U> Options
      ( U&options )
      : BOUNDER_TYPE( options.BOUNDER_TYPE ),
        BOUNDER_ORDER( options.BOUNDER_ORDER ),
        PROPAGATE_BNDT( options.PROPAGATE_BNDT ),
        INTER_WITH_BNDT( options.INTER_WITH_BNDT ),
        SCALE_VARIABLES( options.SCALE_VARIABLES ),
        CENTER_REMAINDER( options.CENTER_REMAINDER ),
        REF_MIDPOINT( options.REF_MIDPOINT ),
        REF_POLY(options.REF_POLY),
	BERNSTEIN_USE(options.BERNSTEIN_USE),
	BERNSTEIN_OPT(options.BERNSTEIN_OPT),
        BERNSTEIN_MAXIT(options.BERNSTEIN_MAXIT),
	BERNSTEIN_TOL(options.BERNSTEIN_TOL),
	DISPLAY_DIGITS(options.DISPLAY_DIGITS)
      {}
    //! @brief Assignment of mc::TModel::Options
    template <typename U> Options& operator =
      ( U&options ){
        BOUNDER_TYPE     = options.BOUNDER_TYPE;
        BOUNDER_ORDER    = options.BOUNDER_ORDER;
        PROPAGATE_BNDT   = options.PROPAGATE_BNDT;
        INTER_WITH_BNDT  = options.INTER_WITH_BNDT;
        SCALE_VARIABLES  = options.SCALE_VARIABLES;
        CENTER_REMAINDER = options.CENTER_REMAINDER;
        REF_MIDPOINT     = options.REF_MIDPOINT,
        REF_POLY         = options.REF_POLY,
	BERNSTEIN_USE    = options.BERNSTEIN_USE,
	BERNSTEIN_OPT    = options.BERNSTEIN_OPT,
        BERNSTEIN_MAXIT  = options.BERNSTEIN_MAXIT,
	BERNSTEIN_TOL    = options.BERNSTEIN_TOL,
	DISPLAY_DIGITS   = options.DISPLAY_DIGITS;
        return *this;
      }
    //! @brief Taylor model range bounder option
    enum BOUNDER{
      NAIVE=0,	//!< Naive polynomial range bounder
      LSB,	//!< Lin & Stadtherr range bounder
      EIGEN,	//!< Eigenvalue decomposition-based bounder
      BERNSTEIN,//!< Bernstein range bounder
      HYBRID	//!< Hybrid LSB + EIGEN range bounder
    };
    //! @brief Taylor model range bounder - See \ref sec_TAYLOR_opt
    int BOUNDER_TYPE;
    //! @brief Order of Bernstein polynomial for Taylor model range bounding (no less than Taylor model order!). Only if mc::TModel::options::BOUNDER_TYPE is set to mc::TModel::options::BERNSTEIN.
    unsigned int BOUNDER_ORDER;
    //! @brief Array of Taylor model range bounder names (for display)
    static const std::string BOUNDER_NAME[5];
    //! @brief Whether to propagate bounds in arithmetic of the template parameter along with Taylor model arithmetic
    bool PROPAGATE_BNDT;
    //! @brief Whether to intersect bounds derived from the Taylor model with those from the template parameter arithmetic. Only if mc::TModel::options::PROPAGATE_BNDT is set to <tt>true</tt>.
    bool INTER_WITH_BNDT;
    //! @brief Whether to scale the variable ranges to [-1,1] internally
    bool SCALE_VARIABLES;
    //! @brief Whether to center the remainder term during Taylor model propagation
    bool CENTER_REMAINDER;
    //! @brief Whether to take the midpoint of the inner range as the reference in the outer composition with a univariate function (true), as opposed to taking the constant coefficient of the inner Taylor model (false).
    bool REF_MIDPOINT;
    //! @brief Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull (see \ref sec_TAYLOR_fct). A value of 0. amounts to selecting the polynomial part of the left operand, whereas a value of 1. selects the right operand.
    double REF_POLY;
    //! @brief Order of Bernstein polynomial for range bounding Maximum number of iterations for determination of the exact remainder bounds in a Berstein model of convex/concave univariates exp, log, inv and sqrt.
    unsigned int BERNSTEIN_ORDER;
    //! @brief Whether to compute a Berstein model [Stancu, 1963] of the outer function in a univariate composition and use it instead of the Taylor model when its remainder term is smaller. This enables tighter Taylor models on wider variable ranges, while retaining the high-order convergence order of Taylor models, but incurs a computational overhead.
    bool BERNSTEIN_USE;
    //! @brief Whether to compute exact remainder bounds for Berstein models of convex/concave univariates exp, log, inv and sqrt. This provides tighter estimators, but incurs an extra computational overhead. The exact remainder bounds are computed using the golden section search method.
    bool BERNSTEIN_OPT;
    //! @brief Maximum number of iterations for determination of the exact remainder bounds in a Berstein model of convex/concave univariates exp, log, inv and sqrt.
    unsigned int BERNSTEIN_MAXIT;
    //! @brief Termination tolerance for determination of the exact remainder bounds in a Berstein model of convex/concave univariates exp, log, inv and sqrt.
    double BERNSTEIN_TOL;
    //! @brief Number of digits in output stream for Taylor model coefficients.
    unsigned int DISPLAY_DIGITS;
  } options;
  /** @} */

private:  
  //! @brief Order of Taylor model
  unsigned int _nord;
  //! @brief Number of variables in Taylor model
  unsigned int _nvar;
  //! @brief Total number of monomial terms in Taylor model
  unsigned int _nmon;
  //! @brief Array of size <tt>_nord</tt> with indices of first monomial term of order <tt>iord=1,...,_nord</tt> in Taylor model
  unsigned int *_posord;
  //! @brief Order used to size _posord
  unsigned _posord_size;
  //! @brief Array of size <tt>_nmon*_nvar</tt> with variable exponents in monomial terms. The exponent for variable <tt>ivar</tt> in monomial term <tt>imon</tt> is at location <tt>imon*nvar()+ivar</tt>.
  unsigned int *_expmon;
  //! @brief Number of monomial coefficients used to size _expmon
  unsigned _expmon_size;
  //! @brief Double array of size <tt>(_nmon+1,<=_nmon)</tt> with indices of monomial terms from product of two monomial terms <tt>imon=1,...,_nmon</tt> and <tt>jmon=1,...,_nmon</tt>.
  unsigned int **_prodmon;
  //! @brief Array of size <tt>_nmon</tt> with bounds on monomial terms <tt>imon=1,...,_nmon</tt>
  T *_bndmon;
  //! @brief Have any of the model variables been modified?
  bool _modvar;
  //! @brief Array of <tt>(_nvar+_nord-1)*(_nord+1)</tt> contining binomial coefficients
  TM_size *_binom;
  //! @brief Maximum binomial coefficients in array _binom
  std::pair<unsigned int, unsigned int> _binom_size;
  //! @brief Double array of size <tt>(_nvar,_nord+1)</tt> with bounds on powers of (possibly scaled) variable ranges around reference point
  T **_bndpow;
  //! @brief Array of size <tt>_nvar</tt> with reference points for the variables
  double *_refpoint;
  //! @brief Array of size <tt>_nvar</tt> with scaling for the variables
  double *_scaling; 
  //! @brief Berntein coefficients for univariate terms
  double *_cbern; 

  //! @brief Internal Taylor variable to speed-up computations and reduce dynamic allocation
  TVar<T>* _TV;

  //! @brief Set Taylor model order <tt>nord</tt> and number of variables <tt>nvar</tt>
  void _size
    ( const unsigned int nvar, const unsigned int nord );

  //! @brief Populate array <tt>_bndpow</tt> for variable <tt>ivar</tt> with range <tt>X</tt>, reference <tt>Xref</tt> and scaling <tt>scaling</tt>
  void _set_bndpow
    ( const unsigned int ivar, const T&X, const double Xref,
      const double scaling );

  //! @brief Populate array <tt>_bndmon</tt>
  void _set_bndmon();

  //! @brief Populate array <tt>_posord</tt> up to order <tt>nord</tt>
  void _set_posord
    ( const unsigned int nord );

  //! @brief Extend array <tt>_posord</tt> up to maximum order <tt>maxord</tt>
  void _ext_posord
    ( const unsigned int maxord );

  //! @brief Populate array <tt>_expmon</tt> up to order <tt>nord</tt>
  void _set_expmon
    ( const unsigned int nord );

  //! @brief Extend array <tt>_expmon</tt> up to order <tt>maxord</tt> to accomodate <tt>maxmon</tt> coefficients
  void _ext_expmon
    ( const unsigned int maxord, const bool full=false );
  
  //! @brief Generate variable exponents <tt>iexp</tt> for subsequent monomial order <tt>iord</tt>
  void _next_expmon
    ( unsigned int *iexp, const unsigned int iord ) const;

  //! @brief Populate array _prodmon w/ exponents resulting from the product of two monomial terms 1,...,nmon
  void _set_prodmon();
    
  //! @brief Get index of monomial term with variable exponents <tt>iexp</tt> in <tt>1,...,_nmon</tt>
  unsigned int _loc_expmon
    ( const unsigned int *iexp ) const;
    
  //! @brief Populate array <tt>_binom</tt> with binomial coefficients up to order <tt>nord</tt>
  void _set_binom
    ( const unsigned int nord );
    
  //! @brief Extend array <tt>_binom</tt> with binomial coefficients up to order <tt>nord</tt>
  void _ext_binom
    ( const unsigned int nord );

  //! @brief Get binomial coefficient \f$\left(\stackrel{n}{k}\right)\f$
  TM_size _get_binom
    ( const unsigned int n, const unsigned int k ) const;

  //! @brief Scale a given coefficient array in range [0,1] with reference point 0 for variable <a>ivar</a>
  void _scale
    ( const unsigned int ivar, double*coef ) const;
      
  //! @brief Reset the variable bound arrays
  void _reset();

  //! @brief Clean up the arrays
  void _cleanup();
    
  //! @brief Prototype real-valued function for univariate terms
  typedef double (puniv)
    ( const double x, const double*rusr, const int*iusr );
    
  //! @brief Prototype real-valued function for univariate terms
  typedef double (punivopt)
    ( const double x, const double*rusr, const int*iusr, puniv df,
      const T&I,  const std::pair<unsigned int,const double*>&bern );

  //! @brief Prototype interval-valued function for univariate terms
  typedef T (punivext)
    ( const T&x, const double*rusr, const int*iusr );

  //! @brief Compute Berstein-Taylor model of a univariate term
  TVar<T> _univ_bernstein
    ( const TVar<T>&TV, puniv f, puniv df, punivext If,
      punivext d2If, const double*rusr, const int*iusr );

  //! @brief Compute gap between a univariate term and its Berstein polynomial at point <a>x</a>
  static double _gap_bernstein
    ( const double x, const double*rusr, const int*iusr, puniv df,
      const T&I, const std::pair<unsigned int,const double*>&bern );

  //! @brief Compute gap derivative between a univariate term and its Berstein polynomial at point <a>x</a>
  static double _dgap_bernstein
    ( const double x, const double*rusr, const int*iusr, puniv df,
      const T&I, const std::pair<unsigned int,const double*>&bern );

  //! @brief Golden section search method for root finding 
  double _goldsect
    ( const double xL, const double xU, punivopt fopt, const double*rusr,
      const int*iusr, puniv df,  const T&Ix,
      const std::pair<unsigned int,const double*>&bern );

  //! @brief Golden section search iterations 
  double _goldsect_iter
    ( const bool init, const double a, const double fa, const double b,
      const double fb, const double c, const double fc, punivopt fopt,
      const double*rusr, const int*iusr, puniv df, const T&Ix,
      const std::pair<unsigned int,const double*>&bern );

  //! @brief Recursive calculation of nonnegative integer powers
  TVar<T> _intpow
    ( const TVar<T>&TV, const int n );

  //! @brief Taylor model of inverse univariate
  TVar<T> _inv_taylor
    ( const TVar<T>&TV );

  //! @brief Bernstein model of inverse univariate
  TVar<T> _inv_bernstein
    ( const TVar<T>&TV );

  //! @brief Taylor model of square-root univariate
  TVar<T> _sqrt_taylor
    ( const TVar<T>&TV );

  //! @brief Bernstein model of square-root univariate
  TVar<T> _sqrt_bernstein
    ( const TVar<T>&TV );

  //! @brief Taylor model of exponential univariate
  TVar<T> _exp_taylor
    ( const TVar<T>&TV );

  //! @brief Bernstein model of exponential univariate
  TVar<T> _exp_bernstein
    ( const TVar<T>&TV );

  //! @brief Taylor model of log univariate
  TVar<T> _log_taylor
    ( const TVar<T>&TV );

  //! @brief Bernstein model of log univariate
  TVar<T> _log_bernstein
    ( const TVar<T>&TV );

} ;

template <typename T> const std::string TModel<T>::Options::BOUNDER_NAME[5]
  = { "NAIVE", "LSB", "EIGEN", "BERNSTEIN", "HYBRID" };

//! @brief C++ class for Taylor model computation of factorable function - Taylor model propagation
////////////////////////////////////////////////////////////////////////
//! mc::TVar is a C++ class for propagation of Taylor models through
//! factorable functions. The template parameter corresponds to the
//! type used in computing the remainder bound.
////////////////////////////////////////////////////////////////////////
template <typename T>
class TVar
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class TVar;
  template <typename U> friend class TModel;

  template <typename U> friend TVar<U> operator+
    ( const TVar<U>& );
  template <typename U, typename V> friend TVar<U> operator+
    ( const TVar<U>&, const TVar<V>& );
  template <typename U, typename V> friend TVar<U> operator+
    ( const TVar<U>&, const V& );
  template <typename U, typename V> friend TVar<U> operator+
    ( const V&, const TVar<U>& );
  template <typename U> friend TVar<U> operator+
    ( const double, const TVar<U>& );
  template <typename U> friend TVar<U> operator+
    ( const TVar<U>&, const double );
  template <typename U> friend TVar<U> operator-
    ( const TVar<U>& );
  template <typename U, typename V> friend TVar<U> operator-
    ( const TVar<U>&, const TVar<V>& );
  template <typename U, typename V> friend TVar<U> operator-
    ( const TVar<U>&, const V& );
  template <typename U, typename V> friend TVar<U> operator-
    ( const V&, const TVar<U>& );
  template <typename U> friend TVar<U> operator-
    ( const double, const TVar<U>& );
  template <typename U> friend TVar<U> operator-
    ( const TVar<U>&, const double );
  template <typename U> friend TVar<U> operator*
    ( const TVar<U>&, const TVar<U>& );
  template <typename U> friend TVar<U> operator*
    ( const double, const TVar<U>& );
  template <typename U> friend TVar<U> operator*
    ( const TVar<U>&, const double );
  template <typename U> friend TVar<U> operator*
    ( const U&, const TVar<U>& );
  template <typename U> friend TVar<U> operator*
    ( const TVar<U>&, const U& );
  template <typename U> friend TVar<U> operator/
    ( const TVar<U>&, const TVar<U>& );
  template <typename U> friend TVar<U> operator/
    ( const double, const TVar<U>& );
  template <typename U> friend TVar<U> operator/
    ( const TVar<U>&, const double );
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const TVar<U>& );

  template <typename U> friend TVar<U> inv
    ( const TVar<U>& );
  template <typename U> friend TVar<U> sqr
    ( const TVar<U>& );
  template <typename U> friend TVar<U> sqrt
    ( const TVar<U>& );
  template <typename U> friend TVar<U> exp
    ( const TVar<U>& );
  template <typename U> friend TVar<U> log
    ( const TVar<U>& );
  template <typename U> friend TVar<U> xlog
    ( const TVar<U>& );
  template <typename U> friend TVar<U> pow
    ( const TVar<U>&, const int );
  template <typename U> friend TVar<U> pow
    ( const TVar<U>&, const double );
  template <typename U> friend TVar<U> pow
    ( const double, const TVar<U>& );
  template <typename U> friend TVar<U> pow
    ( const TVar<U>&, const TVar<U>& );
  template <typename U> friend TVar<U> monomial
    ( const unsigned int, const TVar<U>*, const int* );
  template <typename U> friend TVar<U> cos
    ( const TVar<U>& );
  template <typename U> friend TVar<U> sin
    ( const TVar<U>& );
  template <typename U> friend TVar<U> tan
    ( const TVar<U>& );
  template <typename U> friend TVar<U> acos
    ( const TVar<U>& );
  template <typename U> friend TVar<U> asin
    ( const TVar<U>& );
  template <typename U> friend TVar<U> atan
    ( const TVar<U>& );
  template <typename U> friend TVar<U> hull
    ( const TVar<U>&, const TVar<U>& );
  template <typename U> friend bool inter
    ( TVar<U>&, const TVar<U>&, const TVar<U>& );

private:
  //! @brief Pointer to Taylor model environment
  TModel<T> *_TM;
  //! @brief Pointer to internal Taylor variable in Taylor model environment
  TVar<T>* _TV() const
    { return _TM->_TV; };
  //! @brief Order of Taylor model environment
  unsigned int _nord() const
    { return _TM->_nord; };
  //! @brief Number of variables in Taylor model environment
  unsigned int _nvar() const
    { return _TM->_nvar; };
  //! @brief Total number of monomial terms in Taylor variable
  unsigned int _nmon() const
    { return _TM->_nmon; };
  //! @brief Index of first monomial term of order <tt>iord</tt> in Taylor variable
  unsigned int _posord
    ( const unsigned int iord ) const
    { return _TM->_posord[iord]; };
  //! @brief Const pointer to array of size <tt>nvar()</tt> of variable exponents in monomial term <tt>imon</tt>
  const unsigned int* _expmon
    ( const unsigned int imon ) const
    { return _TM->_expmon+imon*_TM->_nvar; };
  //! @brief Index of monomial term from product of two monomial terms <tt>imon</tt> and <tt>jmon</tt>
  unsigned int _prodmon
    ( const unsigned int imon, const unsigned int jmon ) const
    { return _TM->_prodmon[imon][jmon]; };
  //! @brief Bound on monomial term <tt>imon</tt>
  const T& _bndmon
    ( const unsigned int imon ) const
    { return _TM->_bndmon[imon]; };
  //! @brief Reference point for variable <tt>ivar</tt> in Taylor model
  double _refpoint
    ( const unsigned int ivar ) const
    { return _TM->_refpoint[ivar]; };
  //! @brief Scaling for variable <tt>ivar</tt> in Taylor model
  double _scaling
    ( const unsigned int ivar ) const
    { return _TM->_scaling[ivar]; };

public:
  /** @addtogroup TAYLOR Taylor Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Constructor of Taylor variable for a real scalar
  TVar
    ( const double d=0. );
  //! @brief Constructor of Taylor variable for a remainder bound
  TVar
    ( const T&B );
  //! @brief Constructor of Taylor variable with index <a>ix</a> (starting from 0),  bounded by <a>X</a>, and with reference point <a>Xref</a>
  TVar
    ( TModel<T>*TM, const unsigned int ix, const T&X,
      const double Xref );
  //! @brief Constructor of Taylor variable with index <a>ix</a> (starting from 0),  bounded by <a>X</a>, and with reference point at mid-point <a>Op<T>::mid(X)</a>
  TVar
    ( TModel<T>*TM, const unsigned int ix, const T&X );
  //! @brief Copy constructor of Taylor variable
  TVar
    ( const TVar<T>&TV );
  //! @brief Copy constructor of Taylor variable in different Taylor model environment (with implicit type conversion)
  template <typename U> TVar
    ( TModel<T>*&TM, const TVar<U>&TV );
  //! @brief Copy constructor of Taylor variable in different Taylor model environment (with explicit type conversion as given by class member function <a>method</a>)
  template <typename U> TVar
    ( TModel<T>*&TM, const TVar<U>&TV, const T& (U::*method)() const );
  //! @brief Copy constructor of Taylor variable in different Taylor model environment (with explicit type conversion as given by non-class member function <a>method</a>)
  template <typename U> TVar
    ( TModel<T>*&TM, const TVar<U>&TV, T (*method)( const U& ) );

  //! @brief Destructor of Taylor variable
  ~TVar()
    { delete [] _coefmon; delete [] _bndord; }

  //! @brief Set Taylor variable with index <a>ix</a> (starting from 0),  bounded by <a>X</a>, and with reference point at mid-point <a>Op<T>::mid(X)</a>
  TVar<T>& set
    ( TModel<T>*TM, const unsigned int ix, const T&X )
    { *this = TVar( TM, ix, X ); return *this; }

  //! @brief Set Taylor variable with index <tt>ix</tt> (starting from 0),  bounded by <tt>X</tt>, and with reference point <tt>Xref</tt>
  TVar<T>& set
    ( TModel<T>*TM, const unsigned int ix, const T&X,
      const double Xref )
    { *this = TVar( TM, ix, X, Xref ); return *this; }

  //! @brief Set Taylor model environment in Taylor variable to <tt>env</tt>
  TVar<T>& set
    ( TModel<T>*env )
    { *this = TVar( env ); return *this; }

  //! @brief Set multivariate polynomial coefficients in Taylor variable to <tt>coefmon</tt>
  TVar<T>& set
    ( const double*coefmon )
    { for( unsigned int imon=0; imon<(_TM?_nmon():1); imon++ )
        _coefmon[imon] = coefmon[imon];
      _update_bndord(); return *this; }

  //! @brief Set multivariate polynomial coefficients in Taylor variable to <tt>coefmon</tt> - only first <tt>coefmon.first</tt> coefficients are set
  TVar<T>& set
    ( std::pair<unsigned int, const double*> coefmon )
    { for( unsigned int imon=0; imon<(_TM?_nmon():1); imon++ )
        _coefmon[imon] = ( imon<coefmon.first && coefmon.second?
	                   coefmon.second[imon]: 0. );
      _update_bndord(); return *this; }

  //! @brief Set remainder term in Taylor variable to <tt>bndrem</tt>
  TVar<T>& set
    ( const T&bndrem )
    { *_bndrem = bndrem; _update_bndord(); return *this; }

  //! @brief Set multivariate polynomial coefficients and remainder term to those of Taylor variable <tt>TV</tt>, possibly defined in another Taylor model environment with less variables and a lower expansion order. Coefficients involving other variables or higher order are initialized to 0 if <tt>reset=true</tt> (default), otherwise they are left unmodified. Caution has to be taken when the Taylor models are scaled (see option mc::TModel::Options::SCALE_VARIABLES) since the polynomial coefficients are then scaled.
  TVar<T>& set
    ( TVar<T>& TV, const bool reset=true )
    { if( !_TM || !TV._TM || _nvar() < TV._nvar() || _nord() < TV._nord() )
        return *this;
      // Reset coefficients to 0
      for( unsigned int imon=0; reset && imon<_nmon(); imon++ )
        _coefmon[imon] = 0.;
      // Copy coefficients from TV --> *this
      unsigned int*iexp = new unsigned int[_nvar()];
      for( unsigned int jmon=0; jmon<TV._nmon(); jmon++ ){
        for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
	  iexp[ivar] = (ivar<TV._nvar()? TV._TM->_expmon[jmon*TV._nvar()+ivar]: 0);
        _coefmon[_TM->_loc_expmon(iexp)] = TV._coefmon[jmon];
      }
      delete[] iexp;
      *_bndrem = *TV._bndrem; _update_bndord(); return *this; }

  //! @brief Copy multivariate polynomial coefficients from current Taylor variable into Taylor variable <tt>TV</tt>, possibly defined in another Taylor model environment with less variables and a lower expansion order. Copied coefficients are reset to 0 in current Taylor variable if <tt>reset=true</tt>, otherwise they are left unmodified (default). Caution has to be taken when the Taylor models are scaled (see option mc::TModel::Options::SCALE_VARIABLES) since the polynomial coefficients are then scaled.
  TVar<T>& get
    ( TVar<T>&TV, const bool reset=false )
    { if( !_TM ){
        TV._coefmon[0] = _coefmon[0];
        // Reset coefficients to 0
	if( reset ) _coefmon[0] = 0.;
        return *this;
      }
      if( !TV._TM || _nvar() < TV._nvar() || _nord() < TV._nord() )
        return *this;
      // Copy coefficients from *this --> TV
      unsigned int*iexp = new unsigned int[_nvar()];
      for( unsigned int jmon=0; jmon<(TV._TM?TV._nmon():1); jmon++ ){
        for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
	  iexp[ivar] = (ivar<TV._nvar()? TV._TM->_expmon[jmon*TV._nvar()+ivar]: 0);
        TV._coefmon[jmon] = _coefmon[_TM->_loc_expmon(iexp)];
        // Reset coefficients to 0
	if( reset ) _coefmon[_TM->_loc_expmon(iexp)] = 0.;
      }
      delete[] iexp;
      if( reset ) _update_bndord();
      *TV._bndrem = 0.; TV._update_bndord(); return *this; }

  //! @brief Get pointer to Taylor model environment
  TModel<T>* env() const
    { return _TM; }

  //! @brief Compute bound on all terms of (total) order <tt>iord</tt> in Taylor variable
  T bound
    ( const unsigned int iord ) const
    { return (!iord || (_TM && iord<=_nord()))? _bndord[iord]: 0.; }

  //! @brief Compute bound on Taylor variable (see option mc::TModel::Options::BOUNDER_TYPE)
  T bound() const
    { T bndmod; return _bound( bndmod ); }

  //! @brief Return reference to variable bound in template parameter arithmetic
  const T& boundT() const
    { return _bndT; }

  //! @brief Return remainder term of Taylor variable
  T remainder() const
    { return( *_bndrem ); }

  //! @brief Center remainder term of Taylor variable
  TVar<T>& center()
    { _center_TM(); return *this; }

  //! @brief Return new Taylor variable with same multivariate polynomial part and zero remainder
  TVar<T> polynomial() const
  { TVar<T> TV = *this; *(TV._bndrem) = 0.; return TV; }

  //! @brief Evaluate polynomial part at <tt>x</tt>
  double polynomial
    ( const double*x ) const;

  //! @brief Shortcut to mc::TVar::bound
  T B
    ( const unsigned int iord ) const
    { return iord<=_nord()? _bndord[iord]: 0.; }

  //! @brief Shortcut to mc::TVar::bound
  T B() const
    { T bndmod; return _bound( bndmod ); }

  //! @brief Shortcut to mc::TVar::remainder
  T R() const
    { return remainder(); }

  //! @brief Shortcut to mc::TVar::center
  TVar<T>& C()
    { return center(); }

  //! @brief Shortcut to mc::TVar::polynomial
  TVar<T> P() const
    { return polynomial(); }

  //! @brief Shortcut to mc::TVar::polynomial
  double P
    ( const double*x ) const
    { return polynomial( x ); }

  //! @brief Get pointer to array of size <tt>nvar</tt> with references for all variables
  double* reference() const;

  //! @brief Get coefficient of constant term in Taylor variable
  double constant() const;

  //! @brief Get pointer to array of size <tt>nvar</tt> with coefficients of linear term in Taylor variable
  double* linear() const;

  //! @brief Get coefficients of linear term for variable <tt>ivar</tt> in Taylor variable. The value of this coefficient is reset to 0 if <tt>reset=true</tt>, otherwise it is left unmodified (default).
  double  linear
    ( const unsigned int ivar, const bool reset=false );

  //! @brief Get (possibly scaled) coefficient in monomial term with variable exponents as given in <a>iexp</a>
  double coefmon
    ( const unsigned int*iexp ) const;

  //! @brief Get pair of size of, and const pointer to, array of (possibly scaled) monomial coefficients in multivariate polynomial of Taylor variable
  std::pair<unsigned int, const double*> coefmon() const;

  //! @brief Get pair of size of, and const pointer to, array of monomial exponents in multivariate polynomial of Taylor variable
  std::pair<unsigned int, const unsigned int*> expmon() const;
  /** @} */

  TVar<T>& operator =
    ( const double );
  TVar<T>& operator =
    ( const TVar<T>& );
  TVar<T>& operator =
    ( const T& );
  template <typename U> TVar<T>& operator +=
    ( const TVar<U>& );
  template <typename U> TVar<T>& operator +=
    ( const U& );
  TVar<T>& operator +=
    ( const double );
  template <typename U> TVar<T>& operator -=
    ( const TVar<U>& );
  template <typename U> TVar<T>& operator -=
    ( const U& );
  TVar<T>& operator -=
    ( const double );
  TVar<T>& operator *=
    ( const TVar<T>& );
  TVar<T>& operator *=
    ( const double );
  TVar<T>& operator *=
    ( const T& );
  TVar<T>& operator /=
    ( const TVar<T>& );
  TVar<T>& operator /=
    ( const double );

private:

  //! @brief Private constructor for real scalar in Taylor model environment <tt>TM</tt>
  TVar
    ( TModel<T>*TM, const double d=0. );
  //! @brief Private constructor for remainder bound in Taylor model environment <tt>TM</tt>
  TVar
    ( TModel<T>*TM, const T&B );

  //! @brief Array of size <tt>_nmon</tt> with (possibly scaled) monomial coefficients
  double *_coefmon;
  //! @brief Array of size <tt>_nord+2</tt> with bounds for all terms of degrees <tt>iord=0,...,_nord</tt> as well as the remainder bound at position <tt>_nord+1</tt>
  T * _bndord; 
  //! @brief Pointer to the remainder bound
  T * _bndrem;
  //! @brief Bound evaluated in T arithmetic (see option mc::TModel::Options::SCALE_VARIABLES)
  T _bndT;

  //! @brief Update bounds for all terms of degrees <tt>iord=0,...,_nord</tt> in <tt>_bndord</tt>
  void _update_bndord();
  //! @brief Center remainder error term <tt>_bndrem</tt>
  void _center_TM();

  //! @brief Range bounder
  T& _bound
    ( T& bndmod ) const;
  //! @brief Range bounder - naive approach
  T& _bound_naive
    ( T& bndmod ) const;
  //! @brief Range bounder - Lin & Stadtherr approach
  T& _bound_LSB
    ( T& bndmod ) const;
  //! @brief Range bounder - eigenvalue decomposition-based approach
  T& _bound_eigen
    ( T& bndmod ) const;
  //! @brief Range bounder - Bernstein approach
  T& _bound_bernstein
    ( T& bndmod ) const;
  //! @brief Compute Bernstein coefficient for variable with index <tt>jmon</tt> and exponents <tt>jexp</tt>, given coefficients in monomial form <tt>coefmon</tt> and maximum order <tt>maxord</tt>
  double _coef_bernstein
    ( const double*coefmon, const unsigned int*jexp,
      const unsigned int jmon, const unsigned int maxord ) const;

  //! @brief Initialize private members
  void _init();
  //! @brief Reinitialize private members
  void _reinit();
  //! @brief Clean up private members
  void _clean();

};

////////////////////////////////// TModel //////////////////////////////////////

template <typename T> inline void
TModel<T>::_size
( const unsigned int nvar, const unsigned int nord )
{
  if( !nvar ) throw Exceptions( Exceptions::SIZE );

  //_cleanup();
  _nvar = nvar;
  _nord = nord; 
  _binom = new TM_size[(nvar+nord-1)*(nord+1)];
  _binom_size = std::make_pair( nvar+nord-1, nord+1 );
  _set_binom( nord );
  _posord = new unsigned int[nord+2];
  _posord_size = nord;
  _set_posord( nord );
  _nmon = _posord[_nord+1];
  _expmon = new unsigned int[_nmon*nvar];
  _expmon_size = _nmon;
  _set_expmon( nord );
  _set_prodmon();
  _bndpow = new T*[_nvar];
  for( unsigned int i=0; i<_nvar; i++ ) _bndpow[i] = 0;
  _bndmon = new T[_nmon];  
  _refpoint = new double[_nvar];
  _scaling = new double[_nvar];
  _modvar = true;
  _cbern = new double[_nord+1];

  _TV = new TVar<T>( this );
}

template <typename T> inline void
TModel<T>::_set_bndpow
( const unsigned int ivar, const T&X, const double Xref,
  const double scaling )
{
  if( ivar>=_nvar ) throw Exceptions( Exceptions::INIT );

  delete[] _bndpow[ivar];
  _bndpow[ivar] = new T [_nord+1];
  _refpoint[ivar] = Xref/scaling;
  _scaling[ivar] = scaling;
  T Xr = X/scaling - _refpoint[ivar];
  _bndpow[ivar][0] = 1.;
  for( unsigned int i=1; i<=_nord; i++ ){
    _bndpow[ivar][i] = Op<T>::pow(Xr,(int)i);
  }
  _modvar = true;
}

template <typename T> inline void
TModel<T>::_set_bndmon()
{
  if( !_modvar ) return;
  
  _bndmon[0] = 1.;
  for( unsigned int i=1; i<_nmon; i++ ){
    _bndmon[i] = 1.;
    for( unsigned int j=0; j<_nvar; j++)
      if( _bndpow[j] ) _bndmon[i] *= _bndpow[j][_expmon[i*_nvar+j]];
  }
  _modvar = false;

#ifdef MC__TMODEL_DEBUG
  mc::display( 1, _nmon, _bndmon, 1, "_bndmon", std::cout );
#endif
}

template <typename T> inline void
TModel<T>::_set_posord
( const unsigned int nord )
{
  _posord[0] = 0;
  _posord[1] = 1;
  for( unsigned int i=1; i<=nord; i++ ){
    TM_size _posord_next = _posord[i] + _get_binom( _nvar+i-1, i );
    if( _posord_next > UINT_MAX )
      throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::MAXSIZE );
    _posord[i+1] = _posord_next;
  }

#ifdef MC__TMODEL_DEBUG
  mc::display( 1, nord+2, _posord, 1, "_posord", std::cout );
#endif
}
    
template <typename T> inline void
TModel<T>::_ext_posord
( const unsigned int maxord )
{
  if( maxord < _posord_size ) return;
  delete[] _posord;
  _posord = new unsigned int[maxord+2];
  _posord_size = maxord;
  _set_posord( maxord );
}

template <typename T> inline void
TModel<T>::_set_expmon
( const unsigned int nord )
{
  unsigned int *iexp = new unsigned int[_nvar] ;
  for( unsigned int k=0; k<_nvar; k++ ) _expmon[k] = 0;
  for( unsigned int i=1; i<=nord; i++ ){
    for( unsigned int j=0; j<_nvar; j++ ) iexp[j] = 0;
    for( unsigned int j=_posord[i]; j<_posord[i+1]; j++ ){
      _next_expmon( iexp, i );
      for( unsigned int k=0; k<_nvar; k++ )
        _expmon[j*_nvar+k] = iexp[k];
    }
  }
  delete[] iexp;

#ifdef MC__TMODEL_DEBUG
  mc::display( _nvar, _expmon_size, _expmon, _nvar, "_expmon", std::cout );
#endif
}
  
template <typename T> inline void
TModel<T>::_next_expmon
( unsigned int *iexp, const unsigned int iord ) const
{
  unsigned int curord;
  do{
    iexp[_nvar-1] += iord;
    unsigned int j = _nvar;
    while( j > 0 && iexp[j-1] > iord ){
      iexp[j-1] -= iord + 1;
      j-- ;
      iexp[j-1]++;
    }
    curord = 0;
    for( unsigned int i=0; i<_nvar; i++ ) curord += iexp[i];
  } while( curord != iord );
}
    
template <typename T> inline void
TModel<T>::_ext_expmon
( const unsigned int maxord, const bool full )
{
  _ext_binom( maxord );
  _ext_posord( maxord ); 
  delete[] _expmon;
  if( full ) _expmon_size = std::pow(maxord+1,_nvar);
  else       _expmon_size = _posord[maxord+1];
  _expmon = new unsigned int[ _expmon_size*_nvar];
  _set_expmon( maxord );
  if( !full ) return;

  unsigned int *iexp = new unsigned int[_nvar];
  for( unsigned int iord=maxord+1, jmon=_posord[maxord+1];
       jmon<_expmon_size; iord++ ){
    _ext_binom( iord );
    _ext_posord( iord );
    if( _posord[iord] >= _posord[iord+1] )
      throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::MAXSIZE );
    for( unsigned int ivar=0; ivar<_nvar; ivar++ ) iexp[ivar] = 0;
    for( unsigned int kmon=_posord[iord]; kmon<_posord[iord+1]; kmon++ ){
      _next_expmon( iexp, iord );
      bool lt_maxord = true;
      for( unsigned int ivar=0; ivar<_nvar; ivar++ ){
        if( iexp[ivar] > maxord ){ lt_maxord = false; break; }
        _expmon[jmon*_nvar+ivar] = iexp[ivar];
      }
      if( lt_maxord ){
#ifdef MC__TMODEL_DEBUG
	std::cout << jmon << ":";
        for( unsigned int ivar=0; ivar<_nvar; ivar++ )
	  std::cout << "  " << _expmon[jmon*_nvar+ivar];
	std::cout << std::endl;
#endif
        jmon++;
      }
    }
  }
  delete[] iexp;
  return;
}

template <typename T> inline void
TModel<T>::_set_prodmon()
{
  _prodmon = new unsigned int*[_nmon];
  _prodmon[0] = new unsigned int[_nmon+1];
  _prodmon[0][0] = _nmon;
  for( unsigned int i=1; i<=_nmon; i++ ) _prodmon[0][i] = i-1;
#ifdef MC__TMODEL_DEBUG
  std::ostringstream ohead0;
  ohead0 << "_prodmon[" << 0 << "]";
  mc::display( 1, _nmon+1, _prodmon[0], 1, ohead0.str(), std::cout );
#endif

  unsigned int *iexp = new unsigned int[_nvar];
  for( unsigned int i=1; i<_nord; i++ ){    
    for( unsigned int j=_posord[i]; j<_posord[i+1]; j++ ){
      _prodmon[j] = new unsigned int [_posord[_nord+1-i]+1];
      _prodmon[j][0] = _posord[_nord+1-i];
      for( unsigned int k=0; k<_posord[_nord+1-i]; k++ ){
        for( unsigned int in=0; in<_nvar; in++ ) 
          iexp[in] = _expmon[j*_nvar+in] + _expmon[k*_nvar+in] ;
        _prodmon[j][k+1] = _loc_expmon( iexp );
      }
#ifdef MC__TMODEL_DEBUG
      std::ostringstream oheadj;
      oheadj << "_prodmon[" << j << "]";
      mc::display( 1, _posord[_nord+1-i]+1, _prodmon[j], 1, oheadj.str(), std::cout );
#endif
    }
  }
  delete[] iexp;

  for( unsigned int i=_posord[_nord]; i<_nmon; i++ ){
    _prodmon[i] = new unsigned int[2];
    _prodmon[i][0] = 1;
    _prodmon[i][1] = i;
#ifdef MC__TMODEL_DEBUG
    std::ostringstream oheadi;
    oheadi << "_prodmon[" << i << "]";
    mc::display( 1, 2, _prodmon[i], 1, oheadi.str(), std::cout );
#endif
  }
}
    
template <typename T> inline unsigned int
TModel<T>::_loc_expmon
( const unsigned int *iexp ) const
{
  unsigned int ord = 0;
  for( unsigned int i=0; i<_nvar; i++ ) ord += iexp[i];
  assert( ord<_nord+2 );
  unsigned int pos = _posord[ord];
  
  unsigned int p = _nvar ; 
  for( unsigned int i=0; i<_nvar-1; i++ ){
    p--;
    for( unsigned int j=0; j<iexp[i]; j++ )
      pos += _get_binom( p-1+ord-j, ord-j );
    ord -= iexp[i];
  }

  return pos;    
}
    
template <typename T> inline void
TModel<T>::_set_binom
( const unsigned int nord )
{
  TM_size *p;
  unsigned int k;
  for( unsigned int i=0; i<_nvar+nord-1; i++ ){
    p = &_binom[i*(nord+1)];
    *p = 1;
    p++;
    *p = i+1;
    p++;
    k = ( i+1<nord? i+1: nord );
    for( unsigned int j=2; j<=k; j++, p++ ) *p = *(p-1) * (i+2-j)/j;
    for( unsigned int j=k+1; j<=nord; j++, p++ ) *p = 0.;
  }
#ifdef MC__TMODEL_DEBUG
  mc::display( _binom_size.second, _binom_size.first, _binom,
    _binom_size.second, "_binom", std::cout );
#endif
}
    
template <typename T> inline void
TModel<T>::_ext_binom
( const unsigned int maxord )
{
  if( maxord < _binom_size.second ) return;
  delete[] _binom;
  _binom = new TM_size[(_nvar+maxord-1)*(maxord+1)];
  _binom_size = std::make_pair( _nvar+maxord-1, maxord+1 );
  _set_binom( maxord );
}

template <typename T> inline TM_size
TModel<T>::_get_binom
( const unsigned int n, const unsigned int k ) const
{
#ifdef MC__TMODEL_CHECK
  assert( n<=_binom_size.first );
  assert( k<=_binom_size.second );
  assert( k<=n );
#endif
  return( n? _binom[(n-1)*_binom_size.second+k]: 1. );
}

template <typename T> inline void
TModel<T>::_scale
( const unsigned int ivar, double*coef ) const
{
  if( ivar>=_nvar || !_nord ) return;
  double dscale = Op<T>::diam(_bndpow[ivar][1]);
  double fscale = Op<T>::l(_bndpow[ivar][1]) / dscale;
  unsigned int iexp[_nvar];
  for( unsigned int imon=0; imon<_nmon; imon++ ){
    unsigned int iord = 0;
    for( unsigned int jvar=0; jvar<_nvar; jvar++ ){
      iexp[jvar] = _expmon[imon*_nvar+jvar];
      iord += iexp[jvar];
    }
    double coefmod = coef[imon] * std::pow(dscale,iexp[ivar]);
    iexp[ivar]++;
    for( unsigned int k=1; k<=_nord-iord; k++, iexp[ivar]++ ){
      coefmod += coef[_loc_expmon(iexp)]
        * _get_binom(iexp[ivar],k) * std::pow(fscale,k)
        * std::pow(dscale,iexp[ivar]);
    }
#ifdef  MC__TMODEL_DEBUG_SCALE
    std::cout << "    a" << std::left << std::setw(3) << imon << " = "
        << std::right << std::setw(12) << coef[imon]
        << std::right << std::setw(14) << coefmod << std::endl;
#endif
    coef[imon] = coefmod;
  }
  return;
}

template <typename T> inline TVar<T>
TModel<T>::_univ_bernstein
( const TVar<T>&TV, puniv f, puniv df, punivext If,
  punivext d2If, const double*rusr, const int*iusr )
{
  assert( TV._TM == this );
  
  const T B( TV.B() );
  const TVar<T> TVs( (TV-Op<T>::l(B))/Op<T>::diam(B) );
  TVar<T> TV2( this, 0. ), MON( 1. );
  for( unsigned int j=0; j<=_nord; j++ ){
    double sign = ( j%2? -1.: 1. );
    _cbern[j] = sign * _get_binom(_nord,j) * f( Op<T>::l(B), rusr, iusr );
    for( unsigned int i=1; i<=j; i++ ){
      sign *= -1.;
      _cbern[j] += sign * _get_binom(_nord,i) * _get_binom(_nord-i,j-i)
        * f( Op<T>::l(B)+(double)i/(double)_nord*Op<T>::diam(B), rusr, iusr );
    }
    TV2 += _cbern[j] * MON;
    MON *= TVs;
  }

  T R( - d2If( B, rusr, iusr ) * sqr(Op<T>::diam(B)) / 2. / _nord );
  if( Op<T>::l( R ) * Op<T>::u( R ) < 0 )
    return TV2 + Op<T>::zeroone() * R;

  if( !options.BERNSTEIN_OPT ){
    T R0( If(B,rusr,iusr) );
    R = Op<T>::zeroone() * ( Op<T>::u(R) <= 0.?
         -std::min( Op<T>::diam(R0), -Op<T>::l(R) ):
          std::min( Op<T>::diam(R0),  Op<T>::u(R) ) );
    return TV2 + R;
  }

  std::pair<unsigned int,const double*> bern = std::make_pair(_nord+1,_cbern);
  double xopt = _goldsect( 0., 1., _dgap_bernstein,  rusr, iusr, df, B, bern );
  R = Op<T>::zeroone() * _gap_bernstein( xopt, rusr, iusr, f, B, bern );
  return TV2 + R;
}

template <typename T> inline double
TModel<T>::_dgap_bernstein
( const double x, const double*rusr, const int*iusr, puniv df,
  const T&B, const std::pair<unsigned int,const double*>&bern )
{
  double phi = df(Op<T>::l(B)+x*Op<T>::diam(B),rusr,iusr)*Op<T>::diam(B);
  for( unsigned int j=1; j<bern.first; j++ )
    phi -= j * bern.second[j] * std::pow(x,j-1);
  return phi;
}

template <typename T> inline double
TModel<T>::_gap_bernstein
( const double x, const double*rusr, const int*iusr, puniv f,
  const T&B, const std::pair<unsigned int,const double*>&bern )
{
  double phi = f(Op<T>::l(B)+x*Op<T>::diam(B),rusr,iusr);
  for( unsigned int j=0; j<bern.first; j++ )
    phi -= bern.second[j] * std::pow(x,j);
  return phi;
}

template <typename T> inline double
TModel<T>::_goldsect
( const double xL, const double xU, punivopt fopt, const double*rusr,
  const int*iusr, puniv df,  const T&Ix,
  const std::pair<unsigned int,const double*>&bern )
{
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  const double fL = fopt( xL, rusr, iusr, df, Ix, bern ),
               fU = fopt( xU, rusr, iusr, df, Ix, bern );
  if( fL*fU > 0 ) throw Exceptions( Exceptions::BERNSTEIN );
  const double xm = xU-phi*(xU-xL),
               fm = fopt( xm, rusr, iusr, df, Ix, bern );
  return _goldsect_iter( true, xL, fL, xm, fm, xU, fU, fopt, rusr, iusr,
                         df, Ix, bern );
}

template <typename T> inline double
TModel<T>::_goldsect_iter
( const bool init, const double a, const double fa, const double b,
  const double fb, const double c, const double fc, punivopt fopt,
  const double*rusr, const int*iusr, puniv df, const T&Ix,
  const std::pair<unsigned int,const double*>&bern )
// a and c are the current bounds; the minimum is between them.
// b is a center point
{
  static unsigned int iter;
  iter = ( init? 1: iter+1 );
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  bool b_then_x = ( c-b > b-a );
  double x = ( b_then_x? b+phi*(c-b): b-phi*(b-a) );
  if( std::fabs(c-a) < options.BERNSTEIN_TOL*(std::fabs(b)+std::fabs(x)) 
   || iter > options.BERNSTEIN_MAXIT ) return (c+a)/2.;
  double fx = fopt( x, rusr, iusr, df, Ix, bern );
  if( b_then_x )
    return( fa*fx<0?
      _goldsect_iter( false, a, fa, b, fb, x, fx, fopt, rusr, iusr, df, Ix, bern ):
     _goldsect_iter( false, b, fb, x, fx, c, fc, fopt, rusr, iusr, df, Ix, bern ) );
  return( fa*fb<0?
     _goldsect_iter( false, a, fa, x, fx, b, fb, fopt, rusr, iusr, df, Ix, bern ):
     _goldsect_iter( false, x, fx, b, fb, c, fc, fopt, rusr, iusr, df, Ix, bern ) );
}

template <typename T> inline void
TModel<T>::_reset()
{
  for( unsigned int i=0; i<_nvar; i++ ){
    delete[] _bndpow[i];
    _bndpow[i] = 0;
  }
}

template <typename T> inline void
TModel<T>::_cleanup()
{
  for( unsigned int i=0; i<_nmon; i++ ) delete[] _prodmon[i];
  delete[] _prodmon;
  delete[] _expmon;
  delete[] _posord;
  for( unsigned int i=0; i<_nvar; i++ ) delete[] _bndpow[i];
  delete[] _bndpow;
  delete[] _bndmon;
  delete[] _refpoint;
  delete[] _scaling;
  delete[] _cbern;
  delete[] _binom;
  delete _TV;
}

////////////////////////////////// TVar ///////////////////////////////////////

template <typename T> inline
TVar<T>::TVar
( const double d )
: _TM( 0 ), _bndT( d )
{
  _init();
  _coefmon[0] = d;
  _bndord[0] = 0.;
}

template <typename T> inline TVar<T>&
TVar<T>::operator =
( const double d )
{
  if( _TM ){ _TM = 0; _reinit(); }
  _coefmon[0] = d;
  _bndord[0] = 0.;
  return *this;
}

template <typename T> inline
TVar<T>::TVar
( TModel<T>*TM, const double d )
: _TM( TM )
{
  if( !_TM ){
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INIT );
  }
  _init();
  _coefmon[0] = d;
  for( unsigned int i=1; i<_nmon(); i++ ) _coefmon[i] = 0.;
  _bndord[0] = d;
  for( unsigned int i=1; i<_nord()+2; i++) _bndord[i] = 0.;
  if( _TM->options.PROPAGATE_BNDT ) _bndT = d;
}

template <typename T> inline
TVar<T>::TVar
( const T&B )
: _TM( 0 ), _bndT( B )
{
  _init();
  _coefmon[0] = 0.;
  _bndord[0] = B;
}

template <typename T> inline TVar<T>&
TVar<T>::operator =
( const T&B )
{
  if( _TM ){ _TM = 0; _reinit(); }
  _coefmon[0] = 0.;
  _bndord[0] = B;
  return *this;
}

template <typename T> inline
TVar<T>::TVar
( TModel<T>*TM, const T&B )
: _TM( TM )
{
  if( !_TM ) throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INIT );
  _init();
  for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] = 0.;
  for( unsigned int i=0; i<_nord()+1; i++) _bndord[i] = 0.;
  *_bndrem = B;
  if( _TM->options.PROPAGATE_BNDT ) _bndT = B;
  if( _TM->options.CENTER_REMAINDER ) _center_TM();
}

template <typename T> inline
TVar<T>::TVar
( const TVar<T>&TV )
: _TM(0)
{
  _init();
  *this = TV;
}

template <typename T> inline TVar<T>&
TVar<T>::operator =
( const TVar<T>&TV )
{
  // Same TVar
  if( this == &TV ) return *this;

  // Reinitialization needed?
  if( _TM != TV._TM ){ _TM = TV._TM; _reinit(); }

  // Set to TVar not linked to TModel (either scalar or range)
  if( !_TM ){
    _coefmon[0] = TV._coefmon[0];
    _bndord[0] = TV._bndord[0];
    return *this; 
  }
  // Set to TVar linked to TModel
  for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] = TV._coefmon[i];
  for( unsigned int i=0; i<_nord()+2; i++) _bndord[i] = TV._bndord[i];
  if( _TM->options.PROPAGATE_BNDT ) _bndT = TV._bndT;
  return *this;
}

template <typename T> template <typename U> inline
TVar<T>::TVar
( TModel<T>*&TM, const TVar<U>&TV )
: _TM(TM), _coefmon(0), _bndord(0), _bndrem(0)
{
  _init();
  TVar<U> TVtrunc( TV );
  _coefmon[0] = TVtrunc._coefmon[0];
  TVtrunc._coefmon[0] = 0. ;
  for( unsigned int i=1; _TM && i<_nmon(); i++ ){
    if( TVtrunc._TM && i < TVtrunc._nmon() ){
      _coefmon[i] = TVtrunc._coefmon[i];
      TVtrunc._coefmon[i] = 0.;
    }
    else
      _coefmon[i] = 0.;
  }
  TVtrunc._update_bndord();
  *_bndrem = T( TVtrunc.B() );
  if( !_TM ) return;
  _update_bndord();
  if( _TM->options.PROPAGATE_BNDT ) _bndT = T( TV._bndT );
  return;
}

template <typename T> template <typename U> inline
TVar<T>::TVar
( TModel<T>*&TM, const TVar<U>&TV, T (*method)( const U& ) )
: _TM(TM), _coefmon(0), _bndord(0), _bndrem( 0 )
{
  if( !method )
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INIT );
  _init();
  TVar<U> TVtrunc( TV );
  _coefmon[0] = TVtrunc._coefmon[0];
  TVtrunc._coefmon[0] = 0. ;
  for( unsigned int i=1; _TM && i<_nmon(); i++ ){
    if( TVtrunc._TM && i < TVtrunc._nmon() ){
      _coefmon[i] = TVtrunc._coefmon[i];
      TVtrunc._coefmon[i] = 0.;
    }
    else
      _coefmon[i] = 0.;
  }
  TVtrunc._update_bndord();
  *_bndrem = (*method)( TVtrunc.B() );
  if( !_TM ) return;
  _update_bndord();
  if( _TM->options.PROPAGATE_BNDT ) _bndT = (*method)( TV._bndT );
  return;
}

template <typename T> template <typename U> inline
TVar<T>::TVar
( TModel<T>*&TM, const TVar<U>&TV, const T& (U::*method)() const )
: _TM(TM), _coefmon(0), _bndord(0), _bndrem( 0 )
{
  if( !method )
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INIT );
  _init();
  TVar<U> TVtrunc( TV );
  _coefmon[0] = TVtrunc._coefmon[0];
  TVtrunc._coefmon[0] = 0. ;
  for( unsigned int i=1; _TM && i<_nmon(); i++ ){
    if( TVtrunc._TM && i < TVtrunc._nmon() ){
      _coefmon[i] = TVtrunc._coefmon[i];
      TVtrunc._coefmon[i] = 0.;
    }
    else
      _coefmon[i] = 0.;
  }
  TVtrunc._update_bndord();
  *_bndrem = (TVtrunc.B().*method)();
  if( !_TM ) return;
  _update_bndord();
  if( _TM->options.PROPAGATE_BNDT ) _bndT = (TV._bndT.*method)();
  return;
}

template <typename T> inline
TVar<T>::TVar
( TModel<T>*TM, const unsigned int ivar, const T&X )
: _TM( TM )
{
  if( !TM ){
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INIT );
  }

  // Scale variables and keep track of them in TModel
  double scaling = ( _TM->options.SCALE_VARIABLES? Op<T>::diam(X)/2.: 1. );
  if( isequal( scaling, 0. ) ) scaling = 1.;
  _TM->_set_bndpow( ivar, X, Op<T>::mid(X), scaling );
  _TM->_set_bndmon();
  _init();

  // Populate _coefmon w/ TVar coefficients
  _coefmon[0] = Op<T>::mid(X);
  for( unsigned int i=1; i<_nmon(); i++ ) _coefmon[i] = 0.;
  if( _nord() > 0 ) _coefmon[_nvar()-ivar] = scaling;

  // Populate _bndord w/ bounds on TVar terms
  _bndord[0] = _coefmon[0];
  _bndord[1] = X-_coefmon[0];
  for( unsigned int i=2; i<_nord()+2; i++) _bndord[i] = 0.;
  if( _TM->options.PROPAGATE_BNDT ) _bndT = X;
}

template <typename T> inline
TVar<T>::TVar
( TModel<T>*TM, const unsigned int ivar, const T&X, const double Xref )
: _TM( TM )
{
  if( !TM ){
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INIT );
  }

  // Scale variables and keep track of them in TModel
  double scaling = ( _TM->options.SCALE_VARIABLES? Op<T>::diam(X)/2.: 1. );
  if( isequal( scaling, 0. ) ) scaling = 1.;
  _TM->_set_bndpow( ivar, X, Xref, scaling );
  _TM->_set_bndmon();
  _init();

  // Populate _coefmon w/ TVar coefficients
  _coefmon[0] = Xref;
  for( unsigned int i=1; i<_nmon(); i++ ) _coefmon[i] = 0.;
  if( _nord() > 0 ) _coefmon[_nvar()-ivar] = scaling;

  // Populate _bndord w/ bounds on TVar terms
  _bndord[0] = _coefmon[0];
  _bndord[1] = X-_coefmon[0];
  for( unsigned int i=2; i<_nord()+2; i++) _bndord[i] = 0.;
  if( _TM->options.PROPAGATE_BNDT ) _bndT = X;
}

template <typename T> inline void
TVar<T>::_init()
{
  if( !_TM ){
    _coefmon = new double[1];
    _bndord  = new T[1];
    _bndrem  = _bndord;
    return;
  }
  _coefmon = new double[_nmon()];
  _bndord  = new T[_nord()+2];
  _bndrem  = _bndord + _nord()+1;
}

template <typename T> inline void
TVar<T>::_clean()
{
  delete [] _coefmon; delete [] _bndord;
  _coefmon = 0; _bndord = _bndrem = 0;
}

template <typename T> inline void
TVar<T>::_reinit()
{
  _clean(); _init();
}

template <typename T> inline void
TVar<T>::_update_bndord()
{
  if( !_TM ) return;
  _TM->_set_bndmon();
  _bndord[0] = _coefmon[0];
  for( unsigned int i=1; i<=_nord(); i++ ){
    _bndord[i] = 0.; 
    for( unsigned int j=_posord(i); j<_posord(i+1); j++ )
      _bndord[i] += _coefmon[j] * _bndmon(j);
  }
}

template <typename T> inline void
TVar<T>::_center_TM()
{
  const double remmid = Op<T>::mid(*_bndrem);
  _coefmon[0] += remmid;
  if( _TM ) _bndord[0] = _coefmon[0];
  *_bndrem -= remmid;
}

template <typename T> inline T&
TVar<T>::_bound_eigen
( T& bndmod ) const
{
  static const double TOL = 1e-8;

  bndmod = _coefmon[0];
  if( _nord() == 1 ) bndmod += _bndord[1];

  else if( _nord() > 1 ){
    double*U = new double[_nvar()*_nvar()];
    for( unsigned int i=0; i<_nvar(); i++ ){
      for( unsigned int j=0; j<i; j++ ){
        U[_nvar()*(_nvar()-i-1)+_nvar()-j-1] = 0.;
        U[_nvar()*(_nvar()-j-1)+_nvar()-i-1] = _coefmon[_prodmon(i+1,j+2)]/2.;
      }
      U[(_nvar()+1)*(_nvar()-i-1)] = _coefmon[_prodmon(i+1,i+2)];
    }
    double*D = mc::dsyev_wrapper( _nvar(), U, true );
    if( !D ){
      delete[] U;
      throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::EIGEN );
    }

#ifdef MC__TVAR_HYBRID_EIGEN
    T bndtype1(0.);
#endif
    T bndtype2(0.);
    for( unsigned int i=0; i<_nvar(); i++ ){
      double linaux = 0.;
      T bndaux(0.);
      for( unsigned int k=0; k<_nvar(); k++ ){
        linaux += U[i*_nvar()+k] * _coefmon[_nvar()-k];
        bndaux += U[i*_nvar()+k] * _bndmon(_nvar()-k);
      }
#ifdef MC__TVAR_DEBUG_EIGEN
      std::cout << i << ": LINAUX = " << linaux
                << "  BNDAUX = " << bndaux << std::endl;
#endif
#ifdef MC__TVAR_HYBRID_EIGEN
      bndtype1 += _coefmon[i+1] * _bndmon(i+1) + D[i] * Op<T>::sqr( bndaux );
#endif
#ifdef MC__TVAR_DEBUG_EIGEN
      std::cout << std::endl << "BNDTYPE1: " << bndtype1 << std::endl;
#endif
      if( std::fabs(D[i]) > TOL )
        bndtype2 += D[i] * Op<T>::sqr( linaux/D[i]/2. + bndaux )
          - linaux*linaux/D[i]/4.;
      else     
        //bndtype2 += _coefmon[i+1] * _bndmon(i+1) + D[i] * Op<T>::sqr( bndaux );
        bndtype2 += linaux * bndaux + D[i] * Op<T>::sqr( bndaux );
#ifdef MC__TVAR_DEBUG_EIGEN
        std::cout << "BNDTYPE2: " << bndtype2 << std::endl;
#endif
    }
    delete[] U;
    delete[] D;

#ifdef MC__TVAR_HYBRID_EIGEN
    if( !Op<T>::inter( bndtype1, bndtype1, bndtype2 ) ){
      throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INCON );
#ifdef MC__TVAR_DEBUG_EIGEN
      std::cout << "BNDTYPE3: " << bndtype1 << std::endl;
#endif
    }
    bndmod += bndtype1;
#else
    bndmod += bndtype2;
#endif
  }
#ifdef MC__TVAR_DEBUG_EIGEN
      int tmp; std::cin >> tmp;
#endif

  for( unsigned int i=3; i<=_nord(); i++ ) bndmod += _bndord[i];
  bndmod += *_bndrem;

  return bndmod;
}

template <typename T> inline T&
TVar<T>::_bound_LSB
( T& bndmod ) const
{
  static const double TOL = 1e-8;
  bndmod = _coefmon[0];
  if( _nord() == 1 ) bndmod += _bndord[1];
  else if( _nord() > 1 ){
    for( unsigned int i=1; i<=_nvar(); i++ ){
      // linear and diagonal quadratic terms
      unsigned int ii = _prodmon(i,i+1);
      if( std::fabs(_coefmon[ii]) > TOL )
        bndmod += _coefmon[ii] * Op<T>::sqr( _coefmon[i]/_coefmon[ii]/2.
          + _bndmon(i) ) - _coefmon[i]*_coefmon[i]/_coefmon[ii]/4.;
      else
        bndmod += _coefmon[i] * _bndmon(i) + _coefmon[ii] * _bndmon(ii);
      // off-diagonal quadratic terms
      for( unsigned int k=i+1; k<=_nvar(); k++ ){
	unsigned int ik = _prodmon(i,k+1) ;
	bndmod += _coefmon[ik] * _bndmon(ik);
      }
    }
  }
  // higher-order terms
  for( unsigned int i=3; i<=_nord(); i++ ) bndmod += _bndord[i];
  bndmod += *_bndrem;
  return bndmod;
}

template <typename T> inline T&
TVar<T>::_bound_bernstein
( T& bndmod ) const
{
  // Scale Taylor model in unit hypercube with reference point at the origin
  double *coeftrans = new double[_nmon()];
  for( unsigned int imon=0; imon<_nmon(); imon++ )
    coeftrans[imon] = _coefmon[imon]; 
  for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
    _TM->_scale( ivar, coeftrans );
  
  // Expand binomial coefficient and exponent arrays if needed
  const unsigned int maxord = (_TM->options.BOUNDER_ORDER>_nord()? 
    _TM->options.BOUNDER_ORDER: _nord() );
  const unsigned int maxmon = std::pow(maxord+1,_nvar());
  _TM->_ext_expmon( maxord, true );

  // Compute min/max amongst all Bernstein coefficients
  bndmod = coeftrans[0];
#ifdef  MC__TVAR_DEBUG_BERSTEIN
  std::cout << "\n0:  " << bndmod << std::endl;
#endif
  for( unsigned int jmon=1; jmon<maxmon; jmon++ ){
    const unsigned int*jexp = _TM->_expmon + jmon*_nvar();
    const double coefbern = _coef_bernstein( coeftrans, jexp, jmon, maxord );
    bndmod = Op<T>::hull( bndmod, coefbern );
#ifdef  MC__TVAR_DEBUG_BERNSTEIN
    std::cout << jmon << " ["; 
    for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
      std::cout << std::setw(3) << jexp[ivar];
    std::cout << "] : " << coefbern << " : " << bndmod << std::endl;
#endif
  }

  delete[] coeftrans;
  bndmod += *_bndrem;
  return bndmod;
}

template <typename T> inline double
TVar<T>::_coef_bernstein
( const double*coefmon, const unsigned int*jexp,
  const unsigned int jmon, const unsigned int maxord ) const
{
  // Compute bernstein coefficient with variables indices <tt>jexp</tt>
  double coefbern = coefmon[0];
  for( unsigned int imon=1; imon<=std::min(jmon,_nmon()-1); imon++ ){
    const unsigned int*iexp = _TM->_expmon + imon*_nvar();
    // Only append term if monomial degrees are lower
    bool inrange = true;
    for( unsigned int ivar=0; ivar<_nvar() && inrange; ivar++ )
      if( iexp[ivar] > jexp[ivar] ) inrange = false;
    if( !inrange ) continue;
    double termbern = coefmon[imon];
#ifdef  MC__TVAR_DEBUG_BERSTEIN
    std::cout << "  ["; 
    for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
      std::cout << std::setw(3) << iexp[ivar];
    std::cout << "]";
#endif
    for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
      termbern *= (double)_TM->_get_binom(jexp[ivar],iexp[ivar])
                / (double)_TM->_get_binom(maxord,iexp[ivar]);
    coefbern += termbern;
  }
#ifdef  MC__TVAR_DEBUG_BERSTEIN
  std::cout << std::endl;
#endif
  return coefbern;
}

template <typename T> inline T&
TVar<T>::_bound_naive
( T& bndmod ) const
{
  bndmod = _coefmon[0];
  for( unsigned int i=1; i<=_nord()+1; i++ ) bndmod += _bndord[i];
  return bndmod;
}

template <typename T> inline T&
TVar<T>::_bound
( T& bndmod ) const
{
  if( !_TM ){ bndmod = _coefmon[0] + _bndord[0]; return bndmod; }

  switch( _TM->options.BOUNDER_TYPE ){
  case TModel<T>::Options::NAIVE:     bndmod = _bound_naive(bndmod);     break;
  case TModel<T>::Options::LSB:       bndmod = _bound_LSB(bndmod);       break;
  case TModel<T>::Options::EIGEN:     bndmod = _bound_eigen(bndmod);     break;
  case TModel<T>::Options::BERNSTEIN: bndmod = _bound_bernstein(bndmod); break;
  case TModel<T>::Options::HYBRID: default:{
    T bndlsb(0.), bndeig(0.);
    if( !Op<T>::inter( bndmod, _bound_LSB(bndlsb), _bound_eigen(bndeig) ) )
      throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INCON );
   }
  }

  if( _TM->options.PROPAGATE_BNDT && _TM->options.INTER_WITH_BNDT
    && !Op<T>::inter( bndmod, bndmod, _bndT ) )
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::INCON );

  return bndmod;
}

template <typename T> inline double
TVar<T>::polynomial
( const double*x ) const
{
  if( !_TM ) return _coefmon[0];
  double Pval = _coefmon[0];
  for( unsigned int i=1; i<_nmon(); i++ ){
    double valmon = 1.;
    for( unsigned int k=0; k<_nvar(); k++ )
      valmon *= std::pow( x[k]/_scaling(k)-_refpoint(k), _expmon(i)[k] );
    Pval += _coefmon[i] * valmon;
  }
  return Pval;
}

template <typename T> inline double*
TVar<T>::reference() const
{
  if( !_TM ) return 0;
  if( _nvar() < 1 ) return 0;
  double*pref = new double[_nvar()];
  for( unsigned int i=0; i<_nvar(); i++ ) pref[i] = _refpoint(i)*_scaling(i);
  return pref;
}

template <typename T> inline double
TVar<T>::constant() const
{
  return _coefmon[0];
}

template <typename T> inline double*
TVar<T>::linear() const
{
  if( !_TM || !_nvar() || !_nord() ) return 0;

  double*plin = new double[_nvar()];
  for( unsigned int i=0; i<_nvar(); i++ )
    plin[i] = _coefmon[_nvar()-i] / _scaling(i);
  return plin;
}

template <typename T> inline double
TVar<T>::linear
( const unsigned int ivar, const bool reset )
{
  if( !_TM || ivar>=_nvar() || !_nord() ) return 0.;
  const double coeflin = _coefmon[_nvar()-ivar] / _scaling(ivar);
  if( reset ){ _coefmon[_nvar()-ivar] = 0.; _update_bndord(); }
  return coeflin;
}

template <typename T> inline double
TVar<T>::coefmon
( const unsigned int*iexp ) const
{
  if( !_TM ) return 0;
  const unsigned int imon = _TM->_loc_expmon( iexp );
  return( imon<_nmon()? _coefmon[imon]: 0. );
}

template <typename T> inline std::pair<unsigned int, const double*>
TVar<T>::coefmon() const
{
  return std::make_pair( (_TM?_nmon():1), _coefmon );
}

template <typename T> inline std::pair<unsigned int, const unsigned int*>
TVar<T>::expmon() const
{
  return std::make_pair( (_TM?_nmon()*_nvar():1), _TM->_expmon );
}

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const TVar<T>&TV )
{
  out << std::endl
      << std::scientific << std::setprecision(5)
      << std::right;

  // Constant model
  if( !TV._TM ){
    out << "   a0    = " << std::right << std::setw(12) << TV._coefmon[0]
        << "      0  0"
        << std::endl
        << "   R     = " << *(TV._bndrem) << std::endl;
  }

  // Monomial term coefficients and corresponding exponents
  else{
    out << std::setprecision(TV._TM->options.DISPLAY_DIGITS);
    for( unsigned int i=0; i<TV._nmon(); i++ ){
      double scal = 1.;
      for( unsigned int k=0; k<TV._nvar(); k++ )
        scal *= std::pow( TV._scaling(k), TV._expmon(i)[k] );
      out << "   a" << std::left << std::setw(4) << i << " = "
          << std::right << std::setw(TV._TM->options.DISPLAY_DIGITS+7)
	  << TV._coefmon[i]*scal << "   ";
      for( unsigned int k=0; k<TV._nvar(); k++ )
        out << std::setw(3) << TV._expmon(i)[k];
      out << std::endl;
    }
    // Remainder term
    out << std::right << "   R     =  " << *(TV._bndrem)
        << std::endl;
  }

  // Range bounder
  out << std::right << "   B     =  " << TV.B()
      << std::endl;

  return out;
}

template <typename T> inline TVar<T>
operator +
( const TVar<T>&TV )
{
  return TV;
}

template <typename T> template <typename U> inline TVar<T>&
TVar<T>::operator +=
( const TVar<U>&TV )
{
  if( !TV._TM ){
    if( _TM && _TM->options.PROPAGATE_BNDT ) _bndT += TV._bndT;
    _coefmon[0] += TV._coefmon[0];
    *_bndrem += *(TV._bndrem);
  }
  else if( !_TM ){
    TVar<T> TV2(*this);
    *this = TV; *this += TV2;
  }
  else{
    if( _TM != TV._TM )
      throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::TMODEL );
    if( _TM->options.PROPAGATE_BNDT ) _bndT += TV._bndT;
    for( unsigned int i=0; i<_nmon(); i++ )
      _coefmon[i] += TV._coefmon[i];
    *_bndrem += *(TV._bndrem);
    _update_bndord();
  }
  if( _TM && _TM->options.CENTER_REMAINDER ) _center_TM();
  return *this;
}

template <typename T, typename U> inline TVar<T>
operator +
( const TVar<T>&TV1, const TVar<U>&TV2 )
{
  TVar<T> TV3( TV1 );
  TV3 += TV2;
  return TV3;
}

template <typename T> inline TVar<T>&
TVar<T>::operator +=
( const double c )
{
  _coefmon[0] += c;
  if( _TM ) _bndord[0] = _coefmon[0];
  if( _TM && _TM->options.PROPAGATE_BNDT ) _bndT += c;
  return *this;
}

template <typename T> inline TVar<T>
operator +
( const TVar<T>&TV1, const double c )
{
  TVar<T> TV3( TV1 );
  TV3 += c;
  return TV3;
}

template <typename T> inline TVar<T>
operator +
( const double c, const TVar<T>&TV2 )
{
  TVar<T> TV3( TV2 );
  TV3 += c;
  return TV3;
}

template <typename T> template <typename U> inline TVar<T>&
TVar<T>::operator +=
( const U&I )
{
  *_bndrem += I;
  if( _TM && _TM->options.PROPAGATE_BNDT ) _bndT += I;
  if( _TM && _TM->options.CENTER_REMAINDER ) _center_TM();
  return *this;
}

template <typename T, typename U> inline TVar<T>
operator +
( const TVar<T>&TV1, const U&I )
{
  TVar<T> TV3( TV1 );
  TV3 += I;
  return TV3;
}

template <typename T, typename U> inline TVar<T>
operator +
( const U&I, const TVar<T>&TV2 )
{
  TVar<T> TV3( TV2 );
  TV2 += I;
  return TV3;
}

template <typename T> inline TVar<T>
operator -
( const TVar<T>&TV )
{
  if( !TV._TM ){
    TVar<T> TV2;
    TV2._coefmon[0] = -TV._coefmon[0];
    TV2._bndord[0] = -TV._bndord[0];
    return TV2;
  }
  TVar<T>& TV2 = *TV._TV();
  //TVar<T> TV2( TV._TM );
  for( unsigned int i=0; i<TV._nmon(); i++ ) TV2._coefmon[i] = -TV._coefmon[i];
  for( unsigned int i=0; i<TV._nord()+2; i++ ) TV2._bndord[i] = -TV._bndord[i];
  if( TV._TM->options.PROPAGATE_BNDT ) TV2._bndT = -TV._bndT;
  if( TV._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> template <typename U> inline TVar<T>&
TVar<T>::operator -=
( const TVar<U>&TV )
{
  if( !TV._TM ){
    if( _TM && _TM->options.PROPAGATE_BNDT ) _bndT -= TV._bndT;
    _coefmon[0] -= TV._coefmon[0];
    *_bndrem -= *(TV._bndrem);
  }
  else if( !_TM ){
    TVar<T> TV2(*this);
    *this = -TV; *this += TV2;
  }
  else{
    if( _TM != TV._TM )
      throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::TMODEL );
    if( _TM->options.PROPAGATE_BNDT ) _bndT -= TV._bndT;
    for( unsigned int i=0; i<_nmon(); i++ )
      _coefmon[i] -= TV._coefmon[i];
    *_bndrem -= *(TV._bndrem);
    _update_bndord();
  }
  if( _TM && _TM->options.CENTER_REMAINDER ) _center_TM();
  return *this;
}

template <typename T, typename U> inline TVar<T>
operator-
( const TVar<T>&TV1, const TVar<U>&TV2 )
{
  TVar<T> TV3( TV1 );
  TV3 -= TV2;
  return TV3;
}

template <typename T> inline TVar<T>&
TVar<T>::operator -=
( const double c )
{
  _coefmon[0] -= c;
  if( _TM ) _bndord[0] = _coefmon[0];
  if( _TM && _TM->options.PROPAGATE_BNDT ) _bndT -= c;
  return *this;
}

template <typename T> inline TVar<T>
operator -
( const TVar<T>&TV1, const double c )
{
  TVar<T> TV3( TV1 );
  TV3 -= c;
  return TV3;
}

template <typename T> inline TVar<T>
operator -
( const double c, const TVar<T>&TV2 )
{
  TVar<T> TV3( -TV2 );
  TV3 += c;
  return TV3;
}

template <typename T> template <typename U> inline TVar<T>&
TVar<T>::operator -=
( const U&I )
{
  *_bndrem -= I;
  if( _TM && _TM->options.PROPAGATE_BNDT ) _bndT -= I;
  if( _TM && _TM->options.CENTER_REMAINDER ) _center_TM();
  return *this;
}

template <typename T, typename U> inline TVar<T>
operator -
( const TVar<T>&TV1, const U&I )
{
  TVar<T> TV3( TV1 );
  TV3 -= I;
  return TV3;
}

template <typename T, typename U> inline TVar<T>
operator -
( const U&I, const TVar<T>&TV2 )
{
  TVar<T> TV3( -TV2 );
  TV3 += I;
  return TV3;
}

template <typename T> inline TVar<T>&
TVar<T>::operator *=
( const TVar<T>&TV )
{
   TVar<T> TV2( *this );
   *this = TV * TV2;
   return *this;
}

template <typename T> inline TVar<T>
operator *
( const TVar<T>&TV1, const TVar<T>&TV2 )
{
  if( !TV2._TM )      return( TV1 * TV2._coefmon[0] + TV1 * *(TV2._bndrem) );
  else if( !TV1._TM ) return( TV2 * TV1._coefmon[0] + TV2 * *(TV1._bndrem) );

  if( TV1._TM != TV2._TM )
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::TMODEL );
  TVar<T>& TV3 = *TV1._TV();
  for( unsigned int i=0; i<TV3._nmon(); i++ ) TV3._coefmon[i] = 0.;
  //TVar<T> TV3( TV1._TM, 0. );

  // Populate _coefmon for product term
  for( unsigned int i=0; i<TV3._posord(TV3._nord()/2+1); i++){
    TV3._coefmon[TV3._prodmon(i,i+1)] += TV1._coefmon[i] * TV2._coefmon[i];
    for( unsigned int j=i+1; j<TV3._prodmon(i,0); j++ )
      TV3._coefmon[TV3._prodmon(i,j+1)] += TV1._coefmon[i] * TV2._coefmon[j]
                                         + TV1._coefmon[j] * TV2._coefmon[i];
  }
  // Calculate remainder term _bndrem for product term
  T s1 = 0., s2 = 0.;
  for( unsigned int i=0; i<=TV3._nord()+1; i++ ){
    T r1 = 0., r2 = 0.;
    for( unsigned int j=TV3._nord()+1-i; j<=TV3._nord()+1; j++ ){
      r1 += TV1._bndord[j];
      r2 += TV2._bndord[j];
    }
    s1 += TV2._bndord[i] * r1 ;
    s2 += TV1._bndord[i] * r2 ;
  }
  if( !Op<T>::inter( *(TV3._bndrem), s1, s2) ){
    *(TV3._bndrem) = s1;
    //throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::SQUARE );
  }
  // Populate _bndord for product term (except remainder term)
  TV3._update_bndord();
  if( TV3._TM->options.PROPAGATE_BNDT ) TV3._bndT = TV1._bndT * TV2._bndT;
  if( TV3._TM->options.CENTER_REMAINDER ) TV3._center_TM();
  return TV3;
}

template <typename T> inline TVar<T>
sqr
( const TVar<T>&TV )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] *= TV2._coefmon[0];
    *(TV2._bndrem) *= 2. + *(TV2._bndrem);
    return TV2;
 }

  // Populate _coefmon for product term
  TVar<T> TV2( TV._TM, 0. );
  for( unsigned int i=0; i<TV2._posord(TV2._nord()/2+1); i++){
    TV2._coefmon[TV2._prodmon(i,i+1)] += TV._coefmon[i] * TV._coefmon[i];
    for( unsigned int j=i+1; j<TV2._prodmon(i,0); j++ )
      TV2._coefmon[TV2._prodmon(i,j+1)] += TV._coefmon[i] * TV._coefmon[j] * 2.;
  }

  T s = 0.;
  for( unsigned int i=0; i<=TV2._nord()+1; i++ ){
    unsigned int k = std::max(TV2._nord()+1-i, i+1);
    T r = 0.;
    for( unsigned int j=k; j<=TV2._nord()+1; j++ )
      r += TV._bndord[j];
    s += TV._bndord[i] * r;
  }

  T r = 0.;
  for( unsigned int i=TV2._nord()/2+1; i<=TV2._nord()+1; i++ )
    r += Op<T>::sqr(TV._bndord[i]) ;
  *(TV2._bndrem) = 2. * s + r;
  
  // Populate _bndord for product term (except remainder term)
  TV2._update_bndord();
  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::sqr( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>&
TVar<T>::operator *=
( const double c )
{
  if( !_TM ){
    _coefmon[0] *= c;
    *(_bndrem) *= c;
  }
  else{
    for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] *= c;
    for( unsigned int i=0; i<_nord()+2; i++ ) _bndord[i] *= c;
    if( _TM->options.PROPAGATE_BNDT ) _bndT *= c;
  }
  return *this;
}

template <typename T> inline TVar<T>
operator *
( const TVar<T>&TV1, const double c )
{
  TVar<T> TV3( TV1 );
  TV3 *= c;
  return TV3;
}

template <typename T> inline TVar<T>
operator *
( const double c, const TVar<T>&TV2 )
{
  TVar<T> TV3( TV2 );
  TV3 *= c;
  return TV3;
}

template <typename T> inline TVar<T>&
TVar<T>::operator *=
( const T&I )
{
  if( !_TM ){
    *(_bndrem) += _coefmon[0];
    _coefmon[0] = 0.;
    *(_bndrem) *= I;
  }
  else{
    const double Imid = Op<T>::mid(I);
    T Icur = bound();
    for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] *= Imid;
    for( unsigned int i=0; i<_nord()+2; i++ ) _bndord[i] *= Imid;
    *_bndrem += (I-Imid)*Icur;
  }
  if( _TM && _TM->options.CENTER_REMAINDER ) _center_TM();
  if( _TM && _TM->options.PROPAGATE_BNDT ) _bndT *= I;
  return (*this);
}

template <typename T> inline TVar<T>
operator *
( const TVar<T>&TV1, const T&I )
{
  TVar<T> TV3( TV1 );
  TV3 *= I;
  return TV3;
}

template <typename T> inline TVar<T>
operator *
( const T&I, const TVar<T>&TV2 )
{
  TVar<T> TV3( TV2 );
  TV3 *= I;
  return TV3;
}

template <typename T> inline TVar<T>&
TVar<T>::operator /=
( const TVar<T>&TV )
{
   *this *= inv(TV);
   return *this;
}

template <typename T> inline TVar<T>
operator /
( const TVar<T>&TV1, const TVar<T>&TV2 )
{
  return TV1 * inv(TV2);
}

template <typename T> inline TVar<T>&
TVar<T>::operator /=
( const double c )
{
  if ( isequal( c, 0. ))
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::DIV );
   *this *= (1./c);
   return *this;
}

template <typename T> inline TVar<T>
operator /
( const TVar<T>&TV, const double c )
{
  if ( isequal( c, 0. ))
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::DIV );
  return TV * (1./c);
}

template <typename T> inline TVar<T>
operator /
( const double c, const TVar<T>&TV )
{
  return inv(TV) * c;
}

template <typename T> inline TVar<T>
inv
( const TVar<T>&TV )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::inv(TV._coefmon[0] + *(TV._bndrem));
    TV2._update_bndord();
    return TV2;
  }

  TVar<T> TV2 = TV._TM->_inv_taylor( TV );
  if( TV2._TM->options.BERNSTEIN_USE ){
    TVar<T> TV2b = TV._TM->_inv_bernstein( TV );
    if( Op<T>::diam(*(TV2b._bndrem)) < Op<T>::diam(*(TV2._bndrem)) )
      TV2 = TV2b;
  }

  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::inv( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_inv_taylor
( const TVar<T>&TV )
{
  const T I( TV.B() );
  double x0 = ( TV._TM->options.REF_MIDPOINT? Op<T>::mid(I):
                TV._coefmon[0] );
  const TVar<T> TVmx0( TV - x0 );
  const T Imx0( I - x0 );

  TVar<T> TV2( TV._TM, 1. ), MON( 1. );
  for( unsigned int i=1; i<=TV._nord(); i++ ){
    MON *= TVmx0 / (-x0);
    TV2 += MON;
  }
  TV2 /= x0;
  TV2 += Op<T>::pow( -Imx0, (int)TV2._nord()+1 )
       / Op<T>::pow( Op<T>::zeroone()*Imx0+x0, (int)TV2._nord()+2 );
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_inv_bernstein
( const TVar<T>&TV )
{
  struct loc{
    static double inv
      ( const double x, const double*rusr, const int*iusr )
      { return 1./x; }
    static double dinv
      ( const double x, const double*rusr, const int*iusr )
      { return -1./mc::sqr(x); }
    static T Iinv
      ( const T&x, const double*rusr, const int*iusr )
      { return Op<T>::inv(x); }
    static T d2Iinv
      ( const T&x, const double*rusr, const int*iusr )
      { return 2.*Op<T>::pow( Op<T>::inv(x), 3 ); }
  };
  return TV._TM->_univ_bernstein( TV, loc::inv, loc::dinv, loc::Iinv,
    loc::d2Iinv, 0, 0 );
}

template <typename T> inline TVar<T>
sqrt
( const TVar<T>&TV )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::sqrt(TV._coefmon[0] + *(TV._bndrem));
    TV2._update_bndord();
    return TV2;
  }

  TVar<T> TV2 = TV._TM->_sqrt_taylor( TV );
  if( TV2._TM->options.BERNSTEIN_USE ){
    TVar<T> TV2b = TV._TM->_sqrt_bernstein( TV );
    if( Op<T>::diam(*(TV2b._bndrem)) < Op<T>::diam(*(TV2._bndrem)) )
      TV2 = TV2b;
  }

  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::sqrt( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_sqrt_taylor
( const TVar<T>&TV )
{
  const T I( TV.B() );
  double x0 = ( TV._TM->options.REF_MIDPOINT? Op<T>::mid(I):
                TV._coefmon[0] );
  const TVar<T> TVmx0( TV - x0 );
  const T Imx0( I - x0 );

  double s = 0.5;
  TVar<T> TV2( TV._TM, 1. ), MON( 1. );
  for( unsigned int i=1; i<=TV._nord(); i++ ){
    MON *= TVmx0 / x0;
    TV2 += MON * s;
    s *= -(2.*i-1.)/(2.*i+2.);
  }
  TV2 *= std::sqrt(x0);
  TV2 += s * Op<T>::pow( Imx0, (int)TV2._nord()+1 )
           / Op<T>::pow( Op<T>::zeroone()*Imx0+x0, (int)TV2._nord()+1/2 );
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_sqrt_bernstein
( const TVar<T>&TV )
{
  struct loc{
    static double sqrt
      ( const double x, const double*rusr, const int*iusr )
      { return std::sqrt(x); }
    static double dsqrt
      ( const double x, const double*rusr, const int*iusr )
      { return 1./(2.*std::sqrt(x)); }
    static T Isqrt
      ( const T&x, const double*rusr, const int*iusr )
      { return Op<T>::sqrt(x); }
    static T d2Isqrt
      ( const T&x, const double*rusr, const int*iusr )
      { return -Op<T>::inv( 4.*Op<T>::sqrt(Op<T>::pow(x,3)) ); }
  };
  return TV._TM->_univ_bernstein( TV, loc::sqrt, loc::dsqrt, loc::Isqrt,
    loc::d2Isqrt, 0, 0 );
}

template <typename T> inline TVar<T>
exp
( const TVar<T>&TV )
{ 
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::exp(TV._coefmon[0] + *(TV._bndrem));
    TV2._update_bndord();
    return TV2;
  }

  TVar<T> TV2 = TV._TM->_exp_taylor( TV );
  if( TV2._TM->options.BERNSTEIN_USE ){
    TVar<T> TV2b = TV._TM->_exp_bernstein( TV );
    //std::cout << *(TV2b._bndrem) << " <=? " << *(TV2._bndrem) << std::endl;
    if( Op<T>::diam(*(TV2b._bndrem)) < Op<T>::diam(*(TV2._bndrem)) )
      TV2 = TV2b;
  }

  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::exp( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_exp_taylor
( const TVar<T>&TV )
{
  const T I( TV.B() );
  double x0 = ( TV._TM->options.REF_MIDPOINT? Op<T>::mid(I):
                TV._coefmon[0] );
  const TVar<T> TVmx0( TV - x0 );
  const T Imx0( I - x0 );

  double s = 1.;
  TVar<T> TV2( TV._TM, 1. ), MON( 1. );
  for( unsigned int i=1; i<=TV._nord(); i++ ){
    MON *= TVmx0;
    TV2 += MON * s;
    s /= i+1.;
  }
  TV2 += s * Op<T>::pow( Imx0, (int)TV2._nord()+1 )
           * Op<T>::exp( Op<T>::zeroone()*Imx0 );
  TV2 *= std::exp(x0);
  TV2._update_bndord();
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_exp_bernstein
( const TVar<T>&TV )
{
  struct loc{
    static double exp
      ( const double x, const double*rusr, const int*iusr )
      { return std::exp(x); }
    static double dexp
      ( const double x, const double*rusr, const int*iusr )
      { return std::exp(x); }
    static T Iexp
      ( const T&x, const double*rusr, const int*iusr )
      { return Op<T>::exp(x); }
    static T d2Iexp
      ( const T&x, const double*rusr, const int*iusr )
      { return Op<T>::exp(x); }
  };
  return TV._TM->_univ_bernstein( TV, loc::exp, loc::dexp, loc::Iexp,
    loc::d2Iexp, 0, 0 );
}

template <typename T> inline TVar<T>
log
( const TVar<T>&TV )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::log(TV._coefmon[0] + *(TV._bndrem));
    TV2._update_bndord();
    return TV2;
  }

  TVar<T> TV2 = TV._TM->_log_taylor( TV );
  if( TV2._TM->options.BERNSTEIN_USE ){
    TVar<T> TV2b = TV._TM->_log_bernstein( TV );
    if( Op<T>::diam(*(TV2b._bndrem)) < Op<T>::diam(*(TV2._bndrem)) )
      TV2 = TV2b;
  }

  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::log( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_log_taylor
( const TVar<T>&TV )
{
  const T I( TV.B() );
  double x0 = ( TV._TM->options.REF_MIDPOINT? Op<T>::mid(I):
                TV._coefmon[0] );
  const TVar<T> TVmx0( TV - x0 );
  const T Imx0( I - x0 );

  TVar<T> TV2( TV._TM, 0. ), MON( -1. );
  for( unsigned int i=1; i<=TV._nord(); i++ ){
    MON *= TVmx0 / (-x0);
    TV2 += MON / (double)i;
  }
  TV2._coefmon[0] += std::log(x0);
  TV2._update_bndord();
  TV2 -= Op<T>::pow( - Imx0 / ( Op<T>::zeroone()*Imx0+x0 ),
       (int)TV2._nord()+1 ) / ( TV2._nord()+1. );
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_log_bernstein
( const TVar<T>&TV )
{
  struct loc{
    static double log
      ( const double x, const double*rusr, const int*iusr )
      { return std::log(x); }
    static double dlog
      ( const double x, const double*rusr, const int*iusr )
      { return 1./x; }
    static T Ilog
      ( const T&x, const double*rusr, const int*iusr )
      { return Op<T>::log(x); }
    static T d2Ilog
      ( const T&x, const double*rusr, const int*iusr )
      { return -Op<T>::sqr( Op<T>::inv(x) ); }
  };
  return TV._TM->_univ_bernstein( TV, loc::log, loc::dlog, loc::Ilog,
    loc::d2Ilog, 0, 0 );
}

template <typename T> inline TVar<T>
xlog
( const TVar<T>&TV )
{
  return TV * log( TV );
}

template <typename T> inline TVar<T>
pow
( const TVar<T>&TV, const int n )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::pow(TV._coefmon[0] + *(TV._bndrem), n);
    TV2._update_bndord();
    return TV2;
  }

  if( n < 0 ) return pow( inv( TV ), -n );
  TVar<T> TV2( TV._TM->_intpow( TV, n ) );
  if( TV2._TM && TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::pow( TV._bndT, n );
  if( TV2._TM && TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
TModel<T>::_intpow
( const TVar<T>&TV, const int n )
{
  if( n == 0 ) return 1.;
  else if( n == 1 ) return TV;
  return n%2 ? sqr( _intpow( TV, n/2 ) ) * TV : sqr( _intpow( TV, n/2 ) );
}

template <typename T> inline TVar<T>
pow
( const TVar<T> &TV, const double a )
{
  return exp( a * log( TV ) );
}

template <typename T> inline TVar<T>
pow
( const TVar<T> &TV1, const TVar<T> &TV2 )
{
  return exp( TV2 * log( TV1 ) );
}

template <typename T> inline TVar<T>
pow
( const double a, const TVar<T> &TV )
{
  return exp( TV * std::log( a ) );
}

template <typename T> inline TVar<T>
monomial
(const unsigned int n, const TVar<T>*TV, const int*k)
{
  if( n == 0 ){
    return 1.;
  }
  if( n == 1 ){
    return pow( TV[0], k[0] );
  }
  return pow( TV[0], k[0] ) * monomial( n-1, TV+1, k+1 );
}

template <typename T> inline TVar<T>
cos
( const TVar<T> &TV )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::cos(TV._coefmon[0] + *(TV._bndrem));
    TV2._update_bndord();
    return TV2;
  }

  const T I( TV.B() );
  double x0 = ( TV._TM->options.REF_MIDPOINT? Op<T>::mid(I):
                TV._coefmon[0] );
  const TVar<T> TVmx0( TV - x0 );
  const T Imx0( I - x0 );
  double s = 1., c;

  TVar<T> TV2( TV._TM, 0. ), MON( 1. );
  for( unsigned int i=1; i<=TV._nord(); i++ ){
    switch( i%4 ){
    case 0: c =  std::cos(x0); break;
    case 1: c = -std::sin(x0); break;
    case 2: c = -std::cos(x0); break;
    case 3:
    default: c =  std::sin(x0); break;
    }
    MON *= TVmx0;
    TV2 += c * s * MON;
    s /= i+1;
  }

  TV2._coefmon[0] += std::cos(x0);
  TV2._update_bndord();
  switch( (TV2._nord()+1)%4 ){
  case 0: TV2 += s * Op<T>::pow( Imx0, (int)TV2._nord()+1 )
                   * Op<T>::cos( Op<T>::zeroone()*Imx0+x0 ); break;
  case 1: TV2 -= s * Op<T>::pow( Imx0, (int)TV2._nord()+1 )
                   * Op<T>::sin( Op<T>::zeroone()*Imx0+x0 ); break;
  case 2: TV2 -= s * Op<T>::pow( Imx0, (int)TV2._nord()+1 )
                   * Op<T>::cos( Op<T>::zeroone()*Imx0+x0 ); break;
  case 3: TV2 += s * Op<T>::pow( Imx0, (int)TV2._nord()+1 )
                   * Op<T>::sin( Op<T>::zeroone()*Imx0+x0 ); break;
  }

  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::cos( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
sin
( const TVar<T> &TV )
{
  return cos( TV - PI/2. );
}

template <typename T> inline TVar<T>
asin
( const TVar<T> &TV )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::asin(TV._coefmon[0] + *(TV._bndrem));
    TV2._update_bndord();
    return TV2;
  }

  double x0 = ( TV._TM->options.REF_MIDPOINT? Op<T>::mid(TV.B()):
                TV._coefmon[0] );
  double s = 1., t = 1.;
  TVar<T> G = TV * std::sqrt(1-x0*x0) - x0 * sqrt(1.-sqr(TV));
  TVar<T> TV2( G ), MON( G ), sqrG( sqr(G) );
  T IG( G.B() ), IG0( IG*Op<T>::zeroone() ), sqrIG0( Op<T>::sqr(IG0) ) ;
  T ASIN1, ASIN2( 1./Op<T>::sqrt(1-sqrIG0) ),
    ASIN3( IG0*Op<T>::pow(Op<T>::sqrt(1-sqrIG0),-3) );
  for( unsigned int i=1; i<=TV._nord(); i++ ){
    MON *= sqrG;
    s *= (double)(2*i-1)*(double)(2*i-1)/(double)(2*i)/(double)(2*i+1);
    TV2 += MON * s;
    t *= double(i+1);
    ASIN1 = ASIN2; ASIN2 = ASIN3;
    ASIN3 = ((2*i-1)*IG0*ASIN2+(i-1)*(i-1)*ASIN1)/(1-sqrIG0);
  }
  TV2._coefmon[0] += std::asin(x0);
  TV2._update_bndord();
  TV2 += Op<T>::pow( IG, (int)TV2._nord()+1 ) / t * ASIN2;
  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::asin( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
acos
( const TVar<T> &TV )
{
  return PI/2. - asin( TV );
}

template <typename T> inline TVar<T>
tan
( const TVar<T> &TV )
{
  return sin(TV) / cos(TV);
}

template <typename T> inline TVar<T>
atan
( const TVar<T> &TV )
{
  if( !TV._TM ){
    TVar<T> TV2( TV );
    TV2._coefmon[0] = 0.;
    *(TV2._bndrem) = Op<T>::atan(TV._coefmon[0] + *(TV._bndrem));
    TV2._update_bndord();
    return TV2;
  }

  double x0 = ( TV._TM->options.REF_MIDPOINT? Op<T>::mid(TV.B()):
                TV._coefmon[0] );
  TVar<T> G = ( TV - x0 ) / ( 1. + x0 * TV );
  TVar<T> TV2( G ), MON( G ), msqrG( -sqr(G) );
  T IG( G.B() ), IG0( IG*Op<T>::zeroone() );
  for( unsigned int i=1; i<=TV._nord(); i++ ){
    MON *= msqrG;
    TV2 += MON / (2*i+1);
  }
  TV2._coefmon[0] += std::atan(x0);
  TV2._update_bndord();
  TV2 += Op<T>::pow( IG * Op<T>::cos( Op<T>::atan(IG0) ),
         (int)TV2._nord()+1 ) / (double)(TV2._nord()+1)
         * sin( (TV2._nord()+1) * (atan(IG0)+PI/2.) );
  if( TV2._TM->options.PROPAGATE_BNDT ) TV2._bndT = Op<T>::asin( TV._bndT );
  if( TV2._TM->options.CENTER_REMAINDER ) TV2._center_TM();
  return TV2;
}

template <typename T> inline TVar<T>
hull
( const TVar<T>&TV1, const TVar<T>&TV2 )
{
  // Neither operands associated to TModel -- Make intersection in T type     
  if( !TV1._TM && !TV2._TM ){
    T R1 = TV1._coefmon[0] + TV1._bndord[0];
    T R2 = TV2._coefmon[0] + TV2._bndord[0];
    return Op<T>::hull(R1, R2);
  }

  // First operand not associated to TModel
  else if( !TV1._TM )
    return hull( TV2, TV1 );

  // Second operand not associated to TModel
  else if( !TV2._TM )
    return TV1.P() + Op<T>::hull( TV1.R(), TV2._coefmon[0]+TV2._bndord[0]-TV1.B() );

  // TModel for first and second operands are inconsistent
  else if( TV1._TM != TV2._TM )
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::TMODEL );

  // Perform union
  TVar<T> TV1C( TV1 ), TV2C( TV2 );
  const double eta = TV1._TM->options.REF_POLY;
  T R1C = TV1C.C().R(), R2C = TV2C.C().R(); 
  TV1C.set(T(0.));
  TV2C.set(T(0.));
  T BTVD = (TV1C-TV2C).B();
  return (1.-eta)*TV1C + eta*TV2C + Op<T>::hull( R1C+eta*BTVD, R2C+(eta-1.)*BTVD );
}

template <typename T> inline bool
inter
( TVar<T>&TVR, const TVar<T>&TV1, const TVar<T>&TV2 )
{
  // Neither operands associated to TModel -- Make intersection in T type     
  if( !TV1._TM && !TV2._TM ){
    T R1 = TV1._coefmon[0] + TV1._bndord[0];
    T R2 = TV2._coefmon[0] + TV2._bndord[0];
    T RR( 0. );
    bool flag = Op<T>::inter(RR, R1, R2);
    TVR = RR;
    return flag;
  }

  // First operand not associated to TModel
  else if( !TV1._TM )
    return inter( TVR, TV2, TV1 );

  // Second operand not associated to TModel
  else if( !TV2._TM ){
    TVR = TV1.P();
    return( Op<T>::inter(*(TVR._bndrem), TV1.R(),
      TV2._coefmon[0]+TV2._bndord[0]-TV1.B())? true: false );
  }

  // TModel for first and second operands are inconsistent
  else if( TV1._TM != TV2._TM )
    throw typename TModel<T>::Exceptions( TModel<T>::Exceptions::TMODEL );

  // Perform intersection
  TVar<T> TV1C( TV1 ), TV2C( TV2 );
  const double eta = TV1._TM->options.REF_POLY;
  T R1C = TV1C.C().R(), R2C = TV2C.C().R(); 
  TV1C.set(T(0.));
  TV2C.set(T(0.));
  TVR = (1.-eta)*TV1C + eta*TV2C;
  TV1C -= TV2C;
  T BTVD = TV1C.B();
  return( Op<T>::inter( *(TVR._bndrem), R1C+eta*BTVD, R2C+(eta-1.)*BTVD )?
    true: false );
}

} // namespace mc

#include "mcop.hpp"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure to allow usage of the Taylor model type mc::TVar inside other MC++ type, e.g. mc::McCormick
template <> template<typename T> struct Op< mc::TVar<T> >
{
  typedef mc::TVar<T> TV;
  static TV point( const double c ) { return TV(c); }
  static TV zeroone() { return TV( mc::Op<T>::zeroone() ); }
  static void I(TV& x, const TV&y) { x = y; }
  static double l(const TV& x) { return mc::Op<T>::l(x.B()); }
  static double u(const TV& x) { return mc::Op<T>::u(x.B()); }
  static double abs (const TV& x) { return mc::Op<T>::abs(x.B());  }
  static double mid (const TV& x) { return mc::Op<T>::mid(x.B());  }
  static double diam(const TV& x) { return mc::Op<T>::diam(x.B()); }
  static TV inv (const TV& x) { return mc::inv(x);  }
  static TV sqr (const TV& x) { return mc::sqr(x);  }
  static TV sqrt(const TV& x) { return mc::sqrt(x); }
  static TV log (const TV& x) { return mc::log(x);  }
  static TV xlog(const TV& x) { return x*mc::log(x); }
  static TV fabs(const TV& x) { return TV( mc::Op<T>::fabs(x.B()) ); }
  static TV exp (const TV& x) { return mc::exp(x);  }
  static TV sin (const TV& x) { return mc::sin(x);  }
  static TV cos (const TV& x) { return mc::cos(x);  }
  static TV tan (const TV& x) { return mc::tan(x);  }
  static TV asin(const TV& x) { return mc::asin(x); }
  static TV acos(const TV& x) { return mc::acos(x); }
  static TV atan(const TV& x) { return mc::atan(x); }
  static TV erf (const TV& x) { throw typename mc::TModel<T>::Exceptions( TModel<T>::Exceptions::UNDEF ); }
  static TV erfc(const TV& x) { throw typename mc::TModel<T>::Exceptions( TModel<T>::Exceptions::UNDEF ); }
  static TV hull(const TV& x, const TV& y) { return mc::hull(x,y); }
  static TV min (const TV& x, const TV& y) { return mc::Op<T>::min(x.B(),y.B());  }
  static TV max (const TV& x, const TV& y) { return mc::Op<T>::max(x.B(),y.B());  }
  static TV arh (const TV& x, const double k) { return mc::exp(-k/x); }
  template <typename X, typename Y> static TV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static TV monomial (const unsigned int n, const T* x, const int* k) { return mc::monomial(n,x,k); }
  static bool inter(TV& xIy, const TV& x, const TV& y) { return mc::inter(xIy,x,y); }
  static bool eq(const TV& x, const TV& y) { return mc::Op<T>::eq(x.B(),y.B()); }
  static bool ne(const TV& x, const TV& y) { return mc::Op<T>::ne(x.B(),y.B()); }
  static bool lt(const TV& x, const TV& y) { return mc::Op<T>::lt(x.B(),y.B()); }
  static bool le(const TV& x, const TV& y) { return mc::Op<T>::le(x.B(),y.B()); }
  static bool gt(const TV& x, const TV& y) { return mc::Op<T>::gt(x.B(),y.B()); }
  static bool ge(const TV& x, const TV& y) { return mc::Op<T>::ge(x.B(),y.B()); }
};

} // namespace mc

#endif
