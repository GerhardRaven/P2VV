#ifndef DERIVATIVES_NUMERICAL_H
#define DERIVATIVES_NUMERICAL_H

/* Shamelessy stolen with many thanks to the TooN library:
   http://www.edwardrosten.com/cvd/toon.html
*/

#include <limits>
#include <vector>
#include <cmath>

#include <TMatrixTSym.h>
#include <RooAbsReal.h>
#include <RooArgList.h>

using namespace std;

namespace NumDiv {
   ///@internal
   ///@brief Implementation of Ridder's Extrapolation to zero.
   ///The function input starts from 1.
   ///@param f Functor to extrapolate to zero.
   template<class F> std::pair<double, double> extrapolate_to_zero(F& f)  {
			using std::isfinite;
			using std::max;
			using std::isnan;
			//Use points at c^0, c^-1, ... to extrapolate towards zero.
			const double c = 1.1, t = 2;

			static const int Npts = 400;

			/* Neville's table is:	
         x_0      y_0    P_0                                 
         P_{01}                        
         x_1      y_1    P_1            P_{012}              
         P_{12}             P_{0123}   
         x_2      y_2    P_2            P_{123}              
         P_{23}                        
         x_3      y_3    P_3   

         In Matrix form, this is rearranged as:
			
			
         x_0      y_0    P_0                                 
														  
         x_1      y_1    P_1   P_{01}                 
											   
         x_2      y_2    P_2   P_{12}  P_{012}               
																	   
         x_3      y_3    P_3   P_{23}  P_{123} P_{0123} 

         This is rewritten further as:

         |                  0      1      2       3
         _____|___________________________________________________	
         0   | x_0      y_0    P^00                                
         | 	                                          
         1   | x_1      y_1    P^10  P^11                   
         | 	                               
         2   | x_2      y_2    P^10  P^21    P^22                  
         | 	                                                       
         3   | x_3      y_3    P^10  P^31    P^23    P^33     

         So, P^ij == P_{j-i...j}

         Finally, only the current and previous row are required. This reduces
         the memory cost from quadratic to linear.  The naming scheme is:

         P[i][j]    -> P_i[j]
         P[i-1][j]  -> P_i_1[j]
			*/

			vector<double> P_i(1), P_i_1;
			
			double best_err   = numeric_limits<double>::max();
			double best_point = numeric_limits<double>::max();

      //The first tranche of points might be bad.
      //Don't quit while no good points have ever happened.
      bool ever_ok = 0;

      //Compute f(x) as x goes from 1 towards zero and extrapolate to 0
      double x = 1;
      for (int i = 0; i < Npts; i++) {
         swap(P_i, P_i_1);
         P_i.resize(i + 2);

         //P[i][0] = c^ -i;
         P_i[0] = f(x);

         x /= c;
        
         //Compute the extrapolations
         //cj id c^j
         double cj = 1;
         bool better = 0; //Did we get better?
         bool any_ok = 0;
         for (int j = 1; j <= i; j++)	{
            cj *= c;
            //We have (from above) n = i and k = n - j
            //n-k = i - (n - j) = i - i + j = j
            //Therefore c^(n-k) is just c^j

            P_i[j] = (cj * P_i[j-1] - P_i_1[j-1]) / (cj - 1);

            if(any_ok || isfinite(P_i[j]))
               ever_ok = 1;

            //Compute the difference between the current point (high order)
            //and the corresponding lower order point at the current step size
            double err1 = abs(P_i[j] - P_i[j-1]);

            //Compute the difference between two consecutive points at the
            //corresponding lower order point at the larger stepsize
            double err2 = abs(P_i[j] - P_i_1[j-1]);

            //The error is the larger of these.
            double err = max(err1, err2);

            if (err < best_err && isfinite(err)) {
               best_err = err;
               best_point = P_i[j];
               better=1;
            }
         }
				
         using namespace std;
         //If the highest order point got worse, or went off the rails, 
         //and some good points have been seen, then break.
         if(ever_ok && !better && i > 0 && (abs(P_i[i] - P_i_1[i-1]) > t * best_err|| isnan(P_i[i])))
            break;
			}

			return std::make_pair(best_point, best_err);
   }

   ///@internal
   ///@brief Functor wrapper for computing finite differences along an axis.
   template<class Functor> struct CentralDifferenceGradient {
   public:

			CentralDifferenceGradient(const double v_, const Functor& f_)
         : v(v_), f(f_)
			{
      }
			
      virtual ~CentralDifferenceGradient()
      {
      }

			///Compute central difference.
			double operator()(double hh) 
			{
         using std::max;
         using std::abs;

         //Make the step size be on the scale of the value.
         double h = hh * max(abs(v) * 1e-3, 1e-3);

         double f1 = f(v - h);
         double f2 = f(v + h);

         return (f2 - f1) / (2 * h);
			}
      
   private:

      const double v;         ///< Value at witch to get the derivative
			const Functor& f;       ///< Functor to evaluate

   };

   ///@internal
   ///@brief Functor wrapper for computing finite difference second derivatives along an axis.
   template<class Functor> class CentralDifferenceSecond {
   public:
			CentralDifferenceSecond(const double v_, const Functor& f_)
         : v(v_), f(f_), central(f(v))
			{}
			
			///Compute central difference.
			double operator()(double hh) 
			{
         using std::max;
         using std::abs;

         //Make the step size be on the scale of the value.
         double h = hh * max(abs(v) * 1e-3, 1e-3);

         double f1 = f(v - h);
         double f2 = f(v + h);

         double d =  (f2 - 2 * central + f1) / (h * h);
         return d;
			}

   private:

      const double v;                      ///< Point at which to get the derivative
			const Functor&  f;                   ///< Functor to evaluate
			const double central;                ///< Central point.

   };

   ///@internal
   ///@brief Functor wrapper for computing finite difference cross derivatives along a pair of axes.
   template<class Functor> class CentralCrossDifferenceSecond {
   public:

			CentralCrossDifferenceSecond(const double vx_, const double vy_, const Functor& f_)
         : vx(vx_), vy(vy_), f(f_)
			{}
			
			///Compute central difference.
			double operator()(double hh) 
			{
         using std::max;
         using std::abs;

         //Make the step size be on the scale of the value.
         double hx = hh * max(abs(vx) * 1e-3, 1e-3);
         double hy = hh * max(abs(vy) * 1e-3, 1e-3);

         x[0] = vx + hx;
         x[1] = vy + hy;
         double a = f(x);

         x[0] = vx - hx;
         x[1] = vy + hy;
         double b = f(x);

         x[0] = vx + hx;
         x[1] = vy - hy;
         double c = f(x);

         x[0] = vx - hx;
         x[1] = vy - hy;
         double d = f(x);

         return (a - b - c + d) / (4 * hx * hy);
			}

   private:

			double x[2];              ///< Local copy of v
      const double vx;          ///< Point in x at which to calculate the derivative
      const double vy;          ///< Point in y at which to calculate the derivative
			const Functor&  f;        ///< Functor to evaluate

   };
}


/** Extrapolate a derivative to zero using Ridder's Algorithm.

    Ridder's algorithm works by evaluating derivatives with smaller and smaller step
    sizes, fitting a polynomial and extrapolating its value to zero.
	
    This algorithm is generally more accurate and much more reliable, but much slower than
    using simple finite differences. It is robust to awkward functions and does not
    require careful tuning of the step size, furthermore it provides an estimate of the
    errors. This implementation has been tuned for accuracy instead of speed. For an
    estimate of the errors, see also numerical_gradient_with_errors(). In general it is
    useful to know the errors since some functions are remarkably hard to differentiate
    numerically.

    Neville's Algorithm can be used to find a point on a fitted polynomial at
    \f$h\f$. Taking
    some points \f$h_i\f$ and \f$ y_i = f(h_i)\f$, one can define a table of
    points on various polyomials:
    \f[
    \begin{array}{ccccccc}
    h_0    & y_0  & P_0 &        &         &          \\
    &      &     & P_{01} &         &          \\
    h_1    & y_1  & P_1 &        & P_{012} &          \\
    &      &     & P_{12} &         & P_{0123} \\
    h_2    & y_2  & P_2 &        & P_{123} &          \\
    &      &     & P_{23} &         &          \\
    h_3    & y_3  & P_3 &        &         &          \\
    \end{array}                            
    \f]
    where \f$P_{123}\f$ is the value of a polynomial fitted to datapoints 1, 2 and 3
    evaluated at \f$ h\f$. The initial values are simple to evaluate: \f$P_i = y_i =
    f(h_i)\f$. The remaining values are determined by the recurrance:
    \f{equation}
		P_{k\cdots n} = \frac{(h - h_n)P_{k \cdots n-1} + (h_k - h)P_{k+1 \cdots n}}{h_k - h_n}
    \f}

    For Ridder's algorithm, we define the extrapolation  point \f$ h=0 \f$ and the
    sequence of points to evaluate as \f$ h_i = c^{-i} \f$ where \f$ c \f$ is some
    constant a little greater than 1. 
    Substituting in to the above equation gives:
    \f[
		P_{k \cdots n} = \frac{c^{-k} P_{k+1 \cdots n} - P_{k \cdots n-1}}{c^{n - k} - 1}
    \f]

    To measure errors, when computing (for example) \f$P_{234}\f$, we compare the
    value to the lower order with the same step size \f$P_{34}\f$, and the lower
    order with a larger step size \f$P_{23}\f$. The error estimate is the largest of
    these. The  extrapolation with the smallest error is retained.

    Not only that, but since every value of P is used, every single subset of points
    is used for polynomial fitting. So, bad points for large and small \f$h\f$ do
    not spoil the answer.

    It is generally assumed that as \f$h\f$ decreases, the errors decrease, then
    increase again. Therefore, if the current step size did not yield any
    improvements on the best point so far, then we terminate when the highest order
    point is \f$ t \f$ times worse then the previous order point.

    The parameters are:
    - \f$ c = 1.1 \f$
    - \f$ t = 2 \f$
    - Max iterations = 400

    \section rRef  References                         

    - Ridders, C. J. F, 1982, Advances in Engineering Software, 4(2) 75--76
    - Press, Vetterling, Teukolsky, Flannery, Numerical, 1997, Numerical Recipies in C (2nd ed.), Chapter 5.7

    @param f Functor to differentiate
    @param x Point about which to differentiate.
    @ingroup gFunctions
*/

template<class F> vector<double> numerical_gradient(const F& f, const vector<double>& x)
{
   // FIXME, doesn't work as currently implemented
   using namespace NumDiv;
   vector<double> grad(x.size());

   for (int i = 0; i < x.size(); i++) {
      CentralDifferenceGradient<F> d(x[i], f);
			grad[i] = extrapolate_to_zero<CentralDifferenceGradient<F> >(d).first;
   }

   return grad;
}
	
///Compute numerical gradients with errors.
///See numerical_gradient().
///Gradients are returned in the first row of the returned matrix.
///Errors are returned in the second row. The errors have not been scaled, so they are
///in the same range as the gradients.
///@param f Functor to differentiate
///@param x Point about which to differentiate.
///@ingroup gFunctions 
template<class F> vector<pair<double, double> > numerical_gradient_with_errors(const F& f, const vector<double>& x) {
   // FIXME, doesn't work as currently implemented

   using namespace NumDiv;
   vector<pair<double, double> > g(x.size());

   for (int i = 0; i < x.size(); i++) {	
      CentralDifferenceGradient<F> d(x[i], f);

			pair<double, double> r = extrapolate_to_zero<CentralDifferenceGradient<F> >(d);
      g[i] = r;
   }

   return g;
}

	
///Compute the numerical Hessian using central differences and Ridder's method:
///\f[
/// \frac{\partial^2 f}{\partial x^2} \approx \frac{f(x-h) - 2f(x) + f(x+h)}{h^2}
///\f]
///\f[
/// \frac{\partial^2 f}{\partial x\partial y} \approx \frac{f(x+h, y+h) - f(x-h,y+h) - f(x+h, y-h) + f(x-h, y-h)}{4h^2}
///\f]
///See numerical_gradient().
///@param f Functor to double-differentiate
///@param x Point about which to double-differentiate.
///@ingroup gFunctions 
template<class F> pair<TMatrixTSym<double>, TMatrixTSym<double> > numerical_hessian_with_errors
(const F& f, const vector<double>& x)
{
   // FIXME, doesn't work as currently implemented
   using namespace NumDiv;
   TMatrixTSym<double> hess(x.size());
   TMatrixTSym<double> errors(x.size());

   //Perform the cross differencing.
   for (int r = 0; r < x.size(); r++) {
			for (int c = r + 1; c < x.size(); c++) {
         CentralCrossDifferenceSecond<F> cross(x[r], x[c], f);
         pair<double, double> e = extrapolate_to_zero<CentralCrossDifferenceSecond<F> >(cross);
         hess[r][c] = hess[c][r] = e.first;
         errors[r][c] = errors[c][r] = e.second;
			}
   }

   for (int i = 0; i < x.size(); i++) {
      CentralDifferenceSecond<F> curv(x[i], f);
			pair<double, double> e = extrapolate_to_zero<CentralDifferenceSecond<F> >(curv);
			hess[i][i] = e.first;
			errors[i][i] = e.second;
   }

   return make_pair(hess, errors);
}

TMatrixTSym<double> hessian(const RooAbsReal& f, const RooArgList& params);

#endif
