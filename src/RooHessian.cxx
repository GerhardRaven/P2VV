#include <memory>

#include <P2VV/RooHessian.h>
#include <RooFunctor.h>

///Compute the numerical Hessian and errors. The Hessian is returned as the first
///element, and the errors as the second.
///See numerical_hessian().
///@param f Functor to double-differentiate
///@param x Point about which to double-differentiate.
///@ingroup gFunctions 
TMatrixTSym<double> hessian(const RooAbsReal& f, const RooArgList& params) {

   using namespace NumDiv;
   TMatrixTSym<double> hess(params.getSize());

   std::auto_ptr<RooArgSet> func_params(f.getObservables(RooArgSet(params)));

   //Perform the cross differencing.
   for(int r = 0; r < params.getSize(); r++) {
			for(int c = r + 1; c < params.getSize(); c++) {
         RooAbsReal* x = static_cast<RooAbsReal*>(func_params->find(*params.at(r)));
         RooAbsReal* y = static_cast<RooAbsReal*>(func_params->find(*params.at(c)));
         RooArgList args(*x, *y);
         std::auto_ptr<RooFunctor> func(f.functor(RooArgList(), args));
         CentralCrossDifferenceSecond<RooFunctor> cross(x->getVal(), y->getVal(), *func);
         pair<double, double> e = extrapolate_to_zero<CentralCrossDifferenceSecond<RooFunctor> >(cross);
         hess[r][c] = hess[c][r] = e.first;
			}
   }

   for (int i = 0; i < params.getSize(); i++) {
      RooAbsReal* x = static_cast<RooAbsReal*>(func_params->find(*params.at(i)));
      std::auto_ptr<RooFunctor> func(f.functor(RooArgList(), RooArgList(*x)));
      CentralDifferenceSecond<RooFunctor> curv(x->getVal(), *func);
			pair<double, double> e = extrapolate_to_zero<CentralDifferenceSecond<RooFunctor> >(curv);
			hess[i][i] = e.first;
   }

   return hess;
}

