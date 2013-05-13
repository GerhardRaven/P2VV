/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBSpline.cxx 45780 2012-08-31 15:45:27Z moneta $
 * Authors:                                                                  *
 *   Gerhard Raven
 *                                                                           *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// BSpline basis polynomials are positive-definite in the range [0,1].
// In this implementation, we extend [0,1] to be the range of the parameter.
// There are n+1 BSpline basis polynomials of degree n.
// Thus, by providing N coefficients that are positive-definite, there 
// is a natural way to have well bahaved polynomail PDFs.
// For any n, the n+1 basis polynomials 'form a partition of unity', eg.
//  they sum to one for all values of x. See
// END_HTML
//
// STD & STL
#include <cmath>
#include <iterator>
#include <algorithm>
#include <sstream>

// ROOT
#include "TMath.h"
#include "TH1.h"

// RooFit
#include "RooFit.h"
#include "Riostream.h"
#include "RooMath.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooComplex.h"
#include "RooArgList.h"

// P2VV
#include "P2VV/RooCubicSplineFun.h"
#include "P2VV/RooCubicSplineKnot.h"

using namespace std;

ClassImp(RooCubicSplineFun);

//_____________________________________________________________________________
void RooCubicSplineFun::init(const char* name, const std::vector<double>& knots,
                             const std::vector<double>& heights,
                             const std::vector<double>& errors,
                             double smooth, bool constCoeffs) {
   _aux = new RooCubicSplineKnot( knots.begin(), knots.end() );
   std::vector<double> values(heights);
   if ( smooth > 0 ) { 
      assert(errors.size() == knots.size());
      _aux->smooth( values, errors, smooth );
   }
   _aux->computeCoefficients( values );
   for (unsigned int i=0;i<values.size();++i) { 
      if (constCoeffs) {
         _coefList.add( RooFit::RooConst( values[i] ) );
      } else {
         stringstream name_str;
         name_str << name << "_smoothed_bin_" << i;
         string n(name_str.str());
         RooRealVar* coeff = new RooRealVar(n.c_str(), n.c_str(), values[i], 0.0001, 0.9999);
         if (i == 0 || i == values.size() - 1) coeff->setConstant(true);
         _coefList.add(*coeff);
         _ownList.addOwned(*coeff);
      }
   }
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun()
    : _aux(0)
{
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun(const char* name, const char* title, RooRealVar& x,
                                     const std::vector<double>& knots,
                                     const std::vector<double>& values,
                                     const std::vector<double>& errors,
                                     double smooth, bool constCoeffs) :
   RooAbsReal(name, title),
   _x("x", "Dependent", this, x),
   _coefList("coefficients","List of coefficients",this),
   _ownList("ownList", "List of owned RealVars", this),
   _aux(0)
{
   init(name, knots, values, errors, smooth, constCoeffs);
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun(const char* name, const char* title, 
                                     RooRealVar& x, const TH1* hist, double smooth,
                                     bool constCoeffs) :
  RooAbsReal(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this),
  _ownList("ownList", "List of owned RealVars", this),
  _aux(0)
{
    // bin 0 is underflow, and bin nBins + 1 is overflow...
    int nBins = hist->GetNbinsX();
    std::vector<double> centres, values;
    for (int i=0;i<nBins ;++i) {
        centres.push_back(hist->GetBinCenter(1+i));
        values.push_back(hist->GetBinContent(1+i));
    }

    std::vector<double> errs;
    if (smooth>0) for (int i=0;i<nBins ;++i) errs.push_back(hist->GetBinError(1+i));
    
    init(name, centres, values, errs, smooth, constCoeffs);
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun(const char* name, const char* title, 
                                     RooRealVar& x, const char* knotBinningName,
                                     const RooArgList& coefList): 
  RooAbsReal(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients", "List of coefficients", this),
  _ownList("ownList", "List of owned RealVars", this),
  _aux(0)
{
  // TODO: verify coefList is consistent with knots as specified by the knotBinningName binning
  //    should be N+2 coefficients for N bins...
  const RooAbsBinning* binning = x.getBinningPtr(knotBinningName);
  assert( binning!=0);
  assert( coefList.getSize()==2+binning->numBoundaries());
  _coefList.add(coefList);

  Double_t* boundaries = binning->array();
  _aux = new RooCubicSplineKnot( boundaries, boundaries + binning->numBoundaries() );
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun(const char* name, const char* title, 
                                     RooRealVar& x, const std::vector<double>& knots,
                                     const RooArgList& coefList): 
  RooAbsReal(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients", "List of coefficients", this),
  _ownList("ownList", "list of owned RealVars", this),
  _aux(0)
{
   assert(size_t(coefList.getSize()) == 2 + knots.size());
   _coefList.add(coefList);
   _aux = new RooCubicSplineKnot(knots.begin(), knots.end());
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun(const RooCubicSplineFun& other, const char* name) :
  RooAbsReal(other, name), 
  _x("x", this, other._x), 
  _coefList("coefList", "List of owned RealVars", this),
  _ownList("ownList", "List of owned RealVars", this),
  _aux(new RooCubicSplineKnot(*other._aux))
{
   if (_ownList.getSize()) {
      RooAbsCollection* copy = other._ownList.snapshot();
      _ownList.addOwned(*copy);
      _coefList.add(*copy);
   } else {
      _coefList.add(other._coefList);
   }
}

//_____________________________________________________________________________
RooCubicSplineFun::~RooCubicSplineFun()
{
    delete _aux;
    _ownList.removeAll();
}

//_____________________________________________________________________________
Double_t RooCubicSplineFun::evaluate() const 
{
  return _aux->evaluate(_x,_coefList);
}

//_____________________________________________________________________________
Int_t RooCubicSplineFun::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  // No analytical calculation available (yet) of integrals over subranges
  if (rangeName && strlen(rangeName)) return 0 ;
  if (matchArgs(allVars, analVars, _x)) return 1;
  return 0;
}

RooComplex  RooCubicSplineFun::gaussIntegral(int i, const RooCubicSplineGaussModel::M_n& dM,
                                             const RooCubicSplineGaussModel::K_n& K, double offset,
                                             double* sc) const 
{
        RooComplex sum(0,0);
        RooCubicSplineKnot::S_jk S( _aux->S_jk_sum( i, _coefList ), offset );
        for (int j=0;j<4;++j) for (int k=0;k<4-j;++k) {
            sum = sum + dM(j)*S(j,k)*K(k)*sc[j+k];
        }
        return sum;
}

//_____________________________________________________________________________
Double_t RooCubicSplineFun::analyticalIntegral(Int_t code, const char* /* rangeName */) const
{
  assert(code==1) ;
  return _aux->analyticalIntegral(_coefList);
}

