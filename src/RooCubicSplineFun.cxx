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
RooCubicSplineFun::RooCubicSplineFun()
    : _aux(0)
{
}

RooCubicSplineFun::RooCubicSplineFun(const char* name, const char* title, 
                                     RooRealVar& x, const TH1* hist, double smooth,
                                     bool constCoeffs) :
  RooAbsReal(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this),
  _aux(0)
{
    int nBins = hist->GetNbinsX();
    std::vector<double> boundaries,values;
    for (int i=0;i<nBins ;++i) {
        boundaries.push_back(hist->GetBinCenter(1+i));
        values.push_back(hist->GetBinContent(1+i));
    }
    _aux = new RooCubicSplineKnot( boundaries.begin(), boundaries.end() );
    if ( smooth > 0 ) { 
        std::vector<double> errs;
        for (int i=0;i<nBins ;++i) errs.push_back(hist->GetBinError(1+i));
        _aux->smooth( values, errs, smooth );
    }
    _aux->computeCoefficients( values );
    for (unsigned int i=0;i<values.size();++i) { 
        if (constCoeffs) {
            _coefList.add( RooFit::RooConst( values[i] ) );
        } else {
            stringstream name;
            name << name << "_smoothed_bin_" << i;
            const char* n = name.str().c_str();
            _coefList.addOwned(*new RooRealVar(n, n, values[i], 0.001, 0.999));
        }
    }
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun(const char* name, const char* title, 
                           RooRealVar& x, const char* knotBinningName, const RooArgList& coefList): 
  RooAbsReal(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this),
  _aux(0)
{
  // TODO: verify coefList is consistent with knots as specified by the knotBinningName binning
  //    should be N+2 coefficients for N bins...
  const RooAbsBinning* binning = x.getBinningPtr(knotBinningName);
  assert( binning!=0);
  assert( coefList.getSize()==2+binning->numBins());
 
  // Constructor
  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooCubicSplineFun::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
           << " is not of type RooRealVar" << endl ;
      assert(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter ;

  // Switch from bin boundaries to bin centres.
  Double_t* boundaries = binning->array();
  vector<double> centres;
  centres.reserve(binning->numBins());
  for (int i = 0; i < binning->numBins(); ++i) {
     centres.push_back((boundaries[i] + boundaries[i + 1]) / 2.);
  }
  _aux = new RooCubicSplineKnot(centres.begin(), centres.end());
}

//_____________________________________________________________________________
RooCubicSplineFun::RooCubicSplineFun(const RooCubicSplineFun& other, const char* name) :
  RooAbsReal(other, name), 
  _x("x", this, other._x), 
  _coefList("coefList",this,other._coefList),
  _aux(new RooCubicSplineKnot(*other._aux))
{
}

//_____________________________________________________________________________
RooCubicSplineFun::~RooCubicSplineFun()
{
    delete _aux;
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
Double_t RooCubicSplineFun::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code==1) ;
  return _aux->analyticalIntegral(_coefList);
}

