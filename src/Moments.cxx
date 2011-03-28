#include <cmath>
#include "P2VV.h"
#include "Moments.h"
#include "RooArgSet.h"
#include "RooAbsData.h"

IMoment::IMoment(RooAbsReal& basis, double norm, const std::string& name) :
  _basis(basis), _m0(0.), _m1(0.), _n0(0.), _n1(0.), _n2(0.),
  _norm(norm), _name(name.empty() ? _basis.GetName() : name) {}

void IMoment::inc(double weight)
{
  double x = evaluate();

  // TODO: make a histogram of x... (two, one for accept, one for all)
  _m0 += weight;
  _m1 += weight * x;

  // these we need to compute the error using the jackknife method
  _n0 += weight * weight;
  _n1 += weight * weight * x;
  _n2 += weight * weight * x * x;
}

double IMoment::varmu() const {
  // the following formulas follow either from the jackknife method
  // or from error propagation using the following error on weight_j:
  //     sigma^2( weight_j ) = weight_j^2
  // (this is also exactly how it works with s-weight).

  // jackknife: sigma2 = (N - 1)/N *
  // sum_j ( mj1 - m )^2) where .mJ1 is the value of m if you would
  // leave measurement J away

  // we make one approximation: we ignore the contribution of a
  // single weight to the total in a normalization term

  // var(mu) = 1/m0^2 * sum  w_j^2 (x_j - mu)^2
  double mu    = _m1 / _m0;
  double varmu = (_n2 - 2. * _n1 * mu + _n0 * mu * mu) / (_m0 * _m0);

  return varmu ;
}

double IMoment::significance() const
{
  double mu = _m1 / _m0;
  double var = varmu() ;
  return var > 0 ? std::sqrt(mu * mu / var) : 999;
}

ostream& IMoment::print(ostream& os) const
{
  double mu = _m1 / _m0;
  double var = varmu();
  double sig = var > 0. ? std::sqrt(var) : 0.;

  return os << "moment(" << _name << ") = " << mu << " +- " << sig
      << " (significance: " << significance() << ")" << endl;
}

int _computeMoments(RooAbsData& data, IMomentsVector& moments)
{
  typedef IMomentsVector::iterator IMomIter;

  if (moments.empty()) {
    cout << "P2VV - ERROR: computeMoments: moments vector is empty" << endl;
    return -1;
  }

  for (IMomIter mom = moments.begin(); mom != moments.end(); ++mom)
    (*mom)->reset();

  RooArgSet* obs = moments.front()->basis().getObservables(data);
  
  cout << "P2VV - INFO: computeMoments: computing " << moments.size()
      << " moment(s) for data set '" << data.GetName() << "'" << endl;

  int dataIter = 0;
  ProgressDisplay progress(data.numEntries());
  while (dataIter < data.numEntries()) {
    *obs = *data.get(dataIter++);
    for (IMomIter mom = moments.begin(); mom != moments.end(); ++mom)
      (*mom)->inc(data.isWeighted() ? data.weight() : 1.);

    ++progress;
  }

  cout << endl;

  return dataIter;
}

