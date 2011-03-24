/******************************************************************************
** Moments:                                                                  **
** tools for calculation of efficiency and background angular moments        **
**                                                                           **
** authors:                                                                  **
**   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           **
**   WH,  Wouter Hulsbergen,  Nikhef                                         **
**                                                                           **
******************************************************************************/

#include "Moments.h"
#include "RooAbsData.h"

IMoment::IMoment(RooAbsReal& basis, double norm, const std::string& name) :
  _basis(basis), _m0(0.), _m1(0.), _m2(0.), _norm(norm),
  _name(name.empty() ? _basis.GetName() : name) {}

void IMoment::inc(double weight)
{
  double x = evaluate();

  // TODO: make a histogram of x... (two, one for accept, one for all)
  _m0 += weight;
  _m1 += weight * x;
  _m2 += weight * x * x;
}

ostream& IMoment::print(ostream& os) const
{
  double mu = _m1 / _m0;
  double sig2 = _m2 / _m0 - mu * mu;

  return os << "moment(" << _name << ") = " << mu << " +- "
      << sqrt(sig2 / (_m0 - 1.)) << " significance: " << significance()
      << endl;
}

double IMoment::significance() const
{
  double mu = _m1/_m0;
  double sig2 = _m2/_m0 - mu*mu;
  return fabs(mu/sqrt(sig2/(_m0-1)));
}

int _computeMoments(RooAbsData& data, IMomentsVector& moments)
{
   typedef IMomentsVector::iterator IMomIter;

   if (moments.empty()) return -1;

   RooArgSet *obs = moments.front()->basis().getObservables(data);
   int dataIter = 0;
   while (dataIter < data.numEntries()) {
     *obs = *data.get(dataIter++);
     for (IMomIter mom = moments.begin(); mom != moments.end(); ++mom)
       (*mom)->inc(data.isWeighted() ? data.weight() : 1.);
   }

   return dataIter;
}

