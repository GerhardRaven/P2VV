/******************************************************************************
** Moments:                                                                  **
** tools for calculation of efficiency and background angular moments        **
**                                                                           **
** authors:                                                                  **
**   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           **
**   WH,  Wouter Hulsbergen,  Nikhef                                         **
**                                                                           **
******************************************************************************/

#ifndef MOMENTS_H
#define MOMENTS_H
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include <iosfwd>
#include <string>

class IMoment {
public:
  IMoment(RooAbsReal& basis, double norm = 1.,
      const std::string& name = std::string());
  virtual ~IMoment() {};

  virtual void inc(double weight = 1.);
  virtual RooAbsReal& basis() {return _basis;}
  virtual RooArgSet* getObservables(const RooArgSet* set) { return basis().getObservables(set); }
  RooArgSet* getObservables(const RooAbsData& data) { return getObservables(data.get()); }
  virtual double coefficient(bool normalize = true) const;
  virtual double variance(bool normalize = true) const;
  virtual double significance() const;
  virtual double evaluate() {return _basis.getVal();}

  void reset() {_m0 = _m1 = _n0 = _n1 = _n2 = 0.;}

  virtual std::ostream& print(std::ostream& os, bool normalize = true) const;
  void Print(bool normalize = true) const {print(std::cout, normalize);}

protected:
  RooAbsReal& _basis;
  double _m0, _m1;
  double _n0, _n1, _n2;
  double _norm;
  std::string _name;
};

class Moment : public IMoment {
public:
  Moment(RooAbsReal& x, double norm = 1.) : IMoment(x, norm) {}

private:
};

class EffMoment : public IMoment {
public:
  EffMoment(RooAbsReal& x, double norm, const RooAbsPdf& pdf,
      const RooArgSet& nset) :
    IMoment(x, norm, std::string(x.GetName()) + "_" + pdf.GetName()),
    _pdf(pdf), _nset(nset) {}

   double evaluate() {
       //cout << "EffMoment::evaluate():" << endl;
        //_nset.Print("V");
        //cout << _basis.GetName() << " = " << _basis.getVal() << "  ; " << _pdf.GetName() << " = " << _pdf.getVal(&_nset) << endl;
        return _basis.getVal() / _pdf.getVal(&_nset);}
   virtual RooArgSet* getObservables(const RooArgSet* set) { return _pdf.getObservables(set); }

private:
  const RooAbsPdf& _pdf;
  const RooArgSet& _nset;
};

typedef std::vector<IMoment*> IMomentsVector;

int _computeMoments(RooAbsData& data, IMomentsVector& moments, bool resetFirst = false);

#endif
