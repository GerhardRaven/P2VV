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
#include <string>

class IMoment {
public:
  IMoment(RooAbsReal& basis, double norm = 1.,
      const std::string& name = std::string());
  virtual ~IMoment() {};

  virtual void inc(double weight = 1.);
  virtual RooAbsReal& basis() {return _basis;}
  virtual double coefficient() const {return _norm * _m1 / _m0;}
  virtual double varmu() const;
  virtual double significance() const;
  virtual double evaluate() {return _basis.getVal();}

  void reset() {_m0 = _m1 = _n0 = _n1 = _n2 = 0.;}

  virtual ostream& print(ostream& os) const;
  void Print() const {print(std::cout);}

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

    double evaluate() {return _basis.getVal() / _pdf.getVal(&_nset);}

private:
  const RooAbsPdf& _pdf;
  const RooArgSet& _nset;
};

typedef std::vector<IMoment*> IMomentsVector;

int _computeMoments(RooAbsData& data, IMomentsVector& moments);

#endif
