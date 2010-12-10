#ifndef MOMENT_H
#define MOMENT_H
#include "utils.h"
#include "RooAbsPdf.h"

// TODO: move and/or integrate this code into RooP2VVAngleBasis !!!
class IMoment {
    public:
          IMoment(RooAbsReal &basis, double norm=1, const char *name=0) : _basis(basis), _m0(0),_m1(0),_m2(0),_norm(norm), _name(name ? name : _basis.GetName() ) {}
          virtual ~IMoment() {};
          virtual void inc(bool accepted = true) {
                double x = evaluate();
                // TODO: make a histogram of x... (two, one for accept, one for all)
                _m0 += 1.;
                if (accepted) {
                    _m1 += x;
                    _m2 += x*x;
                }
            }
          virtual ostream& print(ostream& os) const {
                double mu = _m1/_m0;
                double sig2 = _m2/_m0 - mu*mu;
                return os << "moment("<< _name << ") = " << mu << " +- " << sqrt(sig2/(_m0-1)) << " significance: " << significance() << endl;
            }
          virtual RooAbsReal& basis() { return _basis; }
          virtual double coefficient() const { return _norm*_m1/_m0; }
          virtual double significance() const { 
                double mu = _m1/_m0;
                double sig2 = _m2/_m0 - mu*mu;
                return fabs(mu/sqrt(sig2/(_m0-1)));
          }
          virtual double evaluate() { return _basis.getVal(); }
    protected:
          RooAbsReal &_basis;
          double _m0,_m1,_m2;
          double _norm;
          const char *_name;
};


class Moment : public IMoment {
public:
    Moment(RooAbsReal& x, double norm=1) : IMoment(x,norm) {}
private:
};

class EffMoment  : public IMoment{
public:
    EffMoment(RooAbsReal& x, double norm, const RooAbsPdf& pdf, const RooArgSet& nset) : IMoment(x,norm,Format("%s_%s",x.GetName(),pdf.GetName())),_pdf(pdf), _nset(nset)  {}

    double evaluate() { return _basis.getVal()/_pdf.getVal(&_nset); }
private:
    const RooAbsPdf& _pdf;
    const RooArgSet&  _nset;
};

#include "RooAbsData.h"
typedef std::vector<IMoment*> IMomentsVector;
template class std::vector<IMoment*>;
int _computeMoments(RooAbsData& data, IMomentsVector& moments) {
   typedef std::vector<IMoment*>::iterator iter;
   if (moments.empty()) return -1; 
   RooArgSet *obs = moments.front()->basis().getObservables(data);
   int i=0;
   while (i<data.numEntries()) {
       *obs = *data.get(i++);
       for ( iter m = moments.begin(); m!=moments.end(); ++m) (*m)->inc();
   }
   return i;
}

#endif
