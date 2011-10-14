#ifndef MOMENT_H
#define MOMENT_H
#include "RooAbsPdf.h"
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

// TODO: move and/or integrate this code into RooP2VVAngleBasis !!!
// TODO: pick up code from Wouter which computes the proper errors in case of weighted events...
class IMoment {
    public:
          IMoment(RooAbsReal &basis, double norm=1, const std::string& name = std::string()) : _basis(basis), _m0(0),_m1(0),_m2(0),_norm(norm), _name(name.empty() ? _basis.GetName() : name) {}
          virtual ~IMoment() {};
          virtual void inc(double weight = 1.) {
                double x = evaluate();
                // TODO: make a histogram of x... (two, one for accept, one for all)
                _m0 += weight;
                _m1 += weight*x;
                _m2 += weight*x*x;
		//std::cout << "m0 = " << _m0 << std::endl;
		//std::cout << "m1 = " << _m1 << std::endl;
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
          virtual double evaluate() {
	    //std::cout << "evaluate in IMoment" << std::endl;
	    return _basis.getVal(); 
	  }
    protected:
          RooAbsReal &_basis;
          double _m0,_m1,_m2;
          double _norm;
          std::string _name;
};


class Moment : public IMoment {
public:
    Moment(RooAbsReal& x, double norm=1) : IMoment(x,norm) {}
private:
};

class EffMoment  : public IMoment{
public:
    EffMoment(RooAbsReal& x, double norm, const RooAbsPdf& pdf, const RooArgSet& nset) 
        : IMoment(x,norm,std::string(x.GetName())+"_"+pdf.GetName()),_pdf(pdf), _nset(nset)  
    {}
      double evaluate() { 
	//std::cout << "**************************************" << std::endl;
	//std::cout << "basis = " << _basis.GetName() << std::endl;
	//std::cout << "basis in EffMoment = " << _basis.getVal() << std::endl;
	////std::cout << "norm corr basis in EffMoment = " << _basis.getVal()*(1./(16.*sqrt(M_PI))) << std::endl;
	//std::cout << "denominator in EffMoment = " << _pdf.getVal(&_nset) << std::endl;
	////std::cout << "norm corr PDF in EffMoment = " << _pdf.getVal(&_nset)*(1./(16.*sqrt(M_PI))) << std::endl;
	//std::cout << "ratio in EffMoment = " << _basis.getVal()/_pdf.getVal(&_nset) << std::endl;
	return _basis.getVal()/_pdf.getVal(&_nset); 
      }
private:
    const RooAbsPdf& _pdf;
    const RooArgSet&  _nset;
};

#include "RooAbsData.h"
typedef std::vector<IMoment*> IMomentsVector;
int _computeMoments(RooAbsData& data, RooAbsPdf& pdf, IMomentsVector& moments) {
   typedef IMomentsVector::iterator iter;
   if (moments.empty()) return -1; 
   //RooArgSet *obs = moments.front()->basis().getObservables(data);
   RooArgSet *obs = pdf.getObservables(data);
   //std::cout << "observables" << *obs << std::endl;
   int i=0;
   while (i<data.numEntries()) {
       *obs = *data.get(i++);\
       //std::cout << "########################################" << std::endl;
       //std::cout << "!!!!!!!!!!!!!!!!!! i = " << i << std::endl;
       //std::cout << "########################################" << std::endl;
       for ( iter m = moments.begin(); m!=moments.end(); ++m) (*m)->inc( data.isWeighted() ? data.weight() : 1.0) ;
   }
   return i;
}

#endif
