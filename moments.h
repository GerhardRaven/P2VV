
#ifndef MOMENT_H
#define MOMENT_H
class IMoment {
    public:
          IMoment(RooAbsReal &basis, const char *name=0) : _basis(basis), _m0(0),_m1(0),_m2(0), _name(name ? name : _basis.GetName() ) {}
          virtual ~IMoment() {};
          virtual void inc(bool accepted = true) {
                double x = evaluate();
                // TODO: make a histogram of x... (two, one for accept, one for all)
                ++_m0;
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
          virtual double coefficient() const = 0;
          virtual double significance() const { 
                double mu = _m1/_m0;
                double sig2 = _m2/_m0 - mu*mu;
                return fabs(mu/sqrt(sig2/(_m0-1)));
          }
          virtual double evaluate(  ) = 0;
    protected:
          RooAbsReal &_basis;
          double _m0,_m1,_m2;
          const char *_name;
};


class Moment : public IMoment {
public:
    Moment(RooAbsReal& x, double norm=1) : IMoment(x), _norm(norm) {}
    double evaluate() { return _basis.getVal(); }
    double coefficient() const { return _norm*_m1/_m0; }
private:
    double _norm;
};


class EffMoment  : public IMoment{
public:
    EffMoment(RooAbsReal& x, const RooAbsPdf& pdf, const RooArgSet& nset) : IMoment(x,Format("%s_%s",x.GetName(),pdf.GetName())),_pdf(pdf), _nset(nset)  {}

    double evaluate() { return _basis.getVal()/_pdf.getVal(&_nset); }
    double coefficient() const { return _m1/_m0; }
private:
    const RooAbsPdf& _pdf;
    const RooArgSet&  _nset;
};

#endif
