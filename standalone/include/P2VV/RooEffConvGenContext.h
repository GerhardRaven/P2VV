#ifndef ROO_EFF_CONV_GEN_CONTEXT
#define ROO_EFF_CONV_GEN_CONTEXT

#include "RooConvGenContext.h"
#include "RooArgSet.h"

class RooAbsAnaConvPdf;
class RooDataSet;
class RooAbsReal;

class RooEffConvGenContext : public RooConvGenContext {
public:
   RooEffConvGenContext(const RooFFTConvPdf &model, const RooArgSet &vars, const RooDataSet *prototype= 0,
                        const RooArgSet* auxProto=0, Bool_t _verbose= kFALSE);
   RooEffConvGenContext(const RooNumConvPdf &model, const RooArgSet &vars, const RooDataSet *prototype= 0,
                        const RooArgSet* auxProto=0, Bool_t _verbose= kFALSE);
   RooEffConvGenContext(const RooAbsAnaConvPdf &model, const RooArgSet &vars, const RooDataSet *prototype= 0,
                        const RooArgSet* auxProto=0, Bool_t _verbose= kFALSE);
   virtual ~RooEffConvGenContext() { }

   virtual void attach(const RooArgSet& params);

protected:

   virtual void generateEvent(RooArgSet &theEvent, Int_t remaining);

   RooEffConvGenContext(const RooEffConvGenContext& other) ;

private:

   void initEfficiency();
   const RooAbsReal* efficiency();

   double _maxEff;

   ClassDef(RooEffConvGenContext,0) // Context for generating a dataset from a PDF
};

#endif
