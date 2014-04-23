#ifndef ROO_TP_CONV_GEN_CONTEXT
#define ROO_TP_CONV_GEN_CONTEXT

// use C preprocessor to smuggle in a friend class into RooAbsGenContext
//
// *** this is an ugly hack but seems to work well enough for now, and is
// *** better than hacking the header files that come with ROOT for now...
#define RooProdGenContext RooProdGenContext; friend class RooTPConvGenContext;

#include "RooConvGenContext.h"
#include "RooArgSet.h"
#include "RooArgList.h"

class RooAbsAnaConvPdf;
class RooDataSet;
class RooAbsReal;

class RooTPConvGenContext : public RooConvGenContext {
public:
   RooTPConvGenContext(const RooAbsAnaConvPdf &model, const RooArgSet &vars, const RooArgList& tps,
                       const RooDataSet *prototype= 0, const RooArgSet* auxProto=0,
                       Bool_t _verbose= kFALSE);
   virtual ~RooTPConvGenContext() { }

protected:

   virtual void generateEvent(RooArgSet &theEvent, Int_t remaining);

   RooTPConvGenContext(const RooTPConvGenContext& other) ;

private:

   RooArgList _tps;
   const char* _normRange;

   ClassDef(RooTPConvGenContext,0) // Context for generating a dataset from a PDF
};

#endif
