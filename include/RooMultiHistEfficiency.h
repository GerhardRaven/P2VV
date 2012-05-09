#ifndef ROO_MULTIHISTEFFICIENCY
#define ROO_MULTIHISTEFFICIENCY

#include <RooAbsPdf.h>
#include <RooListProxy.h>
#include <TString.h>
#include <TList.h>

class RooRealVar;
class RooArgList ;
class TIterator;
class RooAbsGenContext;

class RooMultiHistEfficiency : public RooAbsPdf {
public:
   // Constructors, assignment etc
   inline RooMultiHistEfficiency() { 
      // Default constructor
   }
   RooMultiHistEfficiency(const char *name, const char *title,
                          const RooArgList& efficiencies, const RooArgList& categories,
                          const TList& sigCatNames);
   RooMultiHistEfficiency(const RooMultiHistEfficiency& other, const char* name=0);

   virtual TObject* clone(const char* newname) const
   {
      return new RooMultiHistEfficiency(*this, newname);
   }
   virtual ~RooMultiHistEfficiency();


   // RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype,
   //                              const RooArgSet* auxProto, Bool_t verbose) const;
   // virtual Int_t getGenerator(const RooArgSet& dv, RooArgSet &gv, Bool_t si) const;
   // virtual void initGenerator(Int_t code);
   // virtual void generateEvent(Int_t code);

   virtual Double_t evaluate() const;

private:

   RooMultiHistEfficiency* operator=(const RooMultiHistEfficiency& other);
   bool allFalse() const;

   RooListProxy _categories;      // Accept/reject categories
   RooListProxy _efficiencies;    // Efficiency modeling functions
   TList _sigCatNames;            // Name of accept state of accept/reject categories
   mutable RooArgSet _pdfGenVars; // Variables to generate for the PDF
   
   typedef std::vector<std::pair<double, TString> > Levels;
   Levels _levels; // 

   TIterator* _nameIter;        //!

   ClassDef(RooMultiHistEfficiency,1) // Generic PDF defined by string expression and list of variables
};

#endif
