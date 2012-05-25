#ifndef ROO_MULTIHISTEFFICIENCY
#define ROO_MULTIHISTEFFICIENCY

#include <map>

#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooListProxy.h>
#include <RooCategoryProxy.h>

#include <TString.h>
#include <TList.h>

#include <RooEffHistProd.h>

class RooAbsGenContext;
class RooArgList ;
class RooEffHistProd;
class RooAbsCategory;
class RooMultiHistEfficiency;

class MultiHistEntry {
public:

   MultiHistEntry();
   MultiHistEntry(const std::map<RooAbsCategory*, std::string>& categories,
                  RooEffHistProd* effProd, RooAbsReal* relative);
   MultiHistEntry(const MultiHistEntry& other, RooMultiHistEfficiency* parent);
   virtual ~MultiHistEntry();

   const RooEffHistProd& effProd() const {
      return const_cast<MultiHistEntry*>(this)->effProd();
      // }
   }

   RooEffHistProd& effProd() {
      if (m_effProd) {
         RooEffHistProd* eff = dynamic_cast<RooEffHistProd*>(m_effProd->absArg());
         return *eff;
      } else {
         return *m_rawEff;
      }
   }
   
   Double_t effVal() {
      return m_effProd ? *m_effProd : m_rawEff->getVal();
   }

   double relative() const {
      return m_relative ? *m_relative : m_rawRel->getVal();
   }

   void setParent(RooMultiHistEfficiency* parent);

   RooArgSet categories() const;

   bool thisEntry() const;

   void setIndex(const Int_t index) {
      m_index = index;
   }
   
   Int_t index() const {
      return m_index;
   }

   void select();

private:

   MultiHistEntry& operator=(const MultiHistEntry& other);

   std::map<RooAbsCategory*, std::string> m_rawCats;
   RooEffHistProd* m_rawEff; //!
   RooAbsReal* m_rawRel; //!

   std::map<RooCategoryProxy*, std::string> m_categories;
   RooRealProxy* m_effProd;
   RooRealProxy* m_relative;
   Int_t m_index;

};

class RooMultiHistEfficiency : public RooAbsPdf {
public:
   // Constructors, assignment etc
   inline RooMultiHistEfficiency() { 
      // Default constructor
   }
   RooMultiHistEfficiency(const char *name, const char *title,
                          std::vector<MultiHistEntry> entries);
   RooMultiHistEfficiency(const RooMultiHistEfficiency& other, const char* name = 0);

   virtual TObject* clone(const char* newname) const
   {
      return new RooMultiHistEfficiency(*this, newname);
   }
   virtual ~RooMultiHistEfficiency();

   std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t min, Double_t max) const;

   RooAbsGenContext* genContext(const RooArgSet &vars, const RooDataSet *prototype,
                                const RooArgSet* auxProto, Bool_t verbose) const;
   virtual Int_t getGenerator(const RooArgSet& dv, RooArgSet &gv, Bool_t si) const;
   virtual void initGenerator(Int_t code);
   virtual void generateEvent(Int_t code);

   virtual Bool_t	forceAnalyticalInt(const RooAbsArg& var) const;
   virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& iset,
                               const char* rangeName) const;
   virtual Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

   virtual Double_t evaluate() const;

private:

   RooMultiHistEfficiency* operator=(const RooMultiHistEfficiency& other);

   // Evaluation
   typedef std::list<double> BinBoundaries;
   BinBoundaries* _binboundaries;
   typedef std::map<Int_t, MultiHistEntry*> HistEntries;
   HistEntries _entries;
      
   // Generation
   mutable RooArgSet _prodGenObs;
   mutable Int_t _prodGenCode;
   RooSuperCategory* _super;
   typedef std::vector<std::pair<double, TString> > Levels;
   Levels _levels;

   // Integration
   mutable Int_t _pdfIntCode;
   mutable RooArgSet _pdfIntObs;

   ClassDef(RooMultiHistEfficiency,1) // Generic PDF defined by string expression and list of variables
};

#endif
