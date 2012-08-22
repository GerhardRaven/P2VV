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
#include <MultiHistEntry.h>

class RooAbsGenContext;
class RooArgList ;
class RooEffHistProd;
class RooAbsCategory;
class RooMultiHistEfficiency;
class RooAddition;

class RooMultiHistEfficiency : public RooAbsPdf {
public:

   typedef MultiHistEntry<RooEffHistProd, RooMultiHistEfficiency> HistEntry;

   // Constructors, assignment etc
   inline RooMultiHistEfficiency() { 
      // Default constructor
   }
   RooMultiHistEfficiency(const char *name, const char *title,
                          std::vector<RooMultiHistEfficiency::HistEntry*> entries);
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
   // virtual Int_t getAnalyticalIntegralWN(RooArgSet& allDeps, RooArgSet& analDeps, 
   //                               const RooArgSet* normSet, const char* rangeName) const;
   virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& iset,
                               const char* rangeName) const;
   virtual Double_t analyticalIntegral(Int_t code, const char* rangeName) const;
   
   // virtual Double_t getValV(const RooArgSet* ns) const;

   const std::map<Int_t, RooMultiHistEfficiency::HistEntry*>& getEntries() const
   {
      return _entries;
   }
   
   const RooSuperCategory* getSuper() const
   {
      return _super;
   }

protected:

   virtual Double_t evaluate() const;

private:

   RooMultiHistEfficiency* operator=(const RooMultiHistEfficiency& other);

   // Evaluation
   typedef std::list<double> BinBoundaries;
   BinBoundaries* _binboundaries;
   typedef std::map<Int_t, RooMultiHistEfficiency::HistEntry*> HistEntries;
   HistEntries _entries;

   typedef std::map<Int_t, double> IntVals;
   mutable IntVals _intVals;
      
   // Generation
   mutable RooArgSet _prodGenObs;
   mutable Int_t _prodGenCode;
   RooSuperCategory* _super;
   typedef std::vector<std::pair<double, TString> > Levels;
   Levels _levels;

   // Integration
   typedef std::vector<std::pair<RooAbsReal*, RooEffHistProd*> > AddEntries;

   class CacheElem : public RooAbsCacheElement {
   public:
      CacheElem(const HistEntries& entries, const RooArgSet& iset, const char* rangeName);
      virtual ~CacheElem();

      virtual Double_t getVal(const Int_t index) const;

      virtual RooArgList containedArgs(Action);

   private:
      // Payload
      typedef std::map<Int_t, RooAbsReal*> Integrals;
      Integrals _I;
      RooArgSet _iset;
   };
   mutable RooObjCacheManager _cacheMgr ; // The cache manager

   ClassDef(RooMultiHistEfficiency,1) // Generic PDF defined by string expression and list of variables
};

#endif
