#include <set>
#include <vector>

#include <RooAddModel.h>

#include <P2VV/RooAbsEffResModel.h>
#include <P2VV/RooEffResAddModel.h>

//_____________________________________________________________________________R
RooEffResAddModel::RooEffResAddModel()
   : RooAddModel(), RooAbsEffResModel()
{

}

//_____________________________________________________________________________
RooEffResAddModel::RooEffResAddModel(const char *name, const char *title, const RooArgList& modelList,
                  const RooArgList& coefList, Bool_t ownPdfList)
   : RooAddModel(name, title, modelList, coefList, ownPdfList),
     RooAbsEffResModel()
{
   RooFIter iter = modelList.fwdIterator();
   while(RooAbsArg* model = iter.next()) {
      RooAbsEffResModel* effModel = dynamic_cast<RooAbsEffResModel*>(model);
      assert(effModel);
   }
}

//_____________________________________________________________________________
RooEffResAddModel::RooEffResAddModel(const RooEffResAddModel& other, const char* name)
   : RooAddModel(other, name),
     RooAbsEffResModel()
{

}
   
//_____________________________________________________________________________
const RooAbsReal* RooEffResAddModel::efficiency() const
{
   std::vector<const RooAbsReal*> effs(efficiencies());
   std::set<const RooAbsReal*> eff_set(effs.begin(), effs.end());
   assert(eff_set.size() == 1);
   return *eff_set.begin();
}

//_____________________________________________________________________________
std::vector<const RooAbsReal*> RooEffResAddModel::efficiencies() const
{
   std::set<const RooAbsReal*> eff_set;
   RooFIter it = pdfList().fwdIterator();
   while(RooAbsArg* arg = it.next()) {
      const RooAbsEffResModel* model = dynamic_cast<const RooAbsEffResModel*>(arg);
      std::vector<const RooAbsReal*> model_effs = model->efficiencies();
      eff_set.insert(model_effs.begin(), model_effs.end());
   }
   return std::vector<const RooAbsReal*>(eff_set.begin(), eff_set.end());
}

//_____________________________________________________________________________
const RooArgSet* RooEffResAddModel::observables() const { 
   // Return pointer to pdf in product
   return new RooArgSet(RooAddModel::convVar());
}
