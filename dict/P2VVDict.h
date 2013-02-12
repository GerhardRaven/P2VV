#ifndef DICT_P2VVDICT_H 
#define DICT_P2VVDICT_H 1

#include "P2VV/RooRealMoments.h"
#include "P2VV/ProgressDisplay.h"
#include "P2VV/RooBTagDecay.h"
#include "P2VV/RooTrivialTagDecay.h"
#include "P2VV/RooMultiCatGenerator.h"
#include "P2VV/RooBinnedPdf.h"
#include "P2VV/RooP2VVAngleBasis.h"
#include "P2VV/RooThresholdPdf.h"
#include "P2VV/RooRelBreitWigner.h"
#include "P2VV/RooTagDecisionWrapper.h"
#include "P2VV/RooRealCategory.h"
#include "P2VV/RooCalibratedDilution.h"
#include "P2VV/RooDataSetToTree.h"
#include "P2VV/Functions.h"
#include "P2VV/RooTransAngle.h"
#include "P2VV/RooCruijff.h"
#include "P2VV/RooEfficiencyBin.h"
#include "P2VV/RooAvEffConstraint.h"
#include "P2VV/RooCorrectedSWeight.h"
#include "P2VV/RooAbsEffResModel.h"
#include "P2VV/RooEffResModel.h"
#include "P2VV/RooMultiEffResModel.h"
#include "P2VV/MultiHistEntry.h"
#include "P2VV/RooComplementCoef.h"
#include "P2VV/RooEffConvGenContext.h"
#include "P2VV/RooBoxPdf.h"
#include "P2VV/RooExplicitNormPdf.h"
#include "P2VV/RooAmoroso.h"
#include "P2VV/RooTPDecay.h"

#include <map>
#include <string>
#include <vector>

template class std::vector<MultiHistEntry*>;
template class std::vector<RooAbsRealMoment*>;

template class std::map<RooAbsCategory*, std::string>;
template class std::map<RooCategoryProxy*, std::string>;
template class std::vector<std::pair<double, TString> >;

template class std::map<Int_t, MultiHistEntry*>;
template class std::pair<Int_t, MultiHistEntry*>;

template class std::pair<RooAbsCategory*, std::string>;

template class std::map<RooRealProxy*, bool>;
template class std::pair<RooRealProxy*,bool>;
template class std::map<RooAbsReal*, bool>;
template class std::pair<RooAbsReal*, bool>;

template class std::pair<Double_t, TString>;
template class std::pair<RooCategoryProxy*, std::string>;

template class std::vector<MultiHistEntry>;

struct Instantiations {
   std::map<int, MultiHistEntry*>::iterator _i00;
};

#endif // DICT_P2VVDICT_H
