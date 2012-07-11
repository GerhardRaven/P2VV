#ifndef P2VVInc_H 
#define P2VVInc_H 1

#include "Moments.h"
#include "RooRealMoments.h"
#include "ProgressDisplay.h"
#include "RooBTagDecay.h"
#include "RooTrivialTagDecay.h"
#include "RooMultiCatGenerator.h"
#include "RooBinnedPdf.h"
#include "RooP2VVAngleBasis.h"
#include "RooThresholdPdf.h"
#include "RooEffHistProd.h"
#include "RooRelBreitWigner.h"
#include "RooTagDecisionWrapper.h"
#include "RooRealCategory.h"
#include "RooCalibratedDilution.h"
#include "RooDataSetToTree.h"
#include "RooTransAngle.h"
#include "RooCruijff.h"
#include "RooMultiHistEfficiency.h"
#include "RooEfficiencyBin.h"
#include "RooAvEffConstraint.h"

#include <map>
#include <string>
#include <vector>

struct Instantiations {

   std::map<RooAbsCategory*, std::string>   _i00;
   std::map<RooCategoryProxy*, std::string> _i01;
   std::vector<std::pair<double, TString> > _i02;
   std::map<Int_t, MultiHistEntry*>         _i03;
   std::pair<Int_t, MultiHistEntry*>        _i04;

   std::pair<RooAbsCategory*, std::string>  _i05;

   std::map<RooRealProxy*, bool> _i06;
   std::pair<RooRealProxy*,bool> _i07;
   std::map<RooAbsReal*, bool>   _i08;
   std::pair<RooAbsReal*, bool>  _i09;
   std::vector<MultiHistEntry> _i10;

   std::map<int, MultiHistEntry*>::iterator _i11;

   std::pair<Double_t, TString> _i12;
   std::pair<RooCategoryProxy*, std::string> _i13;
};
#endif // P2VV_H
