#ifndef DICT_P2VVDICT_H 
#define DICT_P2VVDICT_H 1

#include "TMatrixTSym.h"

#include "P2VV/RooRealMoments.h"
#include "P2VV/RooBTagDecay.h"
#include "P2VV/RooTrivialTagDecay.h"
#include "P2VV/RooMultiCatGenerator.h"
#include "P2VV/RooAbsGaussModelEfficiency.h"
#include "P2VV/RooBinnedPdf.h"
#include "P2VV/RooBinnedFun.h"
#include "P2VV/RooP2VVAngleBasis.h"
#include "P2VV/RooThresholdPdf.h"
#include "P2VV/RooRelBreitWigner.h"
#include "P2VV/RooTagDecisionWrapper.h"
#include "P2VV/RooRealCategory.h"
#include "P2VV/RooCalibratedDilution.h"
#include "P2VV/Functions.h"
#include "P2VV/RooTransAngle.h"
#include "P2VV/RooCruijff.h"
#include "P2VV/RooEfficiencyBin.h"
#include "P2VV/RooAvEffConstraint.h"
#include "P2VV/RooCorrectedWeight.h"
#include "P2VV/RooAbsEffResModel.h"
#include "P2VV/RooEffResModel.h"
#include "P2VV/RooMultiEffResModel.h"
#include "P2VV/MultiHistEntry.h"
#include "P2VV/RooComplementCoef.h"
#include "P2VV/RooEffConvGenContext.h"
#include "P2VV/RooTPConvGenContext.h"
#include "P2VV/RooBoxPdf.h"
#include "P2VV/RooExplicitNormPdf.h"
#include "P2VV/RooCubicSplineKnot.h"
#include "P2VV/RooCubicSplineFun.h"
#include "P2VV/RooGaussEfficiencyModel.h"
#include "P2VV/RooAmoroso.h"
#include "P2VV/RooTPDecay.h"
#include "P2VV/RooTPGen.h"
#include "P2VV/RooMassDependence.h"
#include "P2VV/RooEffResAddModel.h"
#include "P2VV/RooHessian.h"
#include "P2VV/RooCategoryVar.h"
#include "P2VV/RooConvertPolAmp.h"
#include "P2VV/RooIpatia2.h"
#include "P2VV/RooEffConstraint.h"
#include "P2VV/RooCombEffConstraint.h"

#include <map>
#include <string>
#include <vector>
#include <list>

struct Instantiations {

   std::map<RooAbsCategory*, std::string>   _i00;
   std::map<RooCategoryProxy*, std::string> _i01;
   std::vector<std::pair<double, TString> > _i02;

   std::map<Int_t, MultiHistEntry*>         _i031;
   std::pair<Int_t, MultiHistEntry*>        _i041;

   std::pair<RooAbsCategory*, std::string>  _i05;

   std::map<RooRealProxy*, bool> _i06;
   std::pair<RooRealProxy*,bool> _i07;
   std::map<RooAbsReal*, bool>   _i08;
   std::pair<RooAbsReal*, bool>  _i09;

   std::pair<Double_t, TString> _i10;
   std::pair<RooCategoryProxy*, std::string> _i11;

   std::vector<MultiHistEntry> _i16;
   std::map<int, MultiHistEntry*>::iterator _i17;
   std::vector<MultiHistEntry*> _i18;
   std::vector<RooRealMoment*> _i19;

   std::list<RooDataSet*> _i20;
   std::vector<std::pair<double, double> > _i21;
   std::pair<TMatrixDSym, TMatrixDSym> _i22;

   std::map<std::string, RooDataSet*> _i23;
   std::pair<std::string, RooDataSet*> _i24;

};
#endif // DICT_P2VVDICT_H
