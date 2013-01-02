#ifdef __CINT__
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class IMoment;
#pragma link C++ class Moment+;
#pragma link C++ class EffMoment+;
#pragma link C++ function _computeMoments;
#pragma link C++ class vector<IMoment*>;
#pragma link C++ class RooAbsRealMoment+;
#pragma link C++ class RooRealMoment+;
#pragma link C++ class RooRealEffMoment+;
#pragma link C++ function computeRooRealMoments;
#pragma link C++ class vector<RooAbsRealMoment*>+;
#pragma link C++ class ProgressDisplay+;
#pragma link C++ class RooBTagDecay+;
#pragma link C++ class RooTrivialTagDecay+;
#pragma link C++ class RooMultiCatGenerator+;
#pragma link C++ class RooBinnedPdf+;
#pragma link C++ class RooP2VVAngleBasis+;
#pragma link C++ class RooThresholdPdf+;
#pragma link C++ class RooRelBreitWigner+;
#pragma link C++ class RooTagDecisionWrapper+;
#pragma link C++ class RooRealCategory+;
#pragma link C++ class RooCalibratedDilution+;
#pragma link C++ function RooDataSetToTree;
#pragma link C++ class RooTransAngle+;
#pragma link C++ class RooCruijff+;
#pragma link C++ class RooEfficiencyBin+;
#pragma link C++ class RooAvEffConstraint+;
#pragma link C++ class RooCorrectedSWeight+;
#pragma link C++ class RooAbsEffResModel;
#pragma link C++ class RooEffResModel+;
#pragma link C++ class MultiHistEntry+;
#pragma link C++ class RooMultiEffResModel+;
#pragma link C++ class RooComplementCoef+;
#pragma link C++ class RooEffConvGenContext+;
#pragma link C++ function sigmaFromFT;
#pragma link C++ class RooBoxPdf+;
#pragma link C++ class RooExplicitNormPdf+;


#pragma link C++ class std::map<RooAbsCategory*, std::string>;
#pragma link C++ class std::map<RooCategoryProxy*, std::string>;
#pragma link C++ class std::vector<std::pair<double, TString> >;
#pragma link C++ class std::pair<RooAbsCategory*, std::string>;
#pragma link C++ class std::map<RooRealProxy*, bool>;
#pragma link C++ class std::map<RooAbsReal*, bool>;
#pragma link C++ class std::pair<RooAbsReal*, bool>;
#pragma link C++ class std::pair<RooRealProxy*, bool>;
#pragma link C++ class std::pair<double, TString>;
#pragma link C++ class std::pair<RooCategoryProxy*, std::string>;

#pragma link C++ class std::map<Int_t, MultiHistEntry*>;
#pragma link C++ class std::pair<Int_t, MultiHistEntry*>;
#pragma link C++ class std::map<int, MultiHistEntry*>::iterator;
#pragma link C++ class std::vector<MultiHistEntry*>;

#endif
