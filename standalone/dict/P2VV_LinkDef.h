#ifdef __CINT__
#pragma link off all classes;
#pragma link C++ class RooRealMoment+;
#pragma link C++ class RooRealEffMoment+;
#pragma link C++ class vector<RooRealMoment*>+;
#pragma link C++ class RooBTagDecay+;
#pragma link C++ class RooTrivialTagDecay+;
#pragma link C++ class RooMultiCatGenerator+;
#pragma link C++ class RooAbsGaussModelEfficiency+;
#pragma link C++ class RooBinnedPdf+;
#pragma link C++ class RooBinnedFun+;
#pragma link C++ class RooP2VVAngleBasis+;
#pragma link C++ class RooThresholdPdf+;
#pragma link C++ class RooRelBreitWigner+;
#pragma link C++ class RooTagDecisionWrapper+;
#pragma link C++ class RooRealCategory+;
#pragma link C++ class RooCalibratedDilution+;
#pragma link C++ class RooTransAngle+;
#pragma link C++ class RooCruijff+;
#pragma link C++ class RooEfficiencyBin+;
#pragma link C++ class RooAvEffConstraint+;
#pragma link C++ class RooCorrectedWeight+;
#pragma link C++ class RooAbsEffResModel;
#pragma link C++ class RooEffResModel+;
#pragma link C++ class MultiHistEntry+;
#pragma link C++ class RooMultiEffResModel+;
#pragma link C++ class RooComplementCoef+;
#pragma link C++ class RooEffConvGenContext+;
#pragma link C++ class RooCubicSplineKnot::BoundaryConditions+;
#pragma link C++ class RooCubicSplineKnot+;
#pragma link C++ class RooCubicSplineFun+;
#pragma link C++ class RooGaussEfficiencyModel+;
#pragma link C++ class RooBoxPdf+;
#pragma link C++ class RooExplicitNormPdf+;
#pragma link C++ class RooAmoroso+;
#pragma link C++ class RooTPDecay+;
#pragma link C++ class RooMassDependence+;
#pragma link C++ class RooEffResAddModel+;
#pragma link C++ class RooCategoryVar+;
#pragma link C++ class RooConvertPolAmp+;
#pragma link C++ class RooIpatia2+;

#pragma link off all functions;
#pragma link C++ function computeRooRealMoments;
#pragma link C++ function RooDataSetToTree;
#pragma link C++ function TreeToRooDataSet;
#pragma link C++ function addSWeightToTree;
#pragma link C++ function addIntegerToTree;
#pragma link C++ function addCategoryToTree;
#pragma link C++ function addVertexErrors;
#pragma link C++ function sigmaFromFT;
#pragma link C++ function hessian;
#pragma link C++ function hessian_with_errors;
#pragma link C++ function gradient;
#pragma link C++ function gradient_with_errors;
#pragma link C++ function second_gradient;
#pragma link C++ function second_gradient_with_errors;
#pragma link C++ function cross_with_errors;
#pragma link C++ function initial_stepsizes;
#pragma link C++ function HelicityAngles;
#pragma link C++ function GetOwnership;
#pragma link C++ function getRooRealMaxVal;

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
#pragma link C++ class std::list<RooDataSet*>;
#pragma link C++ class std::vector<std::pair<double, double> >;
#pragma link C++ class std::pair<TMatrixDSym, TMatrixDSym>;

#endif
