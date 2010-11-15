
void rootlogon() {
    gSystem->Load("libMathMore.so");
    gROOT->ProcessLine(".L RooLegendre.cxx+");
    gROOT->ProcessLine(".L RooSpHarmonic.cxx+");
    gROOT->ProcessLine(".L RooP2VVAngleBasis.cxx+");
}
