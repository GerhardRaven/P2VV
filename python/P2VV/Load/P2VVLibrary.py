"""load the P2VV library

Load the P2VV library, assuming that the library directory is in $LD_LIBRARYPATH
"""

print "P2VV - INFO: loading Cintex library"
import ROOT
ROOT.gSystem.Load("libCintex.so")
ROOT.Cintex.Enable()
ROOT.gSystem.Load("libP2VVDict.so")

__decorated = False
if not __decorated:
    extra = ['RooAbsEffResModel', 'RooAbsGaussModelEfficiency', 'RooAmoroso', 'RooAvEffConstraint',
             'RooBTagDecay', 'RooBinnedFun', 'RooBinnedPdf', 'RooBoxPdf', 'RooCalibratedDilution',
             'RooCategoryVar', 'RooComplementCoef', 'RooConvertPolAmp', 'RooCorrectedWeight',
             'RooCruijff', 'RooCubicSplineFun', 'RooCubicSplineKnot',
             'RooCubicSplineKnot::BoundaryConditions', 'RooEffConvGenContext', 'RooEffResAddModel',
             'RooEffResModel', 'RooEfficiencyBin', 'RooExplicitNormPdf', 'RooGaussEfficiencyModel',
             'RooIpatia2', 'RooMassDependence', 'RooMultiCatGenerator', 'RooMultiEffResModel',
             'RooP2VVAngleBasis', 'RooRealCategory', 'RooRealEffMoment', 'RooRealMoment',
             'RooRelBreitWigner', 'RooTPDecay', 'RooTagDecisionWrapper', 'RooThresholdPdf',
             'RooTransAngle', 'RooTrivialTagDecay']
    from P2VV.ROOTDecorators import set_classes_creates
    set_classes_creates(dict((cl, []) for cl in extra))
    __decorated = True
