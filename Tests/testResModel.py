import ROOT
from ROOT import (RooDataHist, RooWorkspace,
                  RooBrentRootFinder, RooRealBinding,
                  RooArgSet, RooArgList, RooBinning,
                  RooFit)
from ROOT import (TH1F, TFile, TCanvas, Double)
from RooFitWrappers import *

# Create a HistFunc
obj = RooObject(workspace = 'w')
w = obj.ws()

t = RealVar('t', Title = 'decay time', Unit='ps', Observable = True,  MinMax=(-3,14))


from ROOT import TFile
acceptance_file = TFile.Open('/stuff/PhD/p2vv/data/start_values.root')
if not acceptance_file:
    raise ValueError, "Cannot open histogram file %s" % input_file
hlt1_histogram = acceptance_file.Get('hlt1_shape')
if not hlt1_histogram:
    raise ValueError, 'Cannot get acceptance historgram %s from file' % hlt1_histogram
xaxis = hlt1_histogram.GetXaxis()

from array import array
## bounds = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, hlt1_histogram.GetNbinsX() + 2)))
## heights = [hlt1_histogram.GetBinContent(i) for i in range(1, hlt1_histogram.GetNbinsX() + 1)]

bounds = array('d', [-3, -1, 4, 8, 14])
heights = [0.9, 0.15, 0.45, 0.9]


heights = [RealVar('bin_%03d' % (i + 1), Observable = False, Value = v,
                   MinMax = (0.01, 0.999), Constant = True) for i, v in enumerate(heights)]

binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('acceptance')
t.setBinning(binning, 'acceptance')

## Use EffHistProd to generate events
t_flat = UniformPdf('t_flat', Arguments = [t])
acceptance = BinnedPdf(Name = 'time_acceptance', Observable = t, Coefficients = heights,
                       Binning = binning)
acceptance.setForceUnitIntegral(True)

from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
tres = TimeResolution(time = t)

eff_model = EffResModel(Efficiency = acceptance, ResolutionModel = tres.model())

from P2VVParameterizations.TimePDFs import Single_Exponent_Time
sig_t = Single_Exponent_Time(Name = 'sig_t', time = t, resolutionModel = eff_model)

pdf = sig_t.pdf()
pdf.Print('t')

data = pdf.generate([t], 100000)

result = pdf.fitTo(data, Save = True, Minimizer = 'Minuit2', Verbose = True)

frame = t.frame()
canvas = TCanvas('canvas', 'canvas', 500, 500)
canvas.Divide(1, 1)
canvas.cd(1)
data.plotOn(frame, RooFit.Binning(200))
pdf.plotOn(frame)
frame.Draw()

