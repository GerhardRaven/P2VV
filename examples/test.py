from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))

input_file = '/stuff/PhD/p2vv/data/start_values.root'
histogram = 'hlt2_shape'

from ROOT import TFile
acceptance_file = TFile.Open(input_file)
if not acceptance_file:
    raise ValueError, "Cannot open histogram file %s" % input_file
histogram = acceptance_file.Get(histogram)
if not histogram:
    raise ValueError, 'Cannot get acceptance historgram %s from file' % histogram
xaxis = histogram.GetXaxis()

from array import array
bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX() + 2)))
heights = [histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)]
heights = [RealVar('bin_%03d' % (i + 1), Observable = False, Value = v,
                   MinMax = (0.01, 0.999)) for i, v in enumerate(heights)]

binning_name = 'efficiency_binning'
from ROOT import RooBinning
binning = RooBinning(len(bins) - 1, bins)
t.setBinning(binning, binning_name)

binned_pdf = BinnedPdf(Name = 'efficiency_shape', Observable = t,
                       Binning = binning_name, Coefficients = heights)

data = binned_pdf.generate([t], 10000)

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 500, 500)
frame = t.frame()
data.plotOn(frame)
binned_pdf.plotOn(frame)
frame.Draw()
