import optparse
import sys
import os

parser = optparse.OptionParser(usage = '%prog file year')
(options, args) = parser.parse_args()

if len(args) != 2:
    print parser.usage
    sys.exit(-2)

input_filename = args[0]
data_name = 'wpv_data_%s' % args[1]

from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle

obj = RooObject( workspace = 'w')
w = obj.ws()

t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (0.5, 14))

from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution
tres = Truth_TimeResolution(time = t)

from P2VV.Parameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = tres.model()
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 6.7195e-01)
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.3672, MinMax = (0.1,  10))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 1.3405e-01, MinMax = (0.01, 5))
                         )
psi_t = psi_t.pdf()

from ROOT import TFile
input_file = TFile(input_filename)
sdata = input_file.Get(data_name)
sdata_cut = sdata.reduce('time > %f' % t.getMin())

fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True,
               Verbose = False, Strategy = 1, SumW2Error = False)

result = psi_t.fitTo(sdata_cut, **fitOpts)

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 600, 400)
frame = t.frame()
nBins = 50
sdata_cut.plotOn(frame, Binning = nBins)
psi_t.plotOn(frame)
frame.GetXaxis().SetTitle('decay time [ps]')
frame.GetYaxis().SetTitle('Candidates / (%3.2f ps)' % ((t.getMax() - t.getMin()) / float(nBins)))
frame.Draw()


from P2VV.Utilities import Resolution as ResolutionUtils
result.PrintSpecial(LaTeX = True, ParNames = ResolutionUtils.parNames)
