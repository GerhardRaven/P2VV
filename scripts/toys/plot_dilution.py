import os
import sys
import optparse

parser = optparse.OptionParser(usage = 'usage: %prog file dataset_name')

(options, args) = parser.parse_args()

if len(args) != 2:
    print parser.print_usage()
    sys.exit(-2)

input_file = args[0]
data_name = args[1]
from P2VV.Load import LHCbStyle
from P2VV.RooFitWrappers import *
from ROOT import TFile

f = TFile(input_file)
data = f.Get(data_name)

obj = RooObject( workspace = 'w')
obj.ws().put(data)

da = data.get().find('da')
dft = data.get().find('dft')
diff = data.get().find("diff")
if not diff:
    diff = FormulaVar(Name = 'diff', Formula = '@0 - @1', Arguments = (dft, da), data = data)

mean = RealVar('mean', Value = -0.007, MinMax = (-1, 1))
sigma = RealVar('sigma', Value = 0.003, MinMax = (0.0001, 0.1))
from ROOT import RooGaussian
gauss = Pdf(Name = 'gauss', Type = RooGaussian, Parameters = [diff, mean, sigma])

result = gauss.fitTo(data, Minimizer = 'Minuit2', Save = True, Optimize = 2)

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 600, 400)
frame = diff.frame(Range = (-0.03, 0.03))
data.plotOn(frame, Binning = (150, -0.03, 0.03))
gauss.plotOn(frame)
frame.GetXaxis().SetTitle('D_{FT} - D')
frame.GetYaxis().SetTitle('Candidates / (0.4 fs)')
frame.Draw()
