from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))

widths = [0.04, 0.1, 0.5]
yields = [15000, 7500, 2000]

mean = RealVar(Name = 'mean', Value = 0, MinMax = (-1, 1))
components = []

from ROOT import RooGaussian
for i, w in enumerate(widths):
    sigma = RealVar(Name = 'g%d_sigma' % i, Value = w, MinMax = (-5, 5))
    gauss = Pdf(Name = 'g%d' % i, Type = RooGaussian, Parameters = [t, mean, sigma])
    comp = Component('c%d' % i, (gauss,), Yield = (yields[i], 1, 2 * yields[i]))
    components.append(comp)

pdf = buildPdf(Components = components, Observables = (t, ), Name = 'pdf')
pdf.Print('t')
data = pdf.generate([t], 30000)

result = pdf.fitTo(data, Save = True, Minimizer = 'Minuit2')

## Dilution
from Dilution import dilution
dilution(t, data, result = result, signal = components[:2], subtract = components[-1:])

## Plot
from ROOT import TCanvas
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
canvas = TCanvas('canvas', 'canvas', 500, 500)
from P2VVGeneralUtils import plot
plot(canvas.cd(1), t, pdf = pdf, data = data
         , dataOpts = dict(MarkerSize = 0.8, Binning = 'dilution_binning', MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2)
         , logy = True
         , plotResidHist = False
         , components = { 'g0'   : dict( LineColor = kRed,    LineStyle = kDashed )
                          , 'g1' : dict( LineColor = kOrange, LineStyle = kDashed )
                          , 'g2' : dict( LineColor = kGreen,  LineStyle = kDashed )
                          }
     )
