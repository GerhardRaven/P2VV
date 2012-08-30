from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *
from RooFitWrappers import *
import rootStyle
from ModelBuilders import _buildAngularFunction

myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

from ModelBuilders import *

obj  = RooObject( workspace = 'workspace')

RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Integration))
#RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Tracing))

t = RealVar('time', Title = 'decay time', Unit = 'ps',  Observable = True, MinMax = (0, 10), nBins = 10)
binning = RooBinning(10, 0, 10)
t.setBinning(binning, 'efficiency')
#########################################
### Define variables and simple PDF's ###
#########################################

# Time resolution model
from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
tres = TimeResolution(time = t)

# Signal time pdf
from P2VVParameterizations.TimePDFs import Single_Exponent_Time as TimePDF
pdf = TimePDF(Name = 'pdf', time = t, resolutionModel = tres.model())
pdf = pdf.pdf()

## pdf = UniformPdf(Name = 'pdf', Arguments = [t])
pdf.Print('t')
#####################
### Generate data ###
#####################
data = pdf.generate([t],10000)

##############################################
### Define acceptance function a la Wouter ###
##############################################
nbins = 10
heights = []
for i in range(1,6):
    heights.append(RealVar('bin_%i' % i, Value = 0.1 * i, MinMax = (0.01, 0.9)))
for i in range(6,nbins+1):
    heights.append(RealVar('bin_%i' % i, Value = 0.9, MinMax = (0.1, 0.9)))

binned_pdf = BinnedPdf(Name = 'eff', Observable = t,
                       Binning = 'efficiency', Coefficients = heights)
accpdf = EffProd(Name = 'acc_pdf', Original = pdf, Efficiency = binned_pdf)

##############################
### Generate 'biased' data ###
##############################
accdata1 = accpdf.generate([t],10000)

t.setVal(1)
acc_int = accpdf.createIntegral(RooArgSet(t))
print acc_int.getVal()


lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

#Make Sanity plots
Acc = TCanvas('Acc','Acc')

Acc.cd(1)
tframe = t.frame(Bins = 30)
binned_pdf.plotOn(tframe)
tframe.Draw()

C2 = TCanvas('C2','C2')
C2.Divide(3,2)

C2.cd(1)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(2)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
pdf.plotOn(tframe, LineWidth = 2)
tframe.Draw()

C2.cd(3)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
data.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
pdf.plotOn(tframe, LineWidth = 2)
tframe.Draw()

C2.cd(4)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
accdata1.plotOn(tframe,MarkerSize = 0.5, XErrorSize = 0)
tframe.Draw()

C2.cd(5)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
accpdf.plotOn(tframe,LineWidth = 2)
tframe.Draw()

C2.cd(6)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
accdata1.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
accpdf.plotOn(tframe,LineWidth = 2)
tframe.Draw()
