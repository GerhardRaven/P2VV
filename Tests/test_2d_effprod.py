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

m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5480))
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
time_pdf = TimePDF(Name = 'time_pdf', time = t, resolutionModel = tres.model()).pdf()

mass_mean  = RealVar('mass_mean',   Unit = 'MeV', Value = 5370, MinMax = (5320, 5420))
mass_sigma = RealVar('mass_sigma',  Unit = 'MeV', Value = 13, MinMax = (5, 20))
mass_pdf = Pdf(Name = 'mass_pdf', Type = RooGaussian, Parameters = [m, mass_mean, mass_sigma])

## mass_c = RealVar( 'mass_c',  Unit = '1/MeV', Value = -0.0065, MinMax = (-0.1, -0.0000001))
## bkg_mass = Pdf(Name = 'bkg_mass',  Type = Exponential, Parameters = [m, psi_c])

# Create signal component
signal = Component('signal', (mass_pdf, time_pdf), Yield = (24500,10000,30000))

pdf = buildPdf(Components = (signal,), Observables = (m, t), Name='pdf')

#####################
### Generate data ###
#####################
data = pdf.generate([t, m],10000)

##############################################
### Define acceptance function a la Wouter ###
##############################################
nbins = 10
heights = []
for i in range(1,6):
    heights.append(RealVar('bin_%i' % i, Value = 0.1 * i, MinMax = (0.01, 0.9), Constant = True))
for i in range(6,nbins+1):
    heights.append(RealVar('bin_%i' % i, Value = 0.9, MinMax = (0.1, 0.95), Constant = True))

binned_pdf = BinnedPdf(Name = 'eff', Observable = t,
                       Binning = 'efficiency', Coefficients = heights)
time_acc_pdf = EffProd(Name = 'acc_pdf', Original = time_pdf, Efficiency = binned_pdf)

signal_1 = Component('signal_1', (mass_pdf, time_acc_pdf), Yield = (24500,10000,30000))
pdf_1 = buildPdf(Components = (signal_1,), Observables = (m, t), Name='pdf_1')

pdf_2 = EffProd(Name = 'pdf_2', Original = pdf, Efficiency = binned_pdf)

##############################
### Generate 'biased' data ###
##############################
good_accdata = pdf_1.generate([t, m], 10000)

bad_accdata = pdf_2.generate([t, m], 10000)

pdf_1.fitTo(good_accdata)
pdf_2.fitTo(bad_accdata)

############
### Plot ###
############
#Make Sanity plots
Acc = TCanvas('Acc','Acc')
Acc.Divide(2)

Acc.cd(1)
tframe = t.frame(Bins = 30)
binned_pdf.plotOn(tframe)
tframe.Draw()

Acc.cd(2)
mframe = m.frame(Bins = 30)
binned_pdf.plotOn(mframe)
mframe.Draw()

C2 = TCanvas('C2','C2')
C2.Divide(3,3)

C2.cd(1)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
data.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
tframe.Draw()

C2.cd(2)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
pdf.plotOn(tframe, LineWidth = 2, LineColor = kGreen)
tframe.Draw()

C2.cd(3)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
data.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
pdf.plotOn(tframe, LineWidth = 2, LineColor = kGreen)
tframe.Draw()

C2.cd(4)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
good_accdata.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
tframe.Draw()

C2.cd(5)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
pdf_1.plotOn(tframe, LineWidth = 2, Project = RooArgSet(m))
tframe.Draw()

C2.cd(6)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
good_accdata.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
pdf_1.plotOn(tframe, LineWidth = 2)
tframe.Draw()

C2.cd(7)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
bad_accdata.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
tframe.Draw()

C2.cd(8)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
#pdf_2.plotOn(tframe,LineWidth = 2,RooFit.Project(RooArgSet(ws['m'])))
pdf_2.plotOn(tframe, LineWidth = 2, LineColor = kRed)
tframe.Draw()

C2.cd(9)
#gPad.SetLogy()
tframe = t.frame(Bins = 10)
bad_accdata.plotOn(tframe, MarkerSize = 0.5, XErrorSize = 0)
pdf_2.plotOn(tframe, LineWidth = 2, LineColor = kRed)
tframe.Draw()

#Make Sanity plots
C3 = TCanvas('C3','C3')
C3.Divide(3,3)

C3.cd(1)
mframe = m.frame(Bins = 10)
data.plotOn(mframe, MarkerSize = 0.5, XErrorSize = 0)
mframe.Draw()

C3.cd(2)
mframe = m.frame(Bins = 10)
pdf.plotOn(mframe, LineWidth = 2, LineColor = kGreen)
mframe.Draw()

C3.cd(3)
mframe = m.frame(Bins = 10)
data.plotOn(mframe, MarkerSize = 0.5, XErrorSize = 0)
pdf.plotOn(mframe, LineWidth = 2, LineColor = kGreen)
mframe.Draw()

C3.cd(4)
mframe = m.frame(Bins = 10)
good_accdata.plotOn(mframe, MarkerSize = 0.5, XErrorSize = 0)
mframe.Draw()

C3.cd(5)
mframe = m.frame(Bins = 10)
pdf_1.plotOn(mframe, LineWidth = 2)
mframe.Draw()

C3.cd(6)
mframe = m.frame(Bins = 10)
good_accdata.plotOn(mframe, MarkerSize = 0.5, XErrorSize = 0)
pdf_1.plotOn(mframe, LineWidth = 2)
mframe.Draw()

C3.cd(7)
mframe = m.frame(Bins = 10)
bad_accdata.plotOn(mframe, MarkerSize = 0.5, XErrorSize = 0)
mframe.Draw()

C3.cd(8)
mframe = m.frame(Bins = 10)
pdf_2.plotOn(mframe, LineWidth = 2, LineColor = kRed)
mframe.Draw()

C3.cd(9)
mframe = m.frame(Bins = 10)
bad_accdata.plotOn(mframe, MarkerSize = 0.5, XErrorSize = 0)
pdf_2.plotOn(mframe, LineWidth = 2, LineColor = kRed)
mframe.Draw()
