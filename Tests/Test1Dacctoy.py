from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *
from RooFitWrappers import *
from P2VVLoad import MultiCatGen

from ROOT import MultiHistEntry
category_map = std.map('RooAbsCategory*', 'string')
category_pair = std.pair('RooAbsCategory*', 'string')
var_map = std.map('RooAbsReal*', 'bool')
var_pair = std.pair('RooAbsReal*', 'bool')

RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Integration))
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj  = RooObject( workspace = 'workspace')

t = RealVar('time', Title = 'decay time', Unit = 'ps',  Observable = True, MinMax = (0, 10), nBins = 10)


from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
res_model = TimeResolution(time = t)

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as TimePDF
time_pdf = TimePDF(Name = 'time_pdf', time = t, resolutionModel = res_model.model())
time_pdf = time_pdf.pdf()

## time_pdf = UniformPdf(Name = 'time_pdf', Arguments = [t])


#####################
### Generate data ###
#####################
data = time_pdf.generate([t],10000)

##############################################
### Define acceptance function a la Wouter ###
##############################################
nbins = 10
heights = []
for i in range(nbins):
    name = 'bin_%d' % (i + 1)
    h = 0.1 * (i + 1) if i < 6 else 1.
    heights.append(RooConstVar(name, name, h))
from array import array
bins = array('d', [i for i in range(11)])
from ROOT import RooBinning
binning = RooBinning(len(bins) - 1, bins)
t.setBinning(binning, 'acceptance_binning')

binned_pdf = BinnedPdf(Name = 'acceptance_shape', Observable = t,
                       Binning = binning.GetName(), Coefficients = heights)
pdf = binned_pdf * time_pdf

##############################
### Generate 'biased' data ###
##############################
acc_data = pdf.generate([t],10000)

############
### Plot ###
############

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
data.plotOn(tframe,MarkerSize = 0.5,XErrorSize = 0)
tframe.Draw()

C2.cd(2)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
time_pdf.plotOn(tframe,LineWidth = 2)
tframe.Draw()

C2.cd(3)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
data.plotOn(tframe,MarkerSize = 0.5,XErrorSize = 0)
time_pdf.plotOn(tframe,LineWidth = 2)
tframe.Draw()

C2.cd(4)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
acc_data.plotOn(tframe,MarkerSize = 0.5,XErrorSize = 0)
tframe.Draw()

#effhistprod = ws['accpdf1']
#effhistprod.getAnalyticalIntegralWN(RooArgSet(ws['t']),RooArgSet(ws['t']),RooArgSet(ws['t']))

#print '*******************************'
#fullint = pdf.analyticalIntegralWN(1,0,'')
#print 'full integral pdf = ', fullint
#print '*******************************'

## t.setRange('MyPythonRange',0,1)

## print '*******************************'
## rangeint = pdf.analyticalIntegralWN(1,0,'MyPythonRange')
## print 'range integral pdf = ', rangeint

## print '*******************************'
## normint = pdf.analyticalIntegralWN(1,RooArgSet(ws['t']),'MyPythonRange')
## print 'normintegral pdf = ', normint

#print '*******************************'
#fullint = pdf.analyticalIntegralWN(1,0,'')
#print 'full integral accpdf = ', fullint
#print '*******************************'

## print '*******************************'
## rangeint = pdf.analyticalIntegralWN(1,0,'MyPythonRange')
## print 'range integral accpdf = ', rangeint

## print '*******************************'
## normint = pdf.analyticalIntegralWN(1,RooArgSet(ws['t']),'MyPythonRange')
## print 'normintegral accpdf = ', normint

C2.cd(5)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
pdf.plotOn(tframe,LineWidth = 2)
tframe.Draw()

C2.cd(6)
#gPad.SetLogy()
tframe = t.frame(Bins = 30)
acc_data.plotOn(tframe,MarkerSize = 0.5,XErrorSize = 0)
pdf.plotOn(tframe,LineWidth = 2)
tframe.Draw()

## npoints = 1000
## dt = (t.getMax()-t.getMin())/npoints 
## sum = 0
## for i in range(0,npoints):
##     x = t.getMin() + i*dt 
##     t.setVal( x ) 
##     pdfval = pdf.getVal(RooArgSet(ws['t']) )
##     print "pdfval: ", pdfval
##     sum += dt * pdfval 
## print "Integral: ", sum 
