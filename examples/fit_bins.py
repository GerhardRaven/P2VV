from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *
from RooFitWrappers import *
from P2VVLoad import MultiCatGen

RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Integration))

obj  = RooObject( workspace = 'workspace')

#########################################
### Define variables and simple PDF's ###
#########################################
t = RealVar('time', Title = 'decay time', Unit = 'ps',  Observable = True, MinMax = (0, 10), nBins = 10)

from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
res_model = TimeResolution(time = t)

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as TimePDF
time_pdf = TimePDF(Name = 'time_pdf', time = t, resolutionModel = res_model.model())
time_pdf = time_pdf.pdf()

unbiased = RooCategory('unbiased', 'unbiased')
unbiased.defineType('not_unbiased', 0)
unbiased.defineType('unbiased', 1)

biased = RooCategory('biased', 'biased')
biased.defineType('not_biased', 0)
biased.defineType('biased', 1)

def add_categories(data, categories):
    for category, index in categories.iteritems():
        p = category.getIndex()
        category.setIndex(index)
        data.addColumn(category)
        category.setIndex(p)
    return data

## add_categories(data_unbiased, {biased : 0, unbiased : 1})

# t_sig_tau = dict(Title = 'lifetime', Value = 1.5, Constant = True)
##############################################
### Define acceptance function a la Wouter ###
##############################################

from array import array
boundaries = array('d', [0, 2, 5, 10])
binning = RooBinning(len(boundaries) - 1, boundaries)
t.setBinning(binning, 'default')

bin_heights = [0.001, 0.5, 1]
bin_vars = [RealVar('bin_%03d' % (i + 1), Observable = True, Value = v, MinMax = (0.001, 1)) for i, v in enumerate(bin_heights)]
bin_vars[-1].setConstant()
binned_pdf = BinnedPdf(Name = 'binned_pdf', Observable = t, Binning = 'default',
                       Coefficients = bin_vars)
prescale = RealVar('prescale', Observable = False, Value = 0.2, MinMax = (0.001, 0.999))

efficiency = MultiHistEfficiency(Name = 'efficiency', ConditionalCategories = True,
                                 Efficiencies = {prescale   : (unbiased, 'unbiased'),
                                                 binned_pdf : (biased,   'biased'  )})
pdf = EffProd('pdf', Original = time_pdf, Efficiency = efficiency)
pdf.Print('t')

data = pdf.generate([t, biased, unbiased], 10000)
data.table(biased).Print('v')
data.table(unbiased).Print('v')

## fitOpts = dict(Timer = 1, NumCPU = 1, Save = True, Minimizer = ('Minuit2','minimize'), Verbose = True,
##                Optimize = 1)
## pdf.fitTo(data, **fitOpts)

canvas = TCanvas('canvas', 'canvas', 1000, 500)
canvas.Divide(2, 1)
## canvas.cd(1)
## frame = t.frame()
## data.plotOn(frame, Binning = 100, Cut = "biased == 1 && unbiased == 0")
## pdf.plotOn(frame, Slices = ((biased, '0'), (unbiased, '1')))
## frame.Draw()

canvas.cd(2)
frame = t.frame()
data.plotOn(frame, Binning = 100)
pdf.plotOn(frame)
frame.Draw()

