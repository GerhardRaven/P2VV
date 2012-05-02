from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *
from RooFitWrappers import *
from P2VVLoad import MultiCatGen

RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))

obj  = RooObject( workspace = 'workspace')

#########################################
### Define variables and simple PDF's ###
#########################################
t = RealVar('time', Title = 'decay time', Unit = 'ps',  Observable = True, MinMax = (0, 10), nBins = 10)

from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
res_model = TimeResolution(time = t)

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as TimePDF
time_pdf = TimePDF(Name = 'pdf', time = t, resolutionModel = res_model.model())
time_pdf = time_pdf.pdf()

## data_unbiased = time_pdf.generate([t], 3000)

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
nbins = 10
eff_hist = TH1F("eff_hist","eff_hist", nbins, t.getMin(), t.getMax()) 
for i in range(1, 6):
    eff_hist.SetBinContent(i, 0.1 * i)
for i in range(6, nbins + 1):
    eff_hist.SetBinContent(i,1) 

t.setBins(nbins, 'default')
hist_func = HistFunc(Name = "eff_pdf", Observables = [t], Histogram = eff_hist)
gen_pdf = EffProd(Name = 'gen_pdf', Original = time_pdf, Efficiency = hist_func)

## data_biased = gen_pdf.generate([t], 10000)
## add_categories(data_biased, {biased : 1, unbiased :  0})

## data_both = gen_pdf.generate([t], 2000)
## add_categories(data_both, {biased : 1, unbiased :  1})

## data = data_both.Clone('data_full')
## data.append(data_unbiased)
## data.append(data_biased)

bin_vars = [RealVar('bin_%03d' % (i + 1), Observable = True, Value = 0.5, MinMax = (0.01, 1)) for i in range(nbins)]
bin_vars[-1].setConstant()
binned_pdf = BinnedPdf(Name = 'binned_pdf', Observable = t, Binning = 'default',
                       Coefficients = bin_vars)
prescale = RealVar('prescale', Observable = False, Value = 0.2, MinMax = (0.001, 0.999))

efficiency = MultiEfficiency(Name = 'efficiency', ConditionalCategories = True,
                             Efficiencies = {prescale   : (unbiased, 'unbiased'),
                                             binned_pdf : (biased,   'biased'  )})


fit_pdf = EffProd(Name = 'fit_pdf', Efficiency = efficiency, Original = time_pdf)

fit_pdf.Print('t')

data = fit_pdf.generate([t, biased, unbiased], 20000)
data.table(biased).Print()
## fitOpts = dict(Timer = 1, NumCPU = 1, Save = True, Minimizer = ('Minuit2','minimize'), Verbose = True)
## fit_pdf.fitTo(data, **fitOpts)

## canvas = TCanvas('canvas', 'canvas', 1000, 1000)
## canvas.Divide(2, 2)
## canvas.cd(1)
## frame = t.frame()
## data.plotOn(frame, Binning = 100)
## gen_pdf.plotOn(frame)
## frame.Draw()

## canvas.cd(3)
## frame = t.frame()
## data.plotOn(frame, Binning = 100)
## fit_pdf.plotOn(frame, Slice = (biased, 'not_biased'))
## frame.Draw()

## canvas.cd(4)
## frame = t.frame()
## data.plotOn(frame, Binning = 100)
## fit_pdf.plotOn(frame, Slice = (unbiased, 'not_unbiased'))
## frame.Draw()

