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

obj  = RooObject( workspace = 'workspace')

t = RealVar('time', Title = 'decay time', Unit = 'ps',  Observable = True, MinMax = (0, 10), nBins = 10)

from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
res_model = TimeResolution(time = t)

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as TimePDF
original = TimePDF(Name = 'time_pdf', time = t, resolutionModel = res_model.model())
original = original.pdf()
## original = UniformPdf('uniform', Arguments = [t])

unbiased = Category('unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
biased = Category('biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)

# Binnings
from array import array
biased_bins = array('d', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
biased_heights = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9, 0.9, 0.9, 0.9]

unbiased_bins = array('d', [0, 10])
unbiased_heights = [0.2]

# Spec to build efficiency shapes
spec = {"Bins" : {biased : {'state'   : 'biased',
                            'bounds'  : biased_bins,
                            'heights' : biased_heights},
                  unbiased : {'state' : 'unbiased',
                              'bounds' : unbiased_bins,
                              'heights' : unbiased_heights}
                  },
        "Relative" : {((biased, "biased"),     (unbiased, "unbiased")) : {'Value' : 0.2, 'MinMax' : (0.1, 0.45)},
                      ((biased, "not_biased"), (unbiased, "unbiased")) : {'Value' : 0.3, 'MinMax' : (0.1, 0.45)},
                      ((biased, "biased"),     (unbiased, "not_unbiased")) : None}
        }
mhe = MultiHistEfficiency(Name = "RMHE", Original = original, Observable = t,
                          ConditionalCategories = True, **spec)

entries = mhe.getEntries()
superCat = mhe.getSuper()

data = mhe.generate([t, biased, unbiased], 100000)
data.addColumn(superCat)
## superCat = obj.ws().cat(superCat.GetName())

states = {}
for entry in entries:
    label = superCat.lookupType(entry.first).GetName()
    states[label] = entry.second.effProd()

pdf = SimultaneousPdf(Name = 'pdf', States = states, SplitCategory = superCat)


fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Verbose = True, Minimizer = 'Minuit2', Optimize = 2)
pdf.fitTo(data, **fitOpts)

c = TCanvas("canvas", "canvas", 1000, 1000)
c.Divide(2, 2)
c.cd(1)
f = t.frame()
data.plotOn(f, Binning = 100)
pdf.plotOn(f,  ProjWData = (RooArgSet(biased, unbiased), data))
f.Draw()

c.cd(2)
data_only_unbiased = data.reduce("biased == biased::not_biased && unbiased == unbiased::unbiased")
f = t.frame()
data_only_unbiased.plotOn(f, Binning = 100)
pdf.plotOn(f,  ProjWData = (RooArgSet(biased, unbiased), data_only_unbiased))
f.Draw()

c.cd(3)
data_only_biased = data.reduce("biased == biased::biased && unbiased == unbiased::not_unbiased")
f = t.frame()
data_only_biased.plotOn(f, Binning = 100)
pdf.plotOn(f,  ProjWData = (RooArgSet(biased, unbiased), data_only_biased))
f.Draw()

c.cd(4)
data_both = data.reduce("biased == biased::biased && unbiased == unbiased::unbiased")
f = t.frame()
data_both.plotOn(f, Binning = 100)
pdf.plotOn(f,  ProjWData = (RooArgSet(biased, unbiased), data_both))
f.Draw()
