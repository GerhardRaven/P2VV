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
time_pdf = TimePDF(Name = 'time_pdf', time = t, resolutionModel = res_model.model())
time_pdf = time_pdf.pdf()

unbiased = Category('unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
biased = Category('biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)

# Binnings
from array import array
biased_bins = array('d', [0, 2, 5, 10])
biased_heights = [0.1, 0.5, 0.9]

unbiased_bins = array('d', [0, 10])
unbiased_heights = [0.5]

# Spec to build efficiency shapes
spec = {"Bins" : {biased : {'state'   : 'biased',
                            'bounds'  : biased_bins,
                            'heights' : biased_heights},
                  unbiased : {'state' : 'unbiased',
                              'bounds' : unbiased_bins,
                              'heights' : unbiased_heights}
                  },
        "Relative" : {((biased, "biased"),     (unbiased, "unbiased")) : {'Value' : 0.2, 'MinMax' : (0.1, 0.45)},
                      ((biased, "not_biased"), (unbiased, "unbiased")) : {'Value' : 0.2, 'MinMax' : (0.1, 0.45)},
                      ((biased, "biased"),     (unbiased, "not_unbiased")) : None}
        }
mhe = MultiHistEfficiency(Name = "RMHE", Original = time_pdf, Observable = t, ConditionalCategories = True, **spec)
print mhe.getVal()

data = mhe.generate([t, biased, unbiased], 20000)

fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Verbose = True, Minimizer = 'Minuit2', Optimize = 2)
mhe.fitTo(data, **fitOpts)

f = t.frame()
data.plotOn(f, Binning = 100)
f.Draw()

