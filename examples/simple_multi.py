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

t = RealVar('time', Title = 'decay time', Unit = 'ps',  Observable = True, MinMax = (0, 2), nBins = 10)

original = UniformPdf('uniform', Arguments = [t])

unbiased = Category('unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
biased = Category('biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)

# Binnings
from array import array
biased_bins = array('d', [0, 1, 2])
biased_heights = [0.2, 0.8]

unbiased_bins = array('d', [0, 1, 2])
unbiased_heights = [0.4, 0.2]

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
data = mhe.generate([t, biased, unbiased], 10000)

## fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Verbose = True, Minimizer = 'Minuit2', Optimize = 2)
## mhe.fitTo(data, **fitOpts)

f = t.frame()
data.plotOn(f, Binning = 100)
mhe.plotOn(f,  ProjWData = (RooArgSet(biased, unbiased), data))
f.Draw()

