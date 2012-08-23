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
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj  = RooObject( workspace = 'workspace')

t = RealVar('time', Title = 'Decay time', Unit = 'ps',  Observable = True, MinMax = (0.3, 14))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution
res_model = Moriond2012_TimeResolution(time = t, timeResSFConstraint = True, sigmat = st,
                                       timeResSF =  dict(Value = 1.46, MinMax = ( 0.5, 5. ),
                                                         Constant = True))

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as TimePDF
time_pdf = TimePDF(Name = 'time_pdf', time = t, resolutionModel = res_model.model())
time_pdf = time_pdf.pdf()

## time_pdf = UniformPdf(Name = 'time_pdf', Arguments = [t])
hlt1_excl_biased = Category('hlt1_excl_biased', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)

## Selection
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, st, hlt1_biased, hlt1_unbiased, hlt1_excl_biased,
               hlt2_biased, hlt2_unbiased, selected]

categories = [hlt1_excl_biased, hlt2_biased, hlt2_unbiased, hlt1_biased, hlt1_unbiased]
categories = dict([(c.GetName(), c) for c in categories])

# Binnings
from array import array
biased_bins = array('d', [0, 1, 2, 3.5, 5, 7.5, 10])
biased_heights = [0.1, 0.2, 0.3, 0.4, 0.5, 0.9]

excl_biased_heights = [0.2, 0.3, 0.3, 0.3, 0.8, 0.8]

unbiased_bins = array('d', [0, 10])
unbiased_heights = [0.5]

# Spec to build efficiency shapes
spec = {'Bins' : {hlt1_excl_biased : {'excl_biased' : {'bins'    : biased_bins,
                                                       'heights' : excl_biased_heights,
                                                       'average' : (6.285e-01, 1.633e-02)},
                                      'unbiased'    : {'bins'    : unbiased_bins,
                                                       'heights' : unbiased_heights}
                                      },
                  hlt2_biased      : {'biased'      : {'bins'    : biased_bins,
                                                       'heights' : biased_heights,
                                                       'average' : (6.33e-01, 1.65e-02)}
                                 },
                  hlt2_unbiased    : {'unbiased'    : {'bins'    : unbiased_bins,
                                                       'heights' : unbiased_heights}
                                      }
                  }
        }
rel_spec = {(('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.078,
            (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.027,
            (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : 0.383,
            (('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'not_biased'), ('hlt2_unbiased', 'unbiased')) : 0.01,
            (('hlt1_excl_biased', 'unbiased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'not_unbiased')) : 0.433,
            (('hlt1_excl_biased', 'excl_biased'), ('hlt2_biased', 'biased'), ('hlt2_unbiased', 'unbiased')) : None}
spec['Relative'] = dict([(tuple((categories[c], l) for c, l in k), {'Constant' : True, 'Value' : v} if v else None) for k, v in rel_spec.iteritems()])

## Data
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_2011_biased_unbiased.root'
data = readData(input_file, tree_name, cuts = 'sel == 1 && (hlt1_biased == 1 || hlt1_unbiased == 1) && (hlt2_biased == 1 || hlt2_unbiased == 1)',
                NTuple = True, observables = observables)

## mhe = MultiHistEfficiencyModel(Name = "RMERM", Original = time_pdf, Observable = t,
##                                ConditionalCategories = True, ResolutionModel = res_model.model(),
##                                FitAcceptance = True, UseSingleBinConstraint = False,
##                                ConditionalObservables = res_model.conditionalObservables(),
##                                ExternalConstraints = res_model.externalConstraints(),
##                                **spec)

from P2VVParameterizations.GeneralUtils import valid_combinations
valid_definition = [[(hlt1_excl_biased, 'excl_biased'), (hlt1_excl_biased, 'unbiased')],
                    [(hlt2_biased, 'biased'), (hlt2_unbiased, 'unbiased')]]
valid = valid_combinations(valid_definition)

hists = { hlt1_excl_biased : {  'excl_biased' : { 'histogram' : 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1ExclB_20bins' }
                                , 'unbiased'    : { 'histogram' : 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_20bins'    }
                                }
          }
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance as TimeAcceptance
acceptance = TimeAcceptance(time = time, Input = '/stuff/PhD/p2vv/data/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root',
                            Histogram = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_20bins',
                            Original = time_pdf, ResolutionModel = res_model)
