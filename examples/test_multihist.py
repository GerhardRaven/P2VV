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

unbiased = Category('unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0})
biased = Category('biased', States = {'biased' : 1, 'not_biased' : 0})

counter = 0

# test function to build efficiencies
def build_shapes(binning_name = "efficiency", **kwargs):

    pdf = kwargs.pop('PDF')
    bins = kwargs.pop('Bins')
    relative = kwargs.pop('Relative')

    efficiency_entries = std.vector("MultiHistEntry")()

    coefficients = {}
    base_binning = None
    for (category, entries) in bins.iteritems():
        bounds = entries['bounds']
        heights = entries['heights']
        state = entries['state']
        heights = [RealVar('%s_%s_bin_%03d' % (category.GetName(), state, i + 1),
                           Observable = False, Value = v,
                           MinMax = (0.001, 1)) for i, v in enumerate(heights)]
        coefficients[category] = (bounds, heights)
        if not base_binning or len(bounds) > len(coefficients[base_binning][0]):
            base_binning = category

    # Set the binning on the observable
    base_bounds = coefficients[base_binning][0]
    obs_binning = RooBinning(len(base_bounds) - 1, base_bounds)
    t.setBinning(obs_binning, 'efficiency_binning')

    def make_map(d):
        m = var_map()
        for k, v in d.iteritems():
            p = var_pair(k._target_() if hasattr(k, '_target_') else k, v)
            m.insert(p)
        return m
    
    def find_coefficient(val, bounds, coefficients):
        for i in range(len(bounds) - 1):
            if val > bounds[i] and val < bounds[i + 1]:
                break
        else:
            raise RuntimeError;
        return coefficients[i]

    for categories, relative_efficiency in relative.iteritems():
        heights = []
        bin_vars = [{} for i in range(len(base_bounds) - 1)]
        prefix = '__'.join(['%s_%s' % (c.GetName(), s) for c, s in categories])
        for category, state in categories:
            category_bounds = coefficients[category][0]
            category_heights = coefficients[category][1]
            for i in range(len(base_bounds) - 1):
                val = (base_bounds[i] + base_bounds[i + 1]) / 2
                coefficient = find_coefficient(val, category_bounds, category_heights)
                bin_vars[i][coefficient._target_()] = (state == bins[category]['state'])
        for i, d in enumerate(bin_vars):
            cm = make_map(d)
            name = '%s_%d' % (prefix, i)
            heights.append(RooEfficiencyBin(name, name, cm))

        # Make realvars for relative efficiencies
        efficiency = RealVar('%s_efficiency' % category.GetName(), Observable = False,
                                Value = relative_efficiency, MinMax = (0.001, 0.999))

        binned_pdf = BinnedPdf(Name = '%s_shape' % prefix, Observable = t, Binning = 'efficiency_binning',
                               Coefficients = heights)
        eff_prod = EffProd('%s_efficiency' % prefix, Original = pdf, Efficiency = binned_pdf)

        # MultiHistEntry
        cm = category_map()
        for category, state in categories:
            # cp = category_pair(category._target_(), state)
            cm[category._target_()] = state
        entry = MultiHistEntry(cm, eff_prod._target_(), binned_pdf._target_())
        efficiency_entries.push_back(entry)

    return RooMultiHistEfficiency("RMHE", "RHME", efficiency_entries)


# Binnings
from array import array
biased_bins = array('d', [0, 2, 5, 10])
biased_heights = [0.001, 0.5, 1]

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
        "Relative" : {((biased, "biased"),     (unbiased, "unbiased")) : 0.2,
                      ((biased, "not_biased"), (unbiased, "unbiased")) : 0.2,
                      ((biased, "biased"),     (unbiased, "not_unbiased")) : 0.6}
        }
mhe = MultiHistEfficiency(Name = "RMHE", Original = time_pdf, Observable = t, **spec)
print mhe.getVal()

data = mhe.generate([t, biased, unbiased], 20000)

f = t.frame()
data.plotOn(f)
f.Draw()
