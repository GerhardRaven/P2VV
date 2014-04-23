import os
import sys
import optparse
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog year')
(options, args) = parser.parse_args()

prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
input_data = {'Combined' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats.root'),
              '2012' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_s20r0p1_dv33r6p1_20131107_tupleB_add.root'),
              'MC2012' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_ntupleB_20130904_add.root')}

if len(args) != 1 or args[0] not in input_data.keys() + ['generate']:
    print parser.print_usage()
    print "Possible samples are: %s" % ' '.join(input_data.keys())
    sys.exit(-2)

ntuple_file = None
if args[0] in input_data.keys():
    ntuple_file = input_data[args[0]]

from P2VV.Load import P2VVLibrary
from P2VV.Load import RooFitOutput
from P2VV.RooFitWrappers import *

from itertools import product
from ROOT import RooCBShape as CrystalBall
from P2VV.Parameterizations.GeneralUtils import valid_combinations
#from P2VV.Load import RooFitOutput

from ROOT import RooMsgService
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Integration))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))
m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550))
nPV = RealVar('nPV', Title = 'nPV', Observable = True, MinMax = (0, 15))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

# Categories
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt1_excl_biased_dec = Category('hlt1_excl_biased_dec', States = {'exclB' : 1, 'notExclB' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'B' : 1, 'notB' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'UB' : 1, 'notUB' : 0}, Observable = True)
hlt2_excl_biased = Category('hlt2_excl_biased', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)
run_period = Category('runPeriod', States = {'p2011' : 2011, 'p2012' : 2012}, Observable = True)
## project_vars = [hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]
categories = [hlt1_biased, hlt1_unbiased, hlt1_excl_biased_dec,
              hlt2_biased, hlt2_unbiased, hlt2_excl_biased, run_period]
categories = dict([(c.GetName(), c) for c in categories])

project_vars = [hlt1_excl_biased_dec, hlt1_unbiased, hlt2_biased, hlt2_unbiased, run_period, st]

selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, hlt1_biased, hlt1_unbiased, hlt1_excl_biased_dec,
               hlt2_biased, hlt2_unbiased, hlt2_excl_biased, selected, nPV, run_period]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
tres = TimeResolution(Name = 'resolution', time = t, BiasScaleFactor = False,
                      timeResSigma = dict(Value = 0.045, Constant = True),
                      timeResMu = dict(Value = 0, Constant = True), Cache = False)

# B mass pdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5365, MinMax = (5363,5372)))

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp = dict( Name = 'm_bkg_exp' ) )

# Create components
signal_mass = Component('signal', (sig_m.pdf(),), Yield = (30000,100,100000))
background_mass = Component('background', (bkg_m.pdf(),), Yield = (100000,100,300000) )

## Build mass PDF
mass_pdf = buildPdf(Components = (signal_mass, background_mass), Observables = (m, ), Name='mass_pdf')
mass_pdf.Print("t")

## base_location = '/home/raaij'
base_location = '/stuff/PhD/p2vv'

# Build the acceptance using the histogram as starting values
input_file = os.path.join(base_location, 'data/start_values.root')

## hists = {hlt1_excl_biased_dec : {'excl_biased' : {'histogram' : 'hlt1_shape', 'average' : (6.285e-01, 1.633e-02)},
##                                  'unbiased' : { 'bins' : t.getRange(), 'heights' : [0.5]}}}
## hists = {hlt1_excl_biased_dec : {'exclB' : {'histogram' : 'hlt1_shape'},
##                                  'notExclB' : { 'bins' : t.getRange(), 'heights' : [0.7]}},
##          hlt2_biased : { 'B' : {'histogram' : 'hlt2_shape'}},
##          hlt2_unbiased : { 'UB' : {'histogram' : 'hlt2_shape'}}}
hists = {hlt1_excl_biased_dec : {'exclB' : {'histogram' : 'hlt1_shape'},
                                 'notExclB' : { 'bins' : t.getRange(), 'heights' : [0.7]}},
         hlt2_biased : { 'B' : {'histogram' : 'hlt2_shape'}},
         hlt2_unbiased : { 'UB' : {'bins' : t.getRange(), 'heights' : [0.7]}}}

from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time
pdf_wrapper = Single_Exponent_Time(Name = 'pdf', time = t, resolutionModel = tres.model())
tau = pdf_wrapper._tau
pdf = pdf_wrapper.pdf()
split_cats = [[run_period]]
from itertools import chain
split_pars = [[tau]]
sim_pdf = SimultaneousPdf(pdf.GetName() + '_simul', MasterPdf = pdf,
                          SplitCategories = split_cats, SplitParameters = split_pars)


from ROOT import RooArgSet, RooArgList
split_cat_pars = RooArgSet()
sim_pdf.treeNodeServerList(split_cat_pars)
split_cat = sim_pdf.indexCat()
input_cats = RooArgList(split_cat) if split_cat.isFundamental() else split_cat.inputCatList()
prototype = tres.model()
from P2VV.RooFitWrappers import AddModel
## NOTE, this only works for either a single resolution model or an AddModel which does
## not contain further add models
if isinstance(prototype, AddModel):
    orig_res_params = {}
    for model in prototype.models():
        assert(not isinstance(model, AddModel))
        orig_res_params[model.GetName()] = model['Parameters']
    orig_res_params[prototype.GetName()] = prototype.fractions()
else:
    orig_res_params = {prototype.GetName() : prototype['Parameters']}

split_models = {}
from P2VV.Utilities.Splitting import replacement_model
# Figure out which way the time resolution model is to be split
res_params = set(var.GetName() for var in prototype.getVariables())
split_res_cats = [run_period]
replacements = {( ) : prototype}
for split_state in split_cat:
    comp_pdf = sim_pdf.getPdf(split_state.GetName())
    res_model = comp_pdf.getComponents().find(tres.model().GetName())
    split_models[split_state.GetName()] = replacement_model(prototype, res_model, input_cats,
                                                            {tau._target_() : set([run_period])},
                                                            split_cat_pars, orig_res_params)

## Multiply by acceptances
from P2VV.Parameterizations.TimeAcceptance import Paper2012_csg_TimeAcceptance as TimeAcceptance
acceptances = {}
for s, model in split_models.iteritems():
    acceptance = TimeAcceptance(time = t, ParNamePrefix = '_'.join((run_period.GetName(), s)),
                                ResolutionModel = model, Input = input_file, Histograms = hists,
                                Fit = True, Cache = False)
    acceptances[s] = acceptance
    comp_pdf = sim_pdf.getPdf(s)
    comp_pdf.changeModel(acceptance.acceptance()._target_())
    
## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Optimize = 2,
               Strategy = 2, Offset = True, Minimizer = 'Minuit2')

valid_definition = [[(hlt1_excl_biased_dec, 'exclB'), (hlt1_excl_biased_dec, 'notExclB')], [(hlt2_biased, 'B'), (hlt2_unbiased, 'UB')]]
valid = valid_combinations(valid_definition)

data = None
if ntuple_file:
    from ROOT import TFile
    input_file = TFile(ntuple_file)
    tree_name = 'DecayTree'
    if input_file.FindKey(tree_name):
        input_file.Close()
        from P2VV.Utilities.DataHandling import readData
        data = readData(ntuple_file, tree_name, cuts = 'sel == 1 && (hlt1_biased == 1 || hlt1_unbiased_dec == 1) && (hlt2_biased == 1 || hlt2_unbiased == 1)',
                        NTuple = True, observables = observables)

        for i in range(3):
            mass_result = mass_pdf.fitTo(data, **fitOpts)
            if mass_result.status() == 0:
                break
        assert(mass_result.status() == 0)

        # Plot mass pdf
        from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
        from ROOT import TCanvas
        canvas = TCanvas('mass_canvas', 'mass_canvas', 600, 530)
        obs = [m]
        for (p,o) in zip(canvas.pads(len(obs)), obs):
            from P2VV.Utilities.Plotting import plot
            pdfOpts  = dict()
            plot(p, o, pdf = mass_pdf, data = data
                 , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
                 , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
                 , plotResidHist = True
                 , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed ),
                                  ## 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed ),
                                  'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                  }
                 )
        # Do the sWeights
        # make sweighted dataset. TODO: use mumu mass as well...
        from P2VV.Utilities.SWeights import SData

        for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
        splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
        data = splot.data('signal')
        ## psi_sdata = splot.data('psi_background')
        bkg_sdata = splot.data('background')

        if 'MC' in args[0]:
            import random
            ## Set more events to be unbiased, so we get some HLT2 exclusive biased
            ## sample.
            new_data = RooDataSet("new_data", "new_data", data.get())
            for i, obs in enumerate(data):
                b2 = obs.find('hlt2_biased')
                ub2 = obs.find('hlt2_unbiased')
                eb2 = obs.find('hlt2_excl_biased')
                if b2.getIndex() == 0:
                    pass
                elif random.random() < 0.5:
                    ub2.setIndex(0)
                    eb2.setIndex(1)
                new_data.add(obs)
                if i >= 50000:
                    break
    else:
        dataset_name = 'JpsiKK_sigSWeight'
        data = input_file.Get(dataset_name)
else:
    n_evt = int(5e4)
    from P2VV.Load import MultiCatGen
    
    ## Generate
    data = pdf.generate([t, hlt1_excl_biased_dec, hlt2_unbiased, hlt2_biased], NumEvents = n_evt)
    ## Use the valid combinations to create a cut to remove wrong events
    cut = ' || '.join('(' + ' && '.join('{0} == {0}::{1}'.format(c.GetName(), s) for c, s in comb) + ')' for comb in valid)
    data = data.reduce(Cut = cut, EventRange = (0, 30000))

# Make PDF without acceptance for the constraints
constraints = set()
from copy import copy
obs = copy(categories)
obs['time'] = t
for acceptance in acceptances.itervalues():
    constraints |= acceptance.build_multinomial_constraints(data, obs, extra_cat = run_period)

## Fit
print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
result = sim_pdf.fitTo(data, ExternalConstraints = pdf.ExternalConstraints() | constraints, SumW2Error = False, **fitOpts)
## profiler_stop()

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas, RooBinning
canvas = {}
print 'plotting'

states_signal = set([(state, label) for d in valid_definition for state, label in d])
def sort_combination(combination):
    valid_def = valid_definition[:]
    valid_def.reverse()
    level_left = 0
    n = 0
    c = set(combination)
    for level in valid_def:
        for j, state in enumerate(level):
            n |= int(state in c) << (level_left + j)
        level_left += len(level)
    return n - 1

def make_title(combination):
    title = []
    for level in valid_definition:
        l = level[0][0].GetName()[ : 4]
        level_states = set(level)
        s = [c for c in combination if c in level_states and c in states_signal]
        if len(s) == 1:
            title.append('%s_only_%s' % (l, s[0][1]))
        elif len(s) == 2:
            title.append('%s_both' % l)
    return '_X_'.join(title)
    
# Plot the lifetime shapes
__canvases = []
obs = [t]
for s in run_period:
    cname = '_'.join((run_period.GetName(), s.GetName(), 'canvas'))
    canv = TCanvas(cname, cname, 900, 700)
    __canvases.append(canv)
    for states, (p, o) in zip(sorted(valid, key = sort_combination),
                              (i for i in product(canv.pads(3, 2), obs))):
        name = '__'.join(['%s_%s' % (state.GetName(), label) for state, label in states])
        title = make_title(states)
        cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in states + ((run_period, s.GetName()),)])
        cat_data = data.reduce(cuts)
        project_set = RooArgSet(*project_vars)
        pdfOpts = dict(ProjWData = (project_set, cat_data, True), Slice = (run_period, s.GetName()))
        from P2VV.Utilities.Plotting import plot
        binning = acceptance.shapes()[0].base_binning()
        plot( p, o, cat_data, sim_pdf, plotResidHist = True,
              dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = binning),
              frameOpts = {'Title' : title},
              pdfOpts  = dict(LineWidth = 2, **pdfOpts),
              logy = False, logx = True)
    
# plot the efficiency shapes
__frames = []
def plot_shape(p, o, shape, errorOpts = {}, pdfOpts = {}):
    from operator import itemgetter
    i = shape.createIntegral(RooArgSet(o))
    n = i.getVal()
    p.cd()
    p.SetLogx(True)
    frame = o.frame()
    if errorOpts:
        r = errorOpts.pop('result')
        errorPlots = dict([(x, c) for x, c in errorOpts.iteritems() if type(x) == int])
        for x in errorPlots.keys():
            errorOpts.pop(x)
        entries = sorted(errorPlots.iteritems(), key = itemgetter(0))
        entries.reverse()
        for x, colour in entries:
            shape.plotOn(frame, VisualizeError = (r, x), FillColor = colour, **errorOpts)
    shape.plotOn(frame, **pdfOpts)
    frame.GetXaxis().SetTitle('decay time [ps]')
    n = shape.GetName()
    pos = n.find('hlt')
    title = n[pos : pos + 4]
    frame.GetYaxis().SetTitle(title)
    frame.GetYaxis().SetTitleOffset(1.05)
    frame.Draw()
    __frames.append(frame)

shapes = dict([(s.GetName(), []) for s in run_period])
for s, constraint in product(shapes.iterkeys(), constraints):
    if not s in constraint.GetName():
        continue
    if hasattr(constraint, 'efficiency'):
        shapes.append(s.efficiency())
    elif hasattr(constraint, 'epsB'):
        binning = acceptance.shapes()[0].base_binning()
        from P2VV.RooFitWrappers import BinnedPdf
        for n, l in [('epsA', constraint.epsA()), ('epsB', constraint.epsB())]:
            if l.getSize() <= 1:
                continue
            shapes[s].append(BinnedPdf('_'.join((constraint.GetName(), n, 'shape')), Observable = t,
                                       Binning = binning, Coefficients = [c for c in l]))

eff_canvases = {}
from ROOT import kYellow, kOrange
for period, s in shapes.iteritems():
    n = 'eff_canvas_' + period
    eff_canvas = TCanvas(n, n, 1200, 300)
    eff_canvases[period] = eff_canvas
    for p, shape in zip(eff_canvas.pads(3, 1), sorted(s, key = lambda e: e.GetName())):
        plot_shape(p, t, shape, errorOpts = {'result' : result, 2 : kYellow, 1 : kOrange})
