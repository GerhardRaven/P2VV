import os
import sys
import optparse
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog year')
(options, args) = parser.parse_args()

prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
input_data = {'Combined' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2MassNoMC_6KKMassBins_2TagCats_newTagging_trig.root'),
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
from P2VV.RooFitWrappers import *

from itertools import product
from ROOT import RooCBShape as CrystalBall
from P2VV.Parameterizations.GeneralUtils import valid_combinations
#from P2VV.Load import RooFitOutput

from ROOT import RooMsgService
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.ObjectHandling))
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

## project_vars = [hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]
categories = [hlt1_biased, hlt1_unbiased, hlt1_excl_biased_dec,
              hlt2_biased, hlt2_unbiased, hlt2_excl_biased]
categories = dict([(c.GetName(), c) for c in categories])

project_vars = [hlt1_excl_biased_dec, hlt1_unbiased, hlt2_biased, hlt2_unbiased, st]

selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, hlt1_biased, hlt1_unbiased, hlt1_excl_biased_dec,
               hlt2_biased, hlt2_unbiased, hlt2_excl_biased, selected, nPV]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VV.Parameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
mu = dict(MinMax = (-0.010, 0.010))
mu_values = {'MC11a_incl_Jpsi' : -0.000408, '2011_Reco14' : -0.00259,  
             '2011' : -0.00407301, '2012' : -0.00333,
             'MC2011_Sim08a_incl_Jpsi' : -0.00076}
mu['Value'] = mu_values.get(args[0], 0)
mu['Constant'] = True
from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
tres_args = dict(time = t, sigmat = st, Cache = False,
                 PerEventError = True, Parameterise = 'RMS',
                 TimeResSFParam = 'linear', timeResMu = mu, 
                 ScaleFactors = [(2, 2.00), (1, 1.174)],
                 Fractions = [(2, 0.239)],
                 sf_mean_offset = dict(Value = 1.4887, MinMax = (0.1, 2), Constant = True),
                 sf_mean_slope = dict(Value = -3.88, MinMax = (-5, 5), Constant = True),
                 sf_sigma_offset = dict(Value = 0.4143, MinMax = (0.1, 2), Constant = True),
                 sf_sigma_slope = dict(Value = -2.81, MinMax = (-5, 5), Constant = True),
                 timeResFrac2 = dict(Value = 0.239, MinMax = (0.01, 0.99), Constant = True))
sig_tres = TimeResolution(Name = 'sig_tres', **tres_args)

## from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution
## tres = LP2011_TimeResolution(time = t)
## from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
## tres = TimeResolution(time = t)

# Signal time pdf
from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time
sig_t = Single_Exponent_Time(Name = 'sig_t', time = t, resolutionModel = sig_tres.model())

# B mass pdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5365, MinMax = (5363,5372)))

# J/psi mass pdf
from P2VV.Parameterizations.MassPDFs import DoubleCB_Psi_Mass as PsiMassPdf
psi_m = PsiMassPdf(mpsi, Name = 'psi_m', mpsi_alpha_1 = dict(Value = 1.6832, Constant = True),
                   mpsi_frac = dict(Value = 0.521806), mpsi_mean = dict(Value = 3099.23),
                   mpsi_sigma_1 = dict(Value = 10.3741), mpsi_sigma_sf = dict(Value = 1.60555))

# J/psi background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.00081, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf(Name = 'bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

# Create components
signal_mass = Component('signal', (sig_m.pdf(), psi_m.pdf()), Yield = (30000,100,100000))
psi_background_mass = Component('psi_background', (bkg_m.pdf(), psi_m.pdf()), Yield= (100000,500,200000) )
background_mass = Component('background', (bkg_m.pdf(), bkg_mpsi), Yield = (100000,100,300000) )

## Build mass PDF
mass_pdf = buildPdf(Components = (signal_mass, background_mass), Observables = (m, ), Name='mass_pdf')
mass_pdf.Print("t")

## base_location = '/home/raaij'
base_location = '/stuff/PhD/p2vv'

# Build the acceptance using the histogram as starting values
input_file = os.path.join(base_location, 'data/start_values.root')

## hists = {hlt1_excl_biased_dec : {'excl_biased' : {'histogram' : 'hlt1_shape', 'average' : (6.285e-01, 1.633e-02)},
##                                  'unbiased' : { 'bins' : t.getRange(), 'heights' : [0.5]}}}
hists = {hlt1_excl_biased_dec : {'exclB' : {'histogram' : 'hlt1_shape', 'average' : (6.285e-01, 1.633e-02)},
                                 'notExclB' : { 'bins' : t.getRange(), 'heights' : [0.7]}},
         hlt2_biased : { 'B' : {'histogram' : 'hlt2_shape', 'average' : (6.3290e-01, 1.65e-02)}},
         hlt2_unbiased : { 'UB' : { 'bins' : t.getRange(), 'heights' : [0.5]}}}

from P2VV.Parameterizations.TimeAcceptance import Paper2012_csg_TimeAcceptance as TimeAcceptance
acceptance = TimeAcceptance(time = t, ResolutionModel = sig_tres, Input = input_file,
                            Histograms = hists, Fit = True, Cache = False)
pdf = Single_Exponent_Time(Name = 'pdf', time = t, resolutionModel = acceptance.model())
pdf = pdf.pdf()
pdf.Print('t')

# Read input data
from P2VV.Utilities.DataHandling import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Optimize = 2,
               Strategy = 2, Minimizer = 'Minuit2')

valid_definition = [[(hlt1_excl_biased_dec, 'exclB'), (hlt1_excl_biased_dec, 'notExclB')], [(hlt2_biased, 'B'), (hlt2_unbiased, 'UB')]]
valid = valid_combinations(valid_definition)

data = None
if ntuple_file:
    from ROOT import TFile
    input_file = TFile(ntuple_file)
    if input_file.FindKey(tree_name):
        input_file.Close()
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
        data = data.reduce('runPeriod == runPeriod::p2011')
else:
    ## PDF for sigmat
    from P2VV.Parameterizations.SigmatPDFs import DoubleLogNormal
    dln = DoubleLogNormal(st, frac_ln2 = dict(Value = 0.312508), k1 = dict(Value = 0.801757),
                          k2 = dict(Value = 1.37584), median = dict(Value = 0.0309409))
    ln = dln.pdf()
    gen_pdf = ProdPdf('gen_pdf', [pdf, ln])
    
    ## Generate
    from P2VV.Load import MultiCatGen
    data = gen_pdf.generate([t, st, hlt1_excl_biased_dec, hlt2_unbiased, hlt2_biased], 50000)
    ## Use the valid combinations to create a cut to remove wrong events
    cut = ' || '.join('(' + ' && '.join('{0} == {0}::{1}'.format(c.GetName(), s) for c, s in comb) + ')' for comb in valid)
    data = data.reduce(Cut = cut, EventRange = (0, 30000))

# Make PDF without acceptance for the constraints
constraints = set()
from copy import copy
obs = copy(categories)
obs['time'] = t
constraints |= acceptance.build_multinomial_constraints(data, obs)

pdf.setExternalConstraints(pdf.ExternalConstraints() | constraints)

## Fit
print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
result = pdf.fitTo(data, SumW2Error = False, **fitOpts)
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
canv = TCanvas('canvas', 'canvas', 900, 700)
obs = [t]
for states, (p, o) in zip(sorted(valid, key = sort_combination),
                          (i for i in product(canv.pads(3, 2), obs))):
    name = '__'.join(['%s_%s' % (state.GetName(), label) for state, label in states])
    title = make_title(states)
    cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in states])
    cat_data = data.reduce(cuts)
    project_set = RooArgSet(*project_vars)
    pdfOpts = dict(ProjWData = (project_set, cat_data, True))
    from P2VV.Utilities.Plotting import plot
    binning = acceptance.shapes()[0].base_binning()
    plot( p, o, cat_data, pdf, components = {'sig*' : dict(LineColor = kGreen, LineStyle = kDashed)}
          , plotResidHist = True
          , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = binning)
          , frameOpts = {'Title' : title}
          , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
          , logy = False
          , logx = True
          )
    
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

shapes = []
for s in pdf.ExternalConstraints():
    if hasattr(s, 'efficiency'):
        shapes.append(s.efficiency())
    elif hasattr(s, 'epsB'):
        binning = acceptance.shapes()[0].base_binning()
        from P2VV.RooFitWrappers import BinnedPdf
        shapes.append(BinnedPdf(s.GetName() + '_shape', Observable = t, Binning = binning,
                                Coefficients = (s.epsB() if s.epsB().getSize() > 1 else s.epsA())))

eff_canvases = {}
from ROOT import kYellow, kOrange
for p in ['p2011', 'p2012']:
    n = 'eff_canvas_' + p
    eff_canvas = TCanvas(n, n, 1200, 400)
    eff_canvases[p] = eff_canvas
    for p, shape in zip(eff_canvas.pads(2, 1), sorted(shapes, key = lambda s: s.GetName())):
        plot_shape(p, t, shape, errorOpts = {'result' : result, 2 : kYellow, 1 : kOrange})

output = {'hlt1_shape_10' : 'hlt1_excl_biased_dec_exclB',
          'hlt2_shape_10' : 'hlt2_biased_B'}
output_file = os.path.join(prefix, 'p2vv/data/start_values.root')
if os.path.exists(output_file):
    open_mode = 'update'
else:
    open_mode = 'new'

from ROOT import TFile
output_file = TFile(output_file, open_mode)

allVars = ws.allVars()
from ROOT import TH1D
from itertools import product
from array import array

for (name, pat), period in product(output.items(), ['p2011', 'p2012']):
    binning = time.getBinning(pat + '_binning')
    bins = array('d', [binning.binLow(i) for i in range(binning.numBins())] + [binning.highBound()])
    n = len(bins)
    heights = [v for v in allVars if '_'.join((period, pat)) in v.GetName()]
    heights = sorted(heights, key = lambda v: int(v.GetName().rsplit('_', 1)[-1]))
    v = [(h.getVal(), h.getError()) for h in heights]
    name = '_'.join((name, period))
    hist = TH1D(name, name, n - 1, bins)
    for i in range(1, n):
        hist.SetBinContent(i, v[i - 1][0])
        hist.SetBinError(i, v[i - 1][1])
    output_file.WriteTObject(hist, hist.GetName(), 'overwrite')

output_file.Close()
