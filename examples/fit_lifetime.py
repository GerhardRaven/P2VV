from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.2, 14))
m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550))

# Categories
biased = Category('triggerDecision', States = {'Biased' : 1, 'NotBiased' : 0})
unbiased = Category('triggerDecisionUnbiased', States = {'Unbiased' : 1, 'NotUnbiased' : 0})
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, biased, unbiased, selected]

# Build acceptance
a   = RealVar('a', Title = 'a', Value = 1.45, MinMax = (1, 2))
c   = RealVar('c', Title = 'c', Value = -2.37, MinMax = (-3, -1))
pre = RealVar('eff_pre', Title = 'effective prescale', Value = 0.85, MinMax = (0.5, 0.99), Constant = True)
ub_eff = FormulaVar('pre', '@0', [pre, t])
b_eff  = FormulaVar('det', "(@0 > 0.) ? (1 / (1 + (@1 * @0) ** (@2))) : 0.0001", [t, a, c])

biased_eff   = FormulaVar('biased_eff',   "@0 * (1 - @1)", [b_eff, ub_eff])
unbiased_eff = FormulaVar('unbiased_eff', "@0 * (1 - @1)", [ub_eff, b_eff])
both_eff     = FormulaVar('both_eff',     "@0 * @1",       [ub_eff, b_eff])

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
tres = LP2011_TimeResolution(time = t)['model']

# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, tres, 'SingleSided'])

# Time acceptance
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
time_acceptance = Moriond2012_TimeAcceptance( time = t ).acceptance()
sig_t = time_acceptance * sig_t

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# Create signal component
signal = Component('signal', (sig_m.pdf(), sig_t), Yield = (30000,100,100000))

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = tres
                       , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.3 )
                       , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 1.92, MinMax = (0.5,2.5) )
                       , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.21, MinMax = (0.01,0.5) ) )
background = Component('background', (bkg_m.pdf(), bkg_t.pdf()), Yield = (50000,100,300000) )

# Apply acceptance
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
## input_file = '/stuff/PhD/p2vv/data/B_s0_Output.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120118.root'
data = readData(input_file, tree_name, cuts = '(sel == 1 && triggerDecision == 1)',
                NTuple = True, observables = observables)

## Build PDF
pdf = buildPdf(Components = (signal, background), Observables = (m, t), Name='pdf')
pdf.Print("t")

## from Helpers import Mapping
## mapping = Mapping({m : 'm', t : 'tau'}, data)

## Fit
print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
pdf.fitTo(data, NumCPU = 4 , Timer = 1, Minimizer = 'Minuit2')
## profiler_stop()

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
print 'plotting'
canvas = TCanvas('canvas', 'canvas', 1000, 500)
obs =  [o for o in pdf.Observables() if hasattr(o,'frame')]
for (p,o) in zip(canvas.pads(len(obs)), obs):
    from P2VVGeneralUtils import plot
    plot(p, o, data, pdf, components = {'sig*' : dict(LineColor = kGreen, LineStyle = kDashed),
                                        'bkg*' : dict(LineColor = kBlue,  LineStyle = kDashed)
                                        }
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2)
         , logy = ( o == t )
         )
