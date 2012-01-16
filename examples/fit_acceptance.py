from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))
m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550))
biased = Category('triggerDecision', States = {'True' : 1, 'False' : 0})
unbiased = Category('triggerDecisionUnbiased', States = {'True' : 1, 'False' : 0})
observables = [t, m, biased, unbiased]

# Build acceptance
a = RealVar('a', Title = 'a', Value = 1.45, MinMax = (1, 2))
c = RealVar('c', Title = 'c', Value = -2.37, MinMax = (-3, 2))
eff = FormulaVar('eff_shape', "(@0 > 0.) ? (1 / (1 + (@1 * @0) ** (@2))) : 0.0001", [t, a, c])

from P2VVBinningBuilders import build1DVerticalBinning
binning, eff_func = build1DVerticalBinning('time_binning', eff, t, 0.05, 1.)

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
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t,signal_tau, tres, 'SingleSided'])
acceptance = BinnedPdf(Name = 'time_acceptance', Observables = [t], Function = eff,
                       Binning = binning)
sig_t_acc = acceptance * sig_t

# B mass pdf
m_mean  = RealVar('m_mean',   Unit = 'MeV', Value = 5370, MinMax = (5200, 5500))
m_sigma = RealVar('m_sigma',  Unit = 'MeV', Value = 8, MinMax = (5, 15))
sig_m = Pdf(Name = 'sig_m', Type = Gaussian,  Parameters = (m,m_mean, m_sigma ))

# Create signal component
signal = Component('signal', (sig_m, sig_t), Yield = (30000,100,100000))

# Create combinatorical background component
m_c = RealVar( 'm_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.00001))
bkg_m = Pdf(Name = 'bkg_m', Type = Exponential, Parameters = [m, m_c])
bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1, MinMax = (0.0001, 5))
comb_t = Pdf(Name = 'comb_t', Type = Decay,  Parameters = [t,bkg_tau, tres, 'SingleSided'])
comb_background = Component('comb_background', (bkg_m, comb_t), Yield = (50000,100,300000) )

# Build PDF
spec = {((t,), biased) : {'True' : {signal          : {'PDF'   : [sig_t_acc],
                                                       'Yield' : (20000, 1000, 50000)},
                                    comb_background : {'Yield' : (10000, 500, 50000)}
                                    },
                         'False' : {signal          : {'Yield' : (10000, 500, 50000)},
                                    comb_background : {'Yield' : (200000, 100000, 500000)}
                                    }
                          }
        }

pdf = buildSimultaneousPdf(Components = (signal, comb_background), Observables = (m, t), Spec = spec,
                           Name='pdf')
pdf.Print("t")

# Apply acceptance (dirty way)
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
data = readData(input_file, tree_name, cuts = '(triggerDecisionUnbiased == 1)',
                NTuple = True, observables = [m, t, biased, unbiased])

## from Helpers import Mapping
## mapping = Mapping({m : 'm', t : 'tau'}, data)

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
pdf.fitTo(data, NumCPU = 4 , Timer=1, Minimizer = 'Minuit2')
## profiler_stop()

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
canvas = {}
print 'plotting'
for cat in biased:
    idx, label = cat.getVal(), cat.GetName()
    name = biased.GetName() + '_' + label
    canv = canvas[name] = TCanvas(name, name, 500, 500)
    obs =  [o for o in pdf.Observables() if hasattr(o,'frame')]
    for (p,o) in zip(canv.pads(len(obs)), obs):
        dataOpts = dict(Cut = '{0} == {0}::{1}'.format(biased['Name'], label) )
        pdfOpts  = dict(Slice = (biased, label), ProjWData = (RooArgSet(biased), data))
        from P2VVGeneralUtils import plot
        plot( p, o, data, pdf, components = { 'sig*' : dict(LineColor = kGreen, LineStyle = kDashed)
                                            , 'bkg*' : dict(LineColor = kBlue,  LineStyle = kDashed)
                                              }
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
              , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
              , logy = ( o == t )
              )

## m_frame = m.frame()
## data.plotOn(m_frame, Cut = 'triggerDecision == triggerDecision::True')
## pdf.plotOn(m_frame, LineWidth = 2, Slice = (biased, 'True'),
##            ProjWData = (RooArgSet(biased), data))
## pdf.plotOn(m_frame, Components = ('sig_m'), LineStyle = kDashed, LineWidth = 2,
##            Slice = (biased, 'True'), ProjWData = (RooArgSet(biased), data))
## pdf.plotOn(m_frame, Components = ('bkg_m'), LineStyle = kDashed, LineWidth = 2, LineColor= kRed,
##            Slice = (biased, 'True'), ProjWData = (RooArgSet(biased), data))
