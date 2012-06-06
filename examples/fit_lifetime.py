from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories
biased = Category('triggerDecision', States = {'Biased' : 1, 'NotBiased' : 0})
unbiased = Category('triggerDecisionUnbiased', States = {'Unbiased' : 1, 'NotUnbiased' : 0})
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, biased, unbiased, selected]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as SignalTimeResolution
sig_tres = SignalTimeResolution(Name = 'sig_tres', time = t, timeResSFConstraint = True, sigmat = st,
                                timeResSF = dict( Name = 'timeResSF', Value = 1.45, MinMax = (0.1,5.),
                                                  Constant = False))

# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, sig_tres.model(), 'SingleSided'],
            ConditionalObservables = sig_tres.model().ConditionalObservables(),
            ExternalConstraints = sig_tres.model().ExternalConstraints())

# Time acceptance
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
sig_acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/stuff/PhD/p2vv/data/efficiencies.root', Histogram = 'signal_efficiency_histo_20bins')
sig_t = sig_acceptance * sig_t

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.8, MinMax = (0.5, 3), Constant = True)
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 2, MinMax = (0.1, 4), Constant = True)
psi_m  = Pdf(Name = 'psi_m', Type = CrystalBall, Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# J/psi background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf(Name = 'bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])

# Create signal component
signal = Component('signal', (sig_m.pdf(), psi_m, sig_t), Yield = (21000,10000,30000))

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = sig_tres.model()
                         , bkg_t_fll    = dict( Name = 'bkg_t_fll',    Value = 0.2 )
                         , bkg_t_ll_tau = dict( Name = 'bkg_t_ll_tau', Value = 1.25, MinMax = (0.5,2.5) )
                         , bkg_t_ml_tau = dict( Name = 'bkg_t_ml_tau', Value = 0.16, MinMax = (0.01,0.5) )
                         , ExternalConstraints = sig_tres.model().ExternalConstraints())
bkg_acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/stuff/PhD/p2vv/data/efficiencies.root', Histogram = 'background_efficiency_histo_20bins')
bkg_t = bkg_acceptance * bkg_t.pdf()

# Create psi background component
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = sig_tres.model()
                         , psi_t_fll = dict( Name = 'psi_t_fll',    Value = 0.2 )
                         , psi_t_ll_tau = dict( Name = 'psi_t_ll_tau', Value = 1.25, MinMax = (0.5,2.5) )
                         , psi_t_ml_tau = dict( Name = 'psi_t_ml_tau', Value = 0.16, MinMax = (0.01,0.5) )
                         , ExternalConstraints = sig_tres.model().ExternalConstraints())
psi_acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/stuff/PhD/p2vv/data/efficiencies.root', Histogram = 'psi_background_efficiency_histo_20bins')
psi_t = psi_acceptance * psi_t.pdf()
psi_background = Component('psi_background', (bkg_m.pdf(), psi_m, psi_t), Yield= (10000,500,50000) )

## bkg_t = bkg_t.pdf()
background = Component('background', (bkg_m.pdf(), bkg_mpsi, bkg_t), Yield = (20000,2000,50000) )

# Apply acceptance
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
## input_file = '/stuff/PhD/p2vv/data/B_s0_Output.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
data = readData(input_file, tree_name, cuts = '(sel == 1 && triggerDecision == 1)',
                NTuple = False, observables = observables)
signal_data = data.reduce(CutRange = 'signal')
bkg_data    = data.reduce(CutRange = 'leftsideband' )
bkg_data.append(data.reduce(CutRange = 'rightsideband'))

## Build PDF
pdf = buildPdf(Components = (signal, psi_background, background), Observables = (m, mpsi,t), Name='pdf')
pdf.Print("t")

## from Helpers import Mapping
## mapping = Mapping({m : 'm', t : 'tau'}, data)

from ROOT import RooMsgService
## RooMsgService.instance().addStream(RooFit.DEBUG)
## RooMsgService.instance().getStream(1).removeTopic(RooFit.Caching)
## RooMsgService.instance().getStream(1).removeTopic(RooFit.Eval)

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Verbose = True, Minimizer = 'Minuit2', Optimize = 2)
sig_pdf = signal[m, t]
bkg_pdf = background[t,mpsi]
## signal[m, t].fitTo(signal_data, **fitOpts)
result = pdf.fitTo(data, **fitOpts)
## pdf.fitTo(data, **fitOpts)
## profiler_stop()
result.Print('v')

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
print 'plotting'
canvas = TCanvas('canvas', 'canvas', 1500, 500)
obs = [m, mpsi, t]
for (p,o) in zip(canvas.pads(len(obs)), obs):
    from P2VVGeneralUtils import plot
    pdfOpts  = dict(ProjWData = (RooArgSet(st), data, True))
    plot(p, o, pdf = pdf if o != st else None, data = data
         , plotResidHist = True
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
         , logy = ( o == t )
         , components = { 'psi_*'   : dict( LineColor = kGreen, LineStyle = kDashed )
                          , 'bkg_*' : dict( LineColor = kRed,   LineStyle = kDashed )
                          , 'sig_*' : dict( LineColor = kBlue,  LineStyle = kDashed )
                          }
         )
