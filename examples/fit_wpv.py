from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

from ROOT import RooAbsReal
RooAbsReal.evalErrorLoggingMode = RooAbsReal.PrintErrors
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )

mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories needed for selecting events
unbiased = Category('triggerDecisionUnbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'selected' : 1, 'not_selected' : 0})

observables = [t, m, mpsi, st, unbiased, selected]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VV.Parameterizations.TimeResolution import Moriond2012_TimeResolution as TimeResolution
sig_tres = TimeResolution(time = t, timeResSigma = st, timeResSFConstraint = True, Cache = True,
                          timeResSigmaSF = dict( Name = 'timeResSF', Value = 1.46, MinMax = (0.1,5.)))

# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, sig_tres.model(), 'SingleSided'],
            ConditionalObservables = sig_tres.model().ConditionalObservables(),
            ExternalConstraints = sig_tres.model().ExternalConstraints())

# B mass pdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass
from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as Background_BMass
## sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

mass_mean = RealVar( 'm_sig_mean', Title = 'B Mass', Unit = 'MeV/c^2'
                     , Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
mass_sigma = RealVar('m_sig_sigma', Title = 'B Mass resolution', Unit = 'MeV/c^2'
                     , Value = 6.3,   Error = 0.1,  MinMax = ( 0.1, 20. ) )
sig_m = Pdf(Name ='sig_mass', Type = Gaussian, Parameters = (m, mass_mean, mass_sigma))

# J/psi mass pdf
from P2VV.Parameterizations.MassPDFs import Signal_PsiMass as PsiMassPdf
psi_m = PsiMassPdf(mpsi, Name = 'psi_m')

# J/psi background
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf
bkg_mpsi = PsiBkgPdf(mpsi, Name = 'bkg_mpsi')

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

# Create psi background component
from P2VV.Parameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = sig_tres.model()
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 0.788, MinMax = (0.1, 5))
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 0.795, MinMax = (0.5,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 0.16, MinMax = (0.1, 0.5))
                         )

## from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time as Background_Time
## psi_t = Background_Time(Name = 'psi_t', time = t, resolutionModel = sig_tres.model(),
##                              t_sig_tau  = dict(Name = 'psi_tau', Value = 1.5, MinMax = (0.5, 2.5))
##                              )

psi_t = psi_t.pdf()


bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = sig_tres.model()
                         , bkg_t_fml    = dict(Name = 'bkg_t_fml',    Value = 0.4 )
                         , bkg_t_ll_tau = dict(Name = 'bkg_t_ll_tau', Value = 1.29, MinMax = (0.01, 2.5))
                         , bkg_t_ml_tau = dict(Name = 'bkg_t_ml_tau', Value = 0.1,  MinMax = (0.01, 0.5))
                         )
bkg_t = bkg_t.pdf()

# Components
signal = Component('signal', (sig_m, psi_m.pdf(), sig_t), Yield = (2300, 100, 100000))
background = Component('background', (bkg_m.pdf(), bkg_mpsi.pdf(), bkg_t), Yield = (300000,100,500000) )
psi_background = Component('psi_background', (bkg_m.pdf(), psi_m.pdf(), psi_t), Yield= (3125,500,200000) )

# Wrong PV components
from P2VV.Parameterizations.WrongPV import ShapeBuilder
wpv = ShapeBuilder(t, {'B' : m, 'jpsi' : mpsi}, UseKeysPdf = True, Weights = 'both',
                   Draw = True, sigmat = st)
wpv_signal = wpv.shape('B')
signal_wpv = Component('signal_wpv', (wpv_signal,), Yield = (888, 50, 30000))

assert(False)

# Read data
from P2VV.GeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
## input_file = '/stuff/PhD/p2vv/data/B_s0_Output.root'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_prescaled.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
data = readData(input_file, tree_name, cuts = '(sel == 1 && triggerDecisionUnbiased == 1)',
                NTuple = False, observables = observables)

## Build PDF
## pdf = buildPdf(Components = (psi_background, psi_wrong_pv, background, bkg_wrong_pv), Observables = (mpsi,t), Name='pdf')

fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2)
mass_pdf = buildPdf(Components = (signal, psi_background, background), Observables = (m, mpsi), Name='mass_pdf')
mass_result = mass_pdf.fitTo(data, **fitOpts)

# Plot mass pdf
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
obs = [m, mpsi]
mass_canvas = TCanvas('mass_canvas', 'mass_canvas', len(obs) * 500, 650)
for (p,o) in zip(mass_canvas.pads(len(obs)), obs):
    from P2VV.GeneralUtils import plot
    pdfOpts  = dict()
    plot(p, o, pdf = mass_pdf, data = data
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
         , plotResidHist = True
         , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                          , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                          , 'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                          }
         )

from P2VV.GeneralUtils import SData
for p in mass_pdf.Parameters() :
    p.setConstant( not p.getAttribute('Yield') )
splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
signal_sdata = splot.data('signal')
## psi_sdata = splot.data('psi_background')
bkg_sdata = splot.data('background')

time_pdf = buildPdf(Components = (signal, signal_wpv), Observables = (t,), Name='time_pdf')
time_pdf.Print("t")

## values = {'psi_t_fml'        : ( 7.8806e-01, 2.26e-02),
##           'psi_t_ll_tau'     : ( 7.9488e-01, 8.11e-02),
##           'psi_t_ml_tau'     : ( 1.4411e-01, 9.25e-03),
##           'timeResFrac2'     : ( 1.5820e-01, 1.18e-02),
##           'timeResMu'        : (-1.2650e-01, 4.28e-03),
##           'timeResSigmaSF_1' : ( 1.3093e+00, 8.84e-03),
##           'timeResSigmaSF_2' : ( 2.5016e+00, 5.56e-02)}
## variables = time_pdf.getParameters(signal_sdata)
## for v in variables:
##     if v.GetName() in values:
##         val, err = values[v.GetName()]
##         v.setVal(val)
##         v.setError(err)
##         v.setConstant(True)

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
result = time_pdf.fitTo(signal_sdata, SumW2Error = False, **fitOpts)
## profiler_stop()
## result.Print('v')

from array import array
from ROOT import RooBinning
bounds = array('d', [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(21)])
binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('var_binning')
t.setBinning(binning)

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from ROOT import TCanvas
print 'plotting'
obs = [t]
time_canvas = TCanvas('time_canvas', 'time_canvas', len(obs) * 1000, 650)
for (p,o) in zip(time_canvas.pads(len(obs)), obs):
    from P2VV.GeneralUtils import plot
    pdfOpts  = dict(ProjWData = (RooArgSet(st), signal_sdata, True))
    plot(p, o, pdf = time_pdf if o != st else None, data = signal_sdata
         , dataOpts = dict(MarkerSize = 0.8, Binning = binning, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
         , logy = True
         , plotResidHist = True
         , components = { 'psi_*'   : dict( LineColor = kRed,    LineStyle = kDashed )
                          , 'sig_*' : dict( LineColor = kRed,    LineStyle = kDashed )
                          , 'wpv_*' : dict( LineColor = kGreen,  LineStyle = kDashed )
                          }
         )
