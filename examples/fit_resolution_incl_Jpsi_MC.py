from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))
t_true = RealVar('truetime', Title = 'true decay time', Unit='ps', Observable = True, MinMax=(-10000, 14))
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
unbiased = Category('triggerDecisionUnbiasedPrescaled', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'selected' : 1, 'not_selected' : 0})
clean_tail = Category('sel_cleantail', States = {'selected' : 1, 'not_selected' : 0})

observables = [t, t_true, m, mpsi, st, unbiased, selected, clean_tail]

# Read data
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## Data:
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_prescaled.root'
## Incl J/psi MC11a
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_from_incl_Jpsi_MC11a_ntupleB_for_fitting_20121010.root'
## Signal MC
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20120606.root'
data = readData(input_file, tree_name, cuts = '(sel == 1 && triggerDecisionUnbiasedPrescaled == 1)',
                NTuple = False, observables = observables)

# Add time difference (t - t_true) to data
t_diff = FormulaVar('time_diff', '@1 < -100 ? @0 : @0 - @1', [t, t_true], data = data)
t_diff.setMin(-10)
t_diff.setMax(10)
observables.append(t_diff)

clean_data = data.reduce('sel_cleantail == sel_cleantail::selected')

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, PerEventError = True,
                          BiasScaleFactor = False, Cache = False,
                          bias = dict(Value = -0.17, MinMax = (-1, 1)),
                          sigmaSF  = dict(Value = 1.46, MinMax = (0.1, 2)))

## from P2VVParameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
## sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = False, PerEventError = False,
##                           ScaleFactors = [(3, 0.5), (2, 0.08), (1, 0.04)],
##                           Fractions = [(3, 0.1), (2, 0.2)])
## sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True, PerEventError = True,
##                           ScaleFactors = [(2, 4), (1, 1.)],
##                           Fractions = [(2, 0.2)])

# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, sig_tres.model(), 'SingleSided'],
            ConditionalObservables = sig_tres.model().ConditionalObservables(),
            ExternalConstraints = sig_tres.model().ExternalConstraints())

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# J/psi mass pdf
from P2VVParameterizations.MassPDFs import Signal_PsiMass as PsiMassPdf
psi_m = PsiMassPdf(mpsi, Name = 'psi_m')

# J/psi background
from P2VVParameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf
bkg_mpsi = PsiBkgPdf(mpsi, Name = 'bkg_mpsi')

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

# Create psi background component
from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = sig_tres.model()
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 0.8)
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.25, MinMax = (0.5,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 0.16, MinMax = (0.1, 0.5))
                         )

## from P2VVParameterizations.TimePDFs import Single_Exponent_Time as Background_Time
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

signal = Component('signal', (sig_m.pdf(), psi_m.pdf(), sig_t), Yield = (200000, 500, 500000))
psi_background = Component('psi_background', (psi_m.pdf(), bkg_m.pdf(), psi_t), Yield= (3125,100,500000) )

background = Component('background', (bkg_mpsi.pdf(), bkg_m.pdf(), bkg_t), Yield = (19620,100,500000) )

# Prompt component
from P2VVParameterizations.TimePDFs import Prompt_Peak
prompt_pdf = Prompt_Peak(t, sig_tres.model(), Name = 'prompt_pdf')
psi_prompt = Component('prompt', (prompt_pdf.pdf(), ), Yield = (21582, 100, 500000))
                   
# Wrong PV components
from P2VVParameterizations.WrongPV import ShapeBuilder
wpv = ShapeBuilder(t, {'jpsi' : mpsi}, UseKeysPdf = True, Weights = 'jpsi',
                   InputFile = '/stuff/PhD/mixing/Bs2JpsiPhiPrescaledMC.root',
                   WorkSpace = 'Bs2JpsiPhiPrescaledMC_workspace',
                   Draw = True, t_diff = t_diff, sigmat = st)
wpv_psi = wpv.shape('jpsi')
psi_wpv = Component('psi_wpv', (wpv_psi,), Yield = (888, 50, 30000))

## Build PDF
## pdf = buildPdf(Components = (psi_background, psi_wrong_pv, background, bkg_wrong_pv), Observables = (mpsi,t), Name='pdf')

fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 1)
## mass_pdf = buildPdf(Components = (signal, psi_background), Observables = (m,), Name='mass_pdf')
## mass_pdf.fitTo(data, **fitOpts)

## # Plot mass pdf
## from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
## from ROOT import TCanvas
## mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
## obs = [m]
## for (p,o) in zip(mass_canvas.pads(len(obs)), obs):
##     from P2VVGeneralUtils import plot
##     pdfOpts  = dict()
##     plot(p, o, pdf = mass_pdf, data = data
##          , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
##          , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
##          , plotResidHist = True
##          , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
##                           , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
##                           , 'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
##                           }
##          )

## from P2VVGeneralUtils import SData
## for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
## splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
## signal_sdata = splot.data('signal')
## psi_sdata = splot.data('psi_background')
## bkg_sdata = splot.data('background')

time_pdf = buildPdf(Components = (psi_prompt, psi_background, psi_wpv), Observables = (t,), Name='time_pdf')
## time_pdf = buildPdf(Components = (psi_prompt, psi_background), Observables = (t,), Name='time_pdf')
time_pdf.Print("t")

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
result = time_pdf.fitTo(data, **fitOpts)
## profiler_stop()
## result.Print('v')

from array import array
from ROOT import RooBinning
bounds = array('d', [-5 + i * 0.1 for i in range(47)] + [-0.3 + i * 0.05 for i in range(12)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(21)])
binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('var_binning')
t.setBinning(binning)

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from ROOT import TCanvas
print 'plotting'
obs = [t]
plot_data = data
time_canvas = TCanvas('time_canvas', 'time_canvas', len(obs) * 1000, 650)
for (p,o) in zip(time_canvas.pads(len(obs)), obs):
    from P2VVGeneralUtils import plot
    pdfOpts  = dict(ProjWData = (RooArgSet(st), plot_data, True))
    plot(p, o, pdf = time_pdf if o != st else None, data = plot_data
         , frameOpts = dict(Title = "")
         , dataOpts = dict(MarkerSize = 0.8, Binning = binning, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
         , logy = True
         , plotResidHist = True
         , components = { 'psi_*'      : dict( LineColor = kRed,    LineStyle = kDashed )
                          , 'prompt_*' : dict( LineColor = kOrange, LineStyle = kDashed )
                          , 'wpv_*'    : dict( LineColor = kGreen,  LineStyle = kDashed )
                          }
         )
    
from Dilution import dilution
diff_pdf = wpv.diff_shape('jpsi')
dilution(t_diff, data, diff_pdf, result, psi_background, psi_wpv)
