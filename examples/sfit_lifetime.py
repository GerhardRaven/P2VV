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
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.02, 0.12))

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories
excl_biased = Category('triggerDecisionBiasedExcl', States = {'ExclBiased' : 1, 'Unbiased' : 0})
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, excl_biased, hlt1_biased, selected, hlt1_unbiased, hlt2_biased, hlt2_unbiased]
project_vars = [st, excl_biased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.484, MinMax = (1., 2.5))

# Time resolution model
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as SignalTimeResolution
sig_tres = SignalTimeResolution(Name = 'sig_tres', time = t, timeResSFConstraint = True, sigmat = st,
                                timeResSF = dict( Name = 'timeResSF', Value = 1.46, MinMax = (0.1,5.), Constant = True)
                                )
# from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as BackgroundTimeResolution
# bkg_tres = BackgroundTimeResolution(Name = 'bkg_tres', time = t, timeResSFConstraint = True)

# Use a very simple effective time resolution
## timeResMu = ConstVar(Name = 'timeResMu', Title = 'Decay time resolution mean',  Value = 0.)
## timeResSigma = RealVar(Name = 'timeResSigma', Title = 'Decay time resolution width', Value = 0.05, MinMax = (0.001, 0.1), Constant = True)
## from ROOT import RooGaussModel as GaussModel
## tresGauss = ResolutionModel( Name = 'timeResGaussModel' , Type = GaussModel, Parameters  = [t, timeResMu, timeResSigma])

# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, sig_tres.model(), 'SingleSided'],
            ConditionalObservables = sig_tres.model().ConditionalObservables(),
            ExternalConstraints = sig_tres.model().ExternalConstraints())

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.36, MinMax = (0.5, 3))
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 2, MinMax = (0.1, 4), Constant = True)
psi_m  = Pdf(Name = 'psi_m', Type = CrystalBall, Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# J/psi background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf(Name = 'bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )
background = Component('background', (bkg_m.pdf(), bkg_mpsi), Yield = (20000,2000,50000) )

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = sig_tres.model()
                       , t_bkg_fll    = dict( Name = 'bkg_t_fll',    Value = 0.2 )
                       , t_bkg_ll_tau = dict( Name = 'bkg_t_ll_tau', Value = 1.25, MinMax = (0.5,2.5) )
                       , t_bkg_ml_tau = dict( Name = 'bkg_t_ml_tau', Value = 0.16, MinMax = (0.01,0.5) )
                         )
bkg_t = bkg_t.pdf()

# Create psi background component
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = sig_tres.model()
                       , t_bkg_fll    = dict( Name = 'psi_t_fll',    Value = 0.2 )
                       , t_bkg_ll_tau = dict( Name = 'psi_t_ll_tau', Value = 1.25, MinMax = (0.5,2.5) )
                       , t_bkg_ml_tau = dict( Name = 'psi_t_ml_tau', Value = 0.16, MinMax = (0.01,0.5) )
                         )
psi_t = psi_t.pdf()
psi_background = Component('psi_background', (bkg_m.pdf(), psi_m), Yield= (10000,500,50000) )

# Apply acceptance
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
## input_file = '/stuff/PhD/p2vv/data/B_s0_Output.root'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_2011_biased_unbiased.root'
data = readData(input_file, tree_name, cuts = '(sel == 1 && (hlt1_unbiased == 1 || hlt1_biased == 1) && (hlt2_unbiased == 1 || hlt2_biased == 1))',
                NTuple = False, observables = observables)

# Time acceptance
from P2VVParameterizations.TimeAcceptance import Paper2012_TimeAcceptance
sig_acceptance = Paper2012_TimeAcceptance(time = t, Input = '/stuff/PhD/p2vv/data/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root',
                                          Histograms = {excl_biased : { 'ExclBiased' : {'histogram' : 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1ExclB_20bins'},
                                                                        'Unbiased'   : {'histogram' : 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_20bins'}
                                                                        }
                                                        },
                                          Data = data, Fit = False)
sig_t = sig_acceptance * sig_t

# Create signal component
signal = Component('signal', (sig_m.pdf(), psi_m), Yield = (21000,10000,30000))

## Build PDF
mass_pdf = buildPdf(Components = (signal, background), Observables = (m, ), Name='pdf')
mass_pdf.Print("t")

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True,
               Verbose = True, Optimize = 1, Minimizer = 'Minuit2')

# make sweighted dataset. TODO: use mumu mass as well...
from P2VVGeneralUtils import SData, splot

mass_pdf.fitTo(data, **fitOpts)

# Plot mass pdf
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
obs = [m]
for (p,o) in zip(canvas.pads(len(obs)), obs):
    from P2VVGeneralUtils import plot
    pdfOpts  = dict()
    plot(p, o, pdf = mass_pdf, data = data
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
         , plotResidHist = True
         , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                          ##, 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                          , 'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                          }
         )
assert(False)

for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
signal_sdata = splot.data('signal')
## psi_sdata = splot.data('psi_background')
bkg_sdata = splot.data('background')

## Fit
print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")

from P2VVGeneralUtils import plot
print 'plotting'
canvas = TCanvas('canvas', 'canvas', 500, 500)

obs = [t]
pdf = sig_t
sdata = signal_sdata
for p in pdf.getParameters(sdata):
    p.removeMin()
    p.removeMax()
result = pdf.fitTo(sdata, SumW2Error = False, **fitOpts)
result.Print('v')
pdfOpts  = dict(ProjWData = (RooArgSet(*project_vars), sdata, True))
p = canvas.cd(1)
plot(p, t, pdf = pdf, data = sdata
     , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
     , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
     , plotResidHist = True
     , logy = False
     , logx = True
     )
