#!/usr/bin/env python
import optparse
import sys
import os

parser = optparse.OptionParser(usage = '%prog year model')

parser.add_option("--no-pee", dest = "pee", default = True,
                  action = 'store_false', help = 'Do not use per-event proper-time error')

(options, args) = parser.parse_args()

if len(args) != 2:
    print parser.usage
    sys.exit(-2)
elif args[0] not in ['2011', '2012']:
    print parser.usage
    sys.exit(-2)
elif args[1] not in ['single', 'double', 'triple']:
    print parser.usage
    sys.exit(-2)
    
input_data = {}
if args[0] == '2011':
    input_data['data'] = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_prescaled.root'
    input_data['wpv'] = '/stuff/PhD/mixing/Bs2JpsiPhiPrescaled_2011.root'
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2011_workspace'
else:
    input_data['data'] = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20121218.root'
    input_data['wpv'] = '/stuff/PhD/mixing/Bs2JpsiPhiPrescaled_2012.root'
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2012_workspace'

from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from P2VVLoad import LHCbStyle
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

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
unbiased = Category('triggerDecisionUnbiasedPrescaled', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'selected' : 1, 'not_selected' : 0})
nPV = RealVar('nPV', Title = 'Number of PVs', Observable = True, MinMax = (0, 10))

observables = [t, m, mpsi, st, unbiased, selected, nPV]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

# Time resolution model
sig_tres = None
if args[1] == 'single':
    from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, PerEventError = options.pee,
                              BiasScaleFactor = False, Cache = True,
                              bias = dict(Value = -0.17, MinMax = (-1, 1)),
                              sigmaSF  = dict(Value = 1.46, MinMax = (0.1, 2)))
elif args[1] == 'double':
    from P2VVParameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = options.pee,
                              ScaleFactors = [(2, 2.3), (1, 1.2)],
                              Fractions = [(2, 0.2)])
elif args[1] == 'triple':
    from P2VVParameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = options.pee,
                              ScaleFactors = [(3, 0.5), (2, 0.08), (1, 0.04)],
                              Fractions = [(3, 0.1), (2, 0.2)])

# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, sig_tres.model(), 'SingleSided'],
            ConditionalObservables = sig_tres.model().ConditionalObservables(),
            ExternalConstraints = sig_tres.model().ExternalConstraints())

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
## sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))
m_sig_mean  = RealVar('m_sig_mean',   Unit = 'MeV', Value = 5365, MinMax = (5363, 5372))
m_sig_sigma = RealVar('m_sig_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
sig_m   = Pdf(Name = 'sig_m', Type = Gaussian,  Parameters = (m,m_sig_mean, m_sig_sigma ))

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
                         , bkg_t_fml    = dict(Name = 'bkg_t_fml',    Value = 0.76 )
                         , bkg_t_ll_tau = dict(Name = 'bkg_t_ll_tau', Value = 1., MinMax = (0.01, 2.5))
                         , bkg_t_ml_tau = dict(Name = 'bkg_t_ml_tau', Value = 0.1,  MinMax = (0.01, 0.5))
                         )
bkg_t = bkg_t.pdf()

signal = Component('signal', (sig_m, psi_m.pdf(), sig_t), Yield = (200000, 500, 500000))
psi_background = Component('psi_background', (psi_m.pdf(), bkg_m.pdf(), psi_t), Yield= (4000,100,500000) )

background = Component('background', (bkg_mpsi.pdf(), bkg_m.pdf(), bkg_t), Yield = (19620,100,500000) )

# Prompt component
from P2VVParameterizations.TimePDFs import Prompt_Peak
prompt_pdf = Prompt_Peak(t, sig_tres.model(), Name = 'prompt_pdf')
psi_prompt = Component('prompt', (prompt_pdf.pdf(), ), Yield = (77000, 100, 500000))

# Read data
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
data = readData(input_data['data'], tree_name, NTuple = True, observables = observables,
                cuts = '(sel == 1 && triggerDecisionUnbiasedPrescaled == 1)')
data = data.reduce(EventRange = (0, 400000))

fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True)
mass_pdf = buildPdf(Components = (psi_background, background), Observables = (mpsi,), Name='mass_pdf')
mass_pdf.fitTo(data, **fitOpts)

## # Plot mass pdf
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
obs = [mpsi]
for (p,o) in zip(mass_canvas.pads(len(obs)), obs):
    from P2VVGeneralUtils import plot
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

from P2VVGeneralUtils import SData
for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
## signal_sdata = splot.data('signal')
psi_sdata = splot.data('psi_background')
bkg_sdata = splot.data('background')

# Wrong PV components
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

from P2VVParameterizations import WrongPV
reweigh_data = dict(jpsi = psi_sdata, bkg = bkg_sdata)
wpv = WrongPV.ShapeBuilder(t, {'jpsi' : mpsi}, UseKeysPdf = True, Weights = 'jpsi', Draw = True,
                           InputFile = input_data['wpv'], Workspace = input_data['workspace'],
                           Reweigh = dict(Data = reweigh_data, DataVar = nPV, Binning = PV_bounds),
                           sigmat = st)
wpv_psi = wpv.shape('jpsi')
psi_wpv = Component('psi_wpv', (wpv_psi,), Yield = (1000, 50, 30000))

## wpv = WrongPV.ShapeBuilder(t, {'B' : m}, UseKeysPdf = True, Weights = 'B',
##                            Draw = True, sigmat = st)
## wpv_signal = wpv.shape('B')
## signal_wpv = Component('signal_wpv', (wpv_signal,), Yield = (888, 10, 300000))

print wpv.reweigh_weights('jpsi')

time_pdf = buildPdf(Components = (psi_prompt, psi_background, psi_wpv), Observables = (t,), Name='time_pdf')
time_pdf.Print("t")

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
result = time_pdf.fitTo(psi_sdata, SumW2Error = False, **fitOpts)
## profiler_stop()
## result.Print('v')


from ROOT import RooBinning
bounds = array('d', [-5 + i * 0.1 for i in range(47)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(21)])
## bounds = array('d', [-3 + i * 0.1 for i in range(27)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(21)])
binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('var_binning')
t.setBinning(binning)

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from ROOT import TCanvas
import P2VVGeneralUtils

print 'plotting'
obs = [t]
plot_data = psi_sdata
time_canvas = TCanvas('time_canvas', 'time_canvas', len(obs) * 1000, 650)
for (p,o) in zip(time_canvas.pads(len(obs)), obs):
    pdfOpts  = dict(ProjWData = (RooArgSet(st), plot_data, True))
    P2VVGeneralUtils.plot(p, o, pdf = time_pdf if o != st else None, data = plot_data
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

import Dilution
Dilution.dilution(t, data, result = result, sigmat = st, signal = [psi_prompt], subtract = [psi_background, psi_wpv])
