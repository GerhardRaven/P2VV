#!/usr/bin/env python
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog year model')

parser.add_option("--no-pee", dest = "pee", default = True,
                  action = 'store_false', help = 'Do not use per-event proper-time error')
parser.add_option("-w", "--wpv", dest = "wpv", default = False,
                  action = 'store_true', help = 'Add WPV component')
parser.add_option("--wpv-type", dest = "wpv_type", default = 'Gauss', type = 'string',
                  action = 'store', help = 'Add WPV component [Gauss, Mixing]')
parser.add_option('-p', '--parameterisation', dest = 'parameterise', default = False,
                  action = 'store', help = 'Parameterise sigmas [RMS, Comb]')
parser.add_option("--verbose", dest = "verbose", default = False,
                  action = 'store_true', help = 'Verbose fitting')
parser.add_option("--offset", dest = "offset", default = False,
                  action = 'store_true', help = 'Use sigmat offset')
parser.add_option('-s', "--simultaneous", dest = "simultaneous", default = False,
                  action = 'store_true', help = 'Use sigmat offset')
parser.add_option("--plot", dest = "make_plots", default = False,
                  action = 'store_true', help = 'Make plots')

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
    ## input_data['data'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_prescaled.root'
    ## input_data['wpv'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhiPrescaled_2011.root'
    ## input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2011_workspace'
    input_data['data'] = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_20130117_MagDownMagUp.root'
    input_data['wpv'] = '/stuff/PhD/mixing/Bs2JpsiPhiPrescaled_2011.root'
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2011_workspace'
else:
    input_data['data'] = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20121218.root'
    input_data['wpv'] = '/stuff/PhD/mixing/Bs2JpsiPhiPrescaled_2012.root'
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2012_workspace'

if options.wpv and not options.wpv_type in ['Mixing', 'Gauss']:
    print parser.usage
    sys.exit(-2)
    
from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from P2VVLoad import LHCbStyle
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
if options.wpv:
    t_minmax = (-5, 14)
else:
    t_minmax = (-1.5, 8)
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))

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
                              BiasScaleFactor = False, Cache = False,
                              TimeResSFOffset = options.offset,
                              timeResMu = dict(Value = 0.0, MinMax = (-1, 1), Constant = False),
                              sigmaSF  = dict(Value = 1.46, MinMax = (0.1, 2)))
elif args[1] == 'double':
    from P2VVParameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = options.pee, Parameterise = options.parameterise,
                              TimeResSFOffset = options.offset, timeResMu = dict(Value = 0.0, Constant = False),
                              ScaleFactors = [(2, 2.1), (1, 1.26)] if options.pee else [(2, 0.1), (1, 0.06)],
                              Fractions = [(2, 0.2)])
elif args[1] == 'triple':
    from P2VVParameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = [False, options.pee, options.pee],
                              TimeResSFOffset = options.offset, Parameterise = options.parameterise,
                              ScaleFactors = [(3, 1.5), (2, 4), (1, 1.4)],
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
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 0.67)
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.37, MinMax = (0.5,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 0.13, MinMax = (0.1, 0.5))
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
psi_ll = Component('psi_ll', (psi_m.pdf(), bkg_m.pdf(), psi_t), Yield= (30000,100,500000) )

background = Component('background', (bkg_mpsi.pdf(), bkg_m.pdf(), bkg_t), Yield = (19620,100,500000) )

# Prompt component
from P2VVParameterizations.TimePDFs import Prompt_Peak
prompt_pdf = Prompt_Peak(t, sig_tres.model(), Name = 'prompt_pdf')
psi_prompt = Component('prompt', (prompt_pdf.pdf(), ), Yield = (77000, 100, 500000))

# Read data
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
cut = 'nPV == 3 && sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
if not options.wpv:
    cut += ' && sel_cleantail == 1'
data = readData(input_data['data'], tree_name, NTuple = True, observables = observables,
                ntupleCuts = cut)

if options.simultaneous:
    from array import array
    st_bins = array('d', [0.01 + i * 0.012 for i in range(6)])
    from ROOT import RooBinning
    st_binning = RooBinning( len(st_bins) - 1, st_bins, 'st_binning' )
    st.setBinning(st_binning, 'st_binning')
    st_cat = BinningCategory('sigmat_cat', Observable = st, Binning = st_binning,
                             Fundamental = True, Data = data, CatTypeName = 'bin')

## data = data.reduce(EventRange = (0, 50000))

fitOpts = dict(NumCPU = 6, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True,
               Verbose = options.verbose)

## Build a single mass PDF
mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='mass_pdf')

if options.simultaneous:
    from P2VVGeneralUtils import getSplitPar
    # categories for splitting the PDF
    # get mass parameters that are split
    split_cats = [[st_cat]]
    split_pars = [[par for par in mass_pdf.Parameters() if par.getAttribute('Yield')]]
    ## split_cats.append([st_cat])
    ## split_pars.append([par for par in bkg_mpsi.pdf().Parameters() if not par.isConstant()])
    
    # build simultaneous mass PDF
    from RooFitWrappers import SimultaneousPdf
    sWeight_mass_pdf = SimultaneousPdf(mass_pdf.GetName() + '_simul',
                                       MasterPdf       = mass_pdf,
                                       SplitCategories = split_cats,
                                       SplitParameters = split_pars)

    # set yields for categories
    split_cat_pars = sWeight_mass_pdf.getVariables()
    for ct in st_cat:
        psi_yield = getSplitPar('N_psi_ll', ct.GetName(), split_cat_pars)
        bkg_yield = getSplitPar('N_background', ct.GetName(), split_cat_pars)
        
        sel_str = '!(%s-%d)' % (st_cat.GetName(), ct.getVal())
        
        nEv    = data.sumEntries()
        nEvBin = data.sumEntries(sel_str)
        
        psi_yield.setVal( psi_yield.getVal() * nEvBin / nEv )
        psi_yield.setError( sqrt( psi_yield.getVal() ) )
        psi_yield.setMin(0.)
        psi_yield.setMax(nEvBin)
        bkg_yield.setVal( bkg_yield.getVal() * nEvBin / nEv )
        bkg_yield.setError( sqrt( bkg_yield.getVal() ) )
        bkg_yield.setMin(0.)
        bkg_yield.setMax(nEvBin)

    sWeight_mass_result = sWeight_mass_pdf.fitTo(data, **fitOpts)
    from P2VVGeneralUtils import SData
    sData = SData(Pdf = sWeight_mass_pdf, Data = data, Name = 'SimulMassSPlot')
    psi_sdata = sData.data('psi_ll')
    bkg_sdata = sData.data('background')
    sdatas = dict([(ct.GetName(), psi_sdata.reduce(Cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))) for ct in st_cat])
else:
    ## Fit mass pdf
    mass_pdf.fitTo(data, **fitOpts)

    ## Plot mass pdf
    from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
    from ROOT import TCanvas

    mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
    from P2VVGeneralUtils import SData
    from P2VVGeneralUtils import plot
    pdfOpts  = dict()
    plot(mass_canvas.cd(1), mpsi, pdf = mass_pdf, data = data
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
         , plotResidHist = False
         , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                          , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                          }
         )

    from P2VVGeneralUtils import SData
    for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
    sData = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
    psi_sdata = sData.data('psi_ll')
    bkg_sdata = sData.data('background')
    
# Wrong PV components
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

components = [psi_prompt, psi_ll]
if options.wpv and options.wpv_type == 'Mixing':
    from P2VVParameterizations import WrongPV
    reweigh_data = dict(jpsi = psi_sdata, bkg = bkg_sdata)
    wpv = WrongPV.ShapeBuilder(t, {'jpsi' : mpsi}, UseKeysPdf = True, Weights = 'jpsi', Draw = True,
                               InputFile = input_data['wpv'], Workspace = input_data['workspace'],
                               Reweigh = dict(Data = reweigh_data, DataVar = nPV, Binning = PV_bounds),
                               sigmat = st)
    wpv_psi = wpv.shape('jpsi')
    psi_wpv = Component('psi_wpv', (wpv_psi,), Yield = (1000, 50, 30000))
    components += [psi_wpv]
elif options.wpv and options.wpv_type == 'Gauss':
    wpv_mean = sig_tres._timeResMu
    wpv_sigma = RealVar('wpv_sigma', Value = 0.3, MinMax = (0.01, 1000))
    wpv_pdf = Pdf(Name = 'wpv_pdf', Type = Gaussian, Parameters = (t, wpv_mean, wpv_sigma))
    psi_wpv = Component('wpv', (wpv_pdf, ), Yield = (100, 5, 500000))
    components += [psi_wpv]

time_pdf = buildPdf(Components = components, Observables = (t,), Name='time_pdf')
time_pdf.Print("t")

if options.simultaneous:
    from RooFitWrappers import SimultaneousPdf
    split_pars = []
    if args[1] == 'single':
        split_pars.append([sig_tres._sigmaSF])
    else:
        split_pars.append(sig_tres.splitVars())

    ## if options.wpv and options.wpv_type == 'Gauss':    
    ##     split_pars.append([wpv_sigma, psi_wpv.getYield()])
    time_pdf = SimultaneousPdf(  time_pdf.GetName() + '_simul'
                                 , MasterPdf       = time_pdf
                                 , SplitCategories = [[st_cat]]
                                 , SplitParameters = split_pars)
    time_pdf.Print('t')

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
time_result = time_pdf.fitTo(psi_sdata, SumW2Error = False, **fitOpts)
## profiler_stop()
## result.Print('v')

from ROOT import RooBinning
if options.wpv:
    bounds = array('d', [-5 + i * 0.1 for i in range(47)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(21)])
else:
    bounds = array('d', [-1.5 + i * 0.1 for i in range(12)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])

binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('full_binning')

zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])
zoom_binning = RooBinning(len(zoom_bounds) - 1, zoom_bounds)
zoom_binning.SetName('zoom_binning')

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from ROOT import TCanvas
from P2VVGeneralUtils import plot

print 'plotting'

binnings = [binning, zoom_binning]
plotLog = [True, False]
__canvases = []

for i, (bins, pl) in enumerate(zip(binnings, plotLog)):
    if not options.make_plots:
        continue
    projSet = RooArgSet(st)
    if options.simultaneous:
        projSet.add(time_pdf.indexCat())
        r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
        for ct in st_cat:
            name = 'time_canvas_%s_%d' % (ct.GetName(), i)
            canvas = TCanvas(name, name, 600, 400)
            __canvases.append(canvas)
            p = canvas.cd(1)
            pd = sdatas[ct.GetName()]
            pdfOpts  = dict(ProjWData = (projSet, pd, True))
            plot(p, t, pdf = time_pdf, data = pd
                 , frameOpts = dict(Range = r, Title = "")
                 , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                 , pdfOpts  = dict(LineWidth = 4, Slice = (st_cat, ct.GetName()), **pdfOpts)
                 , logy = pl
                 , plotResidHist = False)
    else:
        canvas = TCanvas('time_canvas_%d' % i, 'time_canvas_%d' % i, 600, 400)
        __canvases.append(canvas)
        p = canvas.cd(1)
        r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
        pdfOpts  = dict(ProjWData = (projSet, psi_sdata, True))
        plot(p, t, pdf = time_pdf, data = pd
             , frameOpts = dict(Range = r, Title = "")
             , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
             , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
             , logy = pl
             , plotResidHist = False)

if options.simultaneous:
    split_bounds = array('d', [1000 * v for v in st_bins])
    
    res_canvas = TCanvas('res_canvas', 'res_canvas', 500, 500)
    from ROOT import TH1D
    
    hist_events = TH1D('hist_events', 'hist_events', len(split_bounds) - 1, split_bounds)
    hist_res = TH1D('hist_res', 'hist_res', len(split_bounds) - 1, split_bounds)
    mass_fpf = sWeight_mass_result.floatParsFinal()
    time_fpf = time_result.floatParsFinal()
    
    for index, ct in enumerate(st_cat):
        bin_name = '_'.join((psi_ll.getYield().GetName(), ct.GetName()))
        events = mass_fpf.find(bin_name)
        d = split_bounds[index + 1] - split_bounds[index]
        hist_events.SetBinContent(index + 1, events.getVal() / d)
        hist_events.SetBinError(index + 1, events.getError() / d)
        
        if args[1] == 'double' and options.parameterise == False:
            from PropagateErrors import propagateScaleFactor
            sf, sf_e = propagateScaleFactor(time_result, '_' + ct.GetName())
        elif args[1] == 'double' and options.parameterise == 'Comb':
            sf_comb = time_fpf.find('timeResComb_%s' % ct.GetName())
            sf, sf_e = sf_comb.getVal(), sf_comb.getError()
        elif args[1] == 'single':
            sf_var = time_fpf.find('sigmaSF_%s' % ct.GetName())
            sf, sf_e = sf_var.getVal(), sf_var.getError()
    
        mean = sdatas[ct.GetName()].mean(st._target_())
        res = mean * sf
        res_e = mean * sf_e
        hist_res.SetBinContent(index + 1, 1000 * res)
        hist_res.SetBinError(index + 1, 1000 * res_e)
    
    scale = 100 / hist_events.GetMaximum()
    hist_events.Scale(scale)
    hist_events.GetYaxis().SetRangeUser(0, 110)
    hist_res.GetYaxis().SetRangeUser(0, 110)
    hist_res.Draw('pe')
    hist_events.Draw('hist, same')
    from ROOT import kGray
    hist_events.SetFillColor(kGray + 1)
    hist_res.Draw('pe, same')
    hist_res.GetXaxis().SetTitle('estimated decay time error [fs]')
    hist_res.GetYaxis().SetTitle('decay time resulution [fs]')
    
    from ROOT import TF1
    fit_func = TF1('fit_func', "pol1", st_bins[0], st_bins[-1])
    fit_result = hist_res.Fit(fit_func, "S0")

import Dilution
Dilution.dilution(t, data, result = time_result, sigmat = st, signal = [psi_prompt],
                  subtract = [psi_ll, psi_wpv] if options.wpv else [psi_ll])



