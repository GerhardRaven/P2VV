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
parser.add_option("--fit-mass", dest = "fit_mass", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')
parser.add_option("--force-write", dest = "write_data", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')

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
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
if args[0] == '2011':
    input_data['data'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_20130207.root')
    input_data['wpv'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011.root')
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2011_workspace'
    input_data['results'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled_st_bins.root')
    input_data['cache'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled_st_bins.root')
else:
    input_data['data'] = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20121218.root'
    input_data['wpv'] = '/stuff/PhD/mixing/Bs2JpsiPhiPrescaled_2012.root'
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2012_workspace'
    input_data['cache'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_2012_Prescaled_st_bins.root'

if options.wpv and not options.wpv_type in ['Mixing', 'Gauss']:
    print parser.usage
    sys.exit(-2)
    

from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle
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
    from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, PerEventError = options.pee,
                              BiasScaleFactor = False, Cache = False,
                              TimeResSFOffset = options.offset,
                              timeResMu = dict(Value = 0.0, MinMax = (-1, 1), Constant = False),
                              sigmaSF  = dict(Value = 1.46, MinMax = (0.1, 2)))
elif args[1] == 'double':
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = options.pee, Parameterise = options.parameterise,
                              TimeResSFOffset = options.offset, timeResMu = dict(Value = 0.0, Constant = False),
                              ScaleFactors = [(2, 2.1), (1, 1.26)] if options.pee else [(2, 0.1), (1, 0.06)],
                              Fractions = [(2, 0.2)])
elif args[1] == 'triple':
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
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
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
## sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))
m_sig_mean  = RealVar('m_sig_mean',   Unit = 'MeV', Value = 5365, MinMax = (5363, 5372))
m_sig_sigma = RealVar('m_sig_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
sig_m   = Pdf(Name = 'sig_m', Type = Gaussian,  Parameters = (m,m_sig_mean, m_sig_sigma ))

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
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 0.67)
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.37, MinMax = (0.5,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 0.13, MinMax = (0.1, 0.5))
                         )

## from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time as Background_Time
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
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
prompt_pdf = Prompt_Peak(t, sig_tres.model(), Name = 'prompt_pdf')
psi_prompt = Component('prompt', (prompt_pdf.pdf(), ), Yield = (77000, 100, 500000))

# Read data
fit_mass = options.fit_mass

# Tree and cut
tree_name = 'DecayTree'
cut = 'nPV == 1 && sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
if not options.wpv:
    cut += ' && sel_cleantail == 1'
hd = ('%d' % hash(cut)).replace('-', 'm')

if options.simultaneous:
    from array import array
    split_bins = array('d', [0.01 + i * 0.012 for i in range(6)])
    directory = '%sbins_%4.2ffs_simul/%s' % (len(split_bins) - 1, (1000 * (split_bins[1] - split_bins[0])), hd)
else:
    directory = '1bin_%4.2ffs_simple/%s' % (1000 * (t.getMax() - t.getMin()), hd)

from ROOT import TFile
if os.path.exists(input_data['cache']):
    cache_file = TFile.Open(input_data['cache'], 'update')
else:
    cache_file = TFile.Open(input_data['cache'], 'new')
cache_dir = cache_file.Get(directory)
if not cache_dir:
    cache_file.mkdir(directory)
    cache_dir = cache_file.Get(directory)
    from ROOT import TObjString
    cut_string = TObjString(cut)
    cache_dir.WriteTObject(cut_string, 'cut')
    fit_mass = True

results = []
tree_name = 'DecayTree'
if not fit_mass:
    ## Read sdata
    sds_name = 'sdata'
    sdata_dir = cache_dir.Get('sdata')
    sdatas = {}
    if not sdata_dir:
        fit_mass = True
    else:
        nkeys = rd.ReadKeys()
        if nkeys != (len(split_bins)  - 1):
            fit_mass = True
        else:
            dss = []
            for e in rd.GetListOfKeys():
                if e.GetClassName() == 'RooDataSet':
                    dss.append(os.path.join(sdata_dir.GetName(), e.GetName()))
            for e in dss:
                sdata = chache_dir.Get(e)
                if not sdata:
                    fit_mass = True
                    break
                else:
                    sdatas[bin_data.GetName()] = sdata
            try:
                psi_sdata = sdatas['psi_sdata']
                bgk_sdata = sdatas['bkg_sdata']
            except KeyError:
                fit_mass = True

    # Read data
    data_dir = cache_dir.Get('data')
    if not data_dir:
        fit_mass = True
    else:
        data = data_dir.Get(tree_name)
        if not data:
            fit_mass = True

    # Read results
    rd = cache_dir.Get('results')
    if not rd:
        fit_mass = True
    else:
        mass_result = rd.Get('mass_result')
        if not mass_result:
            fit_mass = True
        else:
            results.append(mass_result)
        if options.simultaneous:
            sWeight_mass_result = rd.Get('sWeight_mass_result')
            if not sWeight_mass_result:
                fit_mass = True
            else:
                results.append(sWeight_mass_result)

## Fitting opts
fitOpts = dict(NumCPU = 6, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True,
               Verbose = options.verbose)

## List of all plots we make
plots = []

## Build simple mass pdf
if fit_mass:
    from P2VV.GeneralUtils import readData
    data = readData(input_data['data'], tree_name, NTuple = True, observables = observables,
                    ntupleCuts = cut)
    data.SetName(tree_name)
    mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='mass_pdf')
    mass_pdf.Print('t')

    ## Fit mass pdf
    mass_result = mass_pdf.fitTo(data, **fitOpts)
    mass_result.SetName('mass_result')
    results.append(mass_result)
    
    ## Plot mass pdf
    from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
    from ROOT import TCanvas

    mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
    from P2VV.GeneralUtils import SData
    from P2VV.GeneralUtils import plot
    pdfOpts  = dict()
    ps = plot(mass_canvas.cd(1), mpsi, pdf = mass_pdf, data = data
              , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
              , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
              , plotResidHist = False
              , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                               , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                               }
              )
    plots.append(ps)
    
    from P2VV.GeneralUtils import SData
    sData = SData(Pdf = mass_pdf, Data = data, Name = 'MassSPlot')
    psi_sdata = sData.data('psi_ll')
    psi_sdata.SetName('psi_sdata')
    bkg_sdata = sData.data('background')
    bkg_sdata.SetName('background_sdata')
    sdatas = {psi_sdata.GetName() : psi_sdata, bkg_sdata.GetName() : bkg_sdata}

if fit_mass and options.simultaneous:
    from ROOT import RooBinning
    st_binning = RooBinning( len(split_bins) - 1, split_bins, 'st_binning' )
    st.setBinning(st_binning, 'st_binning')
    st_cat = BinningCategory('sigmat_cat', Observable = st, Binning = st_binning,
                             Fundamental = True, Data = data, CatTypeName = 'bin')

    from P2VV.GeneralUtils import getSplitPar
    # categories for splitting the PDF
    # get mass parameters that are split
    split_cats = [[st_cat]]
    split_pars = [[par for par in mass_pdf.Parameters() if par.getAttribute('Yield')]]
    ## split_cats.append([st_cat])
    ## split_pars.append([par for par in bkg_mpsi.pdf().Parameters() if not par.isConstant()])
    
    # build simultaneous mass PDF
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
    sWeight_mass_result.SetName('sWeight_mass_result')
    results.append(sWeight_mass_result)
    from P2VV.GeneralUtils import SData
    sData = SData(Pdf = sWeight_mass_pdf, Data = data, Name = 'SimulMassSPlot')
    psi_sdata = sData.data('psi_ll')
    bkg_sdata = sData.data('background')
    for ct in st_cat:
        bin_data = psi_sdata.reduce(Cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))
        bin_data.SetName('psi_sdata_%s' % ct.GetName())
        sdatas[ct.GetName()] = bin_data

# Wrong PV components
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

components = [psi_prompt, psi_ll]
if options.wpv and options.wpv_type == 'Mixing':
    from P2VV.Parameterizations import WrongPV
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
time_result.SetName('time_result')
results.append(time_result)
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
from P2VV.GeneralUtils import plot

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
            ps = plot(p, t, pdf = time_pdf, data = pd
                      , frameOpts = dict(Range = r, Title = "")
                      , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                      , pdfOpts  = dict(LineWidth = 4, Slice = (st_cat, ct.GetName()), **pdfOpts)
                      , logy = pl
                      , plotResidHist = False)
            plots.append(ps)
    else:
        canvas = TCanvas('time_canvas_%d' % i, 'time_canvas_%d' % i, 600, 400)
        __canvases.append(canvas)
        p = canvas.cd(1)
        r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
        pdfOpts  = dict(ProjWData = (projSet, psi_sdata, True))
        ps = plot(p, t, pdf = time_pdf, data = pd
                  , frameOpts = dict(Range = r, Title = "")
                  , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                  , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                  , logy = pl
                  , plotResidHist = False)
        plots.append(ps)

fit_result = None
if options.simultaneous:
    split_bounds = array('d', [1000 * v for v in split_bins])
    
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
    fit_func = TF1('fit_func', "pol1", split_bins[0], split_bins[-1])
    fit_result = hist_res.Fit(fit_func, "S0")
    fr = fit_result.Get()
    fr.SetName('fit_result')
    results.append(fr)

from P2VV import Dilution
Dilution.dilution(t, data, result = time_result, sigmat = st, signal = [psi_prompt],
                  subtract = [psi_ll, psi_wpv] if options.wpv else [psi_ll])


# Write data to cache file
def get_dir(d):
    tmp = cache_dir.Get(d)
    if not tmp:
        cache_dir.mkdir(d)
        tmp = cache_dir.Get(d)
    return tmp

from ROOT import TObject
if options.write_data or fit_mass:
    sdata_dir = get_dir('sdata')
    data_dir = get_dir('data')
    for (dss, d) in [(sdatas, sdata_dir)]:
        for i, ds in enumerate(dss.itervalues()):
            d.WriteTObject(ds, ds.GetName(), "Overwrite")
        d.Write(d.GetName(), TObject.kOverwrite)

# Always write fit results
results_dir = get_dir('results')
for r in results:
    results_dir.WriteTObject(r, r.GetName(), "Overwrite")
results_dir.Write(results_dir.GetName(), TObject.kOverwrite)
# Delete the input TTree which was automatically attached.
cache_file.Delete('%s;*' % tree_name)
