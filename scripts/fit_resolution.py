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
parser.add_option("--model", dest = "model", default = '', type = 'string',
                  action = 'store', help = 'Type of model used to scale sigma_t')
parser.add_option('-s', "--simultaneous", dest = "simultaneous", default = False,
                  action = 'store_true', help = 'Use sigmat offset')
parser.add_option("--plot", dest = "make_plots", default = False,
                  action = 'store_true', help = 'Make plots')
parser.add_option("--fit-mass", dest = "fit_mass", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')
parser.add_option("--force-write", dest = "write_data", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')
parser.add_option("--reduce", dest = "reduce", default = False,
                  action = 'store_true', help = 'Reduce sdata sets to 2000 entries per bin.')
parser.add_option("--extra-cut", dest = "cut", default = '', type = 'string',
                  action = 'store', help = 'extra cut to apply')
parser.add_option("-b", "--batch", dest = "batch", default = False,
                  action = 'store_true', help = 'run ROOT in batch mode')
(options, args) = parser.parse_args()

if len(args) != 2:
    print parser.usage
    sys.exit(-2)
elif args[0] not in ['2011', '2012', 'MC11a']:
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
elif args[0] == '2012':
    input_data['data'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20121218.root')
    input_data['wpv'] = os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_2012.root')
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2012_workspace'
    input_data['cache'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_Prescaled_st_bins.root')
else:
    ## MC11a
    input_data['data'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130222.root')
    input_data['wpv'] = os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC11a.root')
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_MC11a_workspace'
    input_data['cache'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_Prescaled_cache.root')
    
if options.wpv and not options.wpv_type in ['Mixing', 'Gauss']:
    print parser.usage
    sys.exit(-2)

if options.batch:
    from ROOT import gROOT
    gROOT.SetBatch(True)
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
                              TimeResSFModel = options.model,
                              timeResMu = dict(Value = -0.005, MinMax = (-2, 2), Constant = options.model != '' ),
                              sigmaSF  = dict(Value = 1.46, MinMax = (0.1, 2)))
elif args[1] == 'double':
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = options.pee, Parameterise = options.parameterise,
                              TimeResSFOffset = options.model != '',
                              timeResMu = dict(Value = 0, Constant = options.model != ''),
                              ScaleFactors = [(2, 2.1), (1, 1.26)] if options.pee else [(2, 0.1), (1, 0.06)],
                              Fractions = [(2, 0.2)])
elif args[1] == 'triple':
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = [False, options.pee, options.pee],
                              TimeResSFOffset = options.model != '', Parameterise = options.parameterise,
                              ScaleFactors = [(3, 1.5), (2, 4), (1, 1.4)],
                              Fractions = [(3, 0.1), (2, 0.2)])

# J/psi mass pdf
from P2VV.Parameterizations.MassPDFs import Signal_PsiMass as PsiMassPdf
psi_m = PsiMassPdf(mpsi, Name = 'psi_m')

# J/psi background
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf
bkg_mpsi = PsiBkgPdf(mpsi, Name = 'bkg_mpsi')

# Create combinatorical background component
mean = RealVar(Name = 'mean', Value = 5368, MinMax = (5300, 5400))
sigma = RealVar(Name = 'sigma', Value = 50, MinMax = (1, 500))
bkg_m = Pdf(Name = 'gauss', Type = Gaussian, Parameters = (m, mean, sigma))

# Create psi background component
from P2VV.Parameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = sig_tres.model()
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 0.67)
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.37, MinMax = (0.5,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 0.103, MinMax = (0.1, 0.5))
                         )
psi_t = psi_t.pdf()

# J/psi signal component
psi_ll = Component('psi_ll', (psi_m.pdf(), bkg_m, psi_t), Yield= (30000,100,500000) )

# Background component
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = sig_tres.model()
                         , bkg_t_fml    = dict(Name = 'bkg_t_fml',    Value = 0.76 )
                         , bkg_t_ll_tau = dict(Name = 'bkg_t_ll_tau', Value = 1., MinMax = (0.01, 2.5))
                         , bkg_t_ml_tau = dict(Name = 'bkg_t_ml_tau', Value = 0.1,  MinMax = (0.01, 0.5))
                         )
bkg_t = bkg_t.pdf()

background = Component('background', (bkg_mpsi.pdf(), bkg_m, bkg_t), Yield = (19620,100,500000) )

# B signal component
from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time as Signal_Time
sig_t = Signal_Time(time = t, resolutionModel = sig_tres.model(), Name = 'sig_t',
                    tau = dict(Name = 'sig_tau', Value = 1.5, MinMax = (0.5, 2.5)))

from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_Mass
sig_m = Signal_Mass(Name = 'sig_m', mass = m) 
signal = Component('signal', (sig_m.pdf(), sig_t.pdf()), Yield = (150000, 50000, 1000000))

# Prompt component
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
prompt_pdf = Prompt_Peak(t, sig_tres.model(), Name = 'prompt_pdf')
prompt = Component('prompt', (prompt_pdf.pdf(), ), Yield = (77000, 100, 500000))

# Read data
fit_mass = options.fit_mass or args[0] == 'MC11a'

# Tree and cut
tree_name = 'DecayTree'
cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
if not options.wpv or (options.wpv and options.wpv_type == "Gauss"):
    cut += ' && sel_cleantail == 1'
if args[0] == 'MC11a':
    cut += ' && abs(trueid) == 531'
if options.cut:
    cut = options.cut + ' && ' + cut
hd = ('%d' % hash(cut)).replace('-', 'm')

if options.simultaneous:
    from array import array
    split_bins = array('d', [0.01, 0.021, 0.025, 0.028, 0.031, 0.0345, 0.0375, 0.042, 0.048, 0.07])
    ## split_bins = array('d', [0.01, 0.02, 0.025, 0.0275, 0.0325, 0.035, 0.0375, 0.041, 0.045, 0.05, 0.07])
    ## split_bins = array('d', [0.01 + i * 0.01 for i in range(5)] + [0.07])
    directory = '%sbins_%4.2ffs_simul/%s' % (len(split_bins) - 1, (1000 * (split_bins[1] - split_bins[0])), hd)
else:
    directory = '1bin_%4.2ffs_simple/%s' % (1000 * (t.getMax() - t.getMin()), hd)

from ROOT import TFile
if os.path.exists(input_data['cache']):
    cache_file = TFile.Open(input_data['cache'], 'update')
else:
    cache_file = TFile.Open(input_data['cache'], 'new')
if not cache_file:
    cache_file = TFile.Open(input_data['cache'], 'recreate')
    assert(cache_file)
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

## Canvas for correlation histograms
from ROOT import TCanvas
corr_canvas = TCanvas('corr_canvas', 'corr_canvas', 1000, 500)
corr_canvas.Divide(2, 1)
from ROOT import gStyle
gStyle.SetPalette(1)

## Read Cache
if not fit_mass and args[0] != 'MC11a':
    ## Read sdata
    sds_name = 'sdata'
    sdata_dir = cache_dir.Get('sdata')
    sdatas_full = {}
    if not sdata_dir:
        fit_mass = True
    else:
        nkeys = sdata_dir.ReadKeys()
        if options.simultaneous and nkeys - 4 != (len(split_bins)  - 1):
            fit_mass = True
        elif nkeys != 2:
            fit_mass = True            
        else:
            dss = []
            for e in sdata_dir.GetListOfKeys():
                if e.GetClassName() == 'RooDataSet':
                    dss.append(os.path.join(sdata_dir.GetName(), e.GetName()))
            for e in dss:
                sdata = cache_dir.Get(e)
                if not sdata:
                    fit_mass = True
                    break
                else:
                    sdatas_full[e] = sdata
            try:
                single_bin_sig_sdata = sdatas_full['full_sig_sdata']
                single_bin_bkg_sdata = sdatas_full['full_bkg_sdata']
                if options.simultaneous:
                    sig_sdata_full = sdatas_full['sig_sdata']
                    bkg_sdata_full = sdatas_full['bkg_sdata']
            except KeyError:
                fit_mass = True
            ## CHECK ST BINNING
            
    from copy import copy
    sdatas = copy(sdatas_full)
    if options.reduce and not fit_mass:
        ct_names = set([ct.GetName() for ct in st_cat])
        for k, ds in sdatas.items():
            if k not in ct_names:
                continue
            if ds.numEntries() < 2000:
                continue
            sdatas[k] = ds.reduce(EventRange = (0, 2000))
        bin_datas = filter(lambda e: e.GetName().find('bin') != -1, sdatas.values())
        sig_sdata = bin_datas[0].Clone(sig_sdata_full.GetName())
        for bin_data in bin_datas[1:]:
            sig_sdata.append(bin_data)
    elif not fit_mass:
        sig_sdata = sig_sdata_full
        bkg_sdata = bkg_sdata_full
    
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
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True,
               Verbose = options.verbose)

## List of all plots we make
plots = []

## Build simple mass pdf
if fit_mass:
    from P2VV.GeneralUtils import readData
    data = readData(input_data['data'], tree_name, NTuple = True, observables = observables,
                    ntupleCuts = cut)
    data.SetName(tree_name)

    if args[0] == 'MC11a':
        mass_pdf = buildPdf(Components = (signal, background), Observables = (m,), Name = 'mass_pdf')
        mass_pdf.Print('t')
        signal_name = signal.GetName()
    else:
        mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='mass_pdf')
        mass_pdf.Print('t')
        signal_name = psi_ll.GetName()

    ## Fit mass pdf
    mass_result = mass_pdf.fitTo(data, **fitOpts)
    mass_result.SetName('mass_result')
    results.append(mass_result)

    ## Plot correlation histogram
    corr_canvas.cd(1)
    corr_hist_mass = mass_result.correlationHist()
    corr_hist_mass.GetXaxis().SetLabelSize(0.03)
    corr_hist_mass.GetYaxis().SetLabelSize(0.03)
    corr_hist_mass.SetContour(20)
    corr_hist_mass.Draw('colz')

    ## Plot mass pdf
    from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
    from ROOT import TCanvas

    mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
    from P2VV.GeneralUtils import SData
    from P2VV.GeneralUtils import plot
    pdfOpts  = dict()
    if args[0] == 'MC11a':
        mass_obs = m
    else:
        mass_obs = mpsi
    ps = plot(mass_canvas.cd(1), mass_obs, pdf = mass_pdf, data = data
              , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
              , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
              , plotResidHist = False
              , components = { 'sig_*'     : dict( LineColor = kOrange,   LineStyle = kDashed )
                               , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                               , 'bkg_*'  : dict( LineColor = kRed, LineStyle = kDashed )
                               }
              )
    plots.append(ps)
    
    from P2VV.GeneralUtils import SData
    data_clone = data.Clone(data.GetName())
    sData = SData(Pdf = mass_pdf, Data = data_clone, Name = 'MassSPlot')
    single_bin_sig_sdata = sData.data(signal_name)
    single_bin_bkg_sdata = sData.data('background')
    sdatas_full = {'full_sig_sdata' : single_bin_sig_sdata,
                   'full_bkg_sdata' : single_bin_bkg_sdata}

if fit_mass and options.simultaneous:
    from P2VV.GeneralUtils import getSplitPar
    # categories for splitting the PDF
    # get mass parameters that are split
    from ROOT import RooBinning
    st_binning = RooBinning( len(split_bins) - 1, split_bins, 'st_binning' )
    st.setBinning(st_binning, 'st_binning')
    st_cat = BinningCategory('sigmat_cat', Observable = st, Binning = st_binning,
                             Fundamental = True, Data = data, CatTypeName = 'bin')

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
        psi_yield = getSplitPar('N_' + signal_name, ct.GetName(), split_cat_pars)
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

    for i in range(5):
        sWeight_mass_result = sWeight_mass_pdf.fitTo(data, **fitOpts)
        if sWeight_mass_result.status() == 0:
            break
    assert(sWeight_mass_result.status() == 0)
    sWeight_mass_result.SetName('sWeight_mass_result')
    results.append(sWeight_mass_result)

    ## Plot correlation histogram
    corr_canvas.cd(1)
    corr_hist_sWmass = sWeight_mass_result.correlationHist()
    corr_hist_sWmass.GetXaxis().SetLabelSize(0.03)
    corr_hist_sWmass.GetYaxis().SetLabelSize(0.03)
    corr_hist_sWmass.SetContour(20)
    corr_hist_sWmass.Draw('colz')

    from P2VV.GeneralUtils import SData
    sData = SData(Pdf = sWeight_mass_pdf, Data = data, Name = 'SimulMassSPlot')
    sig_sdata_full = sData.data(signal_name)
    bkg_sdata_full = sData.data('background')
    sdatas_full['sig_sdata'] = sig_sdata_full
    sdatas_full['bkg_sdata'] = bkg_sdata_full

    from copy import copy
    sdatas = copy(sdatas_full)
    if options.reduce:
        ct_names = set([ct.GetName() for ct in st_cat])
        from copy import copy
        sdatas = copy(sdatas_full)
        for k, ds in sdatas.items():
            if k not in ct_names:
                continue
            if ds.numEntries() < 2000:
                continue
            sdatas[k] = ds.reduce(EventRange = (0, 2000))
        bin_datas = filter(lambda e: e.GetName().find('bin') != -1, sdatas.values())
        sig_sdata = bin_datas[0].Clone(sig_sdata_full.GetName())
        for bin_data in bin_datas[1:]:
            sig_sdata.append(bin_data)
    else:
        sig_sdata = sig_sdata_full
        bkg_sdata = bkg_sdata_full
elif fit_mass:
    sig_sdata = single_bin_sig_sdata
    bkg_sdata = single_bin_bkg_sdata

corr_canvas.Update()

# Wrong PV components
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

# Define default components
if args[0] == 'MC11a':
    components = [signal]
else:
    components = [prompt, psi_ll]

if options.wpv and options.wpv_type == 'Mixing':
    from P2VV.Parameterizations import WrongPV
    if args[0] == 'MC11a':
        reweigh_data = dict(B = data)
        masses = {'B' : m}
        weights = 'B'
    else:
        reweigh_data = dict(jpsi = single_bin_sig_sdata, bkg = single_bin_bkg_sdata)
        masses = {'jpsi' : mpsi}
        weights = 'jpsi'
    wpv = WrongPV.ShapeBuilder(t, masses, UseKeysPdf = True, Weights = weights, Draw = True,
                               InputFile = input_data['wpv'], Workspace = input_data['workspace'],
                               Reweigh = dict(Data = reweigh_data, DataVar = nPV, Binning = PV_bounds),
                               sigmat = st)
    if args[0] == 'MC11a':
        wpv_signal = wpv.shape('B')
        signal_wpv = Component('signal_wpv', (wpv_signal,), Yield = (888, 50, 300000))
        components += [signal_wpv]
    else:
        wpv_psi = wpv.shape('jpsi')
        psi_wpv = Component('psi_wpv', (wpv_psi,), Yield = (1000, 50, 30000))
        components += [psi_wpv]
elif options.wpv and options.wpv_type == 'Gauss':
    wpv_mean = sig_tres._timeResMu
    wpv_sigma = RealVar('wpv_sigma', Value = 0.3, MinMax = (0.01, 1000))
    wpv_pdf = Pdf(Name = 'wpv_pdf', Type = Gaussian, Parameters = (t, wpv_mean, wpv_sigma))
    psi_wpv = Component('wpv', (wpv_pdf, ), Yield = (100, 1, 500000))
    components += [psi_wpv]

time_pdf = buildPdf(Components = components, Observables = (t,), Name='time_pdf')
time_pdf.Print("t")

if options.simultaneous:
    split_pars = []
    split_pars = [[par for par in time_pdf.Parameters() if par.getAttribute('Yield')]]
    if args[1] == 'single':
        split_pars[0].append(sig_tres._sigmaSF)
    else:
        split_pars[0] += sig_tres.splitVars()

    if args[0] == 'MC11a':
        split_pars[0].append(sig_t._tau)

    if options.wpv and options.wpv_type == 'Gauss':    
        split_pars.append([wpv_sigma, psi_wpv.getYield()])
    time_pdf = SimultaneousPdf(  time_pdf.GetName() + '_simul'
                                 , MasterPdf       = time_pdf
                                 , SplitCategories = [[st_cat]]
                                 , SplitParameters = split_pars)
    time_pdf.Print('t')

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
for i in range(5):
    time_result = time_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)
    if time_result.status() == 0:
        break
extra_name = args[1]
if options.model:
    extra_name += ('_' + options.model)
time_result.SetName('time_result_%s' % extra_name)

## Draw correlation histogram
corr_canvas.cd(2)
## FIXME
corr_hist_time = time_result.correlationHist()
corr_hist_time.GetXaxis().SetLabelSize(0.03)
corr_hist_time.GetYaxis().SetLabelSize(0.03)
corr_hist_time.SetContour(20)
corr_hist_time.Draw('colz')

results.append(time_result)

## profiler_stop()
## result.Print('v')

from ROOT import RooBinning
if options.wpv and options.wpv_type == 'Mixing':
    bounds = array('d', [-5 + i * 0.1 for i in range(47)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(21)])
    zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])
elif args[0] == 'MC11a':
    bounds = array('d', [-1.5 + i * 0.1 for i in range(12)] + [-0.3 + i * 0.05 for i in range(12)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])
    zoom_bounds = array('d', [-0.1 + i * 0.005 for i in range(81)])
else:
    bounds = array('d', [-1.5 + i * 0.1 for i in range(12)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])
    zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])

binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('full')

zoom_binning = RooBinning(len(zoom_bounds) - 1, zoom_bounds)
zoom_binning.SetName('zoom')

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from P2VV.GeneralUtils import plot

print 'plotting'

binnings = [binning, zoom_binning]
plotLog = [True, False]
__canvases = []

from ROOT import SetOwnership

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
            
            opts = dict(Cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))
            pd = sig_sdata_full.reduce(**opts)
            SetOwnership(pd, False)
            pdfOpts  = dict(ProjWData = (projSet, pd, True))
            ps = plot(p, t, pdf = time_pdf, data = pd
                      , frameOpts = dict(Range = r, Title = "")
                      , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                      , pdfOpts  = dict(LineWidth = 4, Slice = (st_cat, ct.GetName()), **pdfOpts)
                      , logy = pl
                      , plotResidHist = False
                      ## , components = { 'wpv_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                      ##                  , 'prompt_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                      ##                  , 'sig_*'     : dict( LineColor = kOrange,  LineStyle = kDashed )
                      ##                  }
                      )
            pd.Delete()
            del pd
            for frame in ps:
                plot_name = '_'.join((t.GetName(), bins.GetName(), ct.GetName(), frame.GetName()))
                frame.SetName(plot_name)
            
            plots.append(ps)
    else:
        canvas = TCanvas('time_canvas_%d' % i, 'time_canvas_%d' % i, 600, 400)
        __canvases.append(canvas)
        p = canvas.cd(1)
        r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
        pdfOpts  = dict(ProjWData = (projSet, sig_sdata, True))
        ps = plot(p, t, pdf = time_pdf, data = sig_sdata
                  , frameOpts = dict(Range = r, Title = "")
                  , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                  , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                  , logy = pl
                  , plotResidHist = False)
        for frame in ps:
            plot_name = '_'.join((t.GetName(), bins.GetName(), frame.GetName()))
            frame.SetName(plot_name)
        
        plots.append(ps)

fit_results = []
from ROOT import TH1D, TGraphErrors
total = sig_sdata.sumEntries()
if options.simultaneous:
    split_bounds = array('d', [1000 * v for v in split_bins])
    
    res_canvas = TCanvas('res_canvas', 'res_canvas', 500, 500)
    
    hist_events = TH1D('hist_events', 'hist_events', len(split_bounds) - 1, split_bounds)
    if args[0] != 'MC11a':
        mass_fpf = sWeight_mass_result.floatParsFinal()
    time_fpf = time_result.floatParsFinal()
    
    res_x = array('d')
    res = array('d')
    res_e = array('d')
    for index, ct in enumerate(st_cat):
        d = split_bounds[index + 1] - split_bounds[index]
        if args[0] == 'MC11a':
            bin_name = '_'.join((signal.getYield().GetName(), ct.GetName()))
            events = sdatas[ct.GetName()].numEntries()
            hist_events.SetBinContent(index + 1, events / d)
            hist_events.SetBinError(index + 1, 0)
        else:
            bin_name = '_'.join((psi_ll.getYield().GetName(), ct.GetName()))
            events = mass_fpf.find(bin_name)
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
        
        range_cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName())
        mean = sig_sdata.mean(st._target_(), range_cut)
        mean *= total / sig_sdata.sumEntries(range_cut)
        res_x.append(1000 * mean)
        res.append(1000 * mean * sf)
        res_e.append(1000 * mean * sf_e)
    
    res_ex = array('d', [0 for i in range(len(res_x))])
    res_graph = TGraphErrors(len(res_x), res_x, res, res_ex, res_e)
    scale = 100 / hist_events.GetMaximum()
    hist_events.Scale(scale)
    hist_events.GetYaxis().SetRangeUser(0, 110)
        
    from ROOT import TF1
    fit_funcs = {'pol1' : 'S0+', 'pol2' : 'S+'}
    for i, (func, opts) in enumerate(fit_funcs.iteritems()):
        fit_func = TF1('fit_func_%s' % func , func, split_bins[0], split_bins[-1])
        fit_result = res_graph.Fit(fit_func, opts, "L")
        fit_results.append(fit_result)
        fr = fit_result.Get()
        fr.SetName('fit_result_%s_%s' % (func, args[1]))
        results.append(fr)
    
    res_graph.GetYaxis().SetRangeUser(0, 110)
    hist_events.Draw('hist')
    from ROOT import kGray
    hist_events.SetFillColor(kGray + 1)
    res_graph.Draw('P')
    res_graph.GetXaxis().SetTitle('estimated decay time error [fs]')
    res_graph.GetYaxis().SetTitle('decay time resulution [fs]')

from P2VV import Dilution
Dilution.dilution(t, data, result = time_result, sigmat = st, signal = [prompt],
                  subtract = [psi_ll, psi_wpv] if options.wpv else [psi_ll], simultaneous = options.simultaneous)

# Write data to cache file
def get_dir(d):
    tmp = cache_dir.Get(d)
    if not tmp:
        cache_dir.mkdir(d)
        tmp = cache_dir.Get(d)
    return tmp

from ROOT import TObject
if (options.write_data or fit_mass):
    sdata_dir = get_dir('sdata')
    data_dir = get_dir('data')
    for name, ds in sdatas_full.iteritems():
        sdata_dir.WriteTObject(ds, name, "Overwrite")

    if options.simultaneous:
        for ct in st_cat:
            opts = dict(Cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))
            bin_data = sig_sdata_full.reduce(**opts)
            bin_data.SetName('sig_sdata_%s' % ct.GetName())
            sdata_dir.WriteTObject(bin_data, bin_data.GetName(), "Overwrite")
            bin_data.Delete()
            del bin_data
    
    sdata_dir.Write(sdata_dir.GetName(), TObject.kOverwrite)

## Write PDFs
pdf_dir = get_dir('PDFs')
pdf_dir.WriteTObject(time_pdf._target_(), 'time_pdf_' + args[1] + \
                     ('_' + options.parameterise) if options.parameterise else '', "Overwrite")
pdf_dir.WriteTObject(mass_pdf._target_(), 'mass_pdf', "Overwrite")
if options.simultaneous:
    pdf_dir.WriteTObject(sWeight_mass_pdf._target_(), 'sWeight_mass_pdf', "Overwrite")
    
if not options.reduce:
    ## Write fit results
    results_dir = get_dir('results')
    for r in results:
        results_dir.WriteTObject(r, r.GetName(), "Overwrite")
    
    results_dir.Write(results_dir.GetName(), TObject.kOverwrite)
    
    ## Write plots
    plots_dir = get_dir('plots/%s' % extra_name)
    for ps in plots:
        for p in ps:
            plots_dir.WriteTObject(p, p.GetName(), "Overwrite")
    
    plots_dir.Write(plots_dir.GetName(), TObject.kOverwrite)

# Delete the input TTree which was automatically attached.
cache_file.Delete('%s;*' % tree_name)
