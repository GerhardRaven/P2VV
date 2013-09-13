#!/usr/bin/env python
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog year model')

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
parser.add_option("--sf-parameterisation", dest = "sf_param", default = '', type = 'string',
                  action = 'store', help = 'Type of model used to scale sigma_t')
parser.add_option('-s', "--simultaneous", dest = "simultaneous", default = False,
                  action = 'store_true', help = 'Use sigmat offset')
parser.add_option("--plot", dest = "make_plots", default = False,
                  action = 'store_true', help = 'Make plots')
parser.add_option("--no-fit", dest = "fit", default = True,
                  action = 'store_false', help = 'Do not do the time fit')
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
parser.add_option("--reweigh-MC", dest = "reweigh", default = False,
                  action = 'store_true', help = 'reweigh MC PV distribution to match data')
parser.add_option("--no-reuse-result", dest = "reuse_result", default = True,
                  action = 'store_false', help = "Don't reuse an old result if it is available.")
parser.add_option("--wpv-gauss-width", dest = "wpv_gauss_width", default = 0, type = 'float',
                  action = 'store', help = 'Width of the Gaussian used for the WPV component; 0 means floating.')
parser.add_option("--peak-only", dest = "peak_only", default = False, action = 'store_true',
                  help = 'Fit only the peak between -0.2 and 0.2 ps')
parser.add_option("--split-mean", dest = "split_mean", default = False, action = 'store_true',
                  help = 'Split the mean of the Gaussians in a simultaneous fit.')
parser.add_option("--split-frac", dest = "split_frac", default = False, action = 'store_true',
                  help = 'Split the fraction of the Gaussians in a simultaneous fit.')
parser.add_option("--split", dest = "split", default = 'sigmat', action = 'store', type = 'string',
                  help = 'Which categories should be used to split, [sigmat, momentum, pt, ppt]')
parser.add_option("--correct-errors", dest = "correct_errors", default = False, action = 'store_true',
                  help = 'Apply the SumW2 error correction')
parser.add_option("--add-background", dest = "add_background", default = False, action = 'store_true',
                  help = 'Add background to the time fit')

(options, args) = parser.parse_args()

if len(args) != 2:
    print parser.print_usage()
    sys.exit(-2)
elif args[1] not in ['single', 'double', 'triple']:
    print parser.print_usage()
    print "Wrong model type; allowed types are: %s." % ' '.join(['single', 'double', 'triple'])
    sys.exit(-2)

signal_MC = args[0] in ['MC11a', 'MC2012', 'MC2011_Sim08a']
prompt_MC = args[0] in ['MC11a_incl_Jpsi', 'MC2011_Sim08a_incl_Jpsi', 'MC2012_incl_Jpsi']

prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
input_data = {'2011' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011_ntupleB_20130722.root'),
                        'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_Mixing_2011_DataSet.root'),
                        'workspace' : 'Bs2JpsiPhiPrescaled_2011_workspace',
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled.root')},
              '2011_Reco14' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011_Reco14_ntupleB_20130906.root'),
                               'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_Mixing_2011_DataSet.root'),
                               'workspace' : 'Bs2JpsiPhiPrescaled_WPV_2011_workspace',
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Reco14_Prescaled.root')},                        
              '2012' : {'data' :os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_2013095.root'),
                        'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_Mixing_2012_DataSet.root'),
                        'workspace' : 'Bs2JpsiPhiPrescaled_2012_workspace',
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_Prescaled.root')},
              'MC11a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130613.root'),
                         'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC11a.root'),
                         'workspace' : 'Bs2JpsiPhiPrescaled_MC11a_workspace',
                         'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_Prescaled.root')},
              'MC11a_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_incl_Jpsi_ntupleB_20130801.root'),
                                   'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC11a.root'),
                                   'workspace' : 'Bs2JpsiPhiPrescaled_MC11a_workspace',
                                   'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_incl_Jpsi_Prescaled.root')},
              'MC2012' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_ntupleB_20130709_MagUpDown.root'),
                          'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC2012.root'),
                          'workspace' : 'Bs2JpsiPhiPrescaled_MC2012_workspace',
                          'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_Prescaled.root')},
              'MC2012_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_incl_Jpsi_ntupleB_20130709_MagUpDown.root'),
                                    'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC2012.root'),
                                    'workspace' : 'Bs2JpsiPhiPrescaled_MC2012_workspace',
                                    'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_incl_Jpsi_Prescaled.root')},
              'MC2011_Sim08a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_ntupleB_20130909.root'),
                                 'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC2011_Sim08a.root'),
                                 'workspace' : 'Bs2JpsiPhiPrescaled_MC2011_Sim08a_workspace',
                                 'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_Prescaled.root')},
              'MC2011_Sim08a_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_incl_Jpsi_ntupleB_20130909.root'),
                                           'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC2011_Sim08a.root'),
                                           'workspace' : 'Bs2JpsiPhiPrescaled_MC2011_Sim08a_workspace',
                                           'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_incl_Jpsi_Prescaled.root')}
              }

if args[0] not in input_data.keys():
    print parser.print_usage()
    print "Possible samples are: %s" % ' '.join(input_data.keys())
    sys.exit(-2)

if options.wpv and not options.wpv_type in ['Mixing', 'Gauss']:
    print parser.print_usage()
    print "Wring mixing type; allowed types: %s" % ' '.join(['Mixing', 'Gauss'])
    sys.exit(-2)

if options.batch:
    from ROOT import gROOT
    gROOT.SetBatch(True)
from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
if options.peak_only:
    t_minmax = (-0.08, 0.08)
elif options.wpv and options.wpv_type == 'Gauss':
    t_minmax = (-1.5, 14)
else:
    t_minmax = (-5, 14)
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3020, 3165))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))

# add 20 bins for caching the normalization integral
st.setBins(500, 'cache')

# Categories needed for selecting events
unbiased = Category('triggerDecisionUnbiasedPrescaled', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'selected' : 1, 'not_selected' : 0})
momentum = RealVar('B_P', Title = 'B momentum', Unit = 'MeV', Observable = True, MinMax = (0, 1e6))
pt = RealVar('B_Pt', Title = 'B transverse momentum', Unit = 'MeV', Observable = True, MinMax = (0, 1e6))
nPV = RealVar('nPV', Title = 'Number of PVs', Observable = True, MinMax = (0, 10))
zerr = RealVar('B_s0_bpv_zerr', Title = 'Best PV Z error', Unit = 'mm', Observable = True, MinMax = (0, 1))

observables = [t, m, mpsi, st, unbiased, selected, nPV, momentum, pt, zerr]
if signal_MC:
    t_true = RealVar('truetime', Title = 'true decay time', Unit='ps', Observable = True, MinMax=(-1100, 14))
    t_diff = RealVar('time_diff', Unit = 'ps', Observable = True, MinMax = (-1, 1))
    observables.extend([t_true, t_diff])

# Define a time_obs symbol to wrap around selection of t or t_diff
time_obs = t_diff if signal_MC else t
    
# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# Time resolution model
sig_tres = None
if args[1] == 'single':
    from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
    tres_args = dict(time = t, sigmat = st, PerEventError = options.pee,
                     BiasScaleFactor = False, Cache = False,
                     TimeResSFParam = options.sf_param, SplitMean = options.split_mean)
    sig_tres = TimeResolution(Name = 'tres', **tres_args)
    bkg_tres = TimeResolution(Name = 'bkg_tres', ParNamePrefix = 'bkg', **tres_args)
elif args[1] == 'double':
    mu = dict(MinMax = (-0.010, 0.010))
    mu['Constant'] = options.simultaneous and not options.split_mean
    mu_values = {'MC11a_incl_Jpsi' : -0.000408, '2011_Reco14' :-0.00468,  
                 '2011' : -0.00407301, '2012' : -0.00365}
    mu['Value'] = mu_values.get(args[0], 0)
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    tres_args = dict(time = time_obs, sigmat = st, Cache = True,
                              PerEventError = options.pee, Parameterise = options.parameterise,
                              TimeResSFParam = options.sf_param, SplitFracs = options.split_frac,
                              timeResMu = mu, GExp = {2 : signal_MC, 1 : False},
                              ScaleFactors = [(2, 2.1), (1, 1.26)] if options.pee else [(2, 0.1), (1, 0.06)],
                              Fractions = [(2, 0.2)], SplitMean = options.split_mean)
    sig_tres = TimeResolution(Name = 'sig_tres', **tres_args)
    bkg_tres = TimeResolution(Name = 'bkg_tres', ParNamePrefix = 'bkg', **tres_args)
elif args[1] == 'triple':
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = time_obs, sigmat = st, Cache = True,
                              PerEventError = options.pee, GExp = {3 : signal_MC, 2: False, 1 : False},
                              TimeResSFParam = options.sf_param, Parameterise = options.parameterise,
                              ScaleFactors = [(3, 1.5), (2, 4), (1, 1.4)], SplitFracs = options.split_frac,
                              Fractions = [(3, 0.1), (2, 0.2)], SplitMean = options.split_mean)

# J/psi mass pdf
from P2VV.Parameterizations.MassPDFs import DoubleCB_Psi_Mass as PsiMassPdf
psi_m = PsiMassPdf(mpsi, Name = 'psi_m')
psi_m = psi_m.pdf()
    
# J/psi background
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf
bkg_mpsi = PsiBkgPdf(mpsi, Name = 'bkg_mpsi')

# Create combinatorical background component
mean = RealVar(Name = 'mean', Value = 5368, MinMax = (5300, 5400))
sigma = RealVar(Name = 'sigma', Value = 50, MinMax = (1, 500))
bkg_m = Pdf(Name = 'gauss', Type = Gaussian, Parameters = (m, mean, sigma))

# Create psi background component
from P2VV.Parameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = time_obs, resolutionModel = sig_tres.model()
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 6.7195e-01)
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.3672, MinMax = (0.1,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 1.3405e-01, MinMax = (0.01, 0.5))
                         )
psi_t = psi_t.pdf()

# J/psi signal component
psi_ll = Component('psi_ll', (psi_m, bkg_m, psi_t), Yield= (8.5575e+03,100,500000) )

# Background component
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
bkg_t = Prompt_Peak(time_obs, bkg_tres.model(), Name = 'bkg_pdf')
background = Component('background', (bkg_mpsi.pdf(), bkg_m, bkg_t.pdf()), Yield = (19620,100,500000) )

# B signal component
sig_t = Prompt_Peak(time_obs, resolutionModel = sig_tres.model(), Name = 'sig_t')

from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_Mass
sig_m = Signal_Mass(Name = 'sig_m', mass = m) 
signal = Component('signal', (sig_m.pdf(), sig_t.pdf()), Yield = (1e6, 10000, 1e7))

# Prompt component
prompt_pdf = Prompt_Peak(time_obs, sig_tres.model(), Name = 'prompt_pdf')
prompt = Component('prompt', (prompt_pdf.pdf(), psi_m), Yield = (160160, 100, 500000))

# Read data
fit_mass = options.fit_mass

# Tree and cut
tree_name = 'DecayTree'
cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
if not options.wpv or (options.wpv and options.wpv_type == "Gauss"):
    cut += ' && sel_cleantail == 1'
if signal_MC:
    cut += ' && abs(trueid) == 531'
if options.peak_only:
    cut += ' && time > %5.3f && time < %5.3f' % t_minmax
if options.cut:
    cut = options.cut + ' && ' + cut
hd = ('%d' % hash(cut)).replace('-', 'm')

if options.simultaneous:
    split_utils = {'sigmat' : ('Sigmat', [st]), 'ppt' : ('PPT', [momentum, pt]),
                   'momentum' : ('Momentum', [momentum]), 'pt' : ('PT', [pt]),
                   'pv_zerr' : ('PVZerr', [zerr])}
    from P2VV.Utilities import Resolution as ResolutionUtils
    split_opts = split_utils[options.split]
    SplitUtil = getattr(ResolutionUtils, 'Split' + split_opts[0])
    split_util = SplitUtil(args[0], *(split_opts[1]))
    directory = split_util.directory(hd)
else:
    directory = '1bin_%4.2ffs_simple/%s' % (1000 * (t.getMax() - t.getMin()), hd)

from P2VV.CacheUtils import CacheFiles
cache_files = CacheFiles(*input_data[args[0]]['cache'].rsplit('/', 1))
cache_dir, cache_file = cache_files.getFromCache(directory)
if not cache_dir:
    fit_mass = True

results = []
tree_name = 'DecayTree'

## Canvas for correlation histograms
from ROOT import TCanvas
corr_canvas = TCanvas('corr_canvas', 'corr_canvas', 1000, 500)
corr_canvas.Divide(2, 1)
from ROOT import gStyle
gStyle.SetPalette(53)

## Extra name for fit result and plots
extra_name = [args[1]]
for a, n in [('parameterise', None), ('wpv', 'wpv_type'), ('sf_param', None),
             ('peak_only', 'peak_only'), ('add_background', 'cfit')]:
    v = getattr(options, a)
    if v:
        if n and n != v and hasattr(options, n):
            n = getattr(options, n)
        extra_name.append(n if n else v)

## Read Cache
if not fit_mass:
    ## Read sdata
    sds_name = 'sdata'
    sdata_dir = cache_dir.Get('sdata')
    sdatas_full = {}
    if not sdata_dir:
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
            single_bin_sig_sdata = sdatas_full['sdata/full_sig_sdata']
            single_bin_bkg_sdata = sdatas_full['sdata/full_bkg_sdata']
            if options.simultaneous:
                sig_sdata_full = sdatas_full['sdata/sig_sdata']
                bkg_sdata_full = sdatas_full['sdata/bkg_sdata']
            else:
                sig_sdata_full = single_bin_sig_sdata
                bkg_sdata_full = single_bin_bkg_sdata
        except KeyError:
            fit_mass = True
        ## CHECK ST BINNING
            
    from copy import copy
    sdatas = copy(sdatas_full)
    if options.reduce and not fit_mass:
        sig_sdata = sig_sdata_full.reduce(EventRange = (0, 40000))
        bkg_sdata = bkg_sdata_full.reduce(EventRange = (0, 40000))        
    elif not fit_mass:
        sig_sdata = sig_sdata_full
        bkg_sdata = bkg_sdata_full

    # Read data
    data_dir = cache_dir.Get('data')
    if data_dir.GetListOfKeys() == 1:
        data = cache_dir.get(os.path.join(data_dir.GetName(), sdata_dir.GetListOfKeys()))
    
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
        tr = rd.Get('_'.join(['time_result'] + extra_name))
        if tr:
            results.append(tr)

    if options.simultaneous:
        split_cats = [split_util.split_cats(sig_sdata)]
        
## Fitting opts
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True,
               Verbose = options.verbose, Strategy = 1)

## List of all plots we make
plots = []

# PV bins
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

## Build simple mass pdf
if fit_mass:
    from P2VV.Utilities.DataHandling import readData
    data = readData(input_data[args[0]]['data'], tree_name, NTuple = True, observables = observables,
                    ntupleCuts = cut)
    data.SetName(tree_name)
    if options.reduce:
        data = data.reduce(EventRange = (0, int(4e4)))
    elif args[0] == '2012':
        data = data.reduce(EventRange = (0, int(6e5)))
        
    # In case of reweighing 
    sig_mass_pdf = buildPdf(Components = (signal, background), Observables = (m,), Name = 'sig_mass_pdf')
    psi_mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='psi_mass_pdf')
    if signal_MC:
        sig_mass_pdf.Print('t')
        signal_name = signal.GetName()
        mass_pdf = sig_mass_pdf
    else:
        psi_mass_pdf.Print('t')
        signal_name = psi_ll.GetName()
        mass_pdf = psi_mass_pdf
        
    ## Fit mass pdf
    for i in range(3):
        mass_result = mass_pdf.fitTo(data, **fitOpts)
        if mass_result.status() == 0:
            break
    assert(mass_result.status() == 0)
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

    mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 600, 530)
    from P2VV.Utilities.SWeights import SData
    from P2VV.Utilities.Plotting import plot
    pdfOpts  = dict()
    if signal_MC:
        mass_obs = m
    else:
        mass_obs = mpsi
    ps = plot(mass_canvas.cd(1), mass_obs, pdf = mass_pdf, data = data
              , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = 50)
              , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
              , plotResidHist = 'BX'
              , xTitle = '#mu^+#mu^- invariant mass [MeV/c^2]'
              , yTitle = 'Candidates / (%f MeV/c^2)' % ((mass_obs.getMax() - mass_obs.getMin()) / float(50))
              , components = { 'sig_*'     : dict( LineColor = kOrange,   LineStyle = kDashed )
                               , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                               , 'bkg_*'  : dict( LineColor = kRed, LineStyle = kDashed )
                               }
              )
    plots.append(ps)
    
    from P2VV.Utilities.SWeights import SData
    data_clone = data.Clone(data.GetName())
    sData = SData(Pdf = mass_pdf, Data = data_clone, Name = 'MassSPlot')
    single_bin_sig_sdata = sData.data(signal_name)
    single_bin_bkg_sdata = sData.data('background')

    sdatas_full = {'full_sig_sdata' : single_bin_sig_sdata,
                   'full_bkg_sdata' : single_bin_bkg_sdata}

if fit_mass and options.simultaneous:
    from P2VV.Utilities.General import getSplitPar
    # categories for splitting the PDF
    # get mass parameters that are split
    split_cats = [split_util.split_cats(data)]
    split_pars = [[par for par in mass_pdf.Parameters() if par.getAttribute('Yield')]]
    
    # build simultaneous mass PDF
    sWeight_mass_pdf = SimultaneousPdf(mass_pdf.GetName() + '_simul',
                                       MasterPdf       = mass_pdf,
                                       SplitCategories = split_cats,
                                       SplitParameters = split_pars)
    
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
    
    from P2VV.Utilities.SWeights import SData
    sData = SData(Pdf = sWeight_mass_pdf, Data = data, Name = 'SimulMassSPlot')
    sig_sdata_full = sData.data(signal_name)
    bkg_sdata_full = sData.data('background')
    
    if (signal_MC or prompt_MC) and options.reweigh:
        if cut.find('trueid') != -1:
            cut = cut.rsplit('&&', 1)[0]
        reco_data = readData(input_data['2011']['data'], tree_name,
                             NTuple = True, observables = observables[:-1],
                             ntupleCuts = cut, Rename = tree_name + '_2011')
        # We build a mass PDF always wrt the J/psi, since there is not enough
        # signal to do it with. Let's hope this is OK...
        reco_mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='reco_mass_pdf')
        reco_mass_result = reco_mass_pdf.fitTo(reco_data, **fitOpts)
        reco_mass_result.SetName('reco_mass_result')
        reco_sData = SData(Pdf = reco_mass_pdf, Data = reco_data, Name = 'RecoMassSPlot')
        reco_sig_sdata = reco_sData.data(psi_ll.GetName())
        reco_bkg_sdata = reco_sData.data('background')
        from P2VV.Reweighing import reweigh
        for target, source, n in [(sig_sdata_full, reco_sig_sdata, 'full_sig_sdata'),
                                  (bkg_sdata_full, reco_bkg_sdata, 'full_bkg_sdata')]:
            ds, weights = reweigh(target, 'nPV', source, 'nPV', binning = PV_bounds)
            sdatas_full[n] = ds
    else:
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
    if (signal_MC or prompt_MC) and options.reweigh:
        reco_data = readData(input_data['2011']['data'], tree_name,
                             NTuple = True, observables = observables[:-1],
                             ntupleCuts = cut, Rename = tree_name + '_2011')
        reco_mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='reco_mass_pdf')
        reco_mass_result = reco_mass_pdf.fitTo(reco_data, **fitOpts)
        reco_mass_result.SetName('reco_mass_result')
        reco_sData = SData(Pdf = reco_mass_pdf, Data = reco_data, Name = 'RecoMassSPlot')
        reco_sig_sdata = reco_sData.data(psi_ll.GetName())
        reco_bkg_sdata = reco_sData.data('background')
        from P2VV.Reweighing import reweigh
        reweigh_sdatas = {}
        for target, source, n in [(single_bin_sig_sdata, reco_sig_sdata, 'sig_sdata'),
                                  (single_bin_bkg_sdata, reco_bkg_sdata, 'bkg_sdata')]:
            ds, weights = reweigh(target, 'nPV', source, 'nPV', binning = PV_bounds)
            reweigh_sdatas[n] = ds
        sig_sdata = reweigh_sdatas['sig_sdata']
        bkg_sdata = reweigh_sdatas['bkg_sdata']
    else:
        sig_sdata = single_bin_sig_sdata
        bkg_sdata = single_bin_bkg_sdata

corr_canvas.Update()

# Define default components
if signal_MC:
    components = [signal]
else:
    components = [prompt, psi_ll]

if options.peak_only:
    components = components[:1]

if options.add_background:
    components += [background]

if options.wpv and options.wpv_type == 'Mixing':
    from P2VV.Parameterizations import WrongPV
    if signal_MC:
        reweigh_data = dict(B = single_bin_sig_sdata, bkg = single_bin_bkg_sdata)
        masses = {'B' : m}
        weights = 'B' 
        extra_args = {'t_diff' : t_diff}
    else:
        reweigh_data = dict(jpsi = single_bin_sig_sdata, bkg = single_bin_bkg_sdata)
        masses = {'jpsi' : mpsi}
        weights = 'jpsi'
        extra_args = {}
    wpv_builder = WrongPV.ShapeBuilder(t, masses, UseKeysPdf = True, Weights = weights, Draw = True,
                                       InputFile = input_data[args[0]]['wpv'], Workspace = input_data[args[0]]['workspace'],
                                       Reweigh = dict(Data = reweigh_data, DataVar = nPV, Binning = PV_bounds),
                                       sigmat = st, MassResult = mass_result, **extra_args)
    if signal_MC:
        wpv_signal = wpv_builder.shape('B')
        sig_wpv = Component('sig_wpv', (wpv_signal, m), Yield = (888, 50, 300000))
    else:
        wpv_psi = wpv_builder.shape('jpsi')
        sig_wpv = Component('sig_wpv', (wpv_psi, psi_m), Yield = (1000, 50, 30000))
    components += [sig_wpv]
    ## if options.add_background:
    ##     wpv_bkg = wpv_builder.shape('bkg')
    ##     bkg_wpv = Component('wpv', (wpv_psi, bkg_m if signal_mc else bkg_mpsi), Yield = (1000, 50, 30000))
    ##     componets += [bkg_wpv]
elif options.wpv and options.wpv_type == 'Gauss':
    ## if not signal_MC:
    ##     wpv_mean = sig_tres._timeResMu
    ## else:
    def make_wpv_pdf(prefix):
        wpv_mean = RealVar('%swpv_mean' % prefix, Value = 0, MinMax = (-1, 1))
        if options.wpv_gauss_width != 0:
            wpv_sigma = RealVar('%swpv_sigma' % prefix, Value = options.wpv_gauss_width, MinMax = (0.01, 1000), Constant = True)
        else:
            wpv_sigma = RealVar('%swpv_sigma' % prefix, Value = 0.305, MinMax = (0.01, 10))
        return Pdf(Name = '%swpv_pdf' % prefix, Type = Gaussian, Parameters = (time_obs, wpv_mean, wpv_sigma))
    sig_wpv_pdf = make_wpv_pdf('sig_')
    sig_wpv = Component('sig_wpv', (sig_wpv_pdf, psi_m), Yield = (552, 1, 50000))
    components += [sig_wpv]
    if options.add_background:
        bkg_wpv_pdf = make_wpv_pdf('bkg_')
        bkg_wpv = Component('bkg_wpv', (bkg_wpv_pdf, bkg_mpsi.pdf()), Yield = (552, 1, 50000))
        components += [bkg_wpv]
if options.add_background:
    pdf_obs = (time_obs, mpsi)
else:
    pdf_obs = (time_obs,)
time_pdf = buildPdf(Components = components, Observables = pdf_obs, Name='time_pdf')

splitLeaves = RooArgSet()

if options.simultaneous:
    split_pars = [[]]
    if options.wpv:
        for par in time_pdf.Parameters():
            if not par.getAttribute('Yield'):
                continue
            if par.GetName().find(components[-1].GetName()) != -1:
                continue
            split_pars[0].append(par)
    else:
        split_pars[0] = [par for par in time_pdf.Parameters() if par.getAttribute('Yield')]
    split_pars[0] += sig_tres.splitVars()
    
    if options.wpv and options.wpv_type == 'Gauss':    
        split_pars[0].extend([sig_wpv.getYield()])
        ## split_pars[0].extend([wpv_sigma])

    if options.sf_param:
        assert(len(split_cats[0]) == 1)
        split_cat = split_cats[0][0]
        if hasattr(split_cat, '_target_'):
            split_cat = split_cat._target_()
        placeholder = sig_tres.sigmatPlaceHolder()
        from ROOT import RooCustomizer
        customizer = RooCustomizer(time_pdf._target_(), split_cat, splitLeaves)
        to_split = RooArgSet(*(split_pars[0] + [placeholder._target_()]))
        customizer.splitArgs(to_split, split_cat)
        states = dict([(s.GetName(), customizer.build(s.GetName())) for s in split_cat])
        time_pdf = SimultaneousPdf(time_pdf.GetName() + '_simul',
                                   SplitCategory = split_cat, ExternalConstraints = time_pdf.ExternalConstraints(),
                                   States = states, ConditionalObservables = time_pdf.ConditionalObservables())
        means = [sig_sdata.mean(st._target_(), '{0} == {0}::{1}'.format(split_cat.GetName(), s.GetName())) for s in split_cat]
        pars = time_pdf.getParameters(RooArgSet())
        st_mean = sig_sdata.mean(st._target_())
        for m, s in zip(means, split_cat):
            p = pars.find(placeholder.GetName() + '_' + s.GetName())
            p.setConstant(True)
            p.setVal(m - st_mean)
            
    else:
        time_pdf = SimultaneousPdf(time_pdf.GetName() + '_simul'
                                   , MasterPdf       = time_pdf
                                   , SplitCategories = split_cats
                                   , SplitParameters = split_pars)
if options.reuse_result:
    # Check if we have a cached time result, if so, use it as initial values for the fit
    time_result = None
    for i, r in enumerate(results):
        if r.GetName() == '_'.join(['time_result'] + extra_name):
            time_result = results.pop(i)
            break
    
    if time_result:
        pdf_vars = time_pdf.getVariables()
        for p in time_result.floatParsFinal():
            pdf_par = pdf_vars.find(p.GetName())
            if pdf_par and not pdf_par.isConstant():
                pdf_par.setVal(p.getVal())
                pdf_par.setError(p.getError())

time_pdf.Print("t")

## Fit
## print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
obs_arg = RooArgSet()
for o in observables:
    obs_arg.add(o._target_())
parameters = dict([(p.GetName(), p) for p in time_pdf.getParameters(obs_arg)])

## from P2VV import Dilution
## comb = parameters['timeResComb']
## frac2 = parameters['timeResFrac2']
## sf2 = parameters['timeResSigmaSF_2']
## sf1 = (comb.getVal() - frac2.getVal() * sf2.getVal()) / (1 - frac2.getVal())
## Dilution.signal_dilution_dg(sig_sdata, st, sf1, frac2.getVal(), sf2.getVal())

## assert(False)

if options.add_background:
    fit_data = data
    fit_data_full = data
else:
    fit_data = sig_sdata
    if options.simultaneous:
        fit_data_full = sig_sdata_full

if options.fit:
    for i in range(3):
        time_result = time_pdf.fitTo(fit_data, SumW2Error = options.correct_errors, **fitOpts)
        if time_result.status() == 0:
            break    
    time_result.SetName('_'.join(['time_result'] + extra_name))
    
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

def cut_binning(t, binning):
    ## find the first boundary which is equal to or larger than t_min
    mm = [None, None]
    for i in range(len(binning)):
        if binning[i] >= t.getMin() and mm[0] == None:
            mm[0] = i
        if binning[i] == t.getMax():
            mm[1] = i
            break
        if binning[i] >= t.getMax():
            mm[1] = i - 1
            break
    else:
        mm[0] = 0
    if mm[1] == None:
        mm[1] = len(binning)
    return binning[mm[0] : mm[1]]

from ROOT import RooBinning

if options.wpv and options.wpv_type == 'Mixing':
    bounds = array('d', [-5 + i * 0.1 for i in range(48)] + [-0.2 + i * 0.01 for i in range(40)] + \
                   [0.2 + i * 0.1 for i in range(58)])
    zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])
elif signal_MC:
    bounds = array('d', [-1.5 + i * 0.1 for i in range(12)] + [-0.3 + i * 0.05 for i in range(12)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])
    zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])
else:
    bounds = array('d', [-1.5 + i * 0.1 for i in range(12)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])
    zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])

bounds = cut_binning(time_obs, bounds)
zoom_bounds = cut_binning(time_obs, zoom_bounds)
    
binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('full')

zoom_binning = RooBinning(len(zoom_bounds) - 1, zoom_bounds)
zoom_binning.SetName('zoom')

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from P2VV.Utilities.Plotting import plot

print 'plotting'

binnings = [binning, zoom_binning]
plotLog = [True, False]
__canvases = []

from ROOT import SetOwnership

for i, (bins, pl) in enumerate(zip(binnings, plotLog)):
    if not options.make_plots or not time_result or options.split != 'sigmat':
        continue
    if options.simultaneous:
        st_cat = split_cats[0][0]
        r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
        for ct in st_cat:
            name = 'time_canvas_%s_%d' % (ct.GetName(), i)
            canvas = TCanvas(name, name, 600, 400)
            __canvases.append(canvas)
            p = canvas.cd(1)
            
            projSet = RooArgSet(st, time_pdf.indexCat())
            pdfOpts  = dict(Slice = (st_cat, ct.GetName()), ProjWData = (projSet, fit_data_full, True))
            ps = plot(p, time_obs, pdf = time_pdf, data = fit_data_full
                      , frameOpts = dict(Range = r, Title = "")
                      , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack, Cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))
                      , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                      , xTitle = 'decay time [ps]'
                      , yTitle = 'Candidates / (XX ps)'
                      , logy = pl
                      , plotResidHist = 'BX'
                      ## , components = { 'wpv_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                      ##                  , 'prompt_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                      ##                  , 'sig_*'     : dict( LineColor = kOrange,  LineStyle = kDashed )
                      ##                  }
                      )
            for frame in ps:
                plot_name = '_'.join((t.GetName(), bins.GetName(), ct.GetName(), frame.GetName()))
                frame.SetName(plot_name)
            
            plots.append(ps)
    else:
        canvas = TCanvas('time_canvas_%d' % i, 'time_canvas_%d' % i, 600, 533)
        __canvases.append(canvas)
        p = canvas.cd(1)
        r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
        projSet = RooArgSet(st)
        pdfOpts  = dict(ProjWData = (projSet, fit_data, True))
        ps = plot(p, time_obs, pdf = time_pdf, data = fit_data
                  , frameOpts = dict(Range = r, Title = "")
                  , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                  , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                  , xTitle = 'decay time [ps]'
                  , yTitle = 'Candidates / (XX ps)'
                  , logy = pl
                  , plotResidHist = 'BX')
        
        for frame in ps:
            plot_name = '_'.join((t.GetName(), bins.GetName(), frame.GetName()))
            frame.SetName(plot_name)
        
        plots.append(ps)

from P2VV.Utilities import Resolution as ResolutionUtils
time_result.PrintSpecial(LaTeX = True, ParNames = ResolutionUtils.parNames)

from P2VV.CacheUtils import WritableCacheFile
with WritableCacheFile(cache_files, directory) as cache_file:
    cache_dir = cache_file.Get(directory)
    from ROOT import TObjString
    cut_string = TObjString(cut)
    cache_dir.WriteTObject(cut_string, 'cut')
    
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
        
        ## if options.simultaneous:
        ##     for ct in st_cat:
        ##         opts = dict(Cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))
        ##         bin_data = sig_sdata_full.reduce(**opts)
        ##         bin_data.SetName('sig_sdata_%s' % ct.GetName())
        ##         sdata_dir.WriteTObject(bin_data, bin_data.GetName(), "Overwrite")
        ##         bin_data.Delete()
        ##         del bin_data
        
        sdata_dir.Write(sdata_dir.GetName(), TObject.kOverwrite)
        data_dir.WriteTObject(data, data.GetName(), "Overwrite")
        data_dir.Write(data_dir.GetName(), TObject.kOverwrite)
    
    ## Write PDFs
    pdf_dir = get_dir('PDFs')
    pdf_dir.WriteTObject(time_pdf._target_(), 'time_pdf_' + args[1] + \
                         ('_' + options.parameterise) if options.parameterise else '', "Overwrite")
    if (options.write_data or fit_mass):
        pdf_dir.WriteTObject(mass_pdf._target_(), 'mass_pdf', "Overwrite")
        if options.simultaneous:
            pdf_dir.WriteTObject(sWeight_mass_pdf._target_(), 'sWeight_mass_pdf', "Overwrite")
    
    if not options.reduce and options.fit:
        ## Write fit results
        results_dir = get_dir('results')
        for r in results:
            results_dir.WriteTObject(r, r.GetName(), "Overwrite")
        
        results_dir.Write(results_dir.GetName(), TObject.kOverwrite)
        
        ## Write plots
        plots_dir = get_dir('plots/%s' % '_'.join(extra_name))
        for ps in plots:
            for p in ps:
                plots_dir.WriteTObject(p, p.GetName(), "Overwrite")
        
        plots_dir.Write(plots_dir.GetName(), TObject.kOverwrite)
    
    # Delete the input TTree which was automatically attached.
    cache_file.Delete('%s;*' % tree_name)

