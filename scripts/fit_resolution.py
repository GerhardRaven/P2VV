#!/usr/bin/env python
import gc
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog year model')

def mb_callback(option, opt_str, value, parser):
    mb = tuple(int(i) for i in value.split(','))
    setattr(parser.values, option.dest, mb)

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
parser.add_option("--mean-parameterisation", dest = "mu_param", default = '', type = 'string',
                  action = 'store', help = 'Type parameterise the mean of the Gaussians as a function of sigmat')
parser.add_option('-s', "--simultaneous", dest = "simultaneous", default = False,
                  action = 'store_true', help = 'Use sigmat offset')
parser.add_option("--plot", dest = "plot", default = False,
                  action = 'store_true', help = 'Make plots')
parser.add_option("--no-fit", dest = "fit", default = True,
                  action = 'store_false', help = 'Do not do the time fit')
parser.add_option("--fit-mass", dest = "fit_mass", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')
parser.add_option("--force-write", dest = "write_data", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')
parser.add_option("--reduce", dest = "reduce", default = 0, type = 'int',
                  action = 'store', help = 'Reduce sdata sets to size N.')
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
parser.add_option("--split-observable", dest = "split", default = 'sigmat', action = 'store', type = 'string',
                  help = 'Which categories should be used to split, [sigmat, momentum, pt, ppt, nPV, pv_zerr]')
parser.add_option("--correct-errors", dest = "correct_errors", default = False, action = 'store_true',
                  help = 'Apply the SumW2 error correction')
parser.add_option("--add-background", dest = "add_background", default = False, action = 'store_true',
                  help = 'Add background to the time fit')
parser.add_option("--use-refit", dest = "use_refit", default = False, action = 'store_true',
                  help = 'Use new PV refitting.')
parser.add_option("--no-cache", dest = "cache", default = True, action = 'store_false',
                  help = 'Use a cache to store results and reuse them.')
parser.add_option("--constrain", dest = "constrain", default = '', action = 'store', type = 'string',
                  help = 'Which parameters to constrain')
parser.add_option("--write-constraints", dest = "write_constraints", default = '',
                  action = 'store', type = 'string', help = 'Write contraints to database.')
parser.add_option("--mass-parameterisation", dest = "mass_param", default = '',
                  action = 'store', type = 'string', help = 'Reparameterise the mass PDF')
parser.add_option("--make-binning", action="callback", callback=mb_callback, type = 'string',
                  dest = "make_binning", default = (0,), help = 'Make binning with n_bins for each ' +
                  'observable; comma separated list')
parser.add_option("--excl-biased", action="store_true", dest = "excl_biased", default = False,
                  help = 'Use excl biased events, only for signal MC')
parser.add_option("--cleantail", action="store_true", dest = "cleantail", default = False,
                  help = 'Always apply cleantail cut.')

(options, args) = parser.parse_args()

if len(args) != 2:
    print parser.print_usage()
    sys.exit(-2)
elif args[1] not in ['single', 'double', 'triple', 'core']:
    print parser.print_usage()
    print "Wrong model type; allowed types are: %s." % ' '.join(['single', 'double', 'triple', 'core'])
    sys.exit(-2)

signal_MC = args[0] in ['MC11a', 'MC2012', 'MC2011_Sim08a']
prompt_MC = args[0] in ['MC11a_incl_Jpsi', 'MC2011_Sim08a_incl_Jpsi', 'MC2012_incl_Jpsi']
pos = args[0].find('201')
year = args[0][pos : pos + 4]

from P2VV.Utilities.Resolution import input_data

if args[0] not in input_data.keys():
    print parser.print_usage()
    print "Possible samples are: %s" % ' '.join(input_data.keys())
    sys.exit(-2)

if options.wpv and not options.wpv_type in ['Mixing', 'Gauss', 'Rest']:
    print parser.print_usage()
    print "Wring mixing type; allowed types: %s" % ' '.join(['Mixing', 'Gauss'])
    sys.exit(-2)

if options.mass_param and not options.mass_param in ['MeanSigma']:
    print parser.print_usage()
    print "Wring mass parameterisation; allowed types: %s" % ' '.join(['MeanSigma'])
    sys.exit(-2)

if options.batch:
    from ROOT import gROOT
    gROOT.SetBatch(True)
from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle

## Extra name for fit result and plots
extra_name = [args[1]]
for a, n in [('parameterise', None), ('wpv', 'wpv_type'), ('sf_param', None),
             ('peak_only', 'peak'), ('add_background', 'cfit'),
             ('use_refit', 'PVRefit'), ('split_mean', 'split_mean')]:
    v = getattr(options, a)
    if v:
        if n and n != a and n != v and hasattr(options, n):
            n = getattr(options, n)
        extra_name.append(n if n else v)

obj = RooObject( workspace = 'w')
w = obj.ws()

## Code for using RooTrace
## from ROOT import RooTrace
## RooTrace.active(True)
## RooTrace.printObjectCounts()

from math import pi
if options.peak_only:
    if signal_MC:
        t_minmax = (-5, 14)
        tdiff_minmax = (-0.1, 0.1)
    else:
        t_minmax = (-0.08, 0.08)
elif options.wpv and options.wpv_type == 'Gauss':
    t_minmax = (-1.5, 14)
else:
    t_minmax = (-5, 14)
t  = RealVar('time' if not options.use_refit else 'time_refit', Title = 't', Unit='ps', Value = 1, Observable = True, MinMax = t_minmax)
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3025, 3165))
st = RealVar('sigmat' if not options.use_refit else 'sigmat_refit',Title = '#sigma_{t}', Unit = 'ps', Observable = True, Value = 0.03, MinMax = (0.0001, 0.12))

# add 20 bins for caching the normalization integral
st.setBins(50, 'cache')

# Categories needed for selecting events
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'selected' : 1, 'not_selected' : 0})
momentum = RealVar('B_P', Title = 'B momentum', Unit = 'MeV', Observable = True, MinMax = (0, 1e6))
pt = RealVar('B_Pt', Title = 'B transverse momentum', Unit = 'MeV', Observable = True, MinMax = (0, 1e6))
nPV = RealVar('nPV', Title = 'Number of PVs', Observable = True, MinMax = (0, 10))
zerr = RealVar('B_s0_bpv_zerr', Title = 'Best PV Z error', Unit = 'mm', Observable = True, MinMax = (0, 1))

t_st = RealVar('time_sigmat' if not options.use_refit else 'time_sigmat_refit', Title = 'time / sigmat', Observable = True, MinMax = (-100, 500))

observables = [t, t_st, m, mpsi, st, hlt1_unbiased, hlt2_unbiased, selected, nPV, momentum, pt, zerr]
if signal_MC:
    t_true = RealVar('truetime', Title = 'true decay time', Unit='ps', Observable = True, MinMax=(-1100, 14))
    t_diff = RealVar('time_diff' if not options.use_refit else 'time_diff_refit', Unit = 'ps', Observable = True, MinMax = (-0.5, 0.5) if not options.peak_only else tdiff_minmax)
    t_diff_st = RealVar('time_diff_sigmat' if not options.use_refit else 'time_diff_sigmat_refit', Title = 'time / sigmat', Observable = True, MinMax = (-30, 40))
    observables.extend([t_true, t_diff, t_diff_st])

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
    tres_args = dict(time = time_obs, sigmat = st, PerEventError = True,
                     BiasScaleFactor = False, Cache = True, MeanParameterisation = options.mu_param,
                     Simultaneous = options.simultaneous,
                     TimeResSFParam = options.sf_param, SplitMean = options.split_mean)
    sig_tres = TimeResolution(Name = 'tres', **tres_args)
    if options.add_background:
        bkg_tres = TimeResolution(Name = 'bkg_tres', ParNamePrefix = 'bkg', **tres_args)
elif args[1] == 'double':
    mu = dict(MinMax = (-0.010, 0.010))
    mu_values = {'MC11a_incl_Jpsi' : -0.000408, '2011_Reco14' : -0.00259,
                 '2011' : -0.00407301, '2012' : -0.00333,
                 'MC2011_Sim08a_incl_Jpsi' : -0.00076}
    mu['Value'] = mu_values.get(args[0], 0)
    # mu['Constant'] = options.simultaneous and not (options.split_mean or options.mu_param)
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sf_pee = options.simultaneous and options.sf_param
    tres_args = dict(time = time_obs, sigmat = st, Cache = True, Parameterise = options.parameterise,
                     TimeResSFParam = options.sf_param, SplitFracs = options.split_frac,
                     timeResMu = mu, Simultaneous = options.simultaneous,
                     ScaleFactors = [(2, 2.00), (1, 1.174)] if sf_pee else [(2, 0.1), (1, 0.06)],
                     Fractions = [(2, 0.143)], SplitMean = options.split_mean,
                     MeanParameterisation = options.mu_param)
    if args[0] == 'MC2012':
        tres_args['timeResFrac2'] = dict(Value = 0.15)
    sig_tres = TimeResolution(Name = 'sig_tres', **tres_args)
    if options.add_background:
        bkg_tres = TimeResolution(Name = 'bkg_tres', ParNamePrefix = 'bkg', **tres_args)
elif args[1] == 'triple':
    from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    sig_tres = TimeResolution(Name = 'tres', time = time_obs, sigmat = st, Cache = True,
                              SplitFracs = options.split_frac,
                              TimeResSFParam = options.sf_param, Parameterise = options.parameterise,
                              ScaleFactors = [(3, 1.5), (2, 4), (1, 1.4)], MeanParameterisation = options.mu_param,
                              Fractions = [(3, 0.1), (2, 0.2)], SplitMean = options.split_mean,
                              Simultaneous = options.simultaneous)

# J/psi mass pdf
from P2VV.Parameterizations.MassPDFs import DoubleCB_Psi_Mass as PsiMassPdf
psi_m = PsiMassPdf(mpsi, Name = 'psi_m', mpsi_alpha_1 = dict(Value = 2, Constant = (year == '2011')),
                   ParameteriseSigma = options.mass_param)
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
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.3672, MinMax = (1,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 1.3405e-01, MinMax = (0.01, 2.5))
                         )
psi_t = psi_t.pdf()

# J/psi signal component
psi_ll = Component('psi_ll', (psi_m, bkg_m, psi_t), Yield= (8.5575e+03,100,500000) )

# Background component
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
background = Component('background', (bkg_mpsi.pdf(), bkg_m), Yield = (19620,100,500000) )
if options.add_background:
    bkg_t = Prompt_Peak(time_obs, bkg_tres.model(), Name = 'bkg_prompt_pdf')
    background[t] = bkg_t.pdf()

# B signal component
sig_t = Prompt_Peak(time_obs, resolutionModel = sig_tres.model(), Name = 'sig_t')

from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_Mass
sig_m = Signal_Mass(Name = 'sig_m', mass = m)
signal = Component('signal', (sig_m.pdf(), sig_t.pdf()), Yield = (3e5, 10000, 1e7))

# Prompt component
prompt_pdf = Prompt_Peak(time_obs, sig_tres.model(), Name = 'prompt_pdf')
prompt = Component('prompt', (prompt_pdf.pdf(), psi_m), Yield = (160160, 100, 500000))

# Read data
fit_mass = (options.fit_mass or not options.cache) or options.reduce

## from profiler import heap_profiler_start, heap_profiler_stop
## heap_profiler_start("profile.log")

# Tree and cut
tree_name = 'DecayTree'
if signal_MC and options.excl_biased:
    hlt1_cut = 'hlt1_excl_biased == 1'
else:
    hlt1_cut = 'hlt1_unbiased == 1'
cut = 'sel == 1 && ' + hlt1_cut + ' && hlt2_unbiased == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
if options.cleantail or not options.wpv or (options.wpv and (options.wpv_type in ["Gauss", 'Rest'])):
    cut += ' && sel_cleantail == 1'
if signal_MC:
    cut += ' && abs(trueid) == 531'
if options.peak_only:
    cut += ' && %s > %5.3f && %s < %5.3f' % (time_obs.GetName(), time_obs.getMin(), time_obs.GetName(), time_obs.getMax())
if options.cut:
    cut = options.cut + ' && ' + cut
hash_str = cut
if options.use_refit:
    hash_str = cut + ' PVRefit'
hd = ('%d' % hash(hash_str)).replace('-', 'm')

if options.simultaneous:
    split_utils = {'sigmat'   : ('Sigmat', [st]),
                   'ppt'      : ('PPT', [momentum, pt]),
                   'momentum' : ('Momentum', [momentum]),
                   'pt'       : ('PT', [pt]),
                   'pv_zerr'  : ('PVZerr', [zerr]),
                   'nPV'      : ('NPV', [nPV])}
    from P2VV.Utilities import Resolution as ResolutionUtils
    assert(options.split in split_utils)
    split_opts = split_utils[options.split]
    SplitUtil = getattr(ResolutionUtils, 'Split' + split_opts[0])
    split_util = SplitUtil(args[0], *(split_opts[1]))
    directory = split_util.directory(hd)
    sub_dir = split_util.sub_dir(options.make_binning)
else:
    directory = '1bin_%4.2ffs_simple/%s' % (1000 * (t.getMax() - t.getMin()), hd)

if options.cache and not options.reduce:
    if options.simultaneous:
        from P2VV.CacheUtils import SimCache as Cache
        cache = Cache(input_data[args[0]]['cache'].rsplit('/', 1), directory, sub_dir)
    else:
        from P2VV.CacheUtils import Cache
        cache = Cache(input_data[args[0]]['cache'].rsplit('/', 1), directory)

results = {}
tree_name = 'DecayTree'

## Read Cache
if not fit_mass and options.cache:
    data, sdatas = cache.read_data()
    data = None
    try:
        if options.simultaneous:
            single_bin_sig_sdata = sdatas['sig_sdata']
            single_bin_bkg_sdata = sdatas['bkg_sdata']
            sig_sdata = sdatas[sub_dir + '/sig_sdata']
            bkg_sdata = sdatas[sub_dir + '/bkg_sdata']
        else:
            sig_sdata = single_bin_sig_sdata = sdatas['sig_sdata']
            bkg_sdata = single_bin_bkg_sdata = sdatas['bkg_sdata']
        results = cache.read_results()
    except KeyError:
        fit_mass = True
    results = cache.read_results()
    try:
        mass_result = results['mass_result']
        if options.simultaneous:
            sWeight_mass_result = results[sub_dir + '/sWeight_mass_result']
    except KeyError:
        sig_sdata = None
        bkg_sdata = None
        single_bin_sig_sdata = None
        single_bin_bkg_sdata = None
        sdatas = {}
        fit_mass = True
    if not fit_mass and options.simultaneous:
        split_cats = [split_util.split_cats(data = sig_sdata, mb = options.make_binning)]

## Fitting opts
fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True,
               Verbose = options.verbose, Strategy = 1)

# PV bins
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

from ROOT import gStyle
gStyle.SetPalette(53)
corr_canvas = None

## List of all plots we make
from collections import defaultdict
plots = defaultdict(list)

## Build simple mass pdf
if fit_mass:
    from P2VV.Utilities.DataHandling import readData
    data = readData(input_data[args[0]]['data'], tree_name, NTuple = True, observables = observables,
                    ntupleCuts = cut, ImportIntoWS = False)
    data.SetName(tree_name)
    if signal_MC and data.numEntries() > 3e5:
        data = data.reduce(EventRange = (int(2e5), int(5e5)))
    elif data.numEntries() > 6e5:
        data = data.reduce(EventRange = (0, int(6e5)))
    if options.reduce:
        data = data.reduce(EventRange = (0, int(options.reduce)))

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
        if mass_result.status() == 0 and abs(mass_result.minNll()) < 5e5:
            break
        
    assert(mass_result.status() == 0)
    mass_result.SetName('mass_result')
    results['mass_result'] = mass_result

    ## Canvas for correlation histograms
    from ROOT import TCanvas
    corr_canvas = TCanvas('mass_corr_canvas', 'mass_corr_canvas', 1000, 500)
    corr_canvas.Divide(2, 1)

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

    if options.plot:
        mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 600, 530)
        from P2VV.Utilities.Plotting import plot
        pdfOpts  = dict()
        if signal_MC:
            mass_obs = m
        else:
            mass_obs = mpsi
        
        ps = plot(mass_canvas.cd(), mass_obs, pdf = mass_pdf, data = data,
                  dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = 50),
                  pdfOpts  = dict(LineWidth = 2, **pdfOpts),
                  plotResidHist = 'BX',
                  xTitle = 'M_{#mu^{+}#mu^{-}} [MeV/c^{2}]',
                  yTitle = 'Candidates / (%3.2f MeV/c^{2})' % ((mass_obs.getMax() - mass_obs.getMin()) / 50),
                  yTitleOffset = 1 / 0.7)
        plots[''].append(ps)
        from P2VV.Utilities.Resolution import plot_dir
        if signal_MC:
            plot_name = 'prescaled_%s_B_mass.pdf'
        else:
            plot_name = 'prescaled_%s_Jpsi_mass.pdf'
        
        mass_canvas.Print(os.path.join(plot_dir, plot_name % args[0]), EmbedFonts = True)
        
    from P2VV.Utilities.SWeights import SData
    from P2VV.Utilities.DataHandling import correctWeights
    data_clone = data.Clone(data.GetName())
    sData = SData(Pdf = mass_pdf, Data = data_clone, Name = 'MassSPlot')
    ## single_bin_sig_sdata = correctWeights(sData.data(signal_name), splitCatNames = None, ImportIntoWS = False)
    ## single_bin_bkg_sdata = correctWeights(sData.data('background'), splitCatNames = None, ImportIntoWS = False)
    single_bin_sig_sdata = sData.data(signal_name)
    single_bin_bkg_sdata = sData.data('background')
    del sData
    sdatas = {'sig_sdata' : single_bin_sig_sdata,
              'bkg_sdata' : single_bin_bkg_sdata}
    ## Mass PDF is still connected to data_clone, redirect and delete data_clone
    mass_pdf.recursiveRedirectServers(data.get())
    del data_clone
    gc.collect()

if fit_mass and options.simultaneous:
    from P2VV.Utilities.General import getSplitPar
    # categories for splitting the PDF
    # get mass parameters that are split
    split_cats = [split_util.split_cats(data, options.make_binning)]
    split_pars = [[par for par in mass_pdf.Parameters() if par.getAttribute('Yield')]]

    # build simultaneous mass PDF
    sWeight_mass_pdf = SimultaneousPdf(mass_pdf.GetName() + '_simul',
                                       MasterPdf       = mass_pdf,
                                       SplitCategories = split_cats,
                                       SplitParameters = split_pars)

    for i in range(5):
        sWeight_mass_result = sWeight_mass_pdf.fitTo(data, **fitOpts)
        if sWeight_mass_result.status() == 0 and abs(sWeight_mass_result.minNll()) < 5e5:
            break

    assert(sWeight_mass_result.status() == 0)
    sWeight_mass_result.SetName('sWeight_mass_result')
    results[sub_dir + '/sWeight_mass_result'] = sWeight_mass_result

    ## Plot correlation histogram
    corr_canvas.cd(2)
    corr_hist_sWmass = sWeight_mass_result.correlationHist()
    corr_hist_sWmass.GetXaxis().SetLabelSize(0.03)
    corr_hist_sWmass.GetYaxis().SetLabelSize(0.03)
    corr_hist_sWmass.SetContour(20)
    corr_hist_sWmass.Draw('colz')

    from P2VV.Utilities.SWeights import SData
    sData = SData(Pdf = sWeight_mass_pdf, Data = data, Name = 'SimulMassSPlot')
    split_cat = sWeight_mass_pdf.indexCat()
    if split_cat.isFundamental():
        split_names = [split_cat.GetName()]
    else:
        split_names = [cat.GetName() for cat in split_cat.inputCatList()]
    ## sig_sdata_full = correctWeights(sData.data(signal_name), splitCatNames = None, ImportIntoWS = False)
    ## bkg_sdata_full = correctWeights(sData.data('background'), splitCatNames = None, ImportIntoWS = False)
    sig_sdata_full = sData.data(signal_name)
    bkg_sdata_full = sData.data('background')
    mass_pdf.recursiveRedirectServers(data.get())
    del sData
    if (signal_MC or prompt_MC) and options.reweigh:
        if cut.find('trueid') != -1:
            cut = cut.rsplit('&&', 1)[0]
        reco_data = readData(input_data[args[0]]['data'], tree_name,
                             NTuple = True, observables = observables[:-1],
                             ntupleCuts = cut, Rename = tree_name + '_' + year)
        # We build a mass PDF always wrt the J/psi, since there is not enough
        # signal to do it with. Let's hope this is OK...
        reco_mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='reco_mass_pdf')
        reco_mass_result = reco_mass_pdf.fitTo(reco_data, **fitOpts)
        reco_mass_result.SetName('reco_mass_result')
        reco_sData = SData(Pdf = reco_mass_pdf, Data = reco_data, Name = 'RecoMassSPlot')
        del reco_sData
        reco_sig_sdata = reco_sData.data(psi_ll.GetName())
        reco_bkg_sdata = reco_sData.data('background')
        from P2VV.Reweighing import reweigh
        for target, source, n in [(sig_sdata_full, reco_sig_sdata, 'full_sig_sdata'),
                                  (bkg_sdata_full, reco_bkg_sdata, 'full_bkg_sdata')]:
            ds, weights = reweigh(target, 'nPV', source, 'nPV', binning = PV_bounds)
            sdatas[n] = ds
    else:
        sdatas[sub_dir + '/sig_sdata'] = sig_sdata_full
        sdatas[sub_dir + '/bkg_sdata'] = bkg_sdata_full

    sig_sdata = sig_sdata_full
    bkg_sdata = bkg_sdata_full
elif fit_mass:
    if (signal_MC or prompt_MC) and options.reweigh:
        reco_data = readData(input_data[args[0]]['data'], tree_name,
                             NTuple = True, observables = observables[:-1],
                             ntupleCuts = cut, Rename = tree_name + '_' + year)
        reco_mass_pdf = buildPdf(Components = (psi_ll, background), Observables = (mpsi,), Name='reco_mass_pdf')
        reco_mass_result = reco_mass_pdf.fitTo(reco_data, **fitOpts)
        reco_mass_result.SetName('reco_mass_result')
        reco_sData = SData(Pdf = reco_mass_pdf, Data = reco_data, Name = 'RecoMassSPlot')
        reco_sig_sdata = reco_sData.data(psi_ll.GetName())
        reco_bkg_sdata = reco_sData.data('background')
        del reco_sData
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

# Write the result of the mass fit already at this point
if fit_mass and options.cache and not options.reduce:
    cache.write_cut(cut)
    cache.write_data(data, sdatas)
    cache.write_results(dict([(k, v) for k, v in results.iteritems() if 'mass' in k]))

gc.collect()

# Define default components
if signal_MC and options.wpv_type == "Rest":
    from P2VV.Parameterizations.TimeResolution import Rest_TimeResolution
    rest_tres = Rest_TimeResolution(Name = 'rest_tres', CoreModel = sig_tres, OwnMu = False)

    rest_t = Prompt_Peak(time_obs, resolutionModel = rest_tres.model(), Name = 'rest_t')

    rest = Component('rest', (rest_t.pdf(),), Yield = (5000, 1, 1e6))

    components = [signal, rest]
elif signal_MC:
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
    wpv_builder = WrongPV.ShapeBuilder(t, masses, UseKeysPdf = True, Weights = weights, Draw = options.plot,
                                       InputFile = input_data[args[0]]['wpv'], Workspace = input_data[args[0]]['workspace'],
                                       Reweigh = dict(Data = reweigh_data, DataVar = nPV, Binning = PV_bounds),
                                       MassResult = mass_result, **extra_args)
    if signal_MC:
        if time_obs == t:
            wpv_signal = wpv_builder.shape('B')
        else:
            wpv_signal = wpv_builder.diff_shape('B')
        sig_wpv = Component('sig_wpv', (wpv_signal, m), Yield = (888, 50, 300000))
    else:
        wpv_psi = wpv_builder.shape('jpsi')
        sig_wpv = Component('sig_wpv', (wpv_psi, psi_m), Yield = (1000, 50, 30000))
    components += [sig_wpv]
    if options.add_background:
        wpv_bkg = wpv_builder.shape('bkg')
        bkg_wpv = Component('wpv', (wpv_psi, bkg_m if signal_mc else bkg_mpsi), Yield = (1000, 50, 30000))
        componets += [bkg_wpv]
elif options.wpv and options.wpv_type == 'Gauss':
    def make_wpv_pdf(prefix):
        wpv_mean = RealVar('%swpv_mean' % prefix, Value = 0, MinMax = (-1, 1), Constant = not signal_MC)
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

## Define this one here so it doesn't get deleted by python. It owns a lot of
## stuff.
splitLeaves = RooArgSet()

# Which data to fit to
if options.add_background:
    fit_data = data
    fit_data_full = data
else:
    fit_data = sig_sdata
    if options.simultaneous:
        fit_data_full = sig_sdata

st_means = {'2011' : 0.350, '2012' : 0.349}
st_means[year] = sig_sdata.mean(st._target_())

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
    elif options.wpv and signal_MC and options.wpv_type == 'Rest':
        split_pars[0].append(rest.getYield())
        ## split_pars[0].append(rest_tres._left_rlifeSF)
        ## split_pars[0].append(rest_tres._frac_left)
    
    ## Calculate the mean in each bin for the splitting observable
    split_obs = split_util.observables()[0]
    split_obs_mean = st_means[year]
    split_cat = split_cats[0][0]
    if hasattr(split_cat, '_target_'):
        split_cat = split_cat._target_()
    means = [sig_sdata.mean(split_obs._target_(), '{0} == {0}::{1}'.format(split_cat.GetName(), s.GetName())) for s in split_cat]

    if options.sf_param or (options.mu_param and 'sigmat' not in options.mu_param):
        assert(len(split_cats[0]) == 1)

        ## The idea is to create a simultaneous PDF, split using the
        ## placeholder which was put in place when the resolution
        ## model was constructed and then in each bin set the created
        ## split parameter to: (mean_in_bin - overall_mean.

        ## Set the value to the mean sigmat already here, even if it is changed
        ## later.
        placeholder = sig_tres.sfPlaceHolder()
        placeholder.setVal(split_obs_mean)

        ## Customize the PDF to be split using the placeholder and any other
        ## pars which need to be split.
        from ROOT import RooCustomizer
        customizer = RooCustomizer(time_pdf._target_(), split_cat, splitLeaves)
        to_split = RooArgSet(*(split_pars[0]))
        ## CAREFULL: when do we need this?
        ## to_split.add(placeholder._target_())
        customizer.splitArgs(to_split, split_cat)
        states = dict([(s.GetName(), customizer.build(s.GetName())) for s in split_cat])
        time_pdf = SimultaneousPdf(time_pdf.GetName() + '_simul',
                                   SplitCategory = split_cat, ExternalConstraints = time_pdf.ExternalConstraints(),
                                   States = states, ConditionalObservables = time_pdf.ConditionalObservables())
        
        ## Set the split parameters to their calculated value and make
        ## them constant.
        st_mean = lambda st: 0.04921 + st * 1.220
        st_sigma = lambda st: 0.012 + st * 0.165

        pars = time_pdf.getParameters(RooArgSet())

        for m, s in zip(means, split_cat):
            for n, f in (('timeResSFMean', st_mean), ('timeResSFSigma', st_sigma)):
               p = pars.find(n + '_' + s.GetName())
               if not p: continue
               p.setVal(f(m - split_obs_mean))
            p = pars.find(placeholder.GetName() + '_' + s.GetName())
            if not p: continue
            p.setConstant(True)
            p.setVal(m - split_obs_mean)
    else:
        time_pdf = SimultaneousPdf(time_pdf.GetName() + '_simul'
                                   , MasterPdf       = time_pdf
                                   , SplitCategories = split_cats
                                   , SplitParameters = split_pars)
        st_mean = lambda st: 0.04921 + st * 1.220
        st_sigma = lambda st: 0.012 + st * 0.165

        pars = time_pdf.getParameters(RooArgSet())
        for m, s in zip(means, split_cat):
            for n, f in (('timeResSFMean', st_mean), ('timeResSFSigma', st_sigma)):
               p = pars.find(n + '_' + s.GetName())
               p.setVal(f(m - split_obs_mean))

if (not options.simultaneous and options.mu_param) or 'sigmat' in options.mu_param:
    placeholder = sig_tres.muPlaceHolder()
    placeholder.setVal(st_means[year])

if options.sf_param:
    placeholder = sig_tres.sfPlaceHolder()
    placeholder.setVal(st_means[year])

if options.plot:
    cats = RooArgSet(st)
    if options.simultaneous:
        cats.add(time_pdf.indexCat())
        st_data = fit_data_full.reduce(cats)
    else:
        st_data = fit_data.reduce(cats)
    
if options.reuse_result and options.cache:
    # Check if we have a cached time result, if so, use it as initial values for the fit
    if options.simultaneous:
        result_name = sub_dir + '/' + '_'.join(['time_result'] + extra_name)
    else:
        result_name = '_'.join(['time_result'] + extra_name)
    time_result = results.get(result_name, None)
    ## Try extra _pee suffix for backwards compatibility
    if not time_result:
        time_result = results.get(result_name + '_pee', None)
    
    if time_result:
        pdf_vars = time_pdf.getVariables()
        for p in time_result.floatParsFinal():
            pdf_par = pdf_vars.find(p.GetName())
            if pdf_par and not pdf_par.isConstant():
                pdf_par.setVal(p.getVal())
                pdf_par.setError(p.getError())

## Constraints
constraint_pars = set(['sf_mean_offset', 'sf_mean_slope', 'sf_sigma_offset', 'sf_sigma_slope',
                       'timeResFrac2', 'timeResMu_offset', 'timeResMu_quad', 'timeResMu_slope',
                       'timeResSFMean', 'timeResSFSigma', 'timeResRestFracLeft', 'timeResRestLTSF',
                       'timeResRestRSSF', 'timeResRestRTSF', 'timeResRestLSSF'])

import shelve
constraints = set()
if options.constrain:
    pdf_vars = time_pdf.getVariables()
    cp = args
    if options.use_refit:
        cp += ['refit']
    if signal_MC and options.excl_biased:
        cp += ['excl_biased']
    dbase = shelve.open('constraints.db')
    cid = ' '.join(args)
    assert(cid in dbase)

    cinfo = dbase[cid]
    from ROOT import TMatrixDSym
    from ROOT import TVectorD
    param_key = {'mu' : options.mu_param.split('_')[0], 'sf' : options.sf_param}
    parameters = cinfo[str(param_key)]
    if options.constrain != 'all':
        cpars = set()
        for exp in options.constrain.split(','):
            import re
            exp = re.compile(exp)
            cpars |= set(filter(exp.match, parameters.keys()))
    parameters = dict([(p, parameters[p]) for p in cpars])

    cv = RooArgList()
    for p in parameters.keys():
        pdf_par = pdf_vars.find(p)
        # Set values in pdf
        if not pdf_par:
            print 'Cannot find parameter %s in PDF' % p
            continue
        pdf_par.setVal(parameters[p][0])
        pdf_par.setError(parameters[p][1])
        pdf_par.setConstant(True)
        assert(pdf_par)
        cv.add(pdf_par)
    dbase.close()

    print 'constrained parameters:'
    for p in sorted(list(cpars)):
        print p, parameters[p]

time_pdf.Print("t")

## time_pdf.recursiveRedirectServers(sig_sdata.get())

## Fit
print 'fitting data'
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

## RooMsgService.instance().getStream(0).removeTopic(RooFit.Integration)
## nll = time_pdf.createNLL(sig_sdata, **fitOpts)
## RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Eval))

## nll.getVal()

## assert(False)

if options.fit:
    for i in range(4):
        time_result = time_pdf.fitTo(fit_data, SumW2Error = options.correct_errors, ExternalConstraints = constraints, **fitOpts)
        if time_result.status() == 0 and abs(time_result.minNll()) < 5e5:
            break
    time_result.SetName('_'.join(['time_result'] + extra_name))

    ## Draw correlation histogram for decay-time fit
    from ROOT import TCanvas
    time_corr_canvas = TCanvas('corr_canvas', 'corr_canvas', 500, 500)
    time_corr_hist_time = time_result.correlationHist()
    time_corr_hist_time.GetXaxis().SetLabelSize(0.03)
    time_corr_hist_time.GetYaxis().SetLabelSize(0.03)
    time_corr_hist_time.SetContour(20)
    time_corr_hist_time.Draw('colz')

    if options.simultaneous:
        results[sub_dir + '/' + time_result.GetName()] = time_result
    else:
        results[time_result.GetName()] = time_result

## profiler_stop()

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
bounds = {}
zoom_bounds = {}
if options.wpv and options.wpv_type == 'Mixing' and not signal_MC:
    bounds = array('d', [-5 + i * 0.1 for i in range(48)] + [-0.2 + i * 0.01 for i in range(40)] + 
                        [0.2 + i * 0.1 for i in range(58)])
    zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])
elif signal_MC:
    bounds = array('d', [-1.5 + i * 0.1 for i in range(10)] + [-0.5 + i * 0.02 for i in range(50)] + [0.5 + i * 0.1 for i in range(55)] + [6 + i * 0.4 for i in range(6)])
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

binnings = {'' : [binning, (None, None)], 'linear' : [zoom_binning, (None, None)]}
if signal_MC and not options.simultaneous:
    if options.use_refit:
        binnings[''][1] = (5, 2e5)
    elif options.wpv_type == 'Mixing':
        binnings[''][1] = (1, 2e5)
    else:
        binnings[''][1] = (1, 1e5)
plotLog = [True, False]
__canvases = []

from ROOT import SetOwnership
from ROOT import TCanvas
from ROOT import TPaveText
year_label = TPaveText(0.71, 0.72, 0.89, 0.85, "NDC")
year_label.SetFillColor(0)
year_label.AddText(args[0].split('_')[0][-4:])
year_label.SetBorderSize(0)

for (suffix, (bins, yRange)), pl in zip(binnings.items(), plotLog):
    if not options.plot or not time_result:
        continue
    if options.simultaneous:
        split_cat = time_pdf.indexCat()
        r = (bins.lowBound(), bins.highBound())
        for ct in split_cat:
            name = 'time_canvas_%s_%s' % (args[0], ct.GetName())
            if suffix:
                name += ('_' + suffix) 
            if options.use_refit:
                name += '_refit'
            canvas = TCanvas(name, name, 600, 400)
            __canvases.append(canvas)
            p = canvas.cd(1)
            
            projSet = RooArgSet(st, time_pdf.indexCat())
            pdfOpts  = dict(Slice = (split_cat, ct.GetName()), ProjWData = (st_data, True))
            if split_cat.isFundamental():
                cut = '{0} == {0}::{1}'.format(split_cat.GetName(), ct.GetName())
            else:
                split_cat.setLabel(ct.GetName())
                cut = ' && '.join(['{0} == {0}::{1}'.format(sc.GetName(), sc.getLabel()) for sc in split_cat.inputCatList()])
            ps = plot(p, time_obs, pdf = time_pdf, data = fit_data_full
                      , frameOpts = dict(Range = r, Title = "")
                      , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack, Cut = cut)
                      , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                      , xTitle = 't [ps]' if not signal_MC else 't - t_{true} [ps]'
                      , yTitle = 'Candidates / % 4.3f' % bins.averageBinWidth()
                      , yTitleOffset = 1 / 0.7
                      , logy = pl, yScale = yRange
                      , plotResidHist = 'BX'
                      ## , components = { 'wpv_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                      ##                  , 'prompt_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                      ##                  , 'sig_*'     : dict( LineColor = kOrange,  LineStyle = kDashed )
                      ##                  }
                      )
            for frame in ps:
                plot_name = '_'.join((t.GetName(), bins.GetName(), ct.GetName(), frame.GetName()))
                frame.SetName(plot_name)
            canvas.cd(1)
            year_label.Draw()
            canvas.Update()
            plots[sub_dir].append(ps)
    else:
        cname = 'time_canvas_' +args[0]
        if options.use_refit:
            cname += '_refit'
        if suffix:
            cname += ('_' + suffix)
        canvas = TCanvas(cname, cname, 600, 533)
        __canvases.append(canvas)
        p = canvas.cd(1)
        r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
        pdfOpts  = dict(ProjWData = (st_data, True))
        ps = plot(p, time_obs, pdf = time_pdf, data = fit_data
                  , frameOpts = dict(Range = r, Title = "")
                  , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                  , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                  , xTitle = 't [ps]' if not signal_MC else 't - t_{true} [ps]'
                  , yTitle = 'Candidates / (%4.3f ps)' % bins.averageBinWidth()
                  , yTitleOffset = 1 / 0.7
                  , logy = pl, yScale = yRange
                  , plotResidHist = 'BX')
        for frame in ps:
            plot_name = '_'.join((t.GetName(), bins.GetName(), frame.GetName()))
            frame.SetName(plot_name)
        canvas.cd(1)
        year_label.Draw()
        canvas.Update()
        plots[''].append(ps)

from P2VV.Utilities import Resolution as ResolutionUtils
if time_result:
    time_result.PrintSpecial(LaTeX = True, ParNames = ResolutionUtils.parNames)

# Write the result of the fit to the cache file
if options.cache:
    if options.fit and not options.reduce:
        cache.write_results(dict([(k, v) for k, v in results.iteritems() if 'time' in k]))

    if options.plot and not options.reduce:
        ## Write plots
        cache.write_plots(plots)

from P2VV.Utilities.Resolution import plot_dir
name_entries = {'model_ll' : args[1][0]+ 'g', 'wpv' : options.wpv_type.lower() + 'wpv'}
if not signal_MC:
    name_entries['model_ll'] += '_de'

plot_name = 'prescaled_%(model_ll)s_%(wpv)s_' % name_entries
if options.wpv_type == 'Mixing' and options.cleantail:
    plot_name += 'ct_'
def print_pdf():
    for canvas in __canvases:
        canvas.Print(os.path.join(plot_dir, plot_name + canvas.GetName()[12:] +'.pdf'), EmbedFonts = True)

print_pdf()
        
def write_constraints(constraints, mu_key = options.mu_param.replace('sigmat', ''),
                      sf_key = options.sf_param):
    ca = RooArgList()
    for p in sorted(list(constraint_pars)):
        rp = time_result.floatParsFinal().find(p)
        if rp:
            ca.add(rp)
    
    dbase = shelve.open('constraints.db')
    key_pars = args
    if options.use_refit:
        key_pars += ['refit']
    if signal_MC and options.excl_biased:
        key_pars += ['excl_biased']
    
    base_key = ' '.join(key_pars)
    param_key = {'mu' : mu_key, 'sf' : sf_key}
    
    cpars = set()
    for exp in constraints:
        import re
        exp = re.compile(exp)
        cpars |= set(filter(exp.match, parameters.keys()))
    print 'writing values for the following parameters:'
    for p in sorted(list(cpars)):
        print p
    
    dv = dict([(p.GetName(), (p.getVal(), p.getError())) for p in ca if p.GetName() in cpars])
    if base_key in dbase:
        d = dbase[base_key]
        if str(param_key) in d:
            pv = d[str(param_key)]
            pv.update(dv)
        else:
            pv = dv
        d.update({str(param_key) : pv})
        dbase[base_key] = d
    else:
        dbase[base_key] = {str(param_key) : dv}
    dbase.close()

if options.fit and options.write_constraints:
    write_constraints(options.write_constraints)

if False:
    from array import array
    tt_bins = array('d', [0.0, 0.09928, 0.2059, 0.3207, 0.447, 0.5844, 0.7365, 0.9077, 1.101, 1.325, 1.592, 1.916, 2.332, 2.916, 3.92, 14.0])
    
    from ROOT import RooBinning
    tt_binning = RooBinning(len(tt_bins) - 1, tt_bins)
    tt_binning.SetName('tt_binning')
    t_true.setBinning(tt_binning, 'tt_binning')
    
    from P2VV.RooFitWrappers import BinningCategory
    tt_cat = BinningCategory(t_true.GetName() + '_cat', Observable = t_true, Binning = tt_binning,
                             CatTypeName = 'tt_bin_', Data = sig_sdata, Fundamental = True)
    
    bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])
    binning = RooBinning(len(bounds) - 1, bounds)
    binning.SetName('binning')
    
    plot_datas = sorted(sig_sdata.split(tt_cat._target_()), key = lambda x: int(x.GetName().split('_')[-1]))
    
    __labels = []
    for i, (ct, ds) in enumerate(zip(tt_cat, plot_datas)):
        time_pdf.recursiveRedirectServers(ds.get())
        name = 'time_canvas_%s' % ct.GetName()
        canvas = TCanvas(name, name, 600, 533)
        __canvases.append(canvas)
        p = canvas.cd(1)
        projSet = RooArgSet(ds.get().find(st.GetName()))
        pdfOpts  = dict(ProjWData = (projSet, ds, True))
        ps = plot(p, time_obs, pdf = time_pdf, data = ds
                  , frameOpts = dict(Range = (bounds[0], bounds[-1]), Title = "")
                  , dataOpts = dict(MarkerSize = 0.8, Binning = binning, MarkerColor = kBlack)
                  , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                  , xTitle = 't - t_{true} [ps]'
                  , yTitle = 'Candidates / (%4.3f ps)' % bins.averageBinWidth()
                  , yTitleOffset = 1 / 0.7
                  , logy = False, yScale = (0, None)
                  , plotResidHist = 'BX')
        
        canvas.cd(1)
        tt_label = TPaveText(0.60, 0.75, 0.92, 0.83, "NDC")
        tt_label.SetFillColor(0)
        tt_label.SetBorderSize(0)
        __labels.append(tt_label)
        tt_label.AddText('%3.2f < t_{true} < %3.2f' % (tt_bins[i], tt_bins[i + 1]))
        tt_label.Draw()
        canvas.Update()
