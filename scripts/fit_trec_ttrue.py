#!/usr/bin/env python
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog data')

parser.add_option('-s', "--simultaneous", dest = "simultaneous", default = False,
                  action = 'store_true', help = 'Use sigmat offset')
parser.add_option("-b", "--batch", dest = "batch", default = False,
                  action = 'store_true', help = 'run ROOT in batch mode')

(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.usage
    sys.exit(-2)
elif args[0] not in ['generate', 'MC11a']:
    print parser.usage
    sys.exit(-2)

if options.batch:
    from ROOT import gROOT
    gROOT.SetBatch(True)

from P2VV.RooFitWrappers import *


obj = RooObject( workspace = 'w')
w = obj.ws()

t_minmax = (-1.5, 8)
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
t_true = RealVar('truetime', Title = 'true decay time', Unit='ps', Observable = True, MinMax=(-1100, 14))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))

observables = [m, t, t_true]

cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
cut += ' && sel_cleantail == 1'
cut += ' && abs(trueid) == 531'

sig_sdata = None
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'

if args[0] != 'generate':
    ds_filename = os.path.join(prefix, 'p2vv/data/MC11_signal_data.root')
    ## ds_filename = '/tmp/test.root'
    from ROOT import TFile
    if os.path.exists(ds_filename):
        t_diff = RealVar('time_diff', Unit = 'ps', Observable = True, MinMax = (-1, 1))
        f = TFile.Open(ds_filename)
        sig_sdata = f.Get("sig_sdata")
        sig_sdata = sig_sdata.reduce(Cut = 'time_diff > %f && time_diff < %f' % t_diff.getRange(), EventRange = (0, 50000))
    else:
        observables.append(st)
        filename = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130222.root')
        from P2VV.GeneralUtils import readData
        tree_name = 'DecayTree'
        data = readData(filename, tree_name, NTuple = True, observables = observables, ntupleCuts = cut)

        t_diff = FormulaVar('time_diff', '@1 > -900 ? @0 - @1 : @0', [t, t_true], data = data)
        t_diff.setMin(-1)
        t_diff.setMax(1)
else:
    t_diff = RealVar('time_diff', Unit = 'ps', Observable = True, MinMax = (-1, 1))

observables.append(t_diff)
    
# Create combinatorical background component
mean = RealVar(Name = 'mean', Value = 5368, MinMax = (5300, 5400))
sigma = RealVar(Name = 'sigma', Value = 50, MinMax = (1, 500))
from ROOT import RooGaussian as Gaussian
bkg_m = Pdf(Name = 'gauss', Type = Gaussian, Parameters = (m, mean, sigma))
background = Component('background', (bkg_m,), Yield = (19620,100,500000) )

# Make a double Lognormal PDF to generate sigma_t
median  = RealVar('median',   Unit = 'ps', Value = 0.0327, MinMax = (0.0001, 0.12))
k1 = RealVar('k1',  Unit = '', Value = 1.2169, MinMax = (0.00001, 10))
k2 = RealVar('k2',  Unit = '', Value = 1.356, MinMax = (0.00001, 10))
frac = RealVar('frac_ln2', Value = 0.4728, MinMax = (0.01, 0.99))

ln1 = LognormalPdf('ln1', Observable = st, Median = median, Shape = k1)
ln2 = LognormalPdf('ln2', Observable = st, Median = median, Shape = k2)

# Do our own sum pdf to have a fraction
st_pdf = SumPdf(Name = 'ln', PDFs = [ln1, ln2], Yields = {'ln2' : frac})

# Signal t_diff PDF
res_mean = RealVar("res_mean", Value = 0, MinMax = (-0.5, 0.5))
res_rlife = RealVar("res_rlife", Value = 0.1, MinMax = (0.001, 10))
res_sigma_sf = RealVar("res_sigma_sf", Value = 1., MinMax = (0.001, 10))
res_rlife_sf = RealVar("res_rlife_sf", Value = 1., MinMax = (0.001, 10), Constant = True)
mean_sf = ConstVar(Name = "gauss_mean_sf", Value = 1)

from ROOT import RooGExpModel
params = [t_diff, res_mean, st, res_rlife, mean_sf, res_sigma_sf, res_rlife_sf, 'false', 'Normal']
gexp_model = ResolutionModel(Name = "gexp_model_one", Type = RooGExpModel, Parameters = params,
                             ConditionalObservables = [st])

gauss_sigma_sf_one = RealVar(Name = "gauss_sigma_sf_one", Value = 1., MinMax = (0.5, 100))
from ROOT import RooGaussModel as GaussModel
gauss_model_one = ResolutionModel(Name = "gauss_model_one", Type = GaussModel, ConditionalObservables = [st],
                              Parameters = [t_diff, res_mean, st, mean_sf, gauss_sigma_sf_one])

gexp_frac = RealVar(Name = "gexp_frac", Value = 0.2, MinMax = (0.0001, 0.9999))

add_model = AddModel("add_model", Models = [gexp_model, gauss_model_one],
                     Fractions = [gexp_frac],
                     ConditionalObservables = [st])

from ROOT import RooDecay as Decay
peak_tau = RealVar(Name = "peak_tau", Value = 0, Constant = True)
peak = Pdf(Name = "peak", Type = Decay, Parameters = [t_diff, peak_tau, add_model, 'SingleSided'],
           ConditionalObservables = [st])

# signal component
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_Mass
sig_m = Signal_Mass(Name = 'sig_m', mass = m) 
signal = Component('signal', (sig_m.pdf(), st_pdf, peak), Yield = (150000, 10000, 1000000))

sig_mass_pdf = buildPdf(Components = (signal, background), Observables = (m,), Name = 'sig_mass_pdf')
signal_name = signal.GetName()
mass_pdf = sig_mass_pdf

fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 0, Verbose = False,
               Offset = True)

fit_pdf = buildPdf(Components = (signal,), Observables = (t_diff,), Name = 'fit_pdf')

if args[0] is not "generate" and not sig_sdata:
    ## Fit mass pdf
    for i in range(3):
        mass_result = mass_pdf.fitTo(data, **fitOpts)
        if mass_result.status() == 0:
            break

    assert(mass_result.status() == 0)
    mass_result.SetName('mass_result')

    from P2VV.GeneralUtils import SData
    data_clone = data.Clone(data.GetName())
    sData = SData(Pdf = mass_pdf, Data = data_clone, Name = 'MassSPlot')
    sig_sdata = sData.data(signal_name)
    bkg_sdata = sData.data('background')
elif args[0] == 'generate':
    gen_pdf = buildPdf(Components = (signal,), Observables = (t_diff, st), Name = 'gen_pdf')
    sig_sdata = gen_pdf.generate([t_diff, st], 10000)

if options.simultaneous:
    split_pars = [[par for par in fit_pdf.Parameters() if par.getAttribute('Yield')]]
    split_pars[0] += [res_sigma_sf, gauss_sigma_sf_one]

    from array import array
    st_bins = array('d', [0.01000, 0.02066, 0.02375, 0.02616, 0.02833,
                          0.03047, 0.03269, 0.03520, 0.03837, 0.04343, 0.07000])
    from ROOT import RooBinning
    st_binning = RooBinning( len(st_bins) - 1, st_bins, 'st_binning' )
    st.setBinning(st_binning, 'st_binning')
    st_cat = BinningCategory(st.GetName() + '_cat', Observable = st, Binning = st_binning,
                             Fundamental = True, Data = sig_sdata, CatTypeName = 'bin')
    original_pdf = fit_pdf
    
    fit_pdf = SimultaneousPdf(fit_pdf.GetName() + '_simul'
                               , MasterPdf       = fit_pdf
                               , SplitCategories = [[st_cat]]
                               , SplitParameters = split_pars)

result = fit_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)

diff_frame = t_diff.frame(Range = (-0.5, 0.5))
sig_sdata.plotOn(diff_frame)
proj_set = RooArgSet(st)
if options.simultaneous:
    proj_set.add(st_cat._target_())
fit_pdf.plotOn(diff_frame, ProjWData = (proj_set, sig_sdata, True))

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 500, 500)
diff_frame.Draw()
plots = [diff_frame]

hd = ('%d' % hash(cut)).replace('-', 'm')
directory = '%sbins_simul/%s' % (len(st_bins) - 1, hd)

from P2VV.CacheUtils import CacheFiles
cache_filename = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_Prescaled.root')
cache_files = CacheFiles(*cache_filename.rsplit('/', 1))
cache_dir, cache_file = cache_files.getFromCache(directory)

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
    sdata_dir = get_dir('sdata')
    data_dir = get_dir('data')
    sdata_dir.WriteTObject(sig_sdata, sig_sdata.GetName(), "Overwrite")
        
    sdata_dir.Write(sdata_dir.GetName(), TObject.kOverwrite)
    
    ## Write PDFs
    pdf_dir = get_dir('PDFs')
    pdf_dir.WriteTObject(fit_pdf._target_(), 'fit_pdf', "Overwrite")
    
    pdf_dir.WriteTObject(mass_pdf._target_(), 'mass_pdf', "Overwrite")
    
    ## Write fit results
    results_dir = get_dir('results')
    results_dir.WriteTObject(result, result.GetName(), "Overwrite")
    results_dir.Write(results_dir.GetName(), TObject.kOverwrite)
    
    ## Write plots
    plots_dir = get_dir('plots/MC11a')
    for p in plots:
        plots_dir.WriteTObject(p, p.GetName(), "Overwrite")
    
    plots_dir.Write(plots_dir.GetName(), TObject.kOverwrite)
    
    # Delete the input TTree which was automatically attached.
    cache_file.Delete('%s;*' % tree_name)
