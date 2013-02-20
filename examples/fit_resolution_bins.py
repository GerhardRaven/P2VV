#!/usr/bin/env python
import optparse
import sys
import os

parser = optparse.OptionParser(usage = '%prog year model')

parser.add_option("--no-pee", dest = "pee", default = True,
                  action = 'store_false', help = 'Do not use per-event proper-time error')
parser.add_option('-p', '--parameterisation', dest = 'parameterise', default = False,
                  action = 'store', help = 'Parameterise sigmas [False, RMS, Comb]')
parser.add_option("--no-wpv", dest = "wpv", default = True,
                  action = 'store_false', help = 'Add WPV component')
parser.add_option("--fit-mass", dest = "fit_mass", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')
parser.add_option("--force-write", dest = "write_data", default = False,
                  action = 'store_true', help = 'Fit the mass spectrum even if data is available.')
parser.add_option("--split-var", dest = "split_var", default = 'st', type = 'string',
                  action = 'store', help = 'Split data using this variable.')

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
    prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
    input_data['data'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_20130207.root')
    input_data['wpv'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011.root')
    input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2011_workspace'
    input_data['weighted'] = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled_st_bins.root')
    ## input_data['data'] = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_prescaled.root'
    ## input_data['wpv'] = '/stuff/PhD/mixing/Bs2JpsiPhiPrescaled_2011.root'
    ## input_data['workspace'] = 'Bs2JpsiPhiPrescaled_2011_workspace'
else:
    input_data['data'] = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20121218.root'
    input_data['weighted'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_2012_Prescaled_st_bins.root'
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
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))
zerr = RealVar('B_s0_bpv_zerr', Title = 'Best PV Z error', Unit = 'mm', Observable = True, MinMax = (0, 0.1))

assert(options.split_var in ['st', 'zerr'])

from array import array
if options.split_var == 'st':
    ## split_bins = array('d', [0.01 + i * 0.005 for i in range(13)])
    split_bins = array('d', [0.01 + i * 0.01 for i in range(5)] + [0.07])
    ## split_bins = array('d', [0.01 + i * 0.01 for i in range(7)])
    split_var = st
else:
    split_bins = array('d', [0, 0.021, 0.025, 0.03, 0.04, 0.06, 0.1])
    split_var = zerr


# add 20 bins for caching the normalization integral
for i in [ split_var ] : i.setBins( 20 , 'cache' )

# Categories needed for selecting events
unbiased = Category('triggerDecisionUnbiasedPrescaled', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
nPV = RealVar('nPV', Title = 'Number of PVs', Observable = True, MinMax = (0, 10))

observables = [t, m, mpsi, st, unbiased, nPV, zerr]

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
                              timeResMu = dict(Value = -0.17, MinMax = (-1, 1)),
                              sigmaSF  = dict(Value = 1.46, MinMax = (0.1, 5)))
elif args[1] == 'double':
    from P2VVParameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
    scaleFactors = [(2, 2.3), (1, 1.2)]
    if not options.pee:
        scaleFactors = [(n, 0.032 * v) for n, v in scaleFactors]
    sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, Cache = True,
                              PerEventError = options.pee, Parameterise = options.parameterise,
                              ScaleFactors = scaleFactors,
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
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.42, MinMax = (0.5,  2.5), Constant = True)
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 0.145, MinMax = (0.01, 0.5), Constant = True)
                         )

## from P2VVParameterizations.TimePDFs import Single_Exponent_Time as Background_Time
## psi_t = Background_Time(Name = 'psi_t', time = t, resolutionModel = sig_tres.model(),
##                              t_sig_tau  = dict(Name = 'tau', Value = 1.5, MinMax = (0.001, 2.5))
##                              )
psi_t = psi_t.pdf()


bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = sig_tres.model()
                         , bkg_t_fml    = dict(Name = 'bkg_t_fml',    Value = 0.76 )
                         , bkg_t_ll_tau = dict(Name = 'bkg_t_ll_tau', Value = 1., MinMax = (0.01, 2.5))
                         , bkg_t_ml_tau = dict(Name = 'bkg_t_ml_tau', Value = 0.1,  MinMax = (0.01, 0.5))
                         )
bkg_t = bkg_t.pdf()

signal = Component('signal', (sig_m, psi_m.pdf(), sig_t), Yield = (200000, 500, 500000))
psi_background = Component('psi_background', (psi_m.pdf(), bkg_m.pdf(), psi_t), Yield= (4000,1,500000) )

background = Component('background', (bkg_mpsi.pdf(), bkg_m.pdf(), bkg_t), Yield = (19620,1,500000) )

# Prompt component
from P2VVParameterizations.TimePDFs import Prompt_Peak
prompt_pdf = Prompt_Peak(t, sig_tres.model(), Name = 'prompt_pdf')
psi_prompt = Component('prompt', (prompt_pdf.pdf(), ), Yield = (77000, 1, 500000))

# Read data
datas = []
sdatas = []
from ROOT import TFile

fit_mass = options.fit_mass

cut = 'nPV > 3 && sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
## cut = 'B_s0_bpv_zerr > 0 && B_s0_bpv_zerr < 0.034 && sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
if not options.wpv:
    cut += ' && sel_cleantail == 1'
hd = ('%d' % hash(cut)).replace('-', 'm')

if options.split_var == 'st':
    directory = '%sbins_%4.2ffs/%s' % (len(split_bins) - 1, (1000 * (split_bins[1] - split_bins[0])), hd)
elif options.split_var == 'zerr':
    directory = 'zerr_%sbins/%s' % (len(split_bins) - 1, hd)
if os.path.exists(input_data['weighted']):
    cache_file = TFile.Open(input_data['weighted'], 'update')
else:
    cache_file = TFile.Open(input_data['weighted'], 'new')
cache_dir = cache_file.Get(directory)
if not cache_dir:
    cache_file.mkdir(directory)
    cache_dir = cache_file.Get(directory)
    from ROOT import TObjString
    cut_string = TObjString(cut)
    cache_dir.WriteTObject(cut_string, 'cut')
    fit_mass = True

## Fitting options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True)

## save all fit results
from collections import defaultdict
results = defaultdict(list)

## Common imports
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas

if not fit_mass:
    rs = set()
    for i in range(len(split_bins) - 1):
        sds_name = 'sdata/DecayTree_%02d_weighted_psi_background' % i
        sdata = cache_dir.Get(sds_name)
        if sdata:
            sdatas.append(sdata)
        else:
            fit_mass = True
            break
        ds_name = 'data/DecayTree_%02d' % i
        data = cache_dir.Get(ds_name)
        if data:
            datas.append(data)
        else:
            fit_mass = True
            break
    rd = cache_dir.Get('mass_results')
    if not rd:
        fit_mass = True
    else:
        nkeys = rd.ReadKeys()
        if nkeys != (len(split_bins)  - 1):
            fit_mass = True
        else:
            for e in rd.GetListOfKeys():
                if e.GetClassName() == 'RooFitResult':
                    rs.add(os.path.join(rd.GetName(), e.GetName()))
            for e in rs:
                results['mass'].append(cache_dir.Get(e))
            results['mass'] = sorted(results['mass'], key = lambda r: int(r.GetName().split('_', 1)[0]))

tree_name = 'DecayTree'
if fit_mass:
    from P2VVGeneralUtils import readData
    data = readData(input_data['data'], tree_name, NTuple = True, observables = observables,
                    ntupleCuts = cut)
    datas = [data.reduce(Cut = '{0} > {1} && {0} < {2}'.format(split_var.GetName(), split_bins[i], split_bins[i + 1]))
             for i in range(len(split_bins) - 1)]
    for i, ds in enumerate(datas):
        ds_name = ds.GetName().replace('DecayTree', 'DecayTree_%02d' % i)
        ds.SetName(ds_name)

    mass_pdf = buildPdf(Components = (psi_background, background), Observables = (mpsi,), Name='mass_pdf')
    mass_pdf.Print('t')

    mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 1200, 900)
    pads = mass_canvas.pads(4, 3)

    from P2VVGeneralUtils import SData
    sdatas = []
    for i, (p, ds) in enumerate(zip(pads, datas)):
        result = mass_pdf.fitTo(ds, **fitOpts)
        result.SetName('%d_%s' % (i, result.GetName()))
        results['mass'].append(result)
        splot = SData(Pdf = mass_pdf, Data = ds, Name = 'MassSplot')
        sdatas.append(splot.data('psi_background'))
        from P2VVGeneralUtils import plot
        plot(p, mpsi, pdf = mass_pdf, data = ds
             , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
             , pdfOpts  = dict(LineWidth = 2)
             , plotResidHist = False
             , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                              , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                              }
             )

# Wrong PV components
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])
if options.wpv:
    mass_pdf.fitTo(data, **fitOpts)
    from P2VVGeneralUtils import SData
    for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
    splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
    psi_sdata = splot.data('psi_background')
    bkg_sdata = splot.data('background')

    from P2VVParameterizations import WrongPV
    reweigh_data = dict(jpsi = psi_sdata, bkg = bkg_sdata)
    wpv = WrongPV.ShapeBuilder(t, {'jpsi' : mpsi}, UseKeysPdf = True, Weights = 'jpsi', Draw = True,
                               InputFile = input_data['wpv'], Workspace = input_data['workspace'],
                               Reweigh = dict(Data = reweigh_data, DataVar = nPV, Binning = PV_bounds),
                               sigmat = st)
    wpv_psi = wpv.shape('jpsi')
    psi_wpv = Component('psi_wpv', (wpv_psi,), Yield = (1000, 1, 30000))
else:
    wpv_mean = sig_tres._timeResMu
    wpv_sigma = RealVar('wpv_sigma', Value = 0.3, MinMax = (0.01, 20), Constant = True)
    wpv_pdf = Pdf(Name = 'wpv_pdf', Type = Gaussian, Parameters = (t, wpv_mean, wpv_sigma))
    psi_wpv = Component('wpv', (wpv_pdf, ), Yield = (100, 5, 500000))
components = [psi_prompt, psi_background, psi_wpv]

time_pdf = buildPdf(Components = components, Observables = (t,), Name='time_pdf')
time_pdf.Print("t")

from ROOT import RooBinning
from array import array

if options.wpv:
    bounds = array('d', [-5 + i * 0.1 for i in range(47)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(21)])
else:
    bounds = array('d', [-1.5 + i * 0.1 for i in range(12)] + [-0.3 + i * 0.05 for i in range(12)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])

binning = RooBinning(len(bounds) - 1, bounds)
binning.SetName('var_binning')
t.setBinning(binning)

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from ROOT import TCanvas
import P2VVGeneralUtils

time_canvas = TCanvas('time_canvas', 'time_canvas', 1200, 1050)
pads = time_canvas.pads(4, 3)
for i, (p, ds) in enumerate(zip(pads, sdatas)):
    result = time_pdf.fitTo(ds, SumW2Error = False, **fitOpts)
    j = 0
    while result.status() != 0 and j < 4:
        result = time_pdf.fitTo(ds, SumW2Error = False, **fitOpts)
        j += 1
    result.SetName('%d_%s' % (i, result.GetName()))
    results['time'].append(result)
    pdfOpts  = dict(ProjWData = (RooArgSet(st), ds, True))
    P2VVGeneralUtils.plot(p, t, pdf = time_pdf, data = ds
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

for k, res in results.items():
    results[k] = sorted(res, key = lambda r: int(r.GetName().split('_', 1)[0]))

res_canvas = TCanvas('res_canvas', 'res_canvas', 500, 500)
from ROOT import TH1D
from math import sqrt

if options.split_var == 'st':
    split_bounds = array('d', [1000 * v for v in split_bins])
else:
    split_bounds = array('d', [v for v in split_bins])
    
hist_res = TH1D('hist_res', 'hist_res', len(split_bounds) - 1, split_bounds)
hist_events = TH1D('hist_events', 'hist_events', len(split_bounds) - 1, split_bounds)

resolutions = []

from ROOT import TMatrixT

for index, result in enumerate(results['time']):
    # Fill events histo
    events = results['mass'][index].floatParsFinal().find(psi_background.getYield().GetName())
    hist_events.SetBinContent(index + 1, events.getVal())
    hist_events.SetBinError(index + 1, events.getError())
    
    if result.status() != 0:
        resolutions.append((0, 0))
        hist_res.SetBinContent(index + 1, 0)
        hist_res.SetBinError(index + 1, 0)
        continue
    
    fpf = result.floatParsFinal()
    if args[1] == 'double' and options.parameterise == False:
        from PropagateErrors import propagateScaleFactor
        sf, sf_e = propagateScaleFactor(result)
    if args[1] == 'double' and options.parameterise == 'Comb':
        sf_comb = fpf.find('timeResComb')
        sf, sf_e = sf_comb.getVal(), sf_comb.getError()
    elif args[1] == 'single':
        sf = fpf.find('sigmaSF').getVal()
        sf_e = fpf.find('sigmaSF').getError()
    
    mean = sdatas[index].mean(st._target_())
    res = mean * sf
    res_e = mean * sf_e
    
    resolutions.append((res, res_e))
    hist_res.SetBinContent(index + 1, 1000 * res if options.split_var == 'st' else res)
    hist_res.SetBinError(index + 1, 1000 * res_e if options.split_var == 'st' else res_e)

scale = 100 / hist_events.GetMaximum()
hist_events.Scale(scale)
hist_events.GetYaxis().SetRangeUser(0, 110)
hist_res.GetYaxis().SetRangeUser(0, 110)
hist_res.Draw('pe')
hist_events.Draw('hist, same')
from ROOT import kGray
hist_events.SetFillColor(kGray + 1)
hist_res.Draw('pe, same')
if options.split_var == 'st':
    hist_res.GetXaxis().SetTitle('estimated decay time error [fs]')
else:
    hist_res.GetXaxis().SetTitle('#sigma_{PV,Z} [mm]')

hist_res.GetYaxis().SetTitle('decay time resulution [fs]')

from ROOT import TF1
fit_func = TF1('fit_func', "pol1", split_bins[0], split_bins[-1])
fit_result = hist_res.Fit(fit_func, "S0")

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
    for (dss, d) in [(sdatas, sdata_dir), (datas, data_dir)]:
        for i, ds in enumerate(dss):
            d.WriteTObject(ds, ds.GetName(), "Overwrite")
        d.Write(d.GetName(), TObject.kOverwrite)

# Always write fit results
results_dirs = dict((k, get_dir('%s_results' % k)) for k in results.iterkeys())
for k, rs in results.iteritems():
    rd = results_dirs[k]
    for r in rs:
        rd.WriteTObject(r, r.GetName(), "Overwrite")
    rd.Write(rd.GetName(), TObject.kOverwrite)
# Delete the input TTree which was automatically attached.
cache_file.Delete('%s;*' % tree_name)

## import Dilution
## Dilution.dilution(t, data, result = result, sigmat = st, signal = [psi_prompt], subtract = [psi_background])
