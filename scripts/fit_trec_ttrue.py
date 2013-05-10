import os
from P2VV.RooFitWrappers import *

generate = False

obj = RooObject( workspace = 'w')
w = obj.ws()

t_minmax = (-1.5, 8)
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
t_true = RealVar('truetime', Title = 'true decay time', Unit='ps', Observable = True, MinMax=(-1100, 14))
t_diff = RealVar('time_diff', Unit = 'ps', Observable = True, MinMax = (-1, 1))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))

observables = [m, t, t_true]

cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
cut += ' && sel_cleantail == 1'
cut += ' && abs(trueid) == 531'

sig_sdata = None
if not generate:
    prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
    ds_filename = os.path.join(prefix, 'p2vv/data/MC11_signal_data.root')
    ## ds_filename = '/tmp/test.root'
    from ROOT import TFile
    f = TFile.Open(ds_filename)
    if f.IsOpen() and f.GetKey("sig_sdata"):
        sig_sdata = f.Get("sig_sdata")
        sig_sdata = sig_sdata.reduce(Cut = 'time_diff > %f && time_diff < %f' % t_diff.getRange(), EventRange = (0, 50000))
    else:
        observables.append(st)
        filename = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130222.root')
        from P2VV.GeneralUtils import readData
        data = readData(filename, "DecayTree", NTuple = True, observables = observables, ntupleCuts = cut)

        t_diff = FormulaVar('time_diff', '@1 > -900 ? @0 - @1 : @0', [t, t_true], data = data)
        t_diff.setMin(-5)
        t_diff.setMax(5)

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
res_sigma = RealVar("res_sigma", Value = 0.1, MinMax = (0.001, 2))
res_sigma_sf = RealVar("res_sigma_sf", Value = 1., MinMax = (0.001, 10))
res_rlife_sf = RealVar("res_rlife_sf", Value = 1., MinMax = (0.001, 10), Constant = True)
mean_sf = ConstVar(Name = "gauss_mean_sf", Value = 1)

from ROOT import RooGExpModel
params = [t_diff, res_mean, st, res_rlife, mean_sf, res_sigma_sf, res_rlife_sf, 'false', 'Normal']
gexp_model = ResolutionModel(Name = "gexp_model", Type = RooGExpModel, Parameters = params,
                             ConditionalObservables = [st])

gauss_sigma_sf_one = RealVar(Name = "gauss_sigma_sf_one", Value = 1., MinMax = (0.5, 10))
from ROOT import RooGaussModel as GaussModel
gauss_model_one = ResolutionModel(Name = "gauss_model_one", Type = GaussModel, ConditionalObservables = [st],
                              Parameters = [t_diff, res_mean, st, mean_sf, gauss_sigma_sf_one])

gauss_sigma_sf_two = RealVar(Name = "gauss_sigma_sf_two", Value = 1., MinMax = (0.5, 10))
from ROOT import RooGaussModel as GaussModel
gauss_model_two = ResolutionModel(Name = "gauss_model_two", Type = GaussModel, ConditionalObservables = [st],
                              Parameters = [t_diff, res_mean, st, mean_sf, gauss_sigma_sf_two])

gauss_two_frac = RealVar(Name = "gauss_two_frac", Value = 0.2, MinMax = (0.0001, 0.9999))
gexp_frac = RealVar(Name = "gexp_frac", Value = 0.2, MinMax = (0.0001, 0.9999))
add_model = AddModel("add_model", Models = [gexp_model, gauss_model_two, gauss_model_one], Fractions = [gexp_frac, gauss_two_frac],
                     ConditionalObservables = [st])

from ROOT import RooDecay as Decay
peak_tau = RealVar(Name = "peak_tau", Value = 0, Constant = True)
peak = Pdf(Name = "peak", Type = Decay, Parameters = [t_diff, peak_tau, add_model, 'SingleSided'],
           ConditionalObservables = [st])

## peak2 = Pdf(Name = "peak2", Type = Decay, Parameters = [t_diff, peak_tau, gexp_model, 'SingleSided'], ConditionalObservables = [st])

## peak3 = Pdf(Name = "peak3", Type = Decay, Parameters = [t_diff, peak_tau, gauss_model, 'SingleSided'], ConditionalObservables = [st])

# signal component
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_Mass
sig_m = Signal_Mass(Name = 'sig_m', mass = m) 
signal = Component('signal', (sig_m.pdf(), st_pdf, peak), Yield = (150000, 10000, 1000000))

sig_mass_pdf = buildPdf(Components = (signal, background), Observables = (m,), Name = 'sig_mass_pdf')
signal_name = signal.GetName()
mass_pdf = sig_mass_pdf

fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 0, Verbose = False,
               Offset = True)

from ROOT import RooMsgService
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Integration))
## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

fit_pdf = buildPdf(Components = (signal,), Observables = (t_diff,), Name = 'fit_pdf')

if not generate and not sig_sdata:
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
elif generate:
    gen_pdf = buildPdf(Components = (signal,), Observables = (t_diff, st), Name = 'gen_pdf')
    sig_sdata = gen_pdf.generate([t_diff, st], 10000)

result = fit_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)

diff_frame = t_diff.frame(Range = (-0.5, 0.5))
sig_sdata.plotOn(diff_frame)
fit_pdf.plotOn(diff_frame, ProjWData = (RooArgSet(st), sig_sdata, True))

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 500, 500)
diff_frame.Draw()
