import sys
import os
from P2VV.ToyMCUtils import Toy

toy = Toy()
parser = toy.parser()
(options, args) = toy.configure()

from P2VV.RooFitWrappers import *

generate = True

obj = RooObject( workspace = 'w')
w = obj.ws()

t_minmax = (-1.5, 8)
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))
t_true = RealVar('truetime', Title = 'true decay time', Unit='ps', Observable = True, MinMax=(-1100, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))

observables = [m, t, st, t_true]

cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
cut += ' && sel_cleantail == 1'
cut += ' && abs(trueid) == 531'

if not generate:
    prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
    filename = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130222.root')
    from P2VV.Utilities.DataHandling import readData
    data = readData(filename, "DecayTree", NTuple = True, observables = observables, ntupleCuts = cut)

    t_diff = FormulaVar(Name = 'time_diff', Formula = '@1 > -900 ? @0 - @1 : @0', Arguments = [t, t_true], data = data)
    t_diff.setMin(-10)
    t_diff.setMax(10)
else:
    t_diff = RealVar('time_diff', Unit = 'ps', Observable = True, MinMax = (-5, 5))

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
res_mean = RealVar("res_mean", Value = 0.5, MinMax = (-1, 1))
res_sigma = RealVar("res_sigma", Value = 0.1, MinMax = (0.0001, 1))
res_rlife = RealVar("res_rlife", Value = 0.1, MinMax = (0.001, 10))
res_sigma_sf = RealVar("res_sigma_sf", Value = 1., MinMax = (0.5, 10))
res_rlife_sf = RealVar("res_rlife_sf", Value = 1., MinMax = (0.5, 10))

mean_sf = ConstVar(Name = "gauss_mean_sf", Value = 1)

from ROOT import RooGExpModel
params = [t_diff, res_mean, st, res_rlife, mean_sf, res_sigma_sf, res_rlife_sf, 'false', 'Normal']
gexp_model = ResolutionModel(Name = "gexp_model", Type = RooGExpModel, Parameters = params,
                             ConditionalObservables = [st])

gauss_sigma_sf = RealVar(Name = "gauss_sigma_sf", Value = 1., MinMax = (0.5, 10))
from ROOT import RooGaussModel as GaussModel
gauss_model = ResolutionModel(Name = "gauss_model", Type = GaussModel, Parameters = [t_diff, res_mean, st, mean_sf, gauss_sigma_sf], ConditionalObservables = [st])

gexp_frac = RealVar(Name = "gexp_frac", Value = 0.1, MinMax = (0.0001, 0.9999))
add_model = AddModel("add_model", Models = [gexp_model, gauss_model], Fractions = [gexp_frac], ConditionalObservables = [st])

from ROOT import RooDecay as Decay
peak_tau = RealVar(Name = "peak_tau", Value = 0, Constant = True)
peak = Pdf(Name = "peak", Type = Decay, Parameters = [t_diff, peak_tau, add_model, 'SingleSided'], ConditionalObservables = [st])

## peak2 = Pdf(Name = "peak2", Type = Decay, Parameters = [t_diff, peak_tau, gexp_model, 'SingleSided'], ConditionalObservables = [st])

## peak3 = Pdf(Name = "peak3", Type = Decay, Parameters = [t_diff, peak_tau, gauss_model, 'SingleSided'], ConditionalObservables = [st])

# signal component
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_Mass
sig_m = Signal_Mass(Name = 'sig_m', mass = m) 
signal = Component('signal', (sig_m.pdf(), st_pdf, peak), Yield = (150000, 10000, 1000000))

sig_mass_pdf = buildPdf(Components = (signal, background), Observables = (m,), Name = 'sig_mass_pdf')
signal_name = signal.GetName()
mass_pdf = sig_mass_pdf

fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2, Offset = True,
               Verbose = False)

gen_pdf = buildPdf(Components = (signal,), Observables = (t_diff, st), Name = 'gen_pdf')
fit_pdf = buildPdf(Components = (signal,), Observables = (t_diff,), Name = 'fit_pdf')

toy.set_fit_opts(**dict(Verbose = False))
toy.run(Observables = [t_diff, st], Pdf = fit_pdf, GenPdf = gen_pdf)

