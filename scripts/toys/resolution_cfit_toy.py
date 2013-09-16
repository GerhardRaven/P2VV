import sys
import os
try:
    from P2VV.ToyMCUtils import Toy
except ImportError:
    print 'Could not import P2VV version, trying local directory'
    from ToyMCUtils import Toy

toy = Toy()
parser = toy.parser()

parser.add_option("--sfit", dest = "sFit", default = False,
                  action = 'store_true', help = "Generate with cFit, fit with sFit")

(options, args) = toy.configure()

from ROOT import gROOT
#gROOT.SetBatch(True)

from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle
from ROOT import RooMsgService

obj = RooObject( workspace = 'w')
w = obj.ws()

t_minmax = (-5, 14)
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3020, 3165))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))

# add 20 bins for caching the normalization integral
st.setBins(500, 'cache')

# Categories needed for selecting events
unbiased = Category('triggerDecisionUnbiasedPrescaled', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'selected' : 1, 'not_selected' : 0})

observables = [t, mpsi, st, unbiased, selected]
time_obs = t

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# Time resolution model
sig_tres = None
mu = dict(MinMax = (-0.010, 0.010))
mu['Constant'] = False
mu['Value'] = -0.00407301

from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution

## Signal and background time resolutions.
tres_args = dict(time = time_obs, sigmat = st, Cache = True, PerEventError = True,
                 Parameterise = 'Comb', TimeResSFParam = '', timeResMu = mu,
                 ScaleFactors = [(2, 2.1), (1, 1.26)], Fractions = [(2, 0.2)])
sig_tres = TimeResolution(Name = 'sig_tres', **tres_args)

# Resolution models
bkg_tres = TimeResolution(Name = 'bkg_tres', ParNamePrefix = 'bkg', **tres_args)

# J/psi mass pdf
from P2VV.Parameterizations.MassPDFs import DoubleCB_Psi_Mass as PsiMassPdf
psi_m = PsiMassPdf(mpsi, Name = 'psi_m', mpsi_alpha_1 = dict(Name = 'mpsi_alpha_1', Value = 2., Constant = True))
psi_m = psi_m.pdf()
    
# J/psi background
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf
bkg_mpsi = PsiBkgPdf(mpsi, Name = 'bkg_mpsi')

# Longlived time PDF
from P2VV.Parameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = time_obs, resolutionModel = sig_tres.model()
                         , psi_t_fml    = dict(Name = 'psi_t_fml',    Value = 6.7195e-01)
                         , psi_t_ll_tau = dict(Name = 'psi_t_ll_tau', Value = 1.3672, MinMax = (0.5,  2.5))
                         , psi_t_ml_tau = dict(Name = 'psi_t_ml_tau', Value = 1.3405e-01, MinMax = (0.01, 0.5))
                         )
psi_t = psi_t.pdf()

# Sigmat PDFs
from P2VV.Parameterizations.SigmatPDFs import DoubleLogNormal
sig_ln = DoubleLogNormal(st, ParNamePrefix = 'sig')
bkg_ln = DoubleLogNormal(st, ParNamePrefix = 'bkg')

# Background time pdf
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
bkg_t = Prompt_Peak(time_obs, bkg_tres.model(), Name = 'bkg_pdf')

# J/psi signal component
psi_ll = Component('psi_ll', (psi_m, psi_t, sig_ln.pdf()), Yield= (8.5575e+03,100,500000) )

# Background component
background = Component('background', (bkg_mpsi.pdf(), bkg_t.pdf(), bkg_ln.pdf()),
                       Yield = (19620,100,500000) )

# Prompt component
prompt_pdf = Prompt_Peak(time_obs, sig_tres.model(), Name = 'prompt_pdf')
prompt = Component('prompt', (prompt_pdf.pdf(), psi_m, sig_ln.pdf()), Yield = (160160, 100, 500000))

def make_wpv_pdf(prefix):
    wpv_mean = RealVar('%swpv_mean' % prefix, Value = 0, MinMax = (-1, 1), Constant = True)
    wpv_sigma = RealVar('%swpv_sigma' % prefix, Value = 0.305, MinMax = (0.01, 10), Constant = True)
    return Pdf(Name = '%swpv_pdf' % prefix, Type = Gaussian, Parameters = (time_obs, wpv_mean, wpv_sigma))
sig_wpv_pdf = make_wpv_pdf('sig_')

sig_wpv = Component('sig_wpv', (sig_wpv_pdf, psi_m, sig_ln.pdf()), Yield = (552, 1, 50000))

components = [prompt, psi_ll, background, sig_wpv]

## PDF to generate with
gen_obs = (time_obs, mpsi, st)
gen_pdf = buildPdf(Components = components, Observables = gen_obs, Name = 'gen_pdf')

## Load parameter values
from ROOT import TFile
param_file = TFile.Open("gen_params.root")
gen_params = param_file.Get("gen_params")
gen_pdf_pars = gen_pdf.getParameters(RooArgSet(*gen_obs))

for p in gen_pdf_pars:
     gp = gen_params.find(p)
     if not gp and not p.isConstant():
         print 'cannot find %s' % p.GetName()
     elif not p.isConstant():
         p.setVal(gp.getVal())
         p.setError(gp.getError())

## Fit options
fitOpts = dict(  Optimize  = 2
               , Timer     = True
               , Strategy  = 1
               , Save      = True
               , Verbose   = False
               , Minimizer = 'Minuit2')

toy.set_fit_opts(**fitOpts)

if options.sFit:
    # Components for sFit
    components = [prompt, psi_ll, sig_wpv]

    ## PDF for sFit
    fit_pdf = buildPdf(Components = components, Observables = (time_obs,), Name = 'fit_pdf')

    ## SWeights
    from P2VV.ToyMCUtils import SWeightTransform
    mass_pdf = buildPdf(Components = (prompt, background), Observables = (mpsi,), Name = 'mass_pdf')
    toy.set_transform(SWeightTransform(mass_pdf, 'prompt', fitOpts))
else:
    # Components for cFit
    components = [prompt, psi_ll, background, sig_wpv]

    ## PDF for cFit
    fit_pdf = buildPdf(Components = components, Observables = (time_obs, mpsi), Name = 'fit_pdf')

toy.run(Observables = gen_obs, Pdf = fit_pdf, GenPdf = gen_pdf)

toy.write_output()
