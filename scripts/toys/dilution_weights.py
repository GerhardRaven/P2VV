import sys
import os
try:
    from P2VV.ToyMCUtils import Toy
except ImportError:
    print 'Could not import P2VV version, trying local directory'
    from ToyMCUtils import Toy

toy = Toy()
(options, args) = toy.configure()

from P2VV.RooFitWrappers import *
obj = RooObject( workspace = 'w')

mpsi = RealVar('mpsi', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3025, 3165))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (-1, 1))

from P2VV.Parameterizations.SigmatPDFs import DoubleLogNormal
dln = DoubleLogNormal(st, frac_ln2 = dict(Value = 0.312508), k1 = dict(Value = 0.801757),
                      k2 = dict(Value = 1.37584), median = dict(Value = 0.0309409))
ln = dln.pdf()

# Resolution models
from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as Multi_TimeResolution
tres_args = dict(time = t, sigmat = st, PerEventError = True, Cache = True)
tres_1 = Multi_TimeResolution(Name = 'tres', ParNamePrefix = 'one',
                              ScaleFactors = [(2, 2.00), (1, 1.174)],
                              Fractions = [(2, 0.143)], **tres_args)
tres_2 = TimeResolution(Name = 'tres', ParNamePrefix = 'two',
                        timeResSigmaSF = dict(Value = 4.), **tres_args)

# Gaussians for Time
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
g1 = Prompt_Peak(t, tres_1.model(), Name = 'g1')

g2 = Prompt_Peak(t, tres_2.model(), Name = 'g2')

# Mass shapes
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf
bkg_m = PsiBkgPdf(mpsi, Name = 'bkg_mpsi')

from P2VV.Parameterizations.MassPDFs import DoubleCB_Psi_Mass as PsiMassPdf
sig_m = PsiMassPdf(mpsi, Name = 'psi_m', mpsi_alpha_1 = dict(Value = 2, Constant = True),
                   mpsi_sigma_sf = dict(Value = 2.3), mpsi_frac = dict(Value = 0.8))

one = Component('one', [g1.pdf(), sig_m.pdf(), ln], Yield = (25000, 100, 100000))
two = Component('two', [g2.pdf(), bkg_m.pdf(), ln], Yield = (25000, 100, 100000))

pdf = buildPdf(Name = 'pdf', Components = [one, two], Observables = [mpsi, t, st])
mass_pdf = buildPdf(Name = 'mass_pdf', Components = [one, two], Observables = [mpsi])

from P2VV.Utilities.Resolution import SplitSigmat
split = SplitSigmat('', st)
sigmat_cat = split.split_cats()[0]

from

## Fit options
fitOpts = dict(  Optimize  = 2
               , Timer     = True
               , Strategy  = 1
               , Save      = True
               , Verbose   = False
               , Minimizer = 'Minuit2')

toy.set_fit_opts(**fitOpts)

## SWeights
from P2VV.ToyMCUtils import SWeightTransform
toy.set_transform(SWeightTransform(mass_pdf, 'one', fitOpts))

toy.run(Observables = [mpsi, t, st], Pdf = pdf)

toy.write_output()
