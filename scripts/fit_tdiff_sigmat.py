from ROOT import RooRealVar
from ROOT import RooGaussian
from ROOT import RooConstVar
from ROOT import RooGExpModel
from ROOT import RooFormulaVar
from ROOT import RooAddPdf
from ROOT import RooArgList
from ROOT import RooArgSet

t_diff_st = RooRealVar('time_diff_sigmat', 'time_diff_sigmat', -30, 40)
st = RooRealVar("sigmat", "sigmat", 0.0001, 0.12)

from ROOT import TFile
f = TFile("tdiff_sigmat_MC2012.root")
sdata = f.Get("sdata")

sdata_cut = sdata.reduce("sigmat > 0.01 && sigmat < 0.024 && time_diff_sigmat > -6 && time_diff_sigmat < 6")

mean_offset = RooRealVar('mean_offset', 'mean_offset', -0.001723, -0.1, 0.1)
mean_slope = RooRealVar('mean_slope', 'mean_slope', -0.00431, -0.1, 0.1)
mean_quad = RooRealVar('mean_quad', 'mean_quad', -0.00380, -0.1, 0.1)
st_mean = RooConstVar('st_mean', 'st_mean', 0.03276)
formula = '@2 + @3 * (@0 - @1) + @4 * (@0 - @1) * (@0 - @1)'
args = [st, st_mean, mean_offset, mean_slope, mean_quad]
mean = RooFormulaVar("mean_quad", "mean_quad", formula, RooArgList(*args))


# Parameterisation
av_sigma = RooRealVar("av_sigma", "av_sigma", 1.24292, 1, 50)
sigma_sigma = RooRealVar("sigma_sigma", "sigma_sigma", 0.245659, 0.001, 50)
frac_g2 = RooRealVar("frac_g2", "frac_g2", 0.169173, 0.01, 0.99)

# 1st Gauss
sigma1 = RooFormulaVar('sigma1', 'sigma1', '- sqrt(@0 / (1 - @0)) * @1 + @2', RooArgList(frac_g2, sigma_sigma, av_sigma))
g1 = RooGaussian("g1", "g1", t_diff_st, mean, sigma1)

# 2nd Gauss
sigma2 = RooFormulaVar('sigma2', 'sigma2', 'sqrt((1 - @0) / @0) * @1 + @2', RooArgList(frac_g2, sigma_sigma, av_sigma))
g2 = RooGaussian("g2", "g2", t_diff_st, mean, sigma2)

gaussians = RooAddPdf("gaussians", "gaussians", RooArgList(g2, g1), RooArgList(frac_g2))

# 1st GExp
one = RooConstVar("one", "one", 1)
mean_gexp = RooRealVar("mean_gexp", "mean_gexp", -0.0529544, -10, 10)
sigma_gexp = RooRealVar("sigma_gexp", "sigma_gexp", 10, 1, 50)
rlife1 = RooRealVar("rlife1", "rlife1", 3.2416, 0.1, 10)
gexp1 = RooGExpModel("gexp1", "gexp1", t_diff_st, mean_gexp, sigma_gexp, rlife1, one, one, one)

# 2nd GExp
rlife2 = RooRealVar("rlife2", "rlife2", 6.07734, 0.1, 10)
gexp2 = RooGExpModel("gexp2", "gexp2", t_diff_st, mean_gexp, sigma_gexp, rlife2, one, one, one, False, RooGExpModel.Flipped)

frac_gexp2 = RooRealVar("frac_gexp2", "frac_gexp2", 0.184357, 0.01, 0.99)
gexps = RooAddPdf("gexps", "gexps", RooArgList(gexp2, gexp1), RooArgList(frac_gexp2))

frac_gexps = RooRealVar("frac_gexps", "frac_gexps", 0.0118392, 0.001, 0.99)
model = RooAddPdf("model", "model", RooArgList(gexps, gaussians), RooArgList(frac_gexps))
model.setParameterizeIntegral(RooArgSet(st))

from P2VV import RooFitDecorators
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 1, Offset = True,
               Strategy = 1)
result = model.fitTo(sdata, SumW2Error = False, **fitOpts)

from ROOT import TCanvas
from ROOT import kGreen, kDashed
from P2VV.Utilities.Plotting import plot
from P2VV.Load import LHCbStyle
canvas = TCanvas("canvas", "canvas", 600, 400)
plot(canvas, t_diff_st, pdf = model, data = sdata, logy = True,
     frameOpts = dict(Range = (-20, 20)),     
     yTitle = 'Candidates / (0.5)', dataOpts = dict(Binning = 80),
     xTitle = '(t_{rec} - t_{true}) / #sigma_{t}',
     pdfOpts = dict(ProjWData = (RooArgSet(st), sdata, True)),
     components = {'gexps' : dict(LineColor = kGreen, LineStyle = kDashed)})
