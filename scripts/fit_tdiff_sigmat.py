#!/usr/bin/env python
import gc
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog file')

parser.add_option("--refit", dest = "refit", default = False,
                  action = 'store_true', help = 'Use refitted observables.')
parser.add_option("-d", '--dataset-name', dest = "ds_name", default = 'sdata',
                  action = 'store', type = 'string', help = 'Name of dataset in file.')
parser.add_option('--momentum-calib', dest = "momentum_calib", default = False,
                  action = 'store_true', help = 'Get the mean of the distribution in bins of t_true')

(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)
if not os.path.exists(args[0]):
    print 'error: file %s does not exist' % args[0]
    print parser.print_usage()
    sys.exit(-2)
    
from P2VV.Load import LHCbStyle

from ROOT import RooRealVar
from ROOT import RooGaussian
from ROOT import RooConstVar
from ROOT import RooGExpModel
from ROOT import RooFormulaVar
from ROOT import RooAddPdf
from ROOT import RooArgList
from ROOT import RooArgSet

t_true = RooRealVar('truetime', 'truetime', 0, 14)
t_diff_st = RooRealVar('time_diff_sigmat' if not options.refit else 'time_diff_sigmat_refit',
                       'time / sigmat', -30, 40)
st = RooRealVar('sigmat' if not options.refit else 'sigmat_refit', '#sigma_{t}', 0.0001, 0.12)

from ROOT import TFile
f = TFile(args[0])
sdata = f.Get(options.ds_name)

mean_offset = RooRealVar('mean_offset', 'mean_offset', -0.001723, -0.1, 0.1)
mean_slope = RooRealVar('mean_slope', 'mean_slope', -0.00431, -0.1, 0.1)
mean_quad = RooRealVar('mean_quad', 'mean_quad', -0.00380, -0.1, 0.1)
st_mean = RooConstVar('st_mean', 'st_mean', 0.03276)
formula = '@2 + @3 * (@0 - @1) + @4 * (@0 - @1) * (@0 - @1)'
if options.momentum_calib:
    mean = RooRealVar('mean', 'mean', 0., -0.2, 0.1)
else:
    params = [st, st_mean, mean_offset, mean_slope, mean_quad]
    mean = RooFormulaVar("mean_quad", "mean_quad", formula, RooArgList(*params))

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
sigma_gexp = RooRealVar("sigma_gexp", "sigma_gexp", 10, 0.1, 50)
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
from ROOT import RooWorkspace
from P2VV.RooFitWrappers import RooObject
obj = RooObject(workspace = 'w')
w = obj.ws()
w.put(model)
w.put(sdata)

from ROOT import TPaveText
if options.momentum_calib:
    year_label = TPaveText(0.71, 0.32, 0.89, 0.45, "NDC")
else:
    year_label = TPaveText(0.71, 0.72, 0.89, 0.85, "NDC")

year_label.SetFillColor(0)
pos = args[0].find('201')
year_label.AddText(args[0][pos : pos + 4])
year_label.SetBorderSize(0)

fitOpts = dict(NumCPU = 8, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 1,
               Offset = True, Strategy = 1)

pos = args[0].find('MC')
from P2VV.Utilities.Resolution import plot_dir

if options.momentum_calib:
    from array import array
    tt_bins = array('d', [0.0, 0.09928, 0.2059, 0.3207, 0.447, 0.5844, 0.7365, 0.9077, 1.101, 1.325, 1.592, 1.916, 2.332, 2.916, 3.92, 14.0])

    from ROOT import RooBinning
    tt_binning = RooBinning(len(tt_bins) - 1, tt_bins)
    tt_binning.SetName('tt_binning')
    t_true.setBinning(tt_binning, 'tt_binning')

    from P2VV.RooFitWrappers import BinningCategory
    tt_cat = BinningCategory(t_true.GetName() + '_cat', Observable = t_true, Binning = tt_binning,
                             CatTypeName = 'tt_bin_', Data = sdata, Fundamental = True)

    split_pdf = w.factory('SIMCLONE::split_pdf(model,$SplitParam(sigma_sigma,truetime_cat),'
                          + '$SplitParam(av_sigma,truetime_cat),$SplitParam(mean,truetime_cat),'
                          + '$SplitParam(frac_gexps,truetime_cat))')

    fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 1,
                   Offset = True, Strategy = 1)
    for i in range(3):
        result = split_pdf.fitTo(sdata, SumW2Error = False, **fitOpts)
        if result.status() == 0 and abs(result.minNll()) < 5e5:
            break
    assert(result.status() == 0)
    
    fpf = result.floatParsFinal()
    from array import array
    tt_cat = w.cat('truetime_cat')
    x_means = array('d', [])
    x_errs = array('d', [])
    y_means = array('d', [])
    y_errs = array('d', [])
    for s in tt_cat:
        cut = '{0} == {0}::{1}'.format(tt_cat.GetName(), s.GetName())
        par = sdata.meanVar(t_true, cut)
        x_means.append(par.getVal())
        x_errs.append(par.getError())
        par = fpf.find('mean_' + s.GetName())
        y_means.append(par.getVal())
        y_errs.append(par.getError())

    from ROOT import TF1
    l = TF1('pol0', 'pol0', t_true.getMin(), t_true.getMax())

    from ROOT import TGraphErrors
    g = TGraphErrors(len(x_means), x_means, y_means, x_errs, y_errs)
    line_result = g.Fit(l, 'S0+', 'L')

    from ROOT import TCanvas
    cn = 'canvas_' + args[0][pos : -5]
    canvas = TCanvas(cn, cn, 600, 400)
    g.Draw('AP')
    g.GetXaxis().SetTitle("t_{true} [ps]")
    g.GetYaxis().SetTitle("#mu / #sigma_{t}")
    g.GetYaxis().SetRangeUser(-0.12, 0)
    
    from ROOT import kBlue
    l.SetLineColor(kBlue)
    l.Draw('same')
    year_label.Draw()

    plot_name = 'mean_res_st_ttrue_' + args[0][pos : -5] + '.pdf'
    canvas.Print(os.path.join(plot_dir, plot_name), EmbedFonts = True)
else:
    from ROOT import kGreen, kDashed
    from P2VV.Utilities.Plotting import plot
    from ROOT import TCanvas

    for i in range(3):
        result = model.fitTo(sdata, SumW2Error = False, **fitOpts)
        if result.status() == 0 and abs(result.minNll()) < 5e5:
            break

    canvas = TCanvas("canvas", "canvas", 600, 400)
    plot(canvas, t_diff_st, pdf = model, data = sdata, logy = True,
         frameOpts = dict(Range = (-20, 20)),     
         yTitle = 'Candidates / (0.5)', dataOpts = dict(Binning = 80),
         xTitle = '(t - t_{true}) / #sigma_{t}',
         yScale = (1, 400000),
         pdfOpts = dict(ProjWData = (RooArgSet(st), sdata, True)),
         components = {'gexps' : dict(LineColor = kGreen, LineStyle = kDashed)})
    year_label.Draw()
    plot_name = 'tdiff_sigmat_' + args[0][pos : -5] + '.pdf'
    canvas.Print(os.path.join(plot_dir, plot_name), EmbedFonts = True)
