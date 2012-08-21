from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))

from ROOT import TFile
fit_input = TFile.Open('fit_result.root')
fit_result = fit_input.Get('fitresult_signal_Paper2012_FitAcceptance_X_sig_t_angles_x_Eff_X_sig_tagCat_iTag_JpsiphiData_weighted_sigMass')
final_pars = fit_result.floatParsFinal()
hlt1_heights = []
for p in final_pars:
    if p.GetName().find('hlt1') == -1:
        continue
    hlt1_heights.append((p.getVal(), p.getError()))

input_file = '/stuff/PhD/p2vv/data/start_values.root'
histogram = 'hlt1_shape'
from ROOT import TFile
acceptance_file = TFile.Open(input_file)
if not acceptance_file:
    raise ValueError, "Cannot open histogram file %s" % input_file

from array import array
from ROOT import TH1D

fixed_file = TFile.Open('/stuff/PhD/p2vv/data/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root')
fixed_hlt1_shape = fixed_file.Get('Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1B_10bins')

hlt1_shape = acceptance_file.Get(histogram)
xaxis = hlt1_shape.GetXaxis()

bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, hlt1_shape.GetNbinsX() + 2)))
hlt1_hist = TH1D('hlt1_shape', 'hlt1_shape', len(bins) - 1, bins)
for i in range(1, len(bins)):
    hlt1_hist.SetBinContent(i, hlt1_heights[i - 1][0])
    hlt1_hist.SetBinError(i, hlt1_heights[i - 1][1])
hlt1_hist.GetXaxis().SetTitle('proper time (ps)')
hlt1_hist.GetYaxis().SetTitle('#epsilon (au)')
hlt1_hist.GetYaxis().SetTitleOffset(1.2)

from ROOT import TCanvas, gStyle
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
from ROOT import kYellow, kOrange
canvas = TCanvas('canvas', 'canvas', 1000, 500)
canvas.Divide(2, 1)
p = canvas.cd(1)
fixed_hlt1_shape.Draw()
