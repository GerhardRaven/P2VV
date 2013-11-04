import readline
from P2VV.RooFitWrappers import *
from P2VV.Load import LHCbStyle
from ROOT import *
from array import array

f = TFile.Open("/bfys/raaij/p2vv/data/Bs2JpsiPhi_2011_Prescaled_cache.root")
d = f.Get('1bin_19000.00fs_simple/m934737057402830078')
sdata_dir = d.Get('sdata')
sdata = sdata_dir.Get('full_sig_sdata')

pdf_dir = d.Get('PDFs')
pdf = pdf_dir.Get('time_pdf_double_Comb')

t = sdata.get().find('time')
t.setBins(1900)
st = sdata.get().find('sigmat')

projSet = RooArgSet(st)

bounds = array('d', [-1.5 + i * 0.1 for i in range(13)] + [-0.2 + i * 0.01 for i in range(40)] \
               + [0.2 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])
binning = RooBinning(len(bounds) - 1, bounds)

zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])
zoom_binning = RooBinning(len(zoom_bounds) - 1, zoom_bounds)


canvas = TCanvas('canvas', 'canvas', 600, 400)

from ROOT import TPaveText
LHCbLabel = TPaveText(0.9 - gStyle.GetPadRightMargin(), 0.81,
                      0.77 - gStyle.GetPadRightMargin(), 0.89, 'BRNDC')
LHCbLabel.AddText('LHCb')
LHCbLabel.SetFillColor(0)
LHCbLabel.SetTextAlign(12)
LHCbLabel.SetTextSize(0.072)
LHCbLabel.SetBorderSize(0) 

frame = t.frame(RooFit.Range(-1.5, 8))
sdata.plotOn(frame, RooFit.Invisible())
sdata.plotOn(frame, RooFit.Binning(binning))
pdf.plotOn(frame, RooFit.ProjWData(projSet, sdata, True), RooFit.LineWidth(3))
frame.SetMinimum(0.01)
frame.GetXaxis().SetTitle('Decay time [ps]')
frame.GetYaxis().SetTitle('Candidates / (5 fs)')
frame.GetYaxis().SetRangeUser(0.01, 20000)

pad = canvas.cd()
pad.SetLogy()
frame.Draw()
LHCbLabel.Draw()

canvas_2 = TCanvas('canvas_2', 'canvas_2', 600, 400)
frame = t.frame(RooFit.Range(-0.2, 0.2))
sdata.plotOn(frame, RooFit.Binning(zoom_binning))
pdf.plotOn(frame, RooFit.ProjWData(projSet, sdata, True), RooFit.LineWidth(3))
frame.GetXaxis().SetTitle('Decay time [ps]')
frame.GetYaxis().SetTitle('Candidates / (5 fs)')
frame.Draw()
LHCbLabel.Draw()
