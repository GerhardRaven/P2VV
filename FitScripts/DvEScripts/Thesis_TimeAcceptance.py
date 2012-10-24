from ROOT import TFile, TH1F, TCanvas, TLine
from array import array
import sys

import RootStyle
from ROOT import (gROOT,gStyle,TStyle)
MyStyle = RootStyle.MyStyle()
gROOT.SetStyle(MyStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

tfile = TFile("/data/bfys/dveijk/DataJpsiPhi/2012/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root")

fakehisto = TH1F('fake','fake',1,0.,14.)
fakehisto.SetTitle("")
fakehisto.GetXaxis().SetTitle("Decay time (ps)")
fakehisto.GetYaxis().SetTitle("Acceptance")

fakehisto2 = TH1F('fake2','fake2',1,0.,5.5)
fakehisto2.SetTitle("")
fakehisto2.GetXaxis().SetTitle("Decay time (ps)")
fakehisto2.GetYaxis().SetTitle("Acceptance")

histo = tfile.Get("BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins")
histo.SetMarkerColor(1)
histo.SetLineColor(1)
from ROOT import gStyle
gStyle.SetOptStat(0)

line = TLine(0.3,0,0.3,1.1)
line.SetLineStyle(2)

line2 = TLine(0.3,0.9,0.3,1)
line2.SetLineStyle(2)

#Check:
from ROOT import TCanvas
canvas = TCanvas("canvas","canvas",972,600)
canvas.SetFixedAspectRatio(True)
fakehisto.GetYaxis().SetRangeUser(0.,1.1)

fakehisto.Draw()
line.Draw('same')
histo.Draw('same')

canvas2 = TCanvas("canvas2","canvas2",972,600)
canvas2.SetFixedAspectRatio(True)
histo2 = TH1F(histo)
fakehisto2.GetYaxis().SetRangeUser(0.9,1.)
fakehisto2.GetXaxis().SetRangeUser(0.,5.5)

fakehisto2.Draw()
line2.Draw('same')
histo2.Draw('same')

#pad = canvas.cd(2)
#pad.SetLeftMargin(0.14)


