from ROOT import *

gSystem.Load("libp2vv")
from math import sqrt,pi
from array import array

from RooFitDecorators import *

acceptancecurvefilename = '/data/bfys/dveijk/DataJpsiPhi/2011/effhistoOverlap.root'

wsfile = TFile('AccWS.root')
ws = wsfile.Get('ws')

timeaccfile = TFile(acceptancecurvefilename)
effhisto = timeaccfile.Get('effhisto')

effdatahist = RooDataHist('eff','eff',RooArgList(ws['t']),effhisto)
#ws.factory("DataHist::eff(t,%s)"%effhisto)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,eff)")

nbins = effhisto.GetNbinsX()

#read and set limits and contents of bins
for i in range(nbins):
    vars()['height%s'%str(i+1)]= RooRealVar('height%s'%str(i+1),'height%s'%str(i+1),effhisto.GetBinContent(i+1),0.,1.)
    vars()['limit%s'%str(i)]= RooRealVar('limit%s'%str(i),'limit%s'%str(i),effhisto.GetBinLowEdge(i+1))

#last bin limit
vars()['limit%s'%nbins] = RooRealVar('limit%s'%nbins,'limit%s'%nbins,effhisto.GetBinLowEdge(nbins)+effhisto.GetBinWidth(nbins))

coefList = RooArgList()
limitList = RooArgList()
for i in range(nbins):
    coefList.add(vars()['height%s'%str(i+1)])
    limitList.add(vars()['limit%s'%str(i)])
#add last bin limit
limitList.add(vars()['limit%s'%nbins])

## # set kTRUE if you want 1st order connected lines instead of 0th order steps
interpolateBool = False

BiasedEff = RooStepFunction("StepEff","StepEff",ws['t'],coefList,limitList,interpolateBool)

canvas = TCanvas()
canvas.Divide(4,1)

canvas.cd(1)
effhisto.Draw()

canvas.cd(2)
tframe = ws['t'].frame()
effdatahist.plotOn(tframe)
tframe.Draw()

canvas.cd(3)
tframe2 = ws['t'].frame()
ws['effpdf'].plotOn(tframe2)
tframe2.Draw()

canvas.cd(4)
tframe3 = ws['t'].frame()
BiasedEff.plotOn(tframe3)
tframe3.Draw()
