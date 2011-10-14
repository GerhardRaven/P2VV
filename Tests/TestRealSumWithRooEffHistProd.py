from ROOT import *
gSystem.Load("libp2vv")
## from math import sqrt,pi

from RooFitDecorators import *
from math import pi
ws = RooWorkspace("ws")

#########################################
### Define variables and simple PDF's ###
#########################################
ws.factory("{psi[0.2,%f,%f]}"%(0,2*pi))

ws.factory("expr::een('a',{a[1],psi})")
ws.factory("expr::cospsi('cos(psi)',psi)")

ceen = RooRealVar('ceen','ceen',2.)
ccospsi = RooRealVar('ccospsi','ccospsi',1.)

pdf = RooRealSumPdf('pdf','pdf',RooArgList(ws['een'],ws['cospsi']),RooArgList(ceen,ccospsi))
ws.put(pdf)
ws.defineSet("observables","psi")

pdf = ws['pdf']

#####################
### Generate data ###
#####################
data = pdf.generate(ws.set('observables'),10000)

##############################################
### Define acceptance function a la Wouter ###
##############################################
nbins = 10
effh1 = TH1F("effh1","effh1",nbins,ws['psi'].getMin(),ws['psi'].getMax()) 
for i in range(0,6):
    effh1.SetBinContent(i,0.1*i)
for i in range(6,nbins+1):
    effh1.SetBinContent(i,1) 
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['psi']),effh1) 
ws.put(effdatahist)

ws.factory("HistPdf::effpdf(psi,effdatahist)")

ws.factory("EffHistProd::accpdf(pdf,effpdf)")

##############################
### Generate 'biased' data ###
##############################
accpdf = ws['accpdf']
accdata = ws['accpdf'].generate(ws.set('observables'),10000)

accpdf.fitTo(accdata)

############
### Plot ###
############
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

#Make Sanity plots
Acc = TCanvas('Acc','Acc')

psiframe = ws['psi'].frame(RooFit.Bins(30))
ws['effpdf'].plotOn(psiframe)
psiframe.Draw()

psicanvas = TCanvas('psicanvas','psicanvas')
psicanvas.Divide(3,2)

psicanvas.cd(1)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
data.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
psiframe.Draw()

psicanvas.cd(2)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
pdf.plotOn(psiframe,lw,RooFit.LineColor(kGreen))
psiframe.Draw()

psicanvas.cd(3)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
data.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(psiframe,lw,RooFit.LineColor(kGreen))
psiframe.Draw()

psicanvas.cd(4)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
accdata.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
psiframe.Draw()

psicanvas.cd(5)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
accpdf.plotOn(psiframe,lw)
psiframe.Draw()

psicanvas.cd(6)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
accdata.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(psiframe,lw)
psiframe.Draw()
