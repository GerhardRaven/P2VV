from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *

ws = RooWorkspace("ws")

#########################################
### Define variables and simple PDF's ###
#########################################
ws.factory("{psi[0.2,%f,%f],m[10,30]}"%(0,2*pi))

ws.factory("expr::een('a',{a[1],psi})")
ws.factory("expr::cospsi('cos(psi)',psi)")

ceen = RooRealVar('ceen','ceen',1.)
ccospsi = RooRealVar('ccospsi','ccospsi',1.)

psi_pdf = RooRealSumPdf('psi_pdf','psi_pdf',RooArgList(ws['een'],ws['cospsi']),RooArgList(ceen,ccospsi))
ws.put(psi_pdf)

ws.factory("Gaussian::m_pdf(m,m_mean[15,13,17],m_sigma[1.,0.,2.])")

ws.defineSet("observables","psi,m")

ws.factory("PROD::pdf( m_pdf, psi_pdf)")

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

#Case1
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

C2 = TCanvas('C2','C2')
C2.Divide(3,2)

C2.cd(1)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
data.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
psiframe.Draw()

C2.cd(2)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
pdf.plotOn(psiframe,lw,RooFit.LineColor(kGreen))
psiframe.Draw()

C2.cd(3)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
data.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(psiframe,lw,RooFit.LineColor(kGreen))
psiframe.Draw()

C2.cd(4)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
accdata.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
psiframe.Draw()

C2.cd(5)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
accpdf.plotOn(psiframe,lw)
psiframe.Draw()

C2.cd(6)
#gPad.SetLogy()
psiframe = ws['psi'].frame(RooFit.Bins(10))
accdata.plotOn(psiframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(psiframe,lw)
psiframe.Draw()

#Make Sanity plots
C3 = TCanvas('C3','C3')
C3.Divide(3,2)

C3.cd(1)
mframe = ws['m'].frame(RooFit.Bins(10))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
mframe.Draw()

C3.cd(2)
mframe = ws['m'].frame(RooFit.Bins(10))
pdf.plotOn(mframe,lw,RooFit.LineColor(kGreen))
mframe.Draw()

C3.cd(3)
mframe = ws['m'].frame(RooFit.Bins(10))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw,RooFit.LineColor(kGreen))
mframe.Draw()

C3.cd(4)
mframe = ws['m'].frame(RooFit.Bins(10))
accdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
mframe.Draw()

C3.cd(5)
mframe = ws['m'].frame(RooFit.Bins(10))
accpdf.plotOn(mframe,lw)
mframe.Draw()

C3.cd(6)
mframe = ws['m'].frame(RooFit.Bins(10))
accdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(mframe,lw)
mframe.Draw()
