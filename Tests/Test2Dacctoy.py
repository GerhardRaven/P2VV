from ROOT import *
gSystem.Load("libP2VV")
## from math import sqrt,pi

from RooFitDecorators import *

ws = RooWorkspace("ws")

#########################################
### Define variables and simple PDF's ###
#########################################
ws.factory("{t[0,10],m[10,30]}")

#ws.factory("RooTruthModel::tres(t)")
ws.factory("GaussModel::tres(t,tres_m[0],tres_s[0.0005])")

ws.factory("Gaussian::m_pdf(m,m_mean[15,13,17],m_sigma[1.,0.,2.])")

#ws.factory("Gaussian::t_pdf(t,t_mean[5.,0.,10.],t_sigma[1.,0.,2.])")
ws.factory("RooDecay::t_pdf(t,tau[2,0.01,4.0],tres,SingleSided)")
#ws.factory("Uniform::t_pdf(t)")

ws.factory("PROD::pdf( m_pdf, t_pdf)")

ws.defineSet("observables","t,m")

pdf = ws['pdf']

#####################
### Generate data ###
#####################
data = ws['pdf'].generate(ws.set('observables'),10000)

##############################################
### Define acceptance function a la Wouter ###
##############################################
nbins = 10
effh1 = TH1F("effh1","effh1",nbins,ws['t'].getMin(),ws['t'].getMax()) 
for i in range(0,6):
    effh1.SetBinContent(i,0.1*i)
for i in range(6,nbins+1):
    effh1.SetBinContent(i,1) 
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effh1) 
ws.put(effdatahist)

ws.factory("HistPdf::effpdf(t,effdatahist)")

#Case1
ws.factory("SingleHistEfficiency::goodaccpdf(pdf,effpdf)")

#Case2
ws.factory("SingleHistEfficiency::taccpdf(t_pdf,effpdf)")
ws.factory("PROD::badaccpdf(taccpdf,m_pdf)")

##############################
### Generate 'biased' data ###
##############################
goodaccpdf = ws['goodaccpdf']
goodaccdata = ws['goodaccpdf'].generate(ws.set('observables'),10000)

badaccpdf = ws['badaccpdf']
badaccdata = ws['badaccpdf'].generate(ws.set('observables'),10000)

goodaccpdf.fitTo(goodaccdata)
badaccpdf.fitTo(badaccdata)

############
### Plot ###
############
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

#Make Sanity plots
Acc = TCanvas('Acc','Acc')
Acc.Divide(2)

Acc.cd(1)
tframe = ws['t'].frame(RooFit.Bins(30))
ws['effpdf'].plotOn(tframe)
tframe.Draw()

Acc.cd(2)
mframe = ws['m'].frame(RooFit.Bins(30))
ws['effpdf'].plotOn(mframe)
mframe.Draw()

C2 = TCanvas('C2','C2')
C2.Divide(3,3)

C2.cd(1)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(2)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
pdf.plotOn(tframe,lw,RooFit.LineColor(kGreen))
tframe.Draw()

C2.cd(3)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw,RooFit.LineColor(kGreen))
tframe.Draw()

C2.cd(4)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
goodaccdata.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(5)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
goodaccpdf.plotOn(tframe,lw,RooFit.Project(RooArgSet(ws['m'])))
tframe.Draw()

C2.cd(6)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
goodaccdata.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
goodaccpdf.plotOn(tframe,lw)
tframe.Draw()

C2.cd(7)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
badaccdata.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(8)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
#badaccpdf.plotOn(tframe,lw,RooFit.Project(RooArgSet(ws['m'])))
badaccpdf.plotOn(tframe,lw,RooFit.LineColor(kRed))
tframe.Draw()

C2.cd(9)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
badaccdata.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
badaccpdf.plotOn(tframe,lw,RooFit.LineColor(kRed))
tframe.Draw()

#Make Sanity plots
C3 = TCanvas('C3','C3')
C3.Divide(3,3)

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
goodaccdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
mframe.Draw()

C3.cd(5)
mframe = ws['m'].frame(RooFit.Bins(10))
goodaccpdf.plotOn(mframe,lw)
mframe.Draw()

C3.cd(6)
mframe = ws['m'].frame(RooFit.Bins(10))
goodaccdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
goodaccpdf.plotOn(mframe,lw)
mframe.Draw()

C3.cd(7)
mframe = ws['m'].frame(RooFit.Bins(10))
badaccdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
mframe.Draw()

C3.cd(8)
mframe = ws['m'].frame(RooFit.Bins(10))
badaccpdf.plotOn(mframe,lw,RooFit.LineColor(kRed))
mframe.Draw()

C3.cd(9)
mframe = ws['m'].frame(RooFit.Bins(10))
badaccdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
badaccpdf.plotOn(mframe,lw,RooFit.LineColor(kRed))
mframe.Draw()
