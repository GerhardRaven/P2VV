from ROOT import *
gSystem.Load("libp2vv")
## from math import sqrt,pi

from RooFitDecorators import *

ws = RooWorkspace("ws")

#########################################
### Define variables and simple PDF's ###
#########################################
ws.factory("{t[0,10]}")

ws.factory("Gaussian::t_sig(t,t_mean[5.,0.,10.],t_sigma[0.5.,0.,2.])")
ws.factory("Exponential::t_bkg(t,t_bkg_exp[-0.2,-0.01,-0.0001])")

ws.factory("SUM::pdf( f_sig[0.71,0.,1.]*t_sig, t_bkg)")

ws.defineSet("observables","t")

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
ws.factory("EffHistProd::accpdf1(pdf,effpdf)")

#Case2
ws.factory("EffHistProd::sigaccpdf(t_sig,effpdf)")
ws.factory("EffHistProd::bkgaccpdf(t_bkg,effpdf)")
ws.factory("SUM::accpdf2( f_sig*sigaccpdf, bkgaccpdf)")

##############################
### Generate 'biased' data ###
##############################
accpdf1 = ws['accpdf1']
accdata1 = ws['accpdf1'].generate(ws.set('observables'),10000)

accpdf2 = ws['accpdf2']
accdata2 = ws['accpdf2'].generate(ws.set('observables'),10000)

print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print 'Fitting eff(t)*[sigpdf*f_sig+bkfpdf]'
print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
accpdf1.fitTo(accdata1)

print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
print 'Fitting [eff(t)*sigpdf]*f_sig+[eff(t)*bkfpdf]'
print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
accpdf2.fitTo(accdata2)

############
### Plot ###
############
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

#Make Sanity plots
Acc = TCanvas('Acc','Acc')

tframe = ws['t'].frame(RooFit.Bins(30))
ws['effpdf'].plotOn(tframe)
tframe.Draw()

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
accdata1.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(5)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
accpdf1.plotOn(tframe,lw)
tframe.Draw()

C2.cd(6)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
accdata1.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
accpdf1.plotOn(tframe,lw)
tframe.Draw()

C2.cd(7)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
accdata2.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(8)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
accpdf2.plotOn(tframe,lw,RooFit.LineColor(kRed))
tframe.Draw()

C2.cd(9)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(10))
accdata2.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
accpdf2.plotOn(tframe,lw,RooFit.LineColor(kRed))
tframe.Draw()
