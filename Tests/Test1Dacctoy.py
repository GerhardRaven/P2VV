from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *
import rootStyle
from ModelBuilders import _buildAngularFunction

myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

from ModelBuilders import *

ws = RooWorkspace("ws")

RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Integration))

#########################################
### Define variables and simple PDF's ###
#########################################
ws.factory("t[0,10]")

#ws.factory("TruthModel::tres(t)")
ws.factory("GaussModel::tres(t,tres_m[0],tres_s[0.0005])")

ws.factory("Decay::pdf(t,tau[2,0.01,4.0],tres,SingleSided)")
#ws.factory("Gaussian::pdf(t,tmean[5.,1.0,3.],tsigma[1.,0.,2.])")
#ws.factory("Uniform::pdf(t)")

ws.defineSet("observables","t")

pdf = ws['pdf']

#####################
### Generate data ###
#####################
data = pdf.generate(ws.set('observables'),10000)

##############################################
### Define acceptance function a la Wouter ###
##############################################
nbins = 10
effh1 = TH1F("effh1","effh1",nbins,ws['t'].getMin(),ws['t'].getMax()) 
for i in range(1,6):
    effh1.SetBinContent(i,0.1*i)
for i in range(6,nbins+1):
    effh1.SetBinContent(i,1) 
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effh1) 
ws.put(effdatahist)

ws.factory("HistFunc::eff(t,effdatahist)")
ws.factory("RooEffHistProd::accpdf1(pdf,eff)")

accpdf1 = ws['accpdf1']

##############################
### Generate 'biased' data ###
##############################
accdata1 = accpdf1.generate(ws.set('observables'),10000)

ws.put(accdata1)
ws.addClassDeclImportDir("..")
ws.addClassImplImportDir("..")
ws.importClassCode()
## accpdf1.fitTo(accdata1, Optimize = 1)

############
### Plot ###
############
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

#Make Sanity plots
Acc = TCanvas('Acc','Acc')

Acc.cd(1)
tframe = ws['t'].frame(RooFit.Bins(30))
ws['eff'].plotOn(tframe)
tframe.Draw()

C2 = TCanvas('C2','C2')
C2.Divide(3,2)

C2.cd(1)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(2)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
pdf.plotOn(tframe,lw)
tframe.Draw()

C2.cd(3)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw)
tframe.Draw()

C2.cd(4)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
accdata1.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

#effhistprod = ws['accpdf1']
#effhistprod.getAnalyticalIntegralWN(RooArgSet(ws['t']),RooArgSet(ws['t']),RooArgSet(ws['t']))

#print '*******************************'
#fullint = pdf.analyticalIntegralWN(1,0,'')
#print 'full integral pdf = ', fullint
#print '*******************************'

## ws['t'].setRange('MyPythonRange',0,1)

## print '*******************************'
## rangeint = pdf.analyticalIntegralWN(1,0,'MyPythonRange')
## print 'range integral pdf = ', rangeint

## print '*******************************'
## normint = pdf.analyticalIntegralWN(1,RooArgSet(ws['t']),'MyPythonRange')
## print 'normintegral pdf = ', normint

#print '*******************************'
#fullint = accpdf1.analyticalIntegralWN(1,0,'')
#print 'full integral accpdf = ', fullint
#print '*******************************'

## print '*******************************'
## rangeint = accpdf1.analyticalIntegralWN(1,0,'MyPythonRange')
## print 'range integral accpdf = ', rangeint

## print '*******************************'
## normint = accpdf1.analyticalIntegralWN(1,RooArgSet(ws['t']),'MyPythonRange')
## print 'normintegral accpdf = ', normint

C2.cd(5)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
accpdf1.plotOn(tframe,lw)
tframe.Draw()

C2.cd(6)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
accdata1.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
accpdf1.plotOn(tframe,lw)
tframe.Draw()

## npoints = 1000
## dt = (ws['t'].getMax()-ws['t'].getMin())/npoints 
## sum = 0
## for i in range(0,npoints):
##     x = ws['t'].getMin() + i*dt 
##     ws['t'].setVal( x ) 
##     pdfval = accpdf1.getVal(RooArgSet(ws['t']) )
##     print "pdfval: ", pdfval
##     sum += dt * pdfval 
## print "Integral: ", sum 

f = TFile.Open("test_EffHistProd.root", "recreate")
f.WriteTObject(ws)
f.Close()
