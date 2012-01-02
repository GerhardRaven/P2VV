from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *

ws = RooWorkspace("ws")

nbins = 1000

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s", "--seed", dest="seed",type="int",help="set random seed for data generation")

(options, args) = parser.parse_args()

print 'Seed for data generation =', options.seed
generator = RooRandom.randomGenerator().SetSeed(int(options.seed))

#########################################
### Define variables and simple PDF's ###
#########################################
ws.factory("t[-0.005,14]")

tau = 2.1
ws.factory("Exponential::t_exp(t,tau[-%s,-10,-0.5])"%(tau))
#ws.factory("TruthModel::tres(t)")
#ws.factory("GaussModel::tres(t,tres_m[0],tres_s[0.0005])")
#ws.factory("Decay::t_exp(t,tau[%s,0.5,10.],tres,SingleSided)"%(tau))

ws.defineSet("observables","t")

pdf = ws['t_exp']

#####################
### Generate data ###
#####################
data = pdf.generate(ws.set('observables'),100000)

ws.factory("expr::eff('1-a*t',t,a[0.0157])")
effhist = ws['eff'].createHistogram('effhist',ws['t'],RooFit.Binning(nbins,ws['t'].getMin(),ws['t'].getMax()))
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effhist)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,effdatahist)")
ws.factory("EffHistProd::accpdf1(t_exp,effpdf)")

accpdf1 = ws['accpdf1']

##############################
### Generate 'biased' data ###
##############################
accdata1 = accpdf1.generate(ws.set('observables'),100000)

ws.put(accdata1)
ws.addClassDeclImportDir("..")
ws.addClassImplImportDir("..")
ws.importClassCode()

###################
### Make FitSet ###
###################

print "Fit AccPDF"
print '##################################################'
accresult = accpdf1.fitTo(accdata1, RooFit.Save(True))
ws.put(accresult)
#acctau = ws['tau'].getVal()
#accgamma = 1./acctau

ws['tau'].setVal(-tau)
print "Fit normal PDF to AccData"
print '##################################################'
noncorrresult = pdf.fitTo(accdata1, RooFit.Save(True))
ws.put(noncorrresult)

from time import sleep
print "#####################"
print "BEGIN SLEEP"
sleep(20)
print "END SLEEP"
print "#####################"
import sys

#noncorrtau = ws['tau'].getVal()
#a = ws['a'].getVal()
#corrtau = ((1.+a*noncorrtau)/(4.*a))-((sqrt(a**2*noncorrtau**2-6*a*noncorrtau+1))/(4*a))
#corrgamma = 1/corrtau

#print "***"
#print "tau = ", tau
#print "gamma = ", (1.)/(tau)
#print "***"
#print "acctau = ", acctau
#print "accgamma = ", accgamma
#print "***"
#print "corrtau = ", corrtau
#print "corrgamma = ", corrgamma
#print "***"

f = TFile.Open("test_EffHistProd.root", "recreate")
f.WriteTObject(ws)
f.Close()

sys.exit(0)

############
### Plot ###
############
## lw = RooCmdArg(RooFit.LineWidth(2))
## xes = RooCmdArg(RooFit.XErrorSize(0))
## err = RooCmdArg(RooFit.DrawOption('E'))
## dashed = RooCmdArg(RooFit.LineStyle(kDashed))

## #Make Sanity plots
## Acc = TCanvas('Acc','Acc')

## Acc.cd(1)
## tframe = ws['t'].frame(RooFit.Bins(30))
## ws['effpdf'].plotOn(tframe)
## tframe.Draw()

## C2 = TCanvas('C2','C2')
## C2.Divide(3,2)

## C2.cd(1)
## #gPad.SetLogy()
## tframe = ws['t'].frame(RooFit.Bins(30))
## data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
## tframe.Draw()

## C2.cd(2)
## #gPad.SetLogy()
## tframe = ws['t'].frame(RooFit.Bins(30))
## pdf.plotOn(tframe,lw)
## tframe.Draw()

## C2.cd(3)
## #gPad.SetLogy()
## tframe = ws['t'].frame(RooFit.Bins(30))
## data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
## pdf.plotOn(tframe,lw)
## tframe.Draw()

## C2.cd(4)
## #gPad.SetLogy()
## tframe = ws['t'].frame(RooFit.Bins(30))
## accdata1.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
## tframe.Draw()

## C2.cd(5)
## #gPad.SetLogy()
## tframe = ws['t'].frame(RooFit.Bins(30))
## accpdf1.plotOn(tframe,lw)
## tframe.Draw()

## C2.cd(6)
## #gPad.SetLogy()
## tframe = ws['t'].frame(RooFit.Bins(30))
## accdata1.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
## accpdf1.plotOn(tframe,lw)
## tframe.Draw()
