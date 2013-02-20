from optparse import OptionParser
import sys
parser = OptionParser(usage = '%prog nbins')

(options, args) = parser.parse_args()
if len(args) != 1:
    print parser.usage
    sys.exit(-1)
nbins = int(args[0])

from ROOT import (gSystem, TH1F, TFile,
                  RooDataHist, RooFit,
                  TCanvas)
gSystem.Load("libP2VV")
## from math import sqrt,pi
from P2VV.RooFitDecorators import *

w = RooWorkspace("w")

#########################################
### Define variables and simple PDF's ###
#########################################
w.factory('t[0,10]')
w.factory('m[10,30]')

#w.factory("RooTruthModel::tres(t)")
w.factory("GaussModel::tres(t,tres_m[0],tres_s[0.0005])")

w.factory("Gaussian::m_pdf(m,m_mean[15,13,17],m_sigma[1.,0.,2.])")

#w.factory("Gaussian::t_pdf(t,t_mean[5.,0.,10.],t_sigma[1.,0.,2.])")
w.factory("RooDecay::t_pdf(t,tau[2,0.01,4.0],tres,SingleSided)")
#w.factory("Uniform::t_pdf(t)")

w.factory("PROD::pdf( m_pdf, t_pdf)")

w.defineSet("observables","t,m")

##############################################
### Define acceptance function a la Wouter ###
##############################################
effh1 = TH1F( "effh1", "effh1", nbins, w.var('t').getMin(), w.var('t').getMax()) 
for i in range(1, int(0.6 * nbins)):
    effh1.SetBinContent(i, 1. / nbins * i)
for i in range(int(0.6 * nbins), nbins + 1):
    effh1.SetBinContent(i, 1) 

## canvas = TCanvas('canvas', 'canvas', 1000, 1000)
## canvas.Divide(2, 2)
## canvas.cd(1)
## effh1.Draw()

effdatahist = RooDataHist("effdatahist", "effdatahist", RooArgList(w.var('t')), effh1) 

_import = getattr(w, 'import')
_import(effdatahist)

w.factory("HistFunc::eff(t, effdatahist)")

case = 1
if case == 1:
    w.factory("EffHistProd::acc_pdf(pdf, eff)")
else:
    w.factory("EffHistProd::taccpdf(t_pdf, eff)")
    w.factory("PROD::acc_pdf(taccpdf, m_pdf)")

w.addClassDeclImportDir("..")
w.addClassImplImportDir("..")
w.importClassCode()

acc_pdf = w.pdf('acc_pdf')

## data = acc_pdf.generate(w.set('observables'), 10000)
## acc_pdf.fitTo(data, RooFit.Minimizer('Minuit2'))

## canvas.cd(3)
## m_frame = w.var('m').frame()
## data.plotOn(m_frame)
## acc_pdf.plotOn(m_frame)
## m_frame.Draw()

## canvas.cd(4)
## t_frame = w.var('t').frame()
## data.plotOn(t_frame)
## acc_pdf.plotOn(t_frame)
## t_frame.Draw()

from ROOT import RooFit, RooGenFitStudy
from ROOT import RooStudyManager

gfs = RooGenFitStudy();
gfs.setGenConfig("acc_pdf", "m,t", RooFit.NumEvents(10000))
gfs.setFitConfig("acc_pdf", "m,t", RooFit.Minimizer("Minuit2"))

mgr = RooStudyManager(w, gfs)

## mgr.run(1)
mgr.prepareBatchInput('data2D_%d' % nbins, 100, True)

## data = gfs.summaryData();
root_file = TFile.Open('workspace2D_%d.root' % nbins, 'recreate')
root_file.WriteTObject(w, 'w')
root_file.Close()
