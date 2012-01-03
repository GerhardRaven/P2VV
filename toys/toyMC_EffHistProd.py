from ROOT import (RooWorkspace, TH1F,
                  gSystem, RooDataHist,
                  TFile)
gSystem.Load('libP2VV')
from RooFitDecorators import *
import rootStyle

myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

w = RooWorkspace("w")

#########################################
### Define variables and simple PDF's ###
#########################################
w.factory("t[0,10]")

#w.factory("TruthModel::tres(t)")
w.factory("GaussModel::tres(t,tres_m[0],tres_s[0.0005])")

w.factory("Decay::pdf(t,tau[2,0.01,4.0],tres,SingleSided)")
w.defineSet("observables","t")

##############################################
### Define acceptance function a la Wouter ###
##############################################
nbins = 100
interval = 10
effh1 = TH1F( "effh1", "effh1", nbins, w.var('t').getMin(), w.var('t').getMax()) 
for i in range(1, int(0.6 * nbins)):
    effh1.SetBinContent(i, 0.1 * i)
for i in range(int(0.6 * nbins), nbins + 1):
    effh1.SetBinContent(i, 1) 
effdatahist = RooDataHist("effdatahist", "effdatahist", RooArgList(w.var('t')), effh1) 

_import = getattr(w, 'import')
_import(effdatahist)

w.factory("HistFunc::eff(t, effdatahist)")
w.factory("EffHistProd::acc_pdf(pdf, eff)")
w.addClassDeclImportDir("..")
w.addClassImplImportDir("..")
w.importClassCode()

acc_pdf = w.pdf('acc_pdf')

from ROOT import RooFit, RooGenFitStudy
from ROOT import RooStudyManager

gfs = RooGenFitStudy();
gfs.storeDetailedOutput(True)
gfs.setGenConfig("acc_pdf", "t", RooFit.NumEvents(10000))
gfs.setFitConfig("acc_pdf", "t", RooFit.PrintLevel(-1), RooFit.Minimizer("Minuit2"))

mgr = RooStudyManager(w, gfs)

mgr.run(100)
## mgr.prepareBatchInput("aap", 100, kTRUE) ;

data = gfs.summaryData()
details = gfs.detailedData()
root_file = TFile.Open('data1D_%d.root' % nbins, 'recreate')
root_file.WriteTObject(data, 'data')
root_file.WriteTObject(details, 'details')
root_file.Close()
