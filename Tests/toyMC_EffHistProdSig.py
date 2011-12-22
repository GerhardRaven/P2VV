import optparse
import sys
parser = optparse.OptionParser(usage = '%prog nbins')

parser.add_option("-f", "--file", dest = "file", default = 'sig.root',
                  type = 'string', help = "set input filename")
parser.add_option("-w", "--workspace", dest = "workspace", default = 'w',
                  type = 'string', help = "set input filename")
parser.add_option("-p", "--pdf", dest = "pdf", default = 'mc_pdf',
                  type = 'string', help = "set pdf name")
parser.add_option("-o", "--observables", dest = "observables", type = 'string',
                  default = 'trphi,trcostheta,trcospsi,t,tagdecision')

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
from RooFitDecorators import *

root_file = TFile.Open(options.file)
if not root_file:
    print 'cannot open root file %s' % options.file
    sys.exit(-2)

w = root_file.Get(options.workspace)
if not w:
    print 'cannot get workspace %s' % options.workspace
    sys.exit(-2)

##############################################
### Define acceptance function a la Wouter ###
##############################################
effh1 = TH1F( "effh1", "effh1", nbins, w.var('t').getMin(), w.var('t').getMax()) 
for i in range(1, int(0.6 * nbins)):
    effh1.SetBinContent(i, 1. / nbins * i)
for i in range(int(0.6 * nbins), nbins + 1):
    effh1.SetBinContent(i, 1) 

effdatahist = RooDataHist("effdatahist", "effdatahist", RooArgList(w.var('t')), effh1) 

_import = getattr(w, 'import')
_import(effdatahist)

w.factory("HistFunc::eff(t, effdatahist)")

w.factory("EffHistProd::acc_pdf(%s, eff)" % options.pdf )

w.addClassDeclImportDir("..")
w.addClassImplImportDir("..")
w.importClassCode()

acc_pdf = w.pdf('acc_pdf')

from ROOT import RooFit, RooGenFitStudy
from ROOT import RooStudyManager

gfs = RooGenFitStudy();
gfs.setGenConfig("acc_pdf", options.observables, RooFit.NumEvents(1000))
gfs.setFitConfig("acc_pdf", options.observables, RooFit.Minimizer("Minuit2"))

mgr = RooStudyManager(w, gfs)

## mgr.run(1)
#mgr.prepareBatchInput('dataSig_%d' % nbins, 100, True)

## data = gfs.summaryData();
## root_file = TFile.Open('workspaceSig_%d.root' % nbins, 'recreate')
## root_file.WriteTObject(w, 'w')
## root_file.Close()
