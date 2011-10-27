from ROOT import *
gSystem.Load("libP2VV")

from RooFitDecorators import *

##############################
### Create ws, observables ###
##############################

ws = RooWorkspace('ws')

ws.factory("t[0.3,14.]")

##############################
### Proper Time Acceptance ###
##############################
ws.factory("expr::effshape('1/(1+(a*t)**(-c))',t,a[1.45],c[2.37])")
effhist = ws['effshape'].createHistogram('effhist',ws['t'],RooFit.Binning(200,ws['t'].getMin(),ws['t'].getMax()))
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effhist)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,effdatahist)")

#############
### PDF's ###
#############
#Resolution Model
ws.factory("GaussModel::tres(t,tres_mean[0.0],tres_sigma[0.05])")

ws.factory("Decay::t_sig(t,tau[1.4],tres,SingleSided)")

#BKG time
#Single no ERROR, Double does have ERROR!
#ws.factory("RooDecay::t_bkg(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")

ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")
ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres,SingleSided)")
ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")

ws.factory("SUM::pdf_ext( Nsig[1186]*t_sig,Nbkg[568]*t_bkg)")

ws.factory("EffHistProd::acc_sig_pdf(t_sig,effpdf)")
ws.factory("EffHistProd::acc_bkg_pdf(t_bkg,effpdf)")
ws.factory("SUM::accpdf_ext( Nsig*acc_sig_pdf,Nbkg*acc_bkg_pdf)")

#########################
### What do you want? ###
#########################
acc= True

if acc:
    pdf = ws.pdf('accpdf_ext')
    signame = 'acc_sig_pdf'
    bkgname = 'acc_bkg_pdf'
else:
    pdf = ws.pdf('pdf_ext')
    signame = 't_sig'
    bkgname = 't_bkg'

data = pdf.generate(RooArgSet(ws.var('t')),0)

###############################################################
### To see the plot problem for the bkg component of the pdf ###
###############################################################

lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))

sigcolor = RooFit.kGreen 
bkgcolor = RooFit.kRed

testCanvas = TCanvas('testCanvas','testCanvas') 
gPad.SetLogy()
_tb = ws.var('t').frame(RooFit.Bins(100))
data.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,lw)
pdf.plotOn(_tb,RooFit.Components(signame),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
pdf.plotOn(_tb,RooFit.Components(bkgname),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
_tb.Draw()
testCanvas.Update()
################################################################

