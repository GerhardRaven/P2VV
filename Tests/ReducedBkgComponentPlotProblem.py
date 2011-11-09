from ROOT import *
gSystem.Load("libP2VV")

from RooFitDecorators import *

#############
### Flags ###
#############
acc= True
#Bkg is double exponential: test = True: (eff*ml+eff*ll), test = False: eff*(ml+ll)
test=False

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

#SigPDF
ws.factory("Decay::t_sig(t,tau[1.4],tres,SingleSided)")

#Biased SigPDF
ws.factory("EffHistProd::t_sig_B(t_sig,effpdf)")

#BKG time
if test:
    ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")
    ws.factory("EffHistProd::ml_B(ml,effpdf)")
    ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres,SingleSided)")
    ws.factory("EffHistProd::ll_B(ll,effpdf)")
    #BkgPdf
    ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")
    #Biased BkgPdf
    ws.factory("SUM::t_bkg_B(t_bkg_fll[0.3,0.,1.]*ll_B,ml_B)")
else:
    ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")
    ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres,SingleSided)")
    #BkgPdf
    ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")
    #Biased BkgPDF
    ws.factory("EffHistProd::t_bkg_B(t_bkg,effpdf)")

ws.factory("SUM::pdf_ext( Nsig[1186]*t_sig,Nbkg[568]*t_bkg)")

ws.factory("SUM::pdf_ext_B( Nsig*t_sig_B,Nbkg*t_bkg_B)")

#########################
### What do you want? ###
#########################

if acc:
    pdf = ws.pdf('pdf_ext_B')
    signame = 't_sig_B'
    bkgname = 't_bkg_B'
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

