from ROOT import *
import GaudiPython
P2VV = GaudiPython.gbl.P2VV
#to load functions (made with namespace function) like makePVVPdf:
GaudiPython.loaddict('P2VVDict')

##############################
### Create ws, observables ###
##############################

wsfile = TFile('WSJpsiPhiPdf.root')

ws = RooWorkspace(wsfile.Get('workspace'))

t = ws.var('t')

m = RooRealVar("m","B mass",5200,5550)

obs = RooArgSet(m)

getattr(ws,'import')(obs)

ws.defineSet("observables","t")

#################
### Load Data ###
#################
file = TFile('duitsedata.root')
NTupletree = file.Get('Bs2JpsiPhi')

data = RooDataSet('data','data',NTupletree,ws.set('observables'),"t==t")

data2 = data.reduce('t>0.3 && t<12')
data2.SetNameTitle('data2','data2')

getattr(ws,'import')(data)
getattr(ws,'import')(data2)

#########################################################################
### Build the PDF's (can only do that after we have sigmat from data) ###
#########################################################################
#Getting the JpsiPhi signal Pdf
p2vv = ws.pdf('myJpsiphiPdf_noEff')

#background t pdf
ws.factory("{t_bkg_ml_tau[0.207,0.1,0.5],t_bkg_ll_tau[1.92,1.,2.5]}")
ws.factory("FormulaVar::cml('-1/@0',{t_bkg_ml_tau})")
ws.factory("FormulaVar::cll('-1/@0',{t_bkg_ll_tau})")

ws.factory("Exponential::t_bkg_ml(t,cml)")
ws.factory("Exponential::t_bkg_ll(t,cll)")

ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*t_bkg_ll,t_bkg_ml)")

ws.factory("SUM::pdf(f_sig[0.6,0,1]*myJpsiphiPdf_noEff,t_bkg)")

#########################
### What do you want? ###
#########################

pdf= ws.pdf("pdf")

t.setRange(0.3,12.)
t.setRange('largeTime',0.3,12.)

###############################################################
### To see the fit problem for the bkg component of the pdf ###
###############################################################

lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))

sigcolor = RooFit.kGreen 
bkgcolor = RooFit.kRed

testCanvas = TCanvas('testCanvas','testCanvas') 
gPad.SetLogy()
_tb = t.frame(RooFit.Bins(100))
data2.plotOn(_tb,RooFit.Invisible())
pdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),lw)
pdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("myJpsiphiPdf_noEff"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
pdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("t_bkg"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
data2.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
_tb.SetMinimum(0.1) 
_tb.SetTitle("")
_tb.Draw()
testCanvas.Update()
################################################################

