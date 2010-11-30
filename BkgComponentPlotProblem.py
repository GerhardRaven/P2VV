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

print 't', t.Print()
#t.setRange(0.3,12.)

ctheta = ws.var('trcostheta')
cpsi = ws.var('trcospsi')
phi = ws.var('trphi')

tagdecision = ws.var('tagdecision')

ws.defineSet("observables","t,trcostheta,trcospsi,trphi,tagdecision")

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
p2vv = ws.pdf('myJpsiphiPdf_noEff')
t_bkg = ws.pdf('t_bkg')

t.setRange('largeTime',0.3,12.)

###############################################################
### To see the fit problem for the bkg component of the pdf ###
###############################################################

lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))

sigcolor = RooFit.kGreen 
bkgcolor = RooFit.kRed

testCanvas2 = TCanvas('testCanvas2','testCanvas2') 
gPad.SetLogy()
_tb = t.frame(RooFit.Bins(100))
data2.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
t_bkg.plotOn(_tb,RooFit.Range('largeTime'),RooFit.LineColor(kRed),lw)
p2vv.plotOn(_tb,RooFit.Range('largeTime'),RooFit.LineColor(kGreen),lw)
_tb.SetMinimum(0.1) 
_tb.SetTitle("")
_tb.Draw()
testCanvas2.Update()

testCanvas = TCanvas('testCanvas','testCanvas') 
gPad.SetLogy()
_tb = t.frame(RooFit.Bins(100))
data2.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,RooFit.Range('largeTime'),lw)
pdf.plotOn(_tb,RooFit.Range('largeTime'),RooFit.Components("myJpsiphiPdf_noEff"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
pdf.plotOn(_tb,RooFit.Range('largeTime'),RooFit.Components("t_bkg"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
_tb.SetMinimum(0.1) 
_tb.SetTitle("")
_tb.Draw()
testCanvas.Update()



################################################################

