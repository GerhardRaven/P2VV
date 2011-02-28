# DvE Fits the signal PDF only to MC data, EvtGen only, to see if we understand the .dec file.

from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi

import rootStyle
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()


#######################
### Plot ICHEP Like ###
#######################
def plothel(ws,data,pdf,title):

    c = TCanvas('MassTime','MassTime',900,700)
    c.Divide(3,2)

    
    lw = RooCmdArg(RooFit.LineWidth(2))
    xes = RooCmdArg(RooFit.XErrorSize(0))
    err = RooCmdArg(RooFit.DrawOption('E'))
    
    sigcolor = RooFit.kGreen 
    bkgcolor = RooFit.kRed
    
    t = ws.var("t") 
    m = ws.var("m") 

    msigmin = 5345.
    msigmax = 5387.
    mmin = 5200.
    mmax = 5550
    m.setRange('sigRegion',msigmin,msigmax)
    m.setRange('leftSideband',mmin,msigmin)
    m.setRange('rightSideband',msigmax,mmax)

    #===========================================================================================================
    c.cd(1)
    myline1=TLine(tmin,msigmin,tmax,msigmin)
    myline1.SetLineColor(3)
    myline1.SetLineWidth(2)
    myline1.SetLineStyle(1)
    myline2=TLine(tmin,msigmax,tmax,msigmax)
    myline2.SetLineColor(3)
    myline2.SetLineWidth(2)
    myline2.SetLineStyle(1)

    hist = data.createHistogram(t,m)
    hist.SetMarkerSize(0.3)
    hist.SetMarkerStyle(20)
    hist.SetStats(kFALSE)
    hist.GetXaxis().SetTitle(str(t.getTitle()))
    hist.GetYaxis().SetTitle(str(m.getTitle()))
    hist.SetTitle('B_{s} mass vs. proper time')
    hist.Draw()
    myline1.Draw('same')
    myline2.Draw('same')
    c.Update()

    #===========================================================================================================
    c.cd(2)
    _m = m.frame(RooFit.Bins(40),RooFit.Title('B_{s} mass, t>%s ps'%(str(tmin))))
    data.plotOn(_m,RooFit.MarkerSize(0.7),xes)
    pdf.plotOn(_m,RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_m,RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_m,lw)
    _m.Draw() 
    c.Update()

    #===========================================================================================================
    c.cd(3)
    gPad.SetLogy()

    _tb = t.frame(-1,13,50)
    data.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,lw)

    #pdf.paramOn(_tb,RooFit.Parameters(parameterprintset))

    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time full mass range")
    _tb.Draw()
    c.Update()
    
    #===========================================================================================================
    c.cd(4)
    gPad.SetLogy()

    _tb = t.frame(-1,13,50)
    data.plotOn(_tb,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,lw,RooFit.ProjectionRange('sigRegion'))

    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time signal region")
    _tb.Draw()
    c.Update()


    #===========================================================================================================
    c.cd(5)
    gPad.SetLogy()

    _tb = t.frame(-1,13,50)
    data.plotOn(_tb,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,lw,RooFit.ProjectionRange('leftSideband'))
    
    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time left sideband")
    _tb.Draw()
    c.Update()
    #===========================================================================================================
    c.cd(6)
    gPad.SetLogy()

    _tb = t.frame(-1,13,50)
    data.plotOn(_tb,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,lw,RooFit.ProjectionRange('rightSideband'))

    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time right sideband")
    _tb.Draw()
    c.Update()
    
    ####################
    ####################

    c2 = TCanvas('Angles','Angles',900,700)
    c2.Divide(3,4)
    #==========================================================================================================
    c2.cd(1)

    _cthetaK = cthetaK.frame(RooFit.Bins(15),RooFit.Title('cthetaK full mass region'))
    data.plotOn(_cthetaK)
    pdf.plotOn(_cthetaK,RooFit.Components('sig_pdf'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cthetaK,RooFit.Components('bkg_pdf'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cthetaK,lw)
    _cthetaK.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(2)

    _cthetaL = cthetaL.frame(RooFit.Bins(15),RooFit.Title('cthetaL full mass region'))
    data.plotOn(_cthetaL)
    pdf.plotOn(_cthetaL,RooFit.Components('sig_pdf'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cthetaL,RooFit.Components('bkg_pdf'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cthetaL,lw)
    _cthetaL.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(3)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi full mass region'))
    data.plotOn(_phi)
    pdf.plotOn(_phi,RooFit.Components('sig_pdf'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_phi,RooFit.Components('bkg_pdf'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_phi,lw)
    _phi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(4)

    _cthetaK = cthetaK.frame(RooFit.Bins(15),RooFit.Title('cthetaK signal region'))
    data.plotOn(_cthetaK,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_cthetaK,RooFit.ProjectionRange('sigRegion'),lw)
    _cthetaK.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(5)

    _cthetaL = cthetaL.frame(RooFit.Bins(15),RooFit.Title('cthetaL signal region'))
    data.plotOn(_cthetaL,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_cthetaL,RooFit.ProjectionRange('sigRegion'),lw)
    _cthetaL.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(6)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi signal region'))
    data.plotOn(_phi,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_phi,RooFit.ProjectionRange('sigRegion'),lw)
    _phi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(7)

    _cthetaK = cthetaK.frame(RooFit.Bins(15),RooFit.Title('cthetaK left sideband'))
    data.plotOn(_cthetaK,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_cthetaK,RooFit.ProjectionRange('leftSideband'),lw)
    _cthetaK.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(8)

    _cthetaL = cthetaL.frame(RooFit.Bins(15),RooFit.Title('cthetaL left sideband'))
    data.plotOn(_cthetaL,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_cthetaL,RooFit.ProjectionRange('leftSideband'),lw)
    _cthetaL.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(9)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi left sideband'))
    data.plotOn(_phi,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_phi,RooFit.ProjectionRange('leftSideband'),lw)
    _phi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(10)

    _cthetaK = cthetaK.frame(RooFit.Bins(15),RooFit.Title('cthetaK right sideband'))
    data.plotOn(_cthetaK,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_cthetaK,RooFit.ProjectionRange('rightSideband'),lw)
    _cthetaK.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(11)

    _cthetaL = cthetaL.frame(RooFit.Bins(15),RooFit.Title('cthetaL right sideband'))
    data.plotOn(_cthetaL,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_cthetaL,RooFit.ProjectionRange('rightSideband'),lw)
    _cthetaL.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(12)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi right sideband'))
    data.plotOn(_phi,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_phi,RooFit.ProjectionRange('rightSideband'),lw)
    _phi.Draw()
    c2.Update()

    return myline1,myline2,c,c2

def plottrans(ws,data,pdf,title):

    c = TCanvas('MassTime','MassTime',900,300)
    c.Divide(4,1)
    
    lw = RooCmdArg(RooFit.LineWidth(2))
    xes = RooCmdArg(RooFit.XErrorSize(0))
    err = RooCmdArg(RooFit.DrawOption('E'))
    
    sigcolor = RooFit.kGreen 
    bkgcolor = RooFit.kRed
    
    t = ws.var("t") 

    #===========================================================================================================
    c.cd(1)
    gPad.SetLogy()

    _tb = t.frame(-1,13,50)
    data.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,lw)

    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time full mass range")
    _tb.Draw()
    c.Update()
    #==========================================================================================================
    c.cd(2)

    _trcospsi = trcospsi.frame(RooFit.Title('trcospsi full mass region'))
    data.plotOn(_trcospsi)
    pdf.plotOn(_trcospsi,lw)
    _trcospsi.Draw()
    c.Update()

    #==========================================================================================================
    c.cd(3)

    _trcostheta = trcostheta.frame(RooFit.Title('trcostheta full mass region'))
    data.plotOn(_trcostheta)
    pdf.plotOn(_trcostheta,lw)
    _trcostheta.Draw()
    c.Update()

    #==========================================================================================================
    c.cd(4)

    _trphi = trphi.frame(RooFit.Title('trphi full mass region'))
    data.plotOn(_trphi)
    pdf.plotOn(_trphi,lw)
    _trphi.Draw()
    c.Update()

    return c

##############################################################################
##########################   There we go!!!!! ################################
##############################################################################
from ModelBuilders import *

ws = RooWorkspace("ws")

useTransversityAngles = True

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))

ws.factory("{ helcosthetaK[-1,1], helcosthetaL[-1,1], helphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))

if useTransversityAngles:
    ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
    angles = ws.set('transversityangles')
else:
    ws.defineSet("helicityangles","helcosthetaK,helcosthetaL,helphi")
    angles = ws.set('helicityangles')

ab = abasis(ws,angles)

#Fill in the values that went into the MC!
#MC:

#rpar = 0.49
#deltapar = 2.50
#rz = 0.775
#deltaz = 0.
#rperp = 0.4
#deltarperp = -0.17

#dG = 0.06852
#Beta_s = 0.02

#This yields
#rperp2 = 0.16
#rz2 = 0.601

ws.factory("{rz2[0.601,0.,1.],rperp2[0.16,0.,1.]")
ws.factory("RooFormulaVar::rpar2('1-@0-@1',{rz2,rperp2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")

ws.factory("{deltaz[0.,%f,%f],deltapar[2.5,%f,%f],deltaperp[-0.17,%f,%f]}"%(-pi,pi,-pi,pi,-pi,pi))

##########################
### physics parameters ###
##########################
#tau = 1.47 defines gamma:

ws.factory("{#Gamma[0.68,0.4,0.9]}")
#ws.factory("{#Gamma[0.68]}")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")

ws.factory("{t_sig_dG[0.06852,-0.3,0.3]}")
#ws.factory("{t_sig_dG[0.06852]}")

ws.factory("{t_sig_dm[17.8,15.,20.]}")
#ws.factory("{t_sig_dm[20.]}")

ws.factory('{phis[-0.04,-0.1,0.1]}')
#ws.factory('{phis[-0.04]}')

ws.factory("{expr::S('-1*sin(phis)',{phis}),expr::D('cos(phis)',{phis}),C[0]}")
ws.factory("{expr::Sold('sin(phis)',{phis}),expr::Dold('cos(phis)',{phis}),Cold[0]}")

###############################
### Experimental parameters ###
###############################

ws.factory("RooTruthModel::tres_sig(t)")

ws.factory("{wtag[0.0]}")

##################################### building the NEW pdf ###################################
# create observables
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

#This one imports automatically in the workspace!
newpdf = buildJpsiphi(ws,'newpdf', useTransversityAngles) 

##############################
#This sets the range for the events in the dataset that you read, and the range in which you will fit!
##############################
t = ws.var('t')

if useTransversityAngles:
    trcostheta = ws.var('trcostheta')
    trcospsi = ws.var('trcospsi')
    trphi = ws.var('trphi')
else:
    cthetaK = ws.var('helcosthetaK')
    cthetaL = ws.var('helcosthetaL')
    phi = ws.var('helphi')

tagdecision = ws.cat('tagdecision')

obs = RooArgSet(t)

getattr(ws,'import')(obs)

if useTransversityAngles:
    ws.defineSet("observables","t,trcospsi,trcostheta,trphi,tagdecision")
else:
    ws.defineSet("observables","t,helcosthetaL,helcosthetaK,helphi,tagdecision")

###################
### Load MC set ###
###################

file = TFile('/tmp/dvaneijk/MC_13144001_1M_EvtGenOnly.root')
#file = TFile('/tmp/dvaneijk/MC_13144001_100k_EvtGenOnly.root')
NTupletree = file.Get('JpsiPhi')

if useTransversityAngles:
    data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && trcospsi==trcospsi && trcostheta==trcostheta && trphi==trphi && tagdecision==tagdecision')
else:
    data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && helcosthetaL==helcosthetaL && helcosthetaK==helcosthetaK && helphi==helphi && tagdecision==tagdecision')
    
data.table(tagdecision).Print('v')

getattr(ws,'import')(data)

#Getting the JpsiPhi signal Pdf
p2vv = ws.pdf('newpdf')

result = p2vv.fitTo(data,RooFit.NumCPU(8),RooFit.Minos(false),RooFit.Save(true))

if useTransversityAngles:
    #plot = plottrans(ws,data,pdf_ext,'transplot')                              
    plot = plottrans(ws,data,p2vv,'transplot')
else:
    #plot = plothel(ws,data,pdf_ext,'transplot')                                
    plot = plothel(ws,data,p2vv,'helplot')

####################
# Latex code of the fitted parameters
paramlist = p2vv.getParameters(data)
paramlist.printLatex(RooFit.Format("NEU",RooFit.AutoPrecision(3),RooFit.VerbatimName()))
####################
