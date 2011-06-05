#######################################
### Author: Daan van Eijk
### Updated on: Jun 5 11
### Description: This script reads the signal JpsiPhi Pdf in transversity angles from a rootfile
###              Adds bkg PDF's, read data and fits for transversity amplitudes, to be compared later with
###              transversity amplitudes fitted with helicity angles in DvEHelFit.py.
########################################
from ROOT import *
gSystem.Load("libp2vv")
from math import pi

import rootStyle
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()


#######################
### Plot ICHEP Like ###
#######################
def plot(ws,data,pdf,title):

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

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi full mass region'))
    data.plotOn(_cpsi)
    pdf.plotOn(_cpsi,RooFit.Components('sig_pdf'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cpsi,RooFit.Components('bkg_pdf'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cpsi,lw)
    _cpsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(2)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta full mass region'))
    data.plotOn(_ctheta)
    pdf.plotOn(_ctheta,RooFit.Components('sig_pdf'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_ctheta,RooFit.Components('bkg_pdf'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_ctheta,lw)
    _ctheta.Draw()
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

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi signal region'))
    data.plotOn(_cpsi,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_cpsi,RooFit.ProjectionRange('sigRegion'),lw)
    _cpsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(5)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta signal region'))
    data.plotOn(_ctheta,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_ctheta,RooFit.ProjectionRange('sigRegion'),lw)
    _ctheta.Draw()
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

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi left sideband'))
    data.plotOn(_cpsi,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_cpsi,RooFit.ProjectionRange('leftSideband'),lw)
    _cpsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(8)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta left sideband'))
    data.plotOn(_ctheta,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_ctheta,RooFit.ProjectionRange('leftSideband'),lw)
    _ctheta.Draw()
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

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi right sideband'))
    data.plotOn(_cpsi,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_cpsi,RooFit.ProjectionRange('rightSideband'),lw)
    _cpsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(11)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta right sideband'))
    data.plotOn(_ctheta,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_ctheta,RooFit.ProjectionRange('rightSideband'),lw)
    _ctheta.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(12)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi right sideband'))
    data.plotOn(_phi,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_phi,RooFit.ProjectionRange('rightSideband'),lw)
    _phi.Draw()
    c2.Update()

    return myline1,myline2,c,c2

##############################################################################
##########################   There we go!!!!! ################################
##############################################################################

##############################
### Read ws, observables   ###
##############################

wsfile = TFile('TransWS.root')

ws = wsfile.Get('workspace')

t = ws.var('t')

##############################
#This sets the range for the events in the dataset that you read, and the range in which you will fit!
##############################
tmin = 0.3
tmax = 12.
t.setRange(tmin,tmax)

cpsi = ws.var('trcospsi')
ctheta = ws.var('trcostheta')
phi = ws.var('trphi')

tagdecision = ws.cat('tagdecision')

m = RooRealVar("m","B_{s} mass",5200,5550)

obs = RooArgSet(m)

getattr(ws,'import')(obs)

ws.defineSet("observables","m,t,trcostheta,trcospsi,trphi,tagdecision")

#################
### Load Data ###
#################

#Using 'common' dataset:
#file = TFile('duitsedata2.root')
#NTupletree = file.Get('Bs2JpsiPhi')

#Using Wouters dataset:
file = TFile('Bs2JpsiPhiTupleReduced.root')
NTupletree = file.Get('dataset')

data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m && trcostheta==trcostheta && trcospsi==trcospsi && trphi==trphi && tagdecision==tagdecision')
data.table(tagdecision).Print('v')

getattr(ws,'import')(data)

#######################
### Build the PDF's ###
#######################

ws.factory("Gaussian::m_sig(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")
    
#Getting the JpsiPhi signal Pdf
p2vv = ws.pdf('newpdf')
#p2vv = ws.pdf('myJpsiphiPdf_withWeights')

#Getting the resolution model from the JpsiPhi pdf in the workspace
res = ws.pdf('tres_sig')

#background B mass pdf
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

#background propertime 
# TODO: split resolution in tres_sig and tres_nonpsi!!
ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.1,0.5],tres_sig,SingleSided)")
ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,1.,2.5],tres_sig,SingleSided)")
ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")

#background angles: 
#ws.factory("Uniform::bkgang({trcostheta,trcospsi,trphi})")

ws.factory("Chebychev::bkg_trcospsi(trcospsi,{c0_trcospsi[-0.13,-1,1],c1_trcospsi[0.23,-1,1],c2_trcospsi[-0.057,-1,1],c3_trcospsi[-0.0058,-1,1],c4_trcospsi[-0.0154,-1,1],c5_trcospsi[-1,1],c6_trcospsi[-1,1]})")
ws.factory("Chebychev::bkg_trcostheta(trcostheta,{c0_trcostheta[0.08,-1,1],c1_trcostheta[-0.22,-1,1],c2_trcostheta[-0.022,-1,1],c3_trcostheta[0.21,-1,1],c4_trcostheta[0.0125,-1,1],c5_trcostheta[-1,1],c6_trcostheta[-1,1]})")
ws.factory("Chebychev::bkg_trphi(trphi,{c0_trphi[0.10,-1,1],c1_trphi[0.328,-1,1],c2_trphi[0.081,-1,1],c3_trphi[0.316,-1,1],c4_trphi[0.044,-1,1],c5_trphi[-1,1],c6_trphi[-1,1]})")

#ws.factory("Chebychev::bkg_trcospsi(trcospsi,{c0_trcospsi[-0.13],c1_trcospsi[0.23],c2_trcospsi[-0.057],c3_trcospsi[-0.0058],c4_trcospsi[-0.0154]})")#,c5_trcospsi[-1,1],c6_trcospsi[-1,1]})")
#ws.factory("Chebychev::bkg_trcostheta(trcostheta,{c0_trcostheta[0.08],c1_trcostheta[-0.22],c2_trcostheta[-0.022],c3_trcostheta[0.21],c4_trcostheta[0.0125]})")#,c5_trcostheta[-1,1],c6_trcostheta[-1,1]})")
#ws.factory("Chebychev::bkg_trphi(trphi,{c0_trphi[0.10],c1_trphi[0.328],c2_trphi[0.081],c3_trphi[0.316],c4_trphi[0.044]})")#,c5_trphi[-1,1],c6_trphi[-1,1]})")

ws.factory("PROD::bkgang(bkg_trcostheta,bkg_trcospsi,bkg_trphi)")

#P2VV fit
ws.factory("PROD::sig_pdf( m_sig, newpdf)")

ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, bkgang)")

ws.factory("SUM::pdf_ext(Nsig[543,0,1000]*sig_pdf,Nbkg[812,0,1000]*bkg_pdf)")
ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sig_pdf,bkg_pdf)")

#also build the full mass-only pdf
ws.factory("SUM::masspdf(Nsig*m_sig,Nbkg*m_bkg)")

#########################
### What do you want? ###
#########################

#get pdf's from the workspace
bkg_pdf = ws.pdf('bkg_pdf')
sig_pdf = ws.pdf('sig_pdf')
pdf= ws.pdf("pdf")
pdf_ext = ws.pdf('pdf_ext')
masspdf = ws.pdf('masspdf')

#Parameters to be printed

## rperp2                   = ws.var('rperp2')
## rz2                      = ws.var('rz2')
## gamma                    = ws.var('gamma')
## Nbkg                     = ws.var('Nbkg')
## Nsig                     = ws.var('Nsig')
## f_sig                    = ws.var('f_sig')
## phis                     = ws.var('phis')
## dG                       = ws.var('dG')
## dm                       = ws.var('dm')
## m_bkg_exp                = ws.var('m_bkg_exp')
## m_sig_mean               = ws.var('m_sig_mean')
## m_sig_sigma_1            = ws.var('m_sig_sigma_1')
## deltapar                 = ws.var('deltapar')
## deltaperp                = ws.var('deltaperp')
## deltaz                   = ws.var('deltaz')
## t_bkg_fll                = ws.var('t_bkg_fll')
## t_bkg_ll_tau             = ws.var('t_bkg_ll_tau')
## t_bkg_ml_tau             = ws.var('t_bkg_ml_tau')

## parameterprintset = RooArgSet(rperp2)
## parameterprintset.add(rz2)
## parameterprintset.add(Nbkg)
## parameterprintset.add(Nsig)
## #parameterprintset.add(f_sig)
## #parameterprintset.add(dm_s)
## parameterprintset.add(m_bkg_exp)
## parameterprintset.add(m_sig_mean)
## parameterprintset.add(m_sig_sigma_1)
## parameterprintset.add(deltapar)
## parameterprintset.add(deltaperp)
## parameterprintset.add(deltaz)
## parameterprintset.add(t_bkg_fll)
## parameterprintset.add(t_bkg_ll_tau)
## parameterprintset.add(t_bkg_ml_tau)

####################
#MASS ONLY
## m= ws.var('m')
## masspdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))

## masscanvas = TCanvas('masscanvas','masscanvas')
## massframe = m.frame(RooFit.Bins(20))
## data.plotOn(massframe,RooLinkedList())

## masspdf.plotOn(massframe)
## masspdf.plotOn(massframe,RooFit.Components('m_sig'),RooFit.LineColor(kGreen),RooFit.LineStyle(kDashed))
## masspdf.plotOn(massframe,RooFit.Components('m_bkg'),RooFit.LineColor(kRed),RooFit.LineStyle(kDashed))

## masspdf.paramOn(massframe)
## massframe.Draw()
####################

####################
# FIT
#result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(false),RooFit.Minos(false),RooFit.Save(true))
#c = plot(ws,data,pdf,'c')
####################

####################
# EXTENDED FIT
result = pdf_ext.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
cext = plot(ws,data,pdf_ext,'cext')
####################

####################
# Latex code of the fitted parameters
paramlist = pdf_ext.getParameters(data)
paramlist.printLatex(RooFit.Format("NEU",RooFit.AutoPrecision(2),RooFit.VerbatimName()))
####################
