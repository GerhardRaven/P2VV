from ROOT import *
import GaudiPython
P2VV = GaudiPython.gbl.P2VV
#to load functions (made with namespace function) like makePVVPdf:
GaudiPython.loaddict('P2VVDict')
gSystem.Load("libp2vv")
from math import pi

import RootStyle
from ROOT import (gROOT,gStyle,TStyle)
myStyle = RootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

def abasis(w,label,i,j,k,l,c):
    ctheta = w.var('trcostheta')
    cpsi = w.var('trcospsi')
    phi = w.var('trphi')
    string = "%s_%d_%d"%(label,i,j)
    if l<0:
        name = "%s_%d_m%d"%(string,k,l)
    else:
        name = "%s_%d_%d"%(string,k,l)
    #b = w.function(name)
    basisfunc = RooP2VVAngleBasis(name,name,cpsi,ctheta,phi,i,j,k,l,c)
    #im = getattr(w,'import')(basisfunc)
    #if not b:
    
    return basisfunc

def moment(data,basis,allObs,factor):
    m0 = 0
    m1 = 0
    m2 = 0
    
    for i in range(0,data.numEntries()):
        allObs.assignValueOnly(data.get(i))
        div = basis.getVal()
        #print 'div = ', div
        m0+=1
        m1+=div
        m2+=div*div

    mu = factor*m1/m0
    sig2 = m2/m0 - mu*mu
    muerr = sqrt(sig2/(m0-1))
    musigni = mu/sqrt(sig2/m0)
    print  "moment(", basis.GetName(),"_",  ") = ",mu," +- ", muerr ," significance: ", musigni
    print 'm0 =',m0
    print 'm1 =',m1
    print 'm2 =',m2
    print 'mu =',mu
    print 'merr =',muerr
    print 'musigni =',musigni

    ret = {'mu':mu,'muerr':muerr,'musigni':musigni}
    return ret

########################
### My Control Plots ###
########################

def myPlot(data,canvasname):
    
    tframe = t.frame()
    mframe = m.frame()
    mdau1frame = mdau1.frame()

    c1 = TCanvas(canvasname,canvasname,1300,900)
    c1.Divide(4,2)

    tpad = c1.cd(1)
    data.plotOn(tframe)
    #pdf.plotOn(tframe,RooFit.ProjWData(st_data))
    tpad.SetLogy()
    tframe.Draw()

    c1.cd(2)
    data.plotOn(mframe)
    #pdf.plotOn(mframe,RooFit.ProjWData(st_data))
    mframe.Draw()

    c1.cd(3)
    data.plotOn(mdau1frame)
    #pdf.plotOn(mdau1frame,RooFit.ProjWData(st_data))
    mdau1frame.Draw()

    cpsiframe = cpsi.frame(RooFit.Bins(18))
    cthetaframe = ctheta.frame(RooFit.Bins(18))
    phiframe = phi.frame(RooFit.Bins(18))

    c1.cd(4)
    data.plotOn(cpsiframe)
    cpsiframe.Draw()

    c1.cd(5)
    data.plotOn(cthetaframe)
    cthetaframe.Draw()

    c1.cd(6)
    data.plotOn(phiframe)
    phiframe.Draw()

    sigmatframe = sigmat.frame()

    c1.cd(7)
    data.plotOn(sigmatframe,RooLinkedList())
    sigmatframe.Draw()

    hist = data.createHistogram(t,m,200,200)
    hist.SetMarkerSize(0.3)
    hist.SetMarkerStyle(20)
    hist.SetStats(kFALSE)
    hist.GetXaxis().SetTitle(str(t.getTitle()))
    hist.GetYaxis().SetTitle(str(m.getTitle()))
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.SetTitle("")

    c1.cd(8)
    hist.Draw('')
        
    return {'c1':c1}

#######################
### Plot ICHEP Like ###
#######################
def plot(ws,data,pdf,title):

    c = TCanvas(title,title,900,700)
    c.Divide(3,2)

    lw = RooCmdArg(RooFit.LineWidth(2))
    xes = RooCmdArg(RooFit.XErrorSize(0))
    
    sigcolor = RooFit.kGreen 
    bkgcolor = RooFit.kRed
    
    t = ws.var("t") 
    st = ws.var("sigmat") 
    m = ws.var("m") 

    #===========================================================================================================
    c.cd(1)
    hist = data.createHistogram(t,m)
    hist.SetMarkerSize(0.3)
    hist.SetMarkerStyle(20)
    hist.SetStats(kFALSE)
    hist.GetXaxis().SetTitle(str(t.getTitle()))
    hist.GetYaxis().SetTitle(str(m.getTitle()))
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.SetTitle("")
    hist.Draw()

    #===========================================================================================================
    c.cd(2)
    gPad.SetLogy()

    _tb = t.frame(RooFit.Bins(100))
    data.plotOn(_tb,RooFit.Invisible())
    pdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),lw)
    #pdf.paramOn(_tb,RooFit.Parameters(parameterprintset))
    data.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
    _tb.SetMinimum(0.1) 
    _tb.SetTitle("")
    _tb.Draw()
    c.Update()


    #===========================================================================================================
    c.cd(3)
    _m = m.frame(RooFit.Bins(50),RooFit.Title('m'))
    data.plotOn(_m,RooFit.Invisible())
    pdf.plotOn(_m,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_m,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_m,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),lw)
    data.plotOn(_m,RooFit.MarkerSize(0.7),xes)
    _m.Draw() 
    c.Update()

    #==========================================================================================================
    c.cd(4)

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi'))
    data.plotOn(_cpsi)
    pdf.plotOn(_cpsi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cpsi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_cpsi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),lw)
    _cpsi.Draw()
    c.Update()

    #==========================================================================================================
    c.cd(5)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta'))
    data.plotOn(_ctheta)
    pdf.plotOn(_ctheta,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_ctheta,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_ctheta,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),lw)
    _ctheta.Draw()
    c.Update()

    #==========================================================================================================
    c.cd(6)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi'))
    data.plotOn(_phi)
    pdf.plotOn(_phi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_phi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_phi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),lw)
    _phi.Draw()
    c.Update()

    #===========================================================================================================
    #c.cd(7)

    #mpsiplot = mpsi.frame(RooFit.Bins(50));
    #data.plotOn(mpsiplot,RooFit.Invisible())
    #pdf.plotOn(mpsiplot,RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(mpsiplot,RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(mpsiplot,lw)
    #data.plotOn(mpsiplot,RooFit.MarkerSize(0.7),xes)
    #mpsiplot.Draw() 
    #c.Update()

    c.Update() 

    return c

#######################
### Plot Sidebands  ###
#######################
def sidebandplot(ws,data,pdf,mmin,mmax,ctitle):

    sidec = TCanvas(ctitle,ctitle,900,700)
    sidec.Divide(3,2)

    lw = RooCmdArg(RooFit.LineWidth(2))
    xes = RooCmdArg(RooFit.XErrorSize(0))
    
    sigcolor = RooFit.kGreen 
    bkgcolor = RooFit.kRed
    
    t = ws.var("t") 
    st = ws.var("sigmat") 
    m = ws.var("m") 
    #mpsi = ws.var("mdau1") 

    #create some ranges
    m.setRange("sideband",mmin,mmax)

    #===========================================================================================================
    sidec.cd(1)
    hist = data.createHistogram(t,m)
    hist.SetMarkerSize(0.3)
    hist.SetMarkerStyle(20)
    hist.SetStats(kFALSE)
    hist.GetXaxis().SetTitle(str(t.getTitle()))
    hist.GetYaxis().SetTitle(str(m.getTitle()))
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.SetTitle("")
    hist.Draw()

    #===========================================================================================================
    sidec.cd(2)
    gPad.SetLogy()

    _tb = t.frame(RooFit.Bins(10))
    data.plotOn(_tb,RooFit.CutRange('sideband'),RooFit.Invisible())
    pdf.plotOn(_tb,RooFit.ProjectionRange('sideband'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_tb,RooFit.ProjectionRange('sideband'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_tb,lw)
    pdf.paramOn(_tb,RooFit.Parameters(parameterprintset))
    data.plotOn(_tb,RooFit.CutRange('sideband'),RooFit.MarkerSize(0.5),xes)
    _tb.SetMinimum(0.1) 
    _tb.SetTitle("")
    _tb.Draw()
    sidec.Update()


    #===========================================================================================================
    sidec.cd(3)

    _m = m.frame(RooFit.Bins(10),RooFit.Title('m'))
    data.plotOn(_m,RooFit.CutRange('sideband'),RooFit.Invisible())
    pdf.plotOn(_m,RooFit.ProjectionRange('sideband'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_m,RooFit.ProjectionRange('sideband'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_m,lw)
    data.plotOn(_m,RooFit.CutRange('sideband'),RooFit.MarkerSize(0.7),xes)
    _m.Draw() 
    sidec.Update() 

    #==========================================================================================================
    sidec.cd(4)

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi'))
    data.plotOn(_cpsi,RooFit.CutRange('sideband'))
    pdf.plotOn(_cpsi,RooFit.ProjectionRange('sideband'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_cpsi,RooFit.ProjectionRange('sideband'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_cpsi,lw)
    _cpsi.Draw()
    sidec.Update()

    #==========================================================================================================
    sidec.cd(5)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta'))
    data.plotOn(_ctheta,RooFit.CutRange('sideband'))
    pdf.plotOn(_ctheta,RooFit.ProjectionRange('sideband'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_ctheta,RooFit.ProjectionRange('sideband'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_ctheta,lw)
    _ctheta.Draw()
    sidec.Update()

    #==========================================================================================================
    sidec.cd(6)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi'))
    data.plotOn(_phi,RooFit.CutRange('sideband'))
    pdf.plotOn(_phi,RooFit.ProjectionRange('sideband'),RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_phi,RooFit.ProjectionRange('sideband'),RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(_phi,lw)
    _phi.Draw()
    sidec.Update()

    #===========================================================================================================
    #sidec.cd(7)

    #mpsiplot = mpsi.frame(RooFit.Bins(50));
    #data.plotOn(mpsiplot,RooFit.Invisible())
    #pdf.plotOn(mpsiplot,RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(mpsiplot,RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    #pdf.plotOn(mpsiplot,lw)
    #data.plotOn(mpsiplot,RooFit.MarkerSize(0.7),xes)
    #mpsiplot.Draw() 
    #sidec.Update()

    sidec.Update() 

    return sidec

########################
### Plot background  ###
########################
def bkgplot(ws,bkgdata,bkgpdf,ctitle):

    bkgc = TCanvas(ctitle,ctitle,900,700)
    bkgc.Divide(3,2)

    lw = RooCmdArg(RooFit.LineWidth(2))
    xes = RooCmdArg(RooFit.XErrorSize(0))
    
    sigcolor = RooFit.kGreen 
    bkgcolor = RooFit.kRed
    
    t = ws.var("t") 
    st = ws.var("sigmat") 
    m = ws.var("m") 
 
    #===========================================================================================================
    bkgc.cd(1)
    hist = bkgdata.createHistogram(t,m)
    hist.SetMarkerSize(0.3)
    hist.SetMarkerStyle(20)
    hist.SetStats(kFALSE)
    hist.GetXaxis().SetTitle(str(t.getTitle()))
    hist.GetYaxis().SetTitle(str(m.getTitle()))
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.SetTitle("")
    hist.Draw()

    #===========================================================================================================
    bkgc.cd(2)
    gPad.SetLogy()

    _tb = t.frame(RooFit.Bins(10))
    bkgdata.plotOn(_tb,RooFit.Invisible())
    bkgpdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    bkgpdf.paramOn(_tb,RooFit.Parameters(parameterprintset))
    bkgdata.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
    _tb.SetMinimum(0.1) 
    _tb.SetTitle("")
    _tb.Draw()
    bkgc.Update()


    #===========================================================================================================
    bkgc.cd(3)

    _m = m.frame(RooFit.Bins(10),RooFit.Title('m'))
    bkgdata.plotOn(_m,RooFit.Invisible())
    bkgpdf.plotOn(_m,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    bkgdata.plotOn(_m,RooFit.MarkerSize(0.7),xes)
    _m.Draw() 
    bkgc.Update() 

    #==========================================================================================================
    bkgc.cd(4)

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi'))
    bkgdata.plotOn(_cpsi)
    bkgpdf.plotOn(_cpsi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    _cpsi.Draw()
    bkgc.Update()

    #==========================================================================================================
    bkgc.cd(5)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta'))
    bkgdata.plotOn(_ctheta)
    bkgpdf.plotOn(_ctheta,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    _ctheta.Draw()
    bkgc.Update()

    #==========================================================================================================
    bkgc.cd(6)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi'))
    bkgdata.plotOn(_phi)
    bkgpdf.plotOn(_phi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    _phi.Draw()
    bkgc.Update()

    return bkgc

########################
### Plot signal      ###
########################
def sigplot(ws,sigdata,sigpdf,ctitle):

    sigc = TCanvas(ctitle,ctitle,900,700)
    sigc.Divide(3,2)

    lw = RooCmdArg(RooFit.LineWidth(2))
    xes = RooCmdArg(RooFit.XErrorSize(0))
    
    sigcolor = RooFit.kGreen 
    bkgcolor = RooFit.kRed
    
    t = ws.var("t") 
    st = ws.var("sigmat") 
    m = ws.var("m") 
 
    #===========================================================================================================
    sigc.cd(1)
    hist = sigdata.createHistogram(t,m)
    hist.SetMarkerSize(0.3)
    hist.SetMarkerStyle(20)
    hist.SetStats(kFALSE)
    hist.GetXaxis().SetTitle(str(t.getTitle()))
    hist.GetYaxis().SetTitle(str(m.getTitle()))
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.SetTitle("")
    hist.Draw()

    #===========================================================================================================
    sigc.cd(2)
    gPad.SetLogy()

    _tb = t.frame(RooFit.Bins(50))
    sigdata.plotOn(_tb,RooFit.Invisible())
    sigpdf.plotOn(_tb,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    sigpdf.paramOn(_tb,RooFit.Parameters(parameterprintset))
    sigdata.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
    _tb.SetMinimum(0.1) 
    _tb.SetTitle("")
    _tb.Draw()
    sigc.Update()


    #===========================================================================================================
    sigc.cd(3)

    _m = m.frame(RooFit.Bins(50),RooFit.Title('m'))
    sigdata.plotOn(_m,RooFit.Invisible())
    sigpdf.plotOn(_m,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    sigdata.plotOn(_m,RooFit.MarkerSize(0.7),xes)
    _m.Draw() 
    sigc.Update() 

    #==========================================================================================================
    sigc.cd(4)

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi'))
    sigdata.plotOn(_cpsi)
    sigpdf.plotOn(_cpsi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    _cpsi.Draw()
    sigc.Update()

    #==========================================================================================================
    sigc.cd(5)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta'))
    sigdata.plotOn(_ctheta)
    sigpdf.plotOn(_ctheta,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    _ctheta.Draw()
    sigc.Update()

    #==========================================================================================================
    sigc.cd(6)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi'))
    sigdata.plotOn(_phi)
    sigpdf.plotOn(_phi,RooFit.NormRange('largeTime'),RooFit.Range('largeTime'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    _phi.Draw()
    sigc.Update()

    return sigc
##############################################################################
##########################   There we go!!!!! ################################
##############################################################################

decay = 'Bs2JpsiPhi'
#selection = 'Bs2JpsiPhiPrescaled'
selection = 'Bs2JpsiPhiUnbiased'
#selection = 'Bs2JpsiPhiDetached'

#decay = 'Bd2JpsiKstar'
#selection = 'Bd2JpsiKstarPrescaled'
#selection = 'Bd2JpsiKstarUnbiased'
#selection = 'Bd2JpsiKstarDetached'
##############################
### Create ws, observables ###
##############################

wsfile = TFile('/stage/Daan/Beta_s/P2VV/Development/GaanWeWeer/WSJpsiPhiPdf.root')

ws = RooWorkspace(wsfile.Get('workspace'))

t = ws.var('t')

cpsi = ws.var('trcospsi')
ctheta = ws.var('trcostheta')
phi = ws.var('trphi')

tagdecision = ws.cat('tagdecision')

m = RooRealVar("m","B mass",5200,5550)

sigmat = RooRealVar("sigmat","B proper time error",0.000001,0.22)

#mdau1 = RooRealVar("mdau1","dau1 mass",3097-60,3097+60)
#mdau2 = RooRealVar("mdau2","dau2 mass",0,2000)

obs = RooArgSet(m)

obs.add(sigmat)

getattr(ws,'import')(obs)

#This is needed for sPlots only
ws.defineSet("observables","m,t,sigmat,trcostheta,trcospsi,trphi,tagdecision")

ws.var("sigmat").setBins(40) 

#################
### Load Data ###
#################
file = TFile('duitsedata.root')
NTupletree = file.Get('Bs2JpsiPhi')

#file = TFile('/stage/Daan/Beta_s/ReaderWriter/Grid/AllStr_10_11/NTuple_'+selection+'_FromGrid.root')
#NTupletree = file.Get('mDSTReaderAlgo/dataset')

data = RooDataSet('data','data',NTupletree,ws.set('observables'),"t==t && sigmat==sigmat && m==m")

data.table(tagdecision).Print('v')

data2 = data.reduce('t>0.3 && t<12')
data2.SetNameTitle('data2','data2')

sigdata = data2.reduce('m>5320 && m<5402')
sigdata.SetNameTitle('sigdata','sigdata')

#MyTFile = TFile.Open('%s_RDS.root'%(selection),'RECREATE')
#data.Write('RDS')
#MyTFile.Close()

#tree = data.tree()

getattr(ws,'import')(data)
getattr(ws,'import')(data2)

# make a background dataset for description of angular distributions and to fit the background distribution
lowbkgdata = data2.reduce('m>5402')
lowbkgdata.SetNameTitle('lowbkgdata','lowbkgdata')
highbkgdata = data2.reduce('m<5320')
highbkgdata.SetNameTitle('highbkgdata','highbkgdata')
bkgdata = lowbkgdata
bkgdata.append(highbkgdata)
bkgdata.SetNameTitle('bkgdata','bkgdata')
#repeat to restore lowbkgdata, WHOO, THAT'S UGLY!!!!!!!!
lowbkgdata = data2.reduce('m>5402')
lowbkgdata.SetNameTitle('lowbkgdata','lowbkgdata')
lowbkgdata.Print()
highbkgdata.Print()
bkgdata.Print()

# Calculate background angular distributions

mu_vars = []
coeffList = RooArgList('coeffList')

muerr_vars = []
coefferrList = RooArgList('coefferrList')

musigni_vars = []
coeffsigniList = RooArgList('coeffsigniList')

basislist = []

allObs = RooArgSet(cpsi,ctheta,phi)

for i in range(0,4):
    for l in range(0,4):
        for m in range(-l,l+1):
            basis = abasis(ws,'mom',i,0,l,m,float(2*i+1)/2. )
            factor = 2./float(2*i+1)

            basislist.append(basis)
            
            momentdict = moment(bkgdata,basis,allObs,factor)
            
            mu_vars.append(RooRealVar('C_%s_0_%s_%s'%(i,l,m),'C_%s_0_%s_%s'%(i,l,m),momentdict['mu']))
            coeffList.add(mu_vars[-1])
            
            muerr_vars.append(RooRealVar('Cerr_%s_0_%s_%s'%(i,l,m),'Cerr_%s_0_%s_%s'%(i,l,m),momentdict['muerr']))
            coefferrList.add(muerr_vars[-1])
            
            musigni_vars.append(RooRealVar('Csigni_%s_0_%s_%s'%(i,l,m),'Csigni_%s_0_%s_%s'%(i,l,m),momentdict['musigni']))
            coeffsigniList.add(musigni_vars[-1])

basisList = RooArgList('basisList')
for i in basislist:
    basisList.add(i)

basisList.Print()
coeffList.Print()

#bkgang = RooRealSumPdf("bkgang","bkgang",basisList,coeffList)
#getattr(ws,'import')(bkgang)
#BKGC = bkgplot(ws,bkgdata,bkgang,'BKGC')

#########################################################################
### Build the PDF's (can only do that after we have sigmat from data) ###
#########################################################################

#define some constants
ws.factory("const_zero[0]") 
ws.factory("const_one[1]") 

#signal B mass pdf
if decay == 'Bd2JpsiKstar':
    ws.factory("Gaussian::m_sig(m,m_sig_mean[5280,5260,5300],m_sig_sigma[10,3,30])")
if decay == 'Bs2JpsiPhi':
    ws.factory("Gaussian::m_sig_1(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,10.])")
    ws.factory("Gaussian::m_sig_2(m,m_sig_mean,m_sig_sigma_2[18.,10.,20.])")
    ws.factory("SUM::m_sig(m_sig_f_1[0.85,0,1]*m_sig_1,m_sig_2)")
    
#Getting the JpsiPhi signal Pdf
p2vv = ws.pdf('myJpsiphiPdf_noEff')
#p2vv = ws.pdf('myJpsiphiPdf_withWeights')


#Getting the resolution model from the JpsiPhi pdf in the workspace
tres = ws.pdf('tres')

#background B mass pdf
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

#background propertime
#Replace the resolution model with the overall resolution model tres
#ws.factory("Decay::t_bkg_ml(t,t_bkg_ml_tau[0.16,0.05,0.3],tres,SingleSided)")
#ws.factory("Decay::t_bkg_ll(t,t_bkg_ll_tau[1.2,1.0,2.5],tres,SingleSided)")

#Make a rootruthmodel:
#truthres = RooTruthModel("truthres","truthres",t)
#getattr(ws,'import')(truthres)
#ws.factory("Decay::t_bkg_ml(t,t_bkg_ml_tau[0.165,0.05,1.0],truthres,SingleSided)")
#ws.factory("Decay::t_bkg_ll(t,t_bkg_ll_tau[1.14,1.,3.],truthres,SingleSided)")

ws.factory("{t_bkg_ml_tau[0.207,0.1,0.5],t_bkg_ll_tau[1.92,1.,2.5]}")
ws.factory("FormulaVar::cml('-1/@0',{t_bkg_ml_tau})")
ws.factory("FormulaVar::cll('-1/@0',{t_bkg_ll_tau})")

ws.factory("Exponential::t_bkg_ml(t,cml)")
ws.factory("Exponential::t_bkg_ll(t,cll)")

ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*t_bkg_ll,t_bkg_ml)")

#background angles: for now assume flat (in t>0.3 ps fit), we know that it is very non-flat and different for Jpsi and non-Jpsi bkg!
ws.factory("Uniform::bkgang({trcostheta,trcospsi,trphi})")

ws.factory("Chebychev::bkg_trcospsi(trcospsi,{c0_trcospsi[-1,1],c1_trcospsi[-1,1],c2_trcospsi[-1,1],c3_trcospsi[-1,1],c4_trcospsi[-1,1]})")#,c5_trcospsi[-1,1],c6_trcospsi[-1,1]})")
ws.factory("Chebychev::bkg_trcostheta(trcostheta,{c0_trcostheta[-1,1],c1_trcostheta[-1,1],c2_trcostheta[-1,1],c3_trcostheta[-1,1],c4_trcostheta[-1,1]})")#,c5_trcostheta[-1,1],c6_trcostheta[-1,1]})")
ws.factory("Chebychev::bkg_trphi(trphi,{c0_trphi[-1,1],c1_trphi[-1,1],c2_trphi[-1,1],c3_trphi[-1,1],c4_trphi[-1,1]})")#,c5_trphi[-1,1],c6_trphi[-1,1]})")

#ws.factory("PROD::bkgang(bkg_trcostheta,bkg_trcospsi,bkg_trphi)")

#define sigmat pdf
#getattr(ws,'import')(RooDataHist("data_sigt","hist Err Per Ev",RooArgSet(ws.var("sigmat")),ws.data("data")))
#ws.factory("HistPdf::pdf_stE(sigmat,data_sigt)")

#P2VV fit
ws.factory("PROD::sig_pdf( m_sig, myJpsiphiPdf_noEff)")
#ws.factory("PROD::sig_pdf( m_sig, myJpsiphiPdf_withWeights)")

ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, bkgang)")
ws.factory("SUM::pdf_ext(Nsig[543,0,1000]*sig_pdf,Nbkg[812,0,1000]*bkg_pdf)")
ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sig_pdf,bkg_pdf)")
#ws.factory("SUM::pdf(f_sig[0.71]*sig_pdf,bkg_pdf)")

#also build the full mass-only pdf
ws.factory("SUM::masspdf(Nsig*m_sig,Nbkg*m_bkg)") 

#ws.defineSet( "yields","Nsig,Nprompt,Nbkg")

#ws.defineSet( "psionlyyields","Nprompt,Nbkg") 

#########################
### What do you want? ###
#########################

#myPlotDict = myPlot(data,'originaldata')

#t.setRange(0.3,12)

#myPlotDict2 = myPlot(data2,'cutdata')

lifetime = 1.5
if ( decay == 'B2JpsiKplus' ):#B+
    BMass  = 5279.17
elif ( decay == 'Bd2JpsiKstar' ):#B0 
    BMass  = 5279.17
elif ( decay == 'Bs2JpsiPhi' ):#Bs
    BMass  = 5366.3

ws.var("m_sig_mean").setRange(BMass-20,BMass+20 ) 
ws.var("m_sig_mean").setVal( BMass ) 

#get pdf's from the workspace
bkg_pdf = ws.pdf('bkg_pdf')
sig_pdf = ws.pdf('sig_pdf')
pdf= ws.pdf("pdf")#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pdf_ext = ws.pdf('pdf_ext')

#nll = pdf_ext.createNLL(data2,RooFit.Extended(true)) 
#print "nll value at starting point: ", nll.getVal()



pars = pdf.getParameters(data2) 
pars.writeToFile( "fittedparameters.txt") 

#Parameters to be printed

Aperp_sq_jpsiphi         = ws.var('Aperp_sq_jpsiphi')
Azero_sq_jpsiphi         = ws.var('Azero_sq_jpsiphi')
G_s                      = ws.var('G_s')
Nbkg                     = ws.var('Nbkg')
Nsig                     = ws.var('Nsig')
f_sig                    = ws.var('f_sig')
Phi_s                    = ws.var('Phi_s')
dG_s                     = ws.var('dG_s')
dm_s                     = ws.var('dm_s')
#m_bkg_exp                = ws.var('m_bkg_exp')
m_sig_mean               = ws.var('m_sig_mean')
m_sig_sigma_1            = ws.var('m_sig_sigma_1')
m_sig_sigma_2            = ws.var('m_sig_sigma_2')
phi_par_jpsiphi          = ws.var('phi_par_jpsiphi')
phi_perp_jpsiphi         = ws.var('phi_perp_jpsiphi')
phi_zero_jpsiphi         = ws.var('phi_zero_jpsiphi')
t_bkg_fll                = ws.var('t_bkg_fll')
t_bkg_ll_tau             = ws.var('t_bkg_ll_tau')
t_bkg_ml_tau             = ws.var('t_bkg_ml_tau')

parameterprintset = RooArgSet(Aperp_sq_jpsiphi)
parameterprintset.add(Azero_sq_jpsiphi)
parameterprintset.add(Nbkg)
parameterprintset.add(Nsig)
parameterprintset.add(f_sig)
#parameterprintset.add(dm_s)
#parameterprintset.add(m_bkg_exp)
parameterprintset.add(m_sig_mean)
parameterprintset.add(m_sig_sigma_1)
parameterprintset.add(m_sig_sigma_2)
parameterprintset.add(phi_par_jpsiphi)
parameterprintset.add(phi_perp_jpsiphi)
parameterprintset.add(phi_zero_jpsiphi)
parameterprintset.add(t_bkg_fll)
parameterprintset.add(t_bkg_ll_tau)
parameterprintset.add(t_bkg_ml_tau)

t.setRange(0.3,12.)
t.setRange('largeTime',0.3,12.)

m= ws.var('m')
masspdf = ws.pdf('masspdf')
masspdf.fitTo(data2,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))

masscanvas = TCanvas('masscanvas','masscanvas')
massframe = m.frame()
data2.plotOn(massframe,RooLinkedList())
masspdf.plotOn(massframe,RooFit.Range('largeTime'),RooFit.NormRange('largeTime'))
masspdf.paramOn(massframe)
massframe.Draw()
###########################################################################################################################
#FINALLY, BKG IN T AND M AND UNIFORM ANGLES IS COOL!!!!
bkg_pdf.fitTo(bkgdata,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Minos(false),RooFit.Save(true))
bkgc = bkgplot(ws,bkgdata,bkg_pdf,'bkgcanvas')
###########################################################################################################################

###########################################################################################################################
#FINALLY, SIG IN T AND M AND ANGLES IS COOL!!!!
#Fit and plot signal
sig_pdf.fitTo(sigdata,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Minos(false),RooFit.Save(true))
sigc = sigplot(ws,sigdata,sig_pdf,'sigcanvas')
###########################################################################################################################

#result = pdf.fitTo(data2,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Extended(false),RooFit.Minos(false),RooFit.Save(true))
#c = plot(ws,data2,pdf,'c')

##########
result = pdf_ext.fitTo(data2,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
cext = plot(ws,data2,pdf_ext,'cext')
##########

## nll = RooNLLVar('nll','nll',pdf_ext,data2,RooFit.Extended(true))
## #nll = RooNLLVar('nll','nll',pdf,data2,RooFit.Extended(false))

## LLCanvas = TCanvas('LLCanvas','LLCanvas')
## LLCanvas.Divide(3,3)

## LLCanvas.cd(1)
## Nbkgframe = Nbkg.frame(RooFit.Title("-log(L) scan vs Nbkg")) 
## nll.plotOn(Nbkgframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## Nbkgframe.Draw()

## LLCanvas.cd(2)
## Nsigframe = Nsig.frame(RooFit.Title("-log(L) scan vs Nsig")) 
## nll.plotOn(Nsigframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## Nsigframe.Draw()

## LLCanvas.cd(3)
## f_sigframe = f_sig.frame(RooFit.Title("-log(L) scan vs Nsig")) 
## nll.plotOn(f_sigframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## f_sigframe.Draw()

## LLCanvas.cd(4)
## dG_sframe = dG_s.frame(RooFit.Title("-log(L) scan vs dG_s")) 
## nll.plotOn(dG_sframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## dG_sframe.Draw()

## LLCanvas.cd(5)
## G_sframe = G_s.frame(RooFit.Title("-log(L) scan vs dG_s")) 
## nll.plotOn(G_sframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## G_sframe.Draw()

## LLCanvas.cd(6)
## Aperp_sq_jpsiphiframe = Aperp_sq_jpsiphi.frame(RooFit.Title("-log(L) scan vs dG_s")) 
## nll.plotOn(Aperp_sq_jpsiphiframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## Aperp_sq_jpsiphiframe.Draw()

## LLCanvas.cd(7)
## Azero_sq_jpsiphiframe = Azero_sq_jpsiphi.frame(RooFit.Title("-log(L) scan vs dG_s")) 
## nll.plotOn(Azero_sq_jpsiphiframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## Azero_sq_jpsiphiframe.Draw()

## LLCanvas.cd(8)
## m_sig_meanframe = m_sig_mean.frame(RooFit.Title("-log(L) scan vs m")) 
## nll.plotOn(m_sig_meanframe,RooFit.PrintEvalErrors(0),RooFit.ShiftToZero(),RooFit.LineColor(kRed),RooFit.Precision(1e-4),RooFit.EvalErrorValue(nll.getVal()+10))
## m_sig_meanframe.Draw()
