from ROOT import *
#to load functions (made with namespace function) like makePVVPdf:
gSystem.Load("libp2vv")
from math import pi

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
    pdf.paramOn(_m,RooFit.Parameters(RooArgSet(Nbkg,Nsig)))
    pdf.paramOn(_m,RooFit.Parameters(numberprintset))    
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

##############################
### Create ws, observables ###
##############################

wsfile = TFile('WSJpsiPhiPdf.root')

#ws = RooWorkspace(wsfile.Get('workspace'))
ws = wsfile.Get('workspace')

t = ws.var('t')

cpsi = ws.var('trcospsi')
ctheta = ws.var('trcostheta')
phi = ws.var('trphi')

tagdecision = ws.cat('tagdecision')

m = RooRealVar("m","B mass",5200,5550)

obs = RooArgSet(m)

getattr(ws,'import')(obs)

#This is needed for sPlots only
ws.defineSet("observables","m,t,trcostheta,trcospsi,trphi,tagdecision")

#################
### Load Data ###
#################
file = TFile('duitsedata.root')
NTupletree = file.Get('Bs2JpsiPhi')

realdata = RooDataSet('realdata','realdata',NTupletree,ws.set('observables'),"t==t && m==m")

realdata.table(tagdecision).Print('v')

realdata2 = realdata.reduce('t>0.3 && t<12')
realdata2.SetNameTitle('realdata2','realdata2')

sigrealdata = realdata2.reduce('m>5320 && m<5402')
sigrealdata.SetNameTitle('sigrealdata','sigrealdata')

getattr(ws,'import')(realdata)
getattr(ws,'import')(realdata2)

# make a background dataset for description of angular distributions and to fit the background distribution
lowbkgrealdata = realdata2.reduce('m>5402')
lowbkgrealdata.SetNameTitle('lowbkgrealdata','lowbkgrealdata')
highbkgrealdata = realdata2.reduce('m<5320')
highbkgrealdata.SetNameTitle('highbkgrealdata','highbkgrealdata')
bkgrealdata = lowbkgrealdata
bkgrealdata.append(highbkgrealdata)
bkgrealdata.SetNameTitle('bkgrealdata','bkgrealdata')
#repeat to restore lowbkgdata, WHOO, THAT'S UGLY!!!!!!!!
lowbkgrealdata = realdata2.reduce('m>5402')
lowbkgrealdata.SetNameTitle('lowbkgrealdata','lowbkgrealdata')
lowbkgrealdata.Print()
highbkgrealdata.Print()
bkgrealdata.Print()

#########################################################################
### Build the PDF's (can only do that after we have sigmat from data) ###
#########################################################################

ws.factory("Gaussian::m_sig_1(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,10.])")
ws.factory("Gaussian::m_sig_2(m,m_sig_mean,m_sig_sigma_2[18.,10.,20.])")
ws.factory("SUM::m_sig(m_sig_f_1[0.85,0,1]*m_sig_1,m_sig_2)")
    
#Getting the JpsiPhi signal Pdf
p2vv = ws.pdf('jpsiphipdf')
#p2vv = ws.pdf('myJpsiphiPdf_withWeights')

#Getting the resolution model from the JpsiPhi pdf in the workspace
tres = ws.pdf('res')

#background B mass pdf
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

#background propertime

ws.factory("{t_bkg_ml_tau[0.207,0.1,0.5],t_bkg_ll_tau[1.92,1.,2.5]}")
ws.factory("FormulaVar::cml('-1/@0',{t_bkg_ml_tau})")
ws.factory("FormulaVar::cll('-1/@0',{t_bkg_ll_tau})")

ws.factory("Exponential::t_bkg_ml(t,cml)")
ws.factory("Exponential::t_bkg_ll(t,cll)")

ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*t_bkg_ll,t_bkg_ml)")

#background angles: 
ws.factory("Uniform::bkgang({trcostheta,trcospsi,trphi})")

ws.factory("Chebychev::bkg_trcospsi(trcospsi,{c0_trcospsi[-1,1],c1_trcospsi[-1,1],c2_trcospsi[-1,1],c3_trcospsi[-1,1],c4_trcospsi[-1,1]})")#,c5_trcospsi[-1,1],c6_trcospsi[-1,1]})")
ws.factory("Chebychev::bkg_trcostheta(trcostheta,{c0_trcostheta[-1,1],c1_trcostheta[-1,1],c2_trcostheta[-1,1],c3_trcostheta[-1,1],c4_trcostheta[-1,1]})")#,c5_trcostheta[-1,1],c6_trcostheta[-1,1]})")
ws.factory("Chebychev::bkg_trphi(trphi,{c0_trphi[-1,1],c1_trphi[-1,1],c2_trphi[-1,1],c3_trphi[-1,1],c4_trphi[-1,1]})")#,c5_trphi[-1,1],c6_trphi[-1,1]})")

#ws.factory("PROD::bkgang(bkg_trcostheta,bkg_trcospsi,bkg_trphi)")

#P2VV fit
ws.factory("PROD::sig_pdf( m_sig, jpsiphipdf)")
#ws.factory("PROD::sig_pdf( m_sig, myJpsiphiPdf_withWeights)")

ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, bkgang)")
ws.factory("SUM::pdf_ext(Nsig[543,0,1000]*sig_pdf,Nbkg[812,0,1000]*bkg_pdf)")
ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sig_pdf,bkg_pdf)")
#ws.factory("SUM::pdf(f_sig[0.71]*sig_pdf,bkg_pdf)")

#also build the full mass-only pdf
ws.factory("SUM::masspdf(Nsig*m_sig,Nbkg*m_bkg)") 

#########################
### What do you want? ###
#########################

BMass  = 5366.3

ws.var("m_sig_mean").setRange(BMass-20,BMass+20 ) 
ws.var("m_sig_mean").setVal( BMass ) 

#get pdf's from the workspace
bkg_pdf = ws.pdf('bkg_pdf')
sig_pdf = ws.pdf('sig_pdf')
pdf= ws.pdf("pdf")
pdf_ext = ws.pdf('pdf_ext')

#Parameters to be printed

rperp                    = ws.var('rperp')
rz                       = ws.var('rz')
gamma                    = ws.var('gamma')
dG                       = ws.var('dG')
dm                       = ws.var('dm')
Nbkg                     = ws.var('Nbkg')
Nsig                     = ws.var('Nsig')
f_sig                    = ws.var('f_sig')
#m_bkg_exp                = ws.var('m_bkg_exp')
m_sig_mean               = ws.var('m_sig_mean')
m_sig_sigma_1            = ws.var('m_sig_sigma_1')
m_sig_sigma_2            = ws.var('m_sig_sigma_2')
t_bkg_fll                = ws.var('t_bkg_fll')
t_bkg_ll_tau             = ws.var('t_bkg_ll_tau')
t_bkg_ml_tau             = ws.var('t_bkg_ml_tau')

parameterprintset = RooArgSet(rperp)
parameterprintset.add(rz)
#parameterprintset.add(gamma)#blinded
#parameterprintset.add(dG)#blinded
#parameterprintset.add(dm)#constant
parameterprintset.add(Nbkg)
parameterprintset.add(Nsig)
parameterprintset.add(f_sig)
#parameterprintset.add(m_bkg_exp)
parameterprintset.add(m_sig_mean)
parameterprintset.add(m_sig_sigma_1)
parameterprintset.add(m_sig_sigma_2)
parameterprintset.add(t_bkg_fll)
parameterprintset.add(t_bkg_ll_tau)
parameterprintset.add(t_bkg_ml_tau)

numberprintset = RooArgSet(Nbkg,Nsig)

t.setRange(0.3,12.)
t.setRange('largeTime',0.3,12.)

m= ws.var('m')
masspdf = ws.pdf('masspdf')
masspdf.fitTo(realdata2,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))

masscanvas = TCanvas('masscanvas','masscanvas')
massframe = m.frame()
realdata2.plotOn(massframe,RooLinkedList())
masspdf.plotOn(massframe,RooFit.Range('largeTime'),RooFit.NormRange('largeTime'))
masspdf.paramOn(massframe)
massframe.Draw()
###########################################################################################################################
#FINALLY, BKG IN T AND M AND UNIFORM ANGLES IS COOL!!!!
#bkg_pdf.fitTo(bkgrealdata,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Minos(false),RooFit.Save(true))
#bkgc = bkgplot(ws,bkgrealdata,bkg_pdf,'bkgcanvas')
###########################################################################################################################

###########################################################################################################################
#FINALLY, SIG IN T AND M AND ANGLES IS COOL!!!!
#Fit and plot signal
#sig_pdf.fitTo(sigrealdata,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Minos(false),RooFit.Save(true))
#sigc = sigplot(ws,sigrealdata,sig_pdf,'sigcanvas')
###########################################################################################################################
##########
result = pdf.fitTo(realdata2,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Extended(false),RooFit.Minos(false),RooFit.Save(true))
c = plot(ws,realdata2,pdf,'c')
##########

##########
#result = pdf_ext.fitTo(realdata2,RooFit.Range('largeTime'),RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
#cext = plot(ws,realdata2,pdf_ext,'cext')
##########
