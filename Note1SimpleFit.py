from ROOT import *
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

    return myline1,myline2,c

##############################################################################
##########################   There we go!!!!! ################################
##############################################################################

########################
### Define workspace ###
########################

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ t[-1,13.]}")
t= ws.var('t')

new = True

if new:
    ws.factory("{m[5200,5550]}")
else :
    ws.factory("{m[5250,5450]}")

m = ws.var('m')

##########################
### physics parameters ###
##########################

ws.factory("{#tau[1.47,1.,2.]}")

ws.factory("RooUnblindUniform::tau_unblind('BsCalvin',0.2,#tau)")

#tau = 1.47 defines gamma:
#ws.factory("{gamma[0.68,0.4,0.9]}")
#ws.factory("expr::tau('1/@0',{gamma})")

##############################
#This sets the range for the events in the dataset that you read, and the range in which you will fit!
##############################
tmin = -1.
tmax = 13.

t.setRange(tmin,tmax)
ws.defineSet("observables","m,t")

if new:
    file = TFile('Bs2JpsiPhiStripping12_JpsiMassCut.root')
    NTupletree = file.Get('MyTree')
else :
    file = TFile('Bs2JpsiPhiTuple.root')
    NTupletree = file.Get('dataset')

data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m ')

print 'There are ', data.numEntries(), 'entries in the dataset'

getattr(ws,'import')(data)

#######################
### Build the PDF's ###
#######################

### Resolution ###

ws.factory("RooGaussModel::tres_1(t,#mu[0],#sigma_1[0.03,0.01,0.1],1)")
ws.factory("RooGaussModel::tres_2(t,#mu,#sigma_2[0.07,0.04,0.1],1)")
ws.factory("RooGaussModel::tres_3(t,#mu,#sigma_3[0.2,0.1,1.],1)")

ws.factory("AddModel::tres({tres_3,tres_2,tres_1},{tres_f3[0.008,0.,0.1],tres_f2[0.3,0.2,0.5]})")

### Signal ###

# t #
# Blinded
ws.factory("RooDecay::t_sig(t,tau_unblind,tres,SingleSided)")

# Unblinded
#ws.factory("RooDecay::t_sig(t,tau,tres,SingleSided)")

# m #
ws.factory("Gaussian::m_sig(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")

### Background ###

# t #
ws.factory("Decay::prompt(t,0,tres,SingleSided)")
ws.factory("RooDecay::ml(t,#tau_ml[0.15,0.1,0.2],tres,SingleSided)")
if new:
    ws.factory("RooDecay::ll(t,#tau_ll[1.4,0.7,1.7],tres,SingleSided)")
else :
    ws.factory("RooDecay::ll(t,#tau_ll[1.4,0.7,2.0],tres,SingleSided)")
ws.factory("SUM::t_bkg(f_ll[0.03,0.,0.1]*ll,f_ml[0.001,0.,0.1]*ml,prompt)")

# m #
ws.factory("Chebychev::m_bkg(m,{c0_m[-0.03,-0.1,-0.01]})")

#full PDF
ws.factory("PROD::sig_pdf( m_sig, t_sig)")

ws.factory("PROD::bkg_pdf( m_bkg, t_bkg)")

if new:
    ws.factory("SUM::pdf_ext(Nsig[700,500,1000]*sig_pdf,Nbkg[40000,35000,45000]*bkg_pdf)")
else :
    ws.factory("SUM::pdf_ext(Nsig[1000,900,1500]*sig_pdf,Nbkg[67000,60000,70000]*bkg_pdf)")

ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sig_pdf,bkg_pdf)")

# Make mass pdf to determine yields

ws.factory("SUM::m_pdf(Nsig*m_sig,Nbkg*m_bkg)")

#########################
### What do you want? ###
#########################

#get pdf's from the workspace
bkg_pdf = ws.pdf('bkg_pdf')
sig_pdf = ws.pdf('sig_pdf')
pdf= ws.pdf("pdf")
pdf_ext = ws.pdf('pdf_ext')
m_pdf = ws.pdf('m_pdf')

######################
### MASS FIT #########
######################
m_pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
c = TCanvas('c','c')
mframe = m.frame()
data.plotOn(mframe,RooLinkedList())
m_pdf.plotOn(mframe)
mframe.Draw()                   

####################
# EXTENDED FIT
result = pdf_ext.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
cext = plot(ws,data,pdf_ext,'cext')
####################

####################
# FIT
#result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Minos(false),RooFit.Save(true))
#cext = plot(ws,data,pdf,'cext')
####################

####################
# Latex code of the fitted parameters
paramlist = pdf_ext.getParameters(data)
paramlist.printLatex(RooFit.Format("NEU",RooFit.AutoPrecision(3),RooFit.VerbatimName()))
####################

