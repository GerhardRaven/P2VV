from ROOT import *
import GaudiPython
P2VV = GaudiPython.gbl.P2VV
#to load functions (made with namespace function) like makePVVPdf:
GaudiPython.loaddict('P2VVDict')
gSystem.Load("libp2vv")
from math import pi

#######################
### Plot ICHEP Like ###
#######################

rootobjects = []

def plot(ws,data,pdf,title):

    c = TCanvas('MassTime_' + title,'MassTime',900,700)
    rootobjects.append(c)
    c.Divide(4,2)

    
    lw = RooCmdArg(RooFit.LineWidth(2))
    xes = RooCmdArg(RooFit.XErrorSize(0))
    err = RooCmdArg(RooFit.DrawOption('E'))
    
    sigcolor = RooFit.kGreen 
    bkgcolor = RooFit.kRed
    nonpsicolor = RooFit.kOrange
    
    t = ws.var("t") 
    m = ws.var("m")
    
    ws.var("sigmat").setBins(40)
    stdata = RooDataHist("data_sigt","hist Err Per Ev",RooArgSet(ws.var("sigmat")),ws.data("data"))
    projst =  RooCmdArg(RooFit.ProjWData(stdata))

    tmin = t.getMin()
    tmax = t.getMax()
    msigmin = 5345.
    msigmax = 5387.
    mmin = 5200.
    mmax = 5550.
    m.setRange('sigRegion',msigmin,msigmax)
    m.setRange('leftSideband',mmin,msigmin)
    m.setRange('rightSideband',msigmax,mmax)
    t.setRange('largeTime',0.3,tmax)
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
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.SetTitle("")
    hist.Draw()
    myline1.Draw('same')
    myline2.Draw('same')
    c.Update()

    #===========================================================================================================
    c.cd(2)
    _m = m.frame(RooFit.Bins(50),RooFit.Title('m'))
    data.plotOn(_m,RooFit.MarkerSize(0.7),xes)
    pdf.plotOn(_m,RooFit.Components("bkg_pdf"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw,projst)
    pdf.plotOn(_m,RooFit.Components("sig_pdf"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw,projst)
    pdf.plotOn(_m,lw,projst)
    _m.Draw() 
    c.Update()

    #===========================================================================================================
    c.cd(3)
    _m = m.frame(RooFit.Bins(50),RooFit.Title('m (t>0.3)'))
    data.plotOn(_m,RooFit.MarkerSize(0.7),xes,RooFit.CutRange("largeTime"))
    pdf.plotOn(_m,RooFit.Components("bkg_pdf"),RooFit.ProjectionRange("largeTime"),
               RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw,projst)
    pdf.plotOn(_m,RooFit.Components("sig_pdf"),RooFit.ProjectionRange("largeTime"),
               RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw,projst)
    pdf.plotOn(_m,lw,RooFit.ProjectionRange("largeTime"),projst)
    _m.Draw() 
    c.Update()

    #===========================================================================================================
    c.cd(5)
    _tb = t.frame(-0.4,0.4,100)
    data.plotOn(_tb,RooFit.MarkerSize(0.5),xes,err)
    pdf.plotOn(_tb,RooFit.Components("sig_pdf"),projst,RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,RooFit.Components("bkg_pdf"),projst,RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,RooFit.Components("nonpsi_pdf"),projst,RooFit.LineColor(nonpsicolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,lw,projst)

    #pdf.paramOn(_tb,RooFit.Parameters(parameterprintset))

    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time full mass range")
    _tb.Draw()
    c.Update()
    
    #===========================================================================================================
    c.cd(6)
    gPad.SetLogy()

    _tb = t.frame(-1,12.,50)
    data.plotOn(_tb,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,RooFit.Components("bkg_pdf"),RooFit.ProjectionRange('sigRegion'),
               projst,RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,lw,RooFit.ProjectionRange('sigRegion'),projst)

    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time signal region")
    _tb.Draw()
    c.Update()


    #===========================================================================================================
    c.cd(7)
    gPad.SetLogy()

    _tb = t.frame(-1,12,50)
    data.plotOn(_tb,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,RooFit.Components("nonpsi_pdf"),RooFit.ProjectionRange('leftSideband'),
               projst,RooFit.LineColor(nonpsicolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,lw,RooFit.ProjectionRange('leftSideband'),projst)
    
    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time left sideband")
    _tb.Draw()
    c.Update()
    #===========================================================================================================
    c.cd(8)
    gPad.SetLogy()

    _tb = t.frame(-1,12,50)
    data.plotOn(_tb,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,RooFit.Components("nonpsi_pdf"),RooFit.ProjectionRange('rightSideband'),
               projst,RooFit.LineColor(nonpsicolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,lw,RooFit.ProjectionRange('rightSideband'),projst)

    _tb.SetMinimum(0.1) 
    _tb.SetTitle("proper time right sideband")
    _tb.Draw()
    c.Update()
    
    #===========================================================================================================
    c.cd(4)
    mpsiplot = ws.var('mdau1').frame(RooFit.Bins(50))
    data.plotOn(mpsiplot,RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(mpsiplot,RooFit.Components("nonpsi_pdf"),projst,RooFit.LineColor(nonpsicolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(mpsiplot,lw,projst)
    mpsiplot.Draw()
    c.Update()

    c.Print("MassTime.eps")
    return
    
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

    _cpsi = cpsi.frame(RooFit.Bins(15),RooFit.Title('cpsi sidebands'))
    data.plotOn(_cpsi,RooFit.CutRange('leftSideband,rightSideband'))
    pdf.plotOn(_cpsi,RooFit.ProjectionRange('leftSideband,rightSideband'),lw)
    _cpsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(8)

    _ctheta = ctheta.frame(RooFit.Bins(15),RooFit.Title('ctheta sidebands'))
    data.plotOn(_ctheta,RooFit.CutRange('leftSideband,rightSideband'))
    pdf.plotOn(_ctheta,RooFit.ProjectionRange('leftSideband,rightSideband'),lw)
    _ctheta.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(9)

    _phi = phi.frame(RooFit.Bins(15),RooFit.Title('phi sidebands'))
    data.plotOn(_phi,RooFit.CutRange('leftSideband,rightSideband'))
    pdf.plotOn(_phi,RooFit.ProjectionRange('leftSideband,rightSideband'),lw)
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

##################################
### Create WS, build the PDF's ###
##################################

from ModelBuilders import *

ws = RooWorkspace('ws')

buildFullJpsiPhiPdf ( ws )

print 'Ready building pdf'

# first test: dump the pdf parameters to a file

pdf_ext = ws.pdf('pdf_ext')
pdf_ext.getParameters( ws.set('observables')).writeToFile( 'initialfitparameters.txt' )


##############################
#This sets the range for the events in the dataset that you read, and the range in which you will fit!
##############################

#t.setRange(0.3,12.)

#ws.defineSet("observables","m,t,trcostheta,trcospsi,trphi,tagdecision")

#################
### Load Data ###
#################

#file = TFile('duitsedata.root')
#NTupletree = file.Get('Bs2JpsiPhi')

file = TFile('/data/bfys/wouterh/JpsiXAnalysis/Bs2JpsiPhiTupleReduced.root')
NTupletree = file.Get('dataset')

data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m==m')
data.table(ws.cat('tagdecision')).Print('v')

getattr(ws,'import')(data)


# if you want the NL fit to behave reasonable, 

readParameters(ws,'JPsiPhiParameters.txt')

#plot(ws,data,pdf_ext,'joho')

# we just want the fit over the lifetimes and masses
allobservables = ws.set('observables')
allobservables.remove(ws.var('t'))
allobservables.remove(ws.var('sigmat'))
allobservables.remove(ws.var('m'))
allobservables.remove(ws.var('mdau1'))
allobservables.Print()
masstimepdf = pdf_ext.createProjection( allobservables )

setConstant(ws,"tcp*")
setConstant(ws,"tct*")
setConstant(ws,"tp*")
setConstant(ws,"delta*")
setConstant(ws,"r*")
setConstant(ws,"dG")
setConstant(ws,"dm")
setConstant(ws,"gamma")
setConstant(ws,"phis")

pdf_ext.getParameters( ws.set('observables')).writeToFile( 'initialfitparameters.txt' )

#result = pdf_ext.fitTo(data,RooFit.ConditionalObservables(ws.set("conditionalobservables")),
#                       RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(false),
#                       RooFit.Save(true),RooFit.Verbose(false))
#
#pdf_ext.getParameters( ws.set('observables')).writeToFile( 'fittedparameters.txt' )

plot(ws,data,pdf_ext,'joho')

