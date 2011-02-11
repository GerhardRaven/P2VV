# This is the implementation of the Note2 fit. Simple resolution model, flat angular background.
# Including angular acceptance correction can be set with the acceptance flag.
# Angular basis can be chosen using the useTransversityAngles flag.

from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi

import RooFitDecorators
import rootStyle
from ModelBuilders import _buildAngularFunction
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
    pdf.plotOn(_m,RooFit.Components("bkg_pdf_eff"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_m,RooFit.Components("sig_pdf_eff"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_m,lw)
    _m.Draw() 
    c.Update()

    #===========================================================================================================
    c.cd(3)
    gPad.SetLogy()

    _tb = t.frame(-1,13,50)
    data.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
    pdf.plotOn(_tb,RooFit.Components("sig_pdf_eff"),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_tb,RooFit.Components("bkg_pdf_eff"),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
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

    _trcospsi = trcospsi.frame(RooFit.Bins(15),RooFit.Title('trcospsi full mass region'))
    data.plotOn(_trcospsi)
    pdf.plotOn(_trcospsi,RooFit.Components('sig_pdf_eff'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_trcospsi,RooFit.Components('bkg_pdf_eff'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_trcospsi,lw)
    _trcospsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(2)

    _trcostheta = trcostheta.frame(RooFit.Bins(15),RooFit.Title('trcostheta full mass region'))
    data.plotOn(_trcostheta)
    pdf.plotOn(_trcostheta,RooFit.Components('sig_pdf_eff'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_trcostheta,RooFit.Components('bkg_pdf_eff'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_trcostheta,lw)
    _trcostheta.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(3)

    _trphi = trphi.frame(RooFit.Bins(15),RooFit.Title('trphi full mass region'))
    data.plotOn(_trphi)
    pdf.plotOn(_trphi,RooFit.Components('sig_pdf_eff'),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_trphi,RooFit.Components('bkg_pdf_eff'),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
    pdf.plotOn(_trphi,lw)
    _trphi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(4)

    _trcospsi = trcospsi.frame(RooFit.Bins(15),RooFit.Title('trcospsi signal region'))
    data.plotOn(_trcospsi,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_trcospsi,RooFit.ProjectionRange('sigRegion'),lw)
    _trcospsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(5)

    _trcostheta = trcostheta.frame(RooFit.Bins(15),RooFit.Title('trcostheta signal region'))
    data.plotOn(_trcostheta,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_trcostheta,RooFit.ProjectionRange('sigRegion'),lw)
    _trcostheta.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(6)

    _trphi = trphi.frame(RooFit.Bins(15),RooFit.Title('trphi signal region'))
    data.plotOn(_trphi,RooFit.CutRange('sigRegion'))
    pdf.plotOn(_trphi,RooFit.ProjectionRange('sigRegion'),lw)
    _trphi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(7)

    _trcospsi = trcospsi.frame(RooFit.Bins(15),RooFit.Title('trcospsi left sideband'))
    data.plotOn(_trcospsi,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_trcospsi,RooFit.ProjectionRange('leftSideband'),lw)
    _trcospsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(8)

    _trcostheta = trcostheta.frame(RooFit.Bins(15),RooFit.Title('trcostheta left sideband'))
    data.plotOn(_trcostheta,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_trcostheta,RooFit.ProjectionRange('leftSideband'),lw)
    _trcostheta.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(9)

    _trphi = trphi.frame(RooFit.Bins(15),RooFit.Title('trphi left sideband'))
    data.plotOn(_trphi,RooFit.CutRange('leftSideband'))
    pdf.plotOn(_trphi,RooFit.ProjectionRange('leftSideband'),lw)
    _trphi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(10)

    _trcospsi = trcospsi.frame(RooFit.Bins(15),RooFit.Title('trcospsi right sideband'))
    data.plotOn(_trcospsi,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_trcospsi,RooFit.ProjectionRange('rightSideband'),lw)
    _trcospsi.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(11)

    _trcostheta = trcostheta.frame(RooFit.Bins(15),RooFit.Title('trcostheta right sideband'))
    data.plotOn(_trcostheta,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_trcostheta,RooFit.ProjectionRange('rightSideband'),lw)
    _trcostheta.Draw()
    c2.Update()

    #==========================================================================================================
    c2.cd(12)

    _trphi = trphi.frame(RooFit.Bins(15),RooFit.Title('trphi right sideband'))
    data.plotOn(_trphi,RooFit.CutRange('rightSideband'))
    pdf.plotOn(_trphi,RooFit.ProjectionRange('rightSideband'),lw)
    _trphi.Draw()
    c2.Update()

    return myline1,myline2,c,c2

##############################################################################
##########################   There we go!!!!! ################################
##############################################################################
from ModelBuilders import *

ws = RooWorkspace("ws")

useTransversityAngles = True

acceptance = False

#######################
### For acceptance: ###
#######################

#xi = { 'parpar' :   1.
#     , '00'     :   1.
#     , 'perpperp' : 1.
#     , 'perp0'  :   0
#     , 'par0'   :   0
#     , 'parperp' :  0 
#     }

xi = { 'parpar' :   0.978745
     , '00'     :   1.02973
     , 'perpperp' : 1.03121
     , 'perp0'  :   0.00130478
     , 'par0'   :   0.
     , 'parperp' :  -0.00211382
     }

def norm_xi( d ) :
    n =  (d['parpar'] + d['00'] + d['perpperp'])/3
    for i in d.iterkeys() : d[i] = d[i]/n

print 'xi =', xi

#coef = [ (0,0,0,   ( xi['parpar']+xi['00']+xi['perpperp'])/3 ) 
#       , (2,0,0,   ( xi['00']     -xi['parpar']   ) *float(5)/3      )
#       , (0,2,0,   ( xi['parpar'] -xi['perpperp'] ) * sqrt(float(20)/9) )
#       , (0,2,-1,  ( xi['parperp'] )*sqrt(float(5)/3) ) 
#       , (1,2,1,   ( xi['perp0']   )*sqrt(float(5)/6)*float(32)/(3*pi) )
#       , (1,2,-2 , ( xi['par0']    )*sqrt(float(5)/6)*float(32)/(3*pi) )
#       ]

#The difference between this one and the one above is exactly the factor 8*sqrt(pi)/9....

#coef = [ (0,0,0,   ( xi['parpar']+xi['00']+xi['perpperp']) * (3./(8.*sqrt(pi)))     ) 
#       , (2,0,0,   ( xi['00']     -xi['parpar']   ) *(15./(8.*sqrt(pi)))            )                     
#       , (0,2,0,   ( xi['parpar'] -xi['perpperp'] ) * (3./4.) * (sqrt(5./pi))       )
#       , (0,2,-1,  ( xi['parperp'] ) * (9./(8.*sqrt(pi))) * (sqrt(5./3.))           ) 
#       , (1,2,1,   ( xi['perp0']   ) * -1. * (12./(pi*sqrt(pi))) * (sqrt(5./6.))     )
#       , (1,2,-2 , ( xi['par0']    ) * (12/(pi*sqrt(pi))) * (sqrt(5./6.))           )
#       ]

#Now start messing around with the reverse enigineering in the last three terms...

#coef = [ (0,0,0,   ( xi['parpar']+xi['00']+xi['perpperp']) * (3./(8.*sqrt(pi)))     ) 
#       , (2,0,0,   ( xi['00']     -xi['parpar']   ) *(15./(8.*sqrt(pi)))            )                     
#       , (0,2,0,   ( xi['parpar'] -xi['perpperp'] ) * (3./4.) * (sqrt(5./pi))       )
#       , (2,2,-1,  ( xi['parperp'] ) * -5. * (9./(8.*sqrt(pi))) * (sqrt(5./3.))           ) 
#       , (3,2,1,   ( xi['perp0']   ) * 1. * (48./(pi*sqrt(pi))) * (sqrt(5./6.))     )
#       , (3,2,-2 , ( xi['par0']    ) * -1.* (48./(pi*sqrt(pi))) * (sqrt(5./6.))           )
#       ]

#Now start messing around with the reverse engineering in the first three terms....

coef = [ (2,2,0,   ( xi['perpperp']) * (45./8.) * sqrt(5./pi)                                              ) 
       , (2,2,2,   ( xi['parpar'] + xi['00'] - 0.5 * xi['perpperp']  ) * (-15./8.) * sqrt(20./(3.*pi))     )                     
       , (0,2,2,   ( xi['parpar'] -0.5 * xi['00'] + xi['perpperp'] ) * (3./4.) * (20./(3.*pi))             )
       , (0,2,-1,  ( xi['parperp'] ) * (9./(8.*sqrt(pi))) * (sqrt(5./3.))                                  )  
       , (1,2,1,   ( xi['perp0']   ) * -1. * (12./(pi*sqrt(pi))) * (sqrt(5./6.))                           )
       , (1,2,-2 , ( xi['par0']    ) * (12/(pi*sqrt(pi))) * (sqrt(5./6.))                                  )
       ]

print 'reverse engineerd c_ijk:'
print coef

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

ws.factory("{rz2[0.460,0.4,0.7],rperp2[0.347,0.2,0.5]")
ws.factory("RooFormulaVar::rpar2('1-@0-@1',{rz2,rperp2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")

ws.factory("{deltaz[0.],deltapar[%f,0.,%f],deltaperp[-0.17]}"%(pi,2*pi))

##########################
### physics parameters ###
##########################
#tau = 1.47 defines gamma:
#ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.05,-0.3,0.3]}")

ws.factory("{#Gamma[0.68,0.4,0.9]}")
ws.factory("RooUnblindUniform::#Gamma_unblind('BsCalvin',0.4,#Gamma)")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma_unblind})")

ws.factory("{dG[0.05,-0.3,0.3]}")
ws.factory("RooUnblindUniform::t_sig_dG('BsHobbes',0.2,dG)")

ws.factory("{t_sig_dm[17.7]}")

ws.factory('{phis[0.]}')

ws.factory("{expr::S('-1*sin(phis)',{phis}),expr::D('cos(phis)',{phis}),C[0]}")
ws.factory("{expr::Sold('sin(phis)',{phis}),expr::Dold('cos(phis)',{phis}),Cold[0]}")

###############################
### Experimental parameters ###
###############################

ws.factory("RooGaussModel::tres_sig(t,mu[0],sigma[0.05])")

ws.factory("{wtag[0.5]}")

##################################### building the NEW pdf ###################################
# create observables
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

#This one imports automatically in the workspace!
if useTransversityAngles:
    newpdf = buildJpsiphi(ws,'newpdf', True) 
else:
    newpdf = buildJpsiphi(ws,'newpdf', False)

#print 'GOING TO BUILD THE ACCEPTANCE CORRECTED FULL PDF!'
#signalcorrpdf = buildEff_x_PDF(ws,'eff',newpdf,[ ( ab.build('mom_eff',c[0],0,c[1],c[2],1.), c[3] ) for c in coef ] )
#ws.put(signalcorrpdf)

##############################
#This sets the range for the events in the dataset that you read, and the range in which you will fit!
##############################
t = ws.var('t')

tmin = 0.3
tmax = 14.
t.setRange(tmin,tmax)

if useTransversityAngles:
    trcostheta = ws.var('trcostheta')
    trcospsi = ws.var('trcospsi')
    trphi = ws.var('trphi')
else:
    cthetaK = ws.var('helcosthetaK')
    cthetaL = ws.var('helcosthetaL')
    phi = ws.var('helphi')

tagdecision = ws.cat('tagdecision')

#wide
m = RooRealVar("m","B_{s} mass",5200,5550)
#narrow
#m = RooRealVar("m","B_{s} mass",5321.67,5411.67)

obs = RooArgSet(m)

getattr(ws,'import')(obs)

if useTransversityAngles:
    ws.defineSet("observables","m,t,trcospsi,trcostheta,trphi,tagdecision")
else:
    ws.defineSet("observables","m,t,helcosthetaL,helcosthetaK,helphi,tagdecision")

#################
### Load Data ###
#################

#Using Edinburgh file
file = TFile('/data/bfys/dveijk/Data/Edinburgh.root')
NTupletree = file.Get('MyTree')

if useTransversityAngles:
    data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m && trcospsi==trcospsi && trcostheta==trcostheta && trphi==trphi && tagdecision==tagdecision')
else:
    data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m && helcosthetaL==helcosthetaL && helcosthetaK==helcosthetaK && helphi==helphi && tagdecision==tagdecision')
    
data.table(tagdecision).Print('v')

getattr(ws,'import')(data)

#######################
### Build the PDF's ###
#######################

ws.factory("Gaussian::m_sig(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")

#background B mass pdf
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

#background propertime 
# TODO: split resolution in tres_sig and tres_nonpsi!!
ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.1,0.5],tres_sig,SingleSided)")
ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres_sig,SingleSided)")
ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")

#background angles: 

#if useTransversityAngles:
#    ws.factory("Uniform::bkgang({trcospsi,trcostheta,trphi})")
#else:
#    ws.factory("Uniform::bkgang({helcosthetaL,helcosthetaK,helphi})")

#if useTransversityAngles:
#    ws.factory("Chebychev::bkg_trcospsi(trcospsi,{c0_trcospsi[-0.13,-1,1]})")
#    ws.factory("Chebychev::bkg_trcostheta(trcostheta,{c0_trcostheta[0.08,-1,1]})")
#    ws.factory("Chebychev::bkg_trphi(trphi,{c0_trphi[0.10,-1,1]})")
#    ws.factory("PROD::bkgang(bkg_trcospsi,bkg_trcostheta,bkg_trphi)")
#else:
#    ws.factory("Chebychev::bkg_helcosthetaK(helcosthetaK,{c0_helcosthetaK[-0.13,-1,1]})")
#    ws.factory("Chebychev::bkg_helcosthetaL(helcosthetaL,{c0_helcosthetaL[0.08,-1,1]})")
#    ws.factory("Chebychev::bkg_helphi(helphi,{c0_helphi[0.10,-1,1]})")
#    ws.factory("PROD::bkgang(bkg_helcosthetaL,bkg_helcosthetaK,bkg_helphi)")

#Build angular background with RooP2VVAngleBasis
_ba = lambda name,comp : _buildAngularFunction(ws,ab,name,comp)

_ba("Bkg_0000",  [ ( 0,0,0,0,1.) ] )

_ba("Bkg_0010",  [ ( 0,0,1,0,1.) ] )
_ba("Bkg_0011",  [ ( 0,0,1,1,1.) ] )
_ba("Bkg_001m1",  [ ( 0,0,1,-1,1.) ] )

_ba("Bkg_1000",  [ ( 1,0,0,0,1.) ] )
_ba("Bkg_2000",  [ ( 2,0,0,0,1.) ] )
_ba("Bkg_3000",  [ ( 3,0,0,0,1.) ] )

c0000 = RooRealVar('c0000','c0000',1.)

c0010 = RooRealVar('c0010','c0010',0.,-1.,1.)
c0011 = RooRealVar('c0011','c0011',0.,-1.,1.)
c001m1 = RooRealVar('c001m1','c001m1',0.,-1.,1.)

c1000 = RooRealVar('c1000','c1000',0.,-1.,1.)
c2000 = RooRealVar('c2000','c2000',0.,-1.,1.)
c3000 = RooRealVar('c3000','c3000',0.,-1.,1.)

#Flat background in all angles
bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis']),RooArgList(c0000))

#Only higher order in cospsi
#bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c1000,c2000,c3000))

#Add higher order in costheta
#bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c1000,c2000,c3000))

#Hiher roder in all angles
#bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_0011_basis'],ws['Bkg_001m1_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c0011,c001m1,c1000,c2000,c3000))

ws.put(bkg)

#ws.factory("PROD::sig_pdf( m_sig, newpdf)")
ws.factory("PROD::sig_pdf( m_sig, newpdf)")

#ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, bkgang)")
ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, bkg)")

ws.factory("SUM::pdf_ext(Nsig[543,0,1000]*sig_pdf,Nbkg[812,0,1000]*bkg_pdf)")
ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sig_pdf,bkg_pdf)")


#########################
### What do you want? ###
#########################

#get pdf's from the workspace
bkg_pdf = ws.pdf('bkg_pdf')
sig_pdf = ws.pdf('sig_pdf')
pdf= ws.pdf('pdf')
pdf_ext = ws.pdf('pdf_ext')

print 'GOING TO BUILD THE ACCEPTANCE CORRECTED FULL PDF!'
fullcorrpdf = buildEff_x_PDF(ws,'eff',pdf_ext,[ ( ab.build('mom_eff',c[0],0,c[1],c[2],1.), c[3] ) for c in coef ] )
ws.put(fullcorrpdf)
assert False
#result = pdf_ext.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
#result = fullcorrpdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))

if useTransversityAngles:
    #plot = plottrans(ws,data,pdf_ext,'transplot')
    plot = plottrans(ws,data,fullcorrpdf,'transplot')
else:
    #plot = plothel(ws,data,pdf_ext,'transplot')
    plot = plothel(ws,data,fullcorrpdf,'helplot')

####################
# Latex code of the fitted parameters
paramlist = fullcorrpdf.getParameters(data)
paramlist.printLatex(RooFit.Format("NEU",RooFit.AutoPrecision(3),RooFit.VerbatimName()))
####################

