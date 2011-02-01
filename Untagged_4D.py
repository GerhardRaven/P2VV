from ROOT import *
gSystem.Load("libp2vv")
from math import pi

import rootStyle
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()


ncpu = RooCmdArg( RooFit.NumCPU(16) )
sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))
lw = RooCmdArg(RooFit.LineWidth(2))
ms = RooCmdArg(RooFit.MarkerSize(0.4))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))
stash = []

#######################
### Plot ICHEP Like ###
#######################
def plotme(ws,data,pdf,title):

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

##############################################################################
##########################   There we go!!!!! ################################
##############################################################################

from ModelBuilders import *
from RooFitDecorators import *

ws = RooWorkspace('ws')
###################
### Observables ###
###################
mode = 'Bs2Jpsiphi'

declareObservables(ws,mode)

files = { 'Bu2JpsiK'  : 'Bu2JpsiKTuple.root'
          , 'Bs2Jpsiphi': 'Bs2JpsiPhiStripping12_NoJpsiMassCut.root'
          , 'Bd2JpsiKstar' : 'Bd2JpsiKstarTuple.root'
        }

tmin = -1
tmax = 14

ws.var('t').setMin(tmin)
ws.var('t').setMax(tmax)

mmin = 5200
mmax = 5550

ws.var('m').setMin(mmin)
ws.var('m').setMax(mmax)

obs = ws.set('observables')

angles = ws.set('helicityangles')
#angles = ws.set('transversityangles')

obs.add( angles )

file = TFile(files[mode])
#data = RooDataSet('data','data',file.Get('dataset'),obs,'t==t && m==m')
data = RooDataSet('data','data',file.Get('MyTree'),obs,'t==t && m==m')

mb = MassPdfBuilder(ws,ws['m'],ws['mdau1'],ws['mdau2'] if mode != 'Bu2JpsiK' else None,mode)

####
massc = TCanvas()
massc.Divide(2,3)

### TODO: move these plots into the mass PDF builder...
if True :
    mdau1pdf = mb.dau1Pdf()
    #mdau1pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mdau1pdf.fitTo(data,ncpu)
    for i in mdau1pdf.getParameters(data): i.setConstant(True)
    comps = { 'mpsi_sig'       : ( sigcolor,dashed )
            , 'mpsi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( massc.cd(1), mb.dau1Obs(), data, mdau1pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True and mode != 'Bu2JpsiK' :
    massc.cd(2)
    mdau2pdf = mb.dau2Pdf()
    #mdau2pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mdau2pdf.fitTo(data,ncpu)
    for i in mdau2pdf.getParameters(data) : i.setConstant(True)
    comps = { 'mphi_sig'       : ( sigcolor,dashed )
            , 'mphi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( massc.cd(2), mb.dau2Obs(), data, mdau2pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True :
    mpdf = mb.Pdf()
    #mpdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mpdf.fitTo(data,ncpu)
    for i in mpdf.getParameters(data) : i.setConstant(True)
    ## TODO: at some point, we need to add the KK mass to the story...
    for i,obs in enumerate( [ mb.Obs(), mb.dau1Obs()])  : #  , mb.dau2Obs() ] ):
        # TODO: plot current in signal of other, sideband of other....
        comps = { 'm_sig'       : ( sigcolor,dashed )
                , 'm_psibkg'    : ( bkgcolor,dashed ) 
                , 'm_nonpsibkg' : ( nonpsicolor,dashed )
                }
        plot( massc.cd(3+i), obs, data, mpdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

massc.Print("InitialMassFits.eps")
#import sys
#sys.exit(0)

### Now that we have the masses fitted, let's make angle SPlots...
## TODO: move this into the mass PDF builder...
for i in mpdf.getParameters(data) : i.setConstant( i not in mb.yields() )
#mpdf.Print("T")
#mpdf.getParameters(data).Print("V")
splot = RooStats.SPlot("splotdata","splotdata",data,mpdf,mb.yields())
wdata = splot.GetSDataSet()

if True :
    c = TCanvas()
    c2 = TCanvas()
    stash.append(c)
    stash.append(c2)
    observables = RooArgSet(angles)
    observables.add(ws['t'])
    observables.add(ws['sigmat'])
    observables.add(ws['tagomega'])
    observables.add(ws['mdau2'])
    c.Divide(len(observables),3)
    c2.Divide(1,3)
    for (i,sample) in enumerate([ j.GetName() for j in mb.yields() ]):
        dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","%s_sw"%sample) # need a dummy cut, as passing a (const char*)0 is kind of difficult...
        # should only make splots for observables not used in fit...
        for (j,observable) in enumerate( observables ) :
            f = observable.frame()
            dataw.plotOn(f, RooFit.DataError(RooAbsData.SumW2) )
            c.cd(i*len(observables)+j+1)
            f.Draw()
        c2.cd(1+i)
        f = ws['t'].frame(-2.,2.,100)
        dataw.plotOn(f)
        f.Draw()

if True :
    # and now, build the angular distributions for psi and non-psi background, using the SWeights we
    # just got...
    ab = abasis(ws,angles)
    x = BkgAnglePdfBuilder( ws, ab, data
                          , { 'psibkg'    : { 'ranges' : (range(3),range(8),range(-3,4))  , 'weight' : 'N_psibkg_sw'    }
                            , 'nonpsibkg' : { 'ranges' : (range(3),range(12),range(-3,4)) , 'weight' : 'N_nonpsibkg_sw' } 
                            } )
    (c1,c2) = x.makeplots()
## next, we pick up the sigmat distributions for our three components...

sigmat = ws['sigmat']
sigmat.setBins(40)
p= sigmat.frame()
stdata = {}
stpdf = {}
c = TCanvas()
for (f,sample) in enumerate([ 'sig','psibkg','nonpsibkg' ]):
      dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","N_%s_sw"%sample) 
      stdata[sample] = RooDataHist("sigmat_%s_data"%sample,"hist Err Per Ev",RooArgSet(sigmat),dataw)
      stpdf[sample] = ws.put(RooHistPdf("sigmat_%s"%sample,"sigmat_%s"%sample,RooArgSet(sigmat),stdata[sample])) # ,2)
      dataw.plotOn(p) # python has trouble finding the right plotOn...
      stpdf[sample].plotOn(p,RooFit.LineColor( [kRed,kBlue,kBlack][f]))
p.Draw()

y = TimeResolutionBuilder(ws, ws['t'],ws['sigmat'])
z = BkgTimePdfBuilder(ws, y, stpdf)

ws.factory("wtag[0.5]")
ws.factory("{S[0],D[1],C[0]}")

ws.factory("{#Gamma[0.68,0.4,0.9]}")
ws.factory("RooUnblindUniform::#Gamma_unblind('BsCalvin',0.4,#Gamma)")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma_unblind})")

ws.factory("{dG[0.05,-0.3,0.3]}")
ws.factory("RooUnblindUniform::t_sig_dG('BsHobbes',0.2,dG)")

ws.factory("{t_sig_dm[17.7]}")
definePolarAngularAmplitudes(ws)
print 'about to build J/psi phi signal'
sigjpsiphi = buildJpsiphi(ws,'jpsiphisignal',False)
print 'about to createProjection of J/psi phi signal over tagdecision'
#t_sig = sigjpsiphi # ws.put(sigjpsiphi.createProjection( ws.argSet('tagdecision') ))
t_sig = ws.put(sigjpsiphi.createProjection( ws.argSet('tagdecision') ))

######################
### Multiplying... ###
######################

m_sig = ws.pdf('m_sig')

sig = RooProdPdf('sig','sig', m_sig,t_sig)
ws.put(sig)

ws.factory("PROD:psibkg(    m_psibkg,    t_psibkg,    %s )" % x.psibkgPdf().GetName())
ws.factory("PROD:nonpsibkg( m_nonpsibkg, t_nonpsibkg, %s )" % x.nonpsibkgPdf().GetName())

ws.factory("SUM::pdf(f_sig[0.1,0,0.4]*sig, SUM(f_psi[0.5,0.1,0.9]*psibkg,nonpsibkg))")

pdf = ws['pdf']
pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Time')
ws['m_sig_sigma'].setConstant(False)
#ws['m_sig_f1'].setConstant(False)
#ws['m_sig_sigma2'].setConstant(False)
result = pdf.fitTo(data,ncpu,RooFit.Save(True),RooFit.Minos(false))
pdf.getParameters(data).writeToFile('fitresult_%s.txt'%mode)

t = ws['t']
m = ws['m']
t.setRange('largeTime',0.3,t.getMax())

c = TCanvas()
c.Divide(5,2)

#===========================================================================================================
bkgcomps = { 'psibkg'    : ( bkgcolor,dashed ) 
           , 'nonpsibkg' : ( nonpsicolor,dashed )
           }
comps = { 'sig'       : ( sigcolor,dashed )
        , 'psibkg'    : ( bkgcolor,dashed ) 
        , 'nonpsibkg' : ( nonpsicolor,dashed )
        }

#===========================================================================================================
_c1 = plot( c.cd(1),ws['mdau1'],data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), )
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c2 = plot( c.cd(2),ws['mdau1'],data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m(#mu#mu) (t>0.3)') )
          , dataOpts = ( ms, xes, RooFit.CutRange('largeTime') )
          , pdfOpts = ( lw, RooFit.ProjectionRange('largeTime') ) 
          )
#===========================================================================================================
_c3 = plot( c.cd(3),m,data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m') )
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c4 = plot( c.cd(4),m,data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m (t>0.3)') )
          , dataOpts = ( ms, xes, RooFit.CutRange('largeTime') )
          , pdfOpts = ( lw, RooFit.ProjectionRange('largeTime') ) 
          )
#===========================================================================================================
_c5 = plot( c.cd(5),sigmat,data,pdf,comps
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c6 = plot( c.cd(6), t, data, pdf, comps
          , frameOpts = ( RooFit.Range(-0.4,0.4), RooFit.Bins(100), RooFit.Title('proper time, full mass range') )
          , dataOpts = ( ms,xes )
          , pdfOpts = ( lw, )
          )
#===========================================================================================================
_c7 = plot( c.cd(7), t, data, pdf, comps
          , frameOpts = ( RooFit.Title('proper time, signal region'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('sigRegion') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('sigRegion') )
          , logy = True
          )
#===========================================================================================================
_c8 = plot( c.cd(8), t, data, pdf, bkgcomps
          , frameOpts = ( RooFit.Title('proper time, lower sideband'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('leftSideband') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('leftSideband') )
          , logy = True
          )
#===========================================================================================================
_c9 = plot( c.cd(9), t, data, pdf, bkgcomps
          , frameOpts = ( RooFit.Title('proper time, upper sideband'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('rightSideband') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('rightSideband') )
          , logy = True
          )


c.Print("FinalPlots.eps")


####################
# Latex code of the fitted parameters
paramlist = pdf.getParameters(data)
paramlist.printLatex(RooFit.Format("NEU",RooFit.AutoPrecision(2),RooFit.VerbatimName()))
####################
