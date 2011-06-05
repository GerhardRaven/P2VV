########################################
### Author: Daan van Eijk
### Updated on: Jun 5 11
### Description: This script is a more complicated version of Note1SimpleFit.py.
###              It uses GR's classes in ModelBuilders to
###              Take jpsi mass into account in signal and bkg mass pdf (using GR's MassPdfBuilder)
###              Split background in jpsi and non-jpsi
###              Make sPlots from the initial mass fits
###              Use per event proper time errors from sPlots
###              Resolution: TimeResolutionBuilder(ws, ws['t'],ws['sigmat'])
###              Bgk: BkgTimePdfBuilder(ws, x, stpdf)
######################################## 

from ROOT import *
from RooFitDecorators import *
gSystem.Load("libp2vv")
from math import pi

#######################
### Plot ICHEP Like ###
#######################

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

##################################
### Create WS, build the PDF's ###
##################################
from ModelBuilders import *
  
new = True

if new:
    files = { 'Bu2JpsiK'  : 'Bu2JpsiKTuple.root'
              , 'Bs2Jpsiphi': 'Bs2JpsiPhiStripping12_NoJpsiMassCut.root'
              , 'Bd2JpsiKstar' : 'Bd2JpsiKstarTuple.root'
              }
        
else:
    files = { 'Bu2JpsiK'  : 'Bu2JpsiKTuple.root'
              , 'Bs2Jpsiphi': 'Bs2JpsiPhiTuple.root'
              , 'Bd2JpsiKstar' : 'Bd2JpsiKstarTuple.root'
              }
    
ws = RooWorkspace('ws')

#mode = 'Bu2JpsiK'
mode = 'Bs2Jpsiphi'

declareObservables(ws,mode)

tmin = -1
tmax = 14

ws.var('t').setMin(tmin)
ws.var('t').setMax(tmax)

mmin = 5200
mmax = 5550

ws.var('m').setMin(mmin)
ws.var('m').setMax(mmax)

obs = ws.set('observables')
if mode != 'Bu2JpsiK' : 
    angles = ws.set('helicityangles')
    obs.add( angles )

file = TFile(files[mode])

if new:
    data = RooDataSet('data','data',file.Get('MyTree'),obs,'t==t && m==m')
else:
    data = RooDataSet('data','data',file.Get('dataset'),obs,'t==t && m==m')
    
mb = MassPdfBuilder(ws,ws['m'],ws['mdau1'],ws['mdau2'] if mode != 'Bu2JpsiK' else None,mode)
                                            
####
c = TCanvas()
c.Divide(2,3)

### TODO: move these plots into the mass PDF builder...
if True :
    mdau1pdf = mb.dau1Pdf()
    mdau1pdf.fitTo(data,ncpu)
    for i in mdau1pdf.getParameters(data): i.setConstant(True)
    comps = { 'mpsi_sig'       : ( sigcolor,dashed )
            , 'mpsi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( c.cd(1), mb.dau1Obs(), data, mdau1pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True and mode != 'Bu2JpsiK' :
    c.cd(2)
    mdau2pdf = mb.dau2Pdf()
    mdau2pdf.fitTo(data,ncpu)
    for i in mdau2pdf.getParameters(data) : i.setConstant(True)
    comps = { 'mphi_sig'       : ( sigcolor,dashed )
            , 'mphi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( c.cd(2), mb.dau2Obs(), data, mdau2pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True :
    mpdf = mb.Pdf()
    #ws['N_sig'].setVal(1000)  # 13K for J/psi K+
    #ws['N_psibkg'].setVal( 0.6*(data.numEntries()- ws['N_sig'].getVal()) )
    #ws['N_nonpsibkg'].setVal( 0.4*(data.numEntries()- ws['N_sig'].getVal()) )
    mpdf.fitTo(data,ncpu)
    for i in mpdf.getParameters(data) : i.setConstant(True)
    ## TODO: at some point, we need to add the KK mass to the story...
    for i,obs in enumerate( [ mb.Obs(), mb.dau1Obs()])  : #  , mb.dau2Obs() ] ):
        # TODO: plot current in signal of other, sideband of other....
        comps = { 'm_sig'       : ( sigcolor,dashed )
                , 'm_psibkg'    : ( bkgcolor,dashed ) 
                , 'm_nonpsibkg' : ( nonpsicolor,dashed )
                }
        plot( c.cd(3+i), obs, data, mpdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

c.Print("InitialMassFits.eps")

### Now that we have the masses fitted, let's make angle SPlots...
## TODO: move this into the mass PDF builder...
for i in mpdf.getParameters(data) : i.setConstant( i not in mb.yields() )
#mpdf.Print("T")
#mpdf.getParameters(data).Print("V")

splot = RooStats.SPlot("splotdata","splotdata",data,mpdf,mb.yields())

wdata = splot.GetSDataSet()

c2 = TCanvas()
c3 = TCanvas()
stash.append(c2)
stash.append(c3)
#observables = RooArgSet(angles)
#observables.add(ws['t'])
observables = RooArgSet(ws['t'])
observables.add(ws['sigmat'])
#observables.add(ws['mdau2'])
c2.Divide(len(observables),3)
c3.Divide(1,3)
for (i,sample) in enumerate([ j.GetName() for j in mb.yields() ]):
    dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","%s_sw"%sample) # need a dummy cut, as passing a (const char*)0 is kind of difficult...
    # should only make splots for observables not used in fit...
    for (j,observable) in enumerate( observables ) :
        f = observable.frame()
        dataw.plotOn(f, RooFit.DataError(RooAbsData.SumW2) )
        c2.cd(i*len(observables)+j+1)
        f.Draw()
    c3.cd(1+i)
    f = ws['t'].frame(-0.5,1.,100)
    dataw.plotOn(f)
    f.Draw()
## next, we pick up the sigmat distributions for our three components...

sigmat = ws['sigmat']
sigmat.setBins(40)
p= sigmat.frame()
stdata = {}
stpdf = {}
c2 = TCanvas()
for (f,sample) in enumerate([ 'sig','psibkg','nonpsibkg' ]):
      dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","N_%s_sw"%sample) 
      stdata[sample] = RooDataHist("sigmat_%s_data"%sample,"hist Err Per Ev",RooArgSet(sigmat),dataw)
      stpdf[sample] = ws.put(RooHistPdf("sigmat_%s"%sample,"sigmat_%s"%sample,RooArgSet(sigmat),stdata[sample])) # ,2)
      dataw.plotOn(p) # python has trouble finding the right plotOn...
      stpdf[sample].plotOn(p,RooFit.LineColor( [kRed,kBlue,kBlack][f]))
p.Draw()

x = TimeResolutionBuilder(ws, ws['t'],ws['sigmat'])
y = BkgTimePdfBuilder(ws, x, stpdf)

#To reset some defaults from ModelBuilders
ws['t_nonpsibkg_fll'].setVal(0.004)
ws['t_nonpsibkg_fll'].setConstant(False)
ws['t_nonpsibkg_ll_tau'].setVal(1.92)
ws['t_nonpsibkg_ll_tau'].setConstant(False)
ws['t_psibkg_fll'].setVal(0.004)
ws['t_psibkg_fll'].setConstant(False)
ws['t_psibkg_ll_tau'].setVal(1.92)
ws['t_psibkg_ll_tau'].setConstant(False)

if True :
    # create a dummy signal PDF...
    ws.factory("{#tau[1.47,1.,2.]}")

    ws.factory("RooUnblindUniform::tau_unblind('BsCalvin',0.2,#tau)")
    
    ws.factory("PROD::t_sig( Decay(t,tau_unblind,tres_sig,SingleSided)|sigmat,sigmat_sig)")

    ws.factory("PROD::sig(       m_sig,       t_sig         )")
    ws.factory("PROD::psibkg(    m_psibkg,    t_psibkg      )")
    ws.factory("PROD::nonpsibkg( m_nonpsibkg, t_nonpsibkg   )")

ws.factory("SUM::pdf(f_sig[0.1,0,0.4]*sig, SUM(f_psi[0.5,0.1,0.9]*psibkg,nonpsibkg))")

pdf = ws['pdf']

if new:
    print 'New, do nothing'
else:
    pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode)

ws['m_sig_sigma'].setConstant(False)
result = pdf.fitTo(data,ncpu,RooFit.Save(True),RooFit.Minos(false))
pdf.getParameters(data).writeToFile('fitresult_%s.txt'%mode)

t = ws['t']
m = ws['m']
t.setRange('largeTime',0.3,t.getMax())

c4 = TCanvas()
c4.Divide(4,2)

#===========================================================================================================
bkgcomps = { 'psibkg'    : ( bkgcolor,dashed ) 
           , 'nonpsibkg' : ( nonpsicolor,dashed )
           }
comps = { 'sig'       : ( sigcolor,dashed )
        , 'psibkg'    : ( bkgcolor,dashed ) 
        , 'nonpsibkg' : ( nonpsicolor,dashed )
        }

#===========================================================================================================
_c1 = plot( c4.cd(1),ws['mdau1'],data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), )
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c2 = plot( c4.cd(2),m,data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m') )
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c3 = plot( c4.cd(3),m,data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m (t>0.3)') )
          , dataOpts = ( ms, xes, RooFit.CutRange('largeTime') )
          , pdfOpts = ( lw, RooFit.ProjectionRange('largeTime') ) 
          )
#===========================================================================================================
_c4 = plot( c4.cd(4),sigmat,data,pdf,comps
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
## #===========================================================================================================
## _c5 = plot( c4.cd(5), t, data, pdf, comps
##           , frameOpts = ( RooFit.Range(-0.4,0.4), RooFit.Bins(100), RooFit.Title('proper time, full mass range') )
##           , dataOpts = ( ms,xes )
##           , pdfOpts = ( lw, )
##           )
## #===========================================================================================================
## _c6 = plot( c4.cd(6), t, data, pdf, comps
##           , frameOpts = ( RooFit.Title('proper time, signal region'), )
##           , dataOpts = ( ms,xes,RooFit.CutRange('sigRegion') )
##           , pdfOpts = ( lw,RooFit.ProjectionRange('sigRegion') )
##           , logy = True
##           )
## #===========================================================================================================
## _c7 = plot( c4.cd(7), t, data, pdf, bkgcomps
##           , frameOpts = ( RooFit.Title('proper time, lower sideband'), )
##           , dataOpts = ( ms,xes,RooFit.CutRange('leftSideband') )
##           , pdfOpts = ( lw,RooFit.ProjectionRange('leftSideband') )
##           , logy = True
##           )
## #===========================================================================================================
## _c8 = plot( c4.cd(8), t, data, pdf, bkgcomps
##           , frameOpts = ( RooFit.Title('proper time, upper sideband'), )
##           , dataOpts = ( ms,xes,RooFit.CutRange('rightSideband') )
##           , pdfOpts = ( lw,RooFit.ProjectionRange('rightSideband') )
##           , logy = True
##           )



