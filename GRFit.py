from ROOT import *
from RooFitDecorators import *
gSystem.Load("libp2vv")
from math import pi

#######################
### Plot ICHEP Like ###
#######################

ncpu = RooCmdArg( RooFit.NumCPU(8) )
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
ws = RooWorkspace('ws')

mode = 'Bu2JpsiK'
#mode = 'Bs2Jpsiphi'

declareObservables(ws,mode)

files = { 'Bu2JpsiK'  : 'Bu2JpsiKTuple.root'
        , 'Bs2Jpsiphi': 'Bs2JpsiPhiTuple.root'
        , 'Bd2JpsiKstar' : 'Bd2JpsiKstarTuple.root'
        }
obs = ws.set('observables')
if mode != 'Bu2JpsiK' : 
    angles = ws.set('helicityangles')
    obs.add( angles )

file = TFile(files[mode])
data = RooDataSet('data','data',file.Get('dataset'),obs,'t==t && m==m && unbiased==1')


mb = MassPdfBuilder(ws,ws['m'],ws['mdau1'],ws['mdau2'] if mode != 'Bu2JpsiK' else None,mode)
                                            
####
c = TCanvas()
c.Divide(2,3)

### TODO: move these plots into the mass PDF builder...
if True :
    mdau1pdf = mb.dau1Pdf()
    mdau1pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mdau1pdf.fitTo(data,ncpu)
    for i in mdau1pdf.getParameters(data): i.setConstant(True)
    comps = { 'mpsi_sig'       : ( sigcolor,dashed )
            , 'mpsi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( c.cd(1), mb.dau1Obs(), data, mdau1pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True and mode != 'Bu2JpsiK' :
    c.cd(2)
    mdau2pdf = mb.dau2Pdf()
    mdau2pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mdau2pdf.fitTo(data,ncpu)
    for i in mdau2pdf.getParameters(data) : i.setConstant(True)
    comps = { 'mphi_sig'       : ( sigcolor,dashed )
            , 'mphi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( c.cd(2), mb.dau2Obs(), data, mdau2pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True :
    mpdf = mb.Pdf()
    mpdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
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
#import sys
#sys.exit(0)

### Now that we have the masses fitted, let's make angle SPlots...
## TODO: move this into the mass PDF builder...
for i in mpdf.getParameters(data) : i.setConstant( i not in mb.yields() )
#mpdf.Print("T")
#mpdf.getParameters(data).Print("V")
splot = RooStats.SPlot("splotdata","splotdata",data,mpdf,mb.yields())
wdata = splot.GetSDataSet()

if False :
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

if False :
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
      
x = TimeResolutionBuilder(ws, ws['t'],ws['sigmat'])
y = BkgTimePdfBuilder(ws, x, stpdf)

tpb = TagPdfBuilder(ws,[0.25,0.35,0.45],ws['tagomega'])

if True :
    #ab = abasis(ws,angles)
    #x = BkgAnglePdfBuilder( ws, ab, data
    #                        , { 'psibkg'    : { 'ranges' : (range(4),range(8),range(-4,5))  , 'weight' : 'N_psibkg_sw'    }
    #                            , 'nonpsibkg' : { 'ranges' : (range(4),range(12),range(-4,5)) , 'weight' : 'N_nonpsibkg_sw' } 
    #                            } )
    #(c1,c2) = x.makeplots()

    # create a dummy signal PDF...
    ws.factory("PROD::t_sig( Decay(t,t_sig_tau[1.5,1.2,1.8],tres_sig,SingleSided)|sigmat,sigmat_sig)")
    # map decision=+1,JpsiK+ to correct, decision=-1,J/psiKmin to incorrect
    #  ws.factory( "MultiCategory({decaytype,tagdecision}")
    #  ws.factory( "MappedCategory::rmc( JpsiKplus, 0)")
    #  rmc = ws['rmc']
    #  rmc.map('{JpsiKplus,bbar},unmix,+1')
    #  rmc.map('{JpsiKmin,b},unmix,+1')
    #  rmc.map('{JpsiKplus,b},mix,-1')
    #  rmc.map('{JpsiKmin,bbar},mix,-1')
    # ws.factory( Efficiency::tag_sig( unmix,w )
    ### and now multiply and sum together...
    # now add the flavour tagging -- pretend we are J/psi K+, and just fit for 
    # tagging efficiency. Nothing more.
    #ws.factory("PROD::sig(       m_sig,       t_sig       , tagcat_sig        )")
    #ws.factory("PROD::psibkg(    m_psibkg,    t_psibkg    , tagcat_psibkg     )")
    #ws.factory("PROD::nonpsibkg( m_nonpsibkg, t_nonpsibkg , tagcat_nonpsibkg  )")

    ws.factory("PROD::sig(       m_sig,       t_sig         )")
    ws.factory("PROD::psibkg(    m_psibkg,    t_psibkg      )")
    ws.factory("PROD::nonpsibkg( m_nonpsibkg, t_nonpsibkg   )")

else :
    ws.factory("wtag[0.5]")
    ws.factory("{S[0],D[1],C[0]}")
    ws.factory("{t_sig_tau[1.4,1.2,1.8],t_sig_dm[17.7],t_sig_dG[0.05,-0.3,0.3]}")
    definePolarAngularAmplitudes(ws)
    print 'about to build J/psi phi signal'
    sigjpsiphi = buildJpsiphi(ws,'jpsiphisignal',False)
    print 'about to createProjection of J/psi phi signal over tagdecision'
    #t_sig = sigjpsiphi # ws.put(sigjpsiphi.createProjection( ws.argSet('tagdecision') ))
    t_sig = ws.put(sigjpsiphi.createProjection( ws.argSet('tagdecision') ))
    print 'about to multiply signal'
    ws.factory("PROD:ta_sig( %s |sigmat,sigmat_sig)" % t_sig.GetName() )
    ### and now multiply and sum together...
    print 'about to multiply background'
    ws.factory("PROD:sig(       m_sig,       ta_sig         )")
    ws.factory("PROD:psibkg(    m_psibkg,    t_psibkg,    %s )" % x.psibkgPdf().GetName())
    ws.factory("PROD:nonpsibkg( m_nonpsibkg, t_nonpsibkg, %s )" % x.nonpsibkgPdf().GetName())

ws.factory("SUM::pdf(f_sig[0.1,0,0.4]*sig, SUM(f_psi[0.5,0.1,0.9]*psibkg,nonpsibkg))")

pdf = ws['pdf']
pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Time')
ws['m_sig_sigma'].setConstant(False)
ws['m_sig_sigma2'].setConstant(False)
ws['m_sig_f1'].setConstant(False)
ws['m_sig_fpsi'].setConstant(False)

### Maybe we should first fit the three time splots with the three components to get them 
### to a reasonable starting point??
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


if False :

    # pick up the psi and nonpsi mass PDFs and multiply them by their angles...

    ws.factory("PROD::psibkg(%s,%s)"%(mb.psibkgPdf().GetName(),abf.psibkgPdf().GetName() ) )
    ws.factory("PROD::nonpsibkg(%s,%s)"%(mb.nonpsibkgPdf().GetName(),abf.nonpsibkgPdf().GetName() ))
    ws.factory("PROD::sig(%s,%s)"%(mb.sigPdf().GetName(),abf.dummysigPdf().GetName() ))
    ws.factory("SUM::bkg(fpsi[0.5,0.01,0.99]*psibkg,nonpsibkg)")
    ws.factory("SUM::pdf(Nsig[12000,0,1000]*sig,Nbkg[81200,0,1000000]*bkg)")

    pdf = ws['pdf']
    pdf.fitTo(data,ncpu)
    plotAngles( data,pdf, angles )

