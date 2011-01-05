from ROOT import *
gSystem.Load("libp2vv")
from RooFitDecorators import *

#RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Plotting))

ncpu = RooCmdArg( RooFit.NumCPU(1) )
sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))
lw = RooCmdArg(RooFit.LineWidth(2))
ms = RooCmdArg(RooFit.MarkerSize(0.4))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

#RooAbsPdf.plotComponenentOn = plotPdfComponents

##################################
### Create WS, build the PDF's ###
##################################

from ModelBuilders import *
ws = RooWorkspace('ws')

sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))

declareObservables(ws)
file = TFile('Bs2JpsiPhiTuple.root')
obs = ws.set('observables')
angles = ws.set('helicityangles')
obs.add( angles )
data = RooDataSet('data','data',file.Get('dataset'),obs,'t==t && m==m && ( m<5340 || m>5390 ) ')


# todo: make this a J/psi phi builder, so that we can also have a J/psi K* one ;-]
#       and we can tweak the initial values of the masses, set the right mass windows...
mb = MassPdfBuilder(ws,ws['m'],ws['mdau1'],ws['mdau2']) 

psipdf = mb.sigDau1Pdf()
bkgpdf = mb.bkgDau1Pdf()
ws.factory("SUM::dau1pdf(N_psibkg*%s,N_nonpsibkg*%s)"% ( mb.sigDau1Pdf().GetName(),mb.bkgDau1Pdf().GetName() ) )

####
c = TCanvas()
c.Divide(2,3)

### TODO: move these plots into the mass PDF builder...
c.cd(1)
mdau1pdf = ws['dau1pdf']
mdau1pdf.fitTo(data,ncpu)
for i in mdau1pdf.getParameters(data): i.setConstant(True)
f_dau1 = mb.dau1Obs().frame()
data.plotOn(f_dau1)
mdau1pdf.plotOn(f_dau1,sigcolor)
mdau1pdf.plotOn(f_dau1,nonpsicolor,RooFit.Components("*bkg*"))
f_dau1.Draw()


### Now that we have the J/psi fitted, let's make angle SPlots...
## TODO: move this into the mass PDF builder...
for i in mdau1pdf.getParameters(data) : i.setConstant( 'N_' not in i.GetName() )
splot = RooStats.SPlot("splotdata","splotdata",data,mdau1pdf,RooArgList(ws['N_psibkg'],ws['N_nonpsibkg']))
wdata = splot.GetSDataSet()

## next, we pick up the sigmat distributions for our three components...

sigmat = ws['sigmat']
sigmat.setBins(40)
p= sigmat.frame()
stdata = {}
stpdf = {}
c = TCanvas()
for (f,sample) in enumerate([ 'psibkg','nonpsibkg' ]):
      dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","N_%s_sw"%sample) 
      stdata[sample] = RooDataHist("sigmat_%s_data"%sample,"hist Err Per Ev",RooArgSet(sigmat),dataw)
      stpdf[sample] = ws.put(RooHistPdf("sigmat_%s"%sample,"sigmat_%s"%sample,RooArgSet(sigmat),stdata[sample])) # ,2)
      #stdata[sample].plotOn(p) # python has trouble finding the right plotOn...
      stpdf[sample].plotOn(p,RooFit.LineColor( [kRed,kBlue,kBlack][f]))
p.Draw()
      
trb = TimeResolutionBuilder(ws, ws['t'],ws['sigmat'])
btb = BkgTimePdfBuilder(ws, trb, stpdf)

# TODO: make a nice pythonic wrapper around ws operators...
ws.factory("PROD:psibkg(    %s, %s )"      % (mb.sigDau1Pdf().GetName(),btb.psibkgPdf().GetName()))
ws.factory("PROD:nonpsibkg( %s, %s )"% (mb.bkgDau1Pdf().GetName(),btb.nonpsibkgPdf().GetName()))
ws.factory("SUM::pdf(N_psibkg*psibkg,N_nonpsibkg*nonpsibkg)")
pdf = ws['pdf']
#pdf.getParameters(data).readFromFile('bkginitialvalues.txt')
### Maybe we should first fit the three time splots with the three components to get them 
### to a reasonable starting point??
#result = pdf.fitTo(data,ncpu,RooFit.Save(True),RooFit.Minos(false))
#pdf.getParameters(data).writeToFile('bkgfitresult.txt')
pdf.getParameters(data).readFromFile('bkgfitresult.txt')

t = ws['t']
mpsi = ws['mdau1']
t.setRange('negativeNonZeroTime',t.getMin(),-0.1)
t.setRange('negativeTime',t.getMin(),-0.3)
t.setRange('nonZeroTime',0.1,t.getMax())
t.setRange('largeTime',0.3,t.getMax())

c = TCanvas()
c.Divide(5,2)
#===========================================================================================================
comps = { 'psibkg'    : ( bkgcolor,dashed ) 
        , 'nonpsibkg' : ( nonpsicolor,dashed )
        }

#===========================================================================================================
_c1 = plot( c.cd(1),mpsi,data,pdf,comps
          , frameOpts = ( RooFit.Bins(30), )
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c2 = plot( c.cd(2), mpsi, data, pdf, comps
          , frameOpts = ( RooFit.Bins(20), RooFit.Title('mpsi (t<-0.3)') )
          , dataOpts  = ( ms, xes, RooFit.CutRange('negativeTime') )
          , pdfOpts   = ( lw, RooFit.ProjectionRange('negativeTime') )
          )
#===========================================================================================================
_c3 = plot( c.cd(3), mpsi, data, pdf, comps
          , frameOpts = ( RooFit.Bins(30), RooFit.Title('mpsi (t<-0.1)') )
          , dataOpts = ( ms,xes,RooFit.CutRange('negativeNonZeroTime') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('negativeNonZeroTime') )
          )
#===========================================================================================================
_c4 = plot( c.cd(4), mpsi, data, pdf, comps
          , frameOpts = ( RooFit.Bins(40), RooFit.Title('mpsi (t>0.1)') )
          , dataOpts = ( ms,xes,RooFit.CutRange('nonZeroTime') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('nonZeroTime') )
          )

#===========================================================================================================
_c5 = plot( c.cd(5), mpsi, data, pdf, comps
          , frameOpts = ( RooFit.Bins(20), RooFit.Title('mpsi (t>0.3)') )
          , dataOpts = ( ms,xes,RooFit.CutRange('largeTime') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('largeTime') )
          )
#===========================================================================================================
_c6 = plot( c.cd(6), sigmat, data, pdf, comps
          , frameOpts = None
          , dataOpts = ( ms,xes )
          , pdfOpts = ( lw, )
          )
#===========================================================================================================
c.cd(7)
_c7 = plot( c.cd(7), t, data, pdf, comps
          , frameOpts = ( RooFit.Range(-0.4,0.4), RooFit.Bins(100), RooFit.Title('proper time, full mass range') )
          , dataOpts = ( ms,xes,err )
          , pdfOpts = ( lw, )
          )
#===========================================================================================================
_c8 = plot( c.cd(8), t, data, pdf, comps
          , frameOpts = ( RooFit.Title('proper time, lower sideband'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('leftSideband') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('leftSideband') )
          , logy = True
          )
#===========================================================================================================
_c9 = plot( c.cd(9), t, data, pdf, comps
          , frameOpts = ( RooFit.Title('proper time, upper sideband'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('rightSideband') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('rightSideband') )
          , logy = True
          )
