from ROOT import *
gSystem.Load("libp2vv")
from RooFitDecorators import *

#RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Plotting))

ncpu = RooCmdArg( RooFit.NumCPU(1) )
sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))
lw = RooCmdArg(RooFit.LineWidth(2))
ms = RooCmdArg(RooFit.MarkerSize(0.7))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))
stash = []


def plotPdfComponents( pdf, obs, *args  ) :
    pdf.plotOn(obs,RooFit.Components("psibkg"), bkgcolor,dashed,lw)
    pdf.plotOn(obs,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
    pdf.plotOn(obs,RooFit.Components("sig"), sigcolor,dashed,lw)
    pdf.plotOn(obs,lw)

RooAbsPdf.plotComponenentOn = plotPdfComponents

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
c.cd(1)
_mpsi= mpsi.frame(RooFit.Bins(60))
data.plotOn(_mpsi,ms,xes)
#pdf.plotComponentsOn(_mpsi)
pdf.plotOn(_mpsi,RooFit.Components("psibkg"),bkgcolor,dashed,lw)
pdf.plotOn(_mpsi,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
pdf.plotOn(_mpsi,lw)
_mpsi.Draw()
c.Update()

#===========================================================================================================
c.cd(2)
_mpsi = mpsi.frame(RooFit.Bins(40),RooFit.Title('mpsi (t<-0.3)'))
data.plotOn(_mpsi,ms,xes,RooFit.CutRange("negativeTime"))
pdf.plotOn(_mpsi,RooFit.ProjectionRange("negativeTime"),lw)
pdf.plotOn(_mpsi,RooFit.Components("psibkg"),RooFit.ProjectionRange("negativeTime"), bkgcolor,dashed,lw)
pdf.plotOn(_mpsi,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange("negativeTime"), nonpsicolor,dashed,lw)
_mpsi.Draw() 
c.Update()

#===========================================================================================================
c.cd(3)
_mpsi = mpsi.frame(RooFit.Bins(40),RooFit.Title('mpsi (t<-0.1)'))
data.plotOn(_mpsi,ms,xes,RooFit.CutRange("negativeNonZeroTime"))
pdf.plotOn(_mpsi,RooFit.ProjectionRange("negativeNonZeroTime"),lw)
pdf.plotOn(_mpsi,RooFit.Components("psibkg"),RooFit.ProjectionRange("negativeNonZeroTime"), bkgcolor,dashed,lw)
pdf.plotOn(_mpsi,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange("negativeNonZeroTime"), nonpsicolor,dashed,lw)
_mpsi.Draw() 
c.Update()

#===========================================================================================================
c.cd(4)
_mpsi = mpsi.frame(RooFit.Bins(40),RooFit.Title('mpsi (t>0.1)'))
data.plotOn(_mpsi,ms,xes,RooFit.CutRange("nonZeroTime"))
pdf.plotOn(_mpsi,RooFit.ProjectionRange("nonZeroTime"),lw)
pdf.plotOn(_mpsi,RooFit.Components("psibkg"),RooFit.ProjectionRange("nonZeroTime"), bkgcolor,dashed,lw)
pdf.plotOn(_mpsi,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange("nonZeroTime"), nonpsicolor,dashed,lw)
_mpsi.Draw() 
c.Update()

#===========================================================================================================
c.cd(5)
_mpsi = mpsi.frame(RooFit.Bins(20),RooFit.Title('mpsi (t>0.3)'))
data.plotOn(_mpsi,RooFit.CutRange("largeTime"),ms,xes)
#pdf.plotComponentsOn(_mpsi,RooFit.ProjectionRange("largeTime"))
pdf.plotOn(_mpsi,RooFit.ProjectionRange("largeTime"),RooFit.Components("psibkg"),    bkgcolor,dashed,lw)
pdf.plotOn(_mpsi,RooFit.ProjectionRange("largeTime"),RooFit.Components("nonpsibkg"), nonpsicolor,dashed,lw)
pdf.plotOn(_mpsi,RooFit.ProjectionRange("largeTime"),                                lw)
_mpsi.Draw() 
c.Update()

#===========================================================================================================
c.cd(6)
_sigmat = sigmat.frame()
data.plotOn(_sigmat,RooFit.MarkerSize(0.7),xes)
pdf.plotOn(_sigmat,RooFit.Components("psibkg"), bkgcolor,dashed,lw)
pdf.plotOn(_sigmat,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
pdf.plotOn(_sigmat,lw)
_sigmat.Draw() 
c.Update()

#===========================================================================================================
c.cd(7)
_tb = t.frame(-0.4,0.4,100)
data.plotOn(_tb,RooFit.MarkerSize(0.7),xes,err)
pdf.plotOn(_tb,RooFit.Components("psibkg"),bkgcolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
pdf.plotOn(_tb,lw)
_tb.SetMinimum(0.1) 
_tb.SetTitle("proper time full mass range")
_tb.Draw()
c.Update()

#===========================================================================================================
c.cd(8)
gPad.SetLogy()

_tb = t.frame()
data.plotOn(_tb,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange('leftSideband'), nonpsicolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("psibkg"),RooFit.ProjectionRange('leftSideband'), bkgcolor,dashed,lw)
pdf.plotOn(_tb,lw,RooFit.ProjectionRange('leftSideband'))

_tb.SetMinimum(0.1) 
_tb.SetTitle("proper time left sideband")
_tb.Draw()
c.Update()
#===========================================================================================================
c.cd(9)
gPad.SetLogy()

_tb = t.frame()
data.plotOn(_tb,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange('rightSideband'), nonpsicolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("psibkg"),RooFit.ProjectionRange('rightSideband'), bkgcolor,dashed,lw)
pdf.plotOn(_tb,lw,RooFit.ProjectionRange('rightSideband'))

_tb.SetMinimum(0.1) 
_tb.SetTitle("proper time right sideband")
_tb.Draw()
c.Update()


