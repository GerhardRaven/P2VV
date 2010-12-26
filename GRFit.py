from ROOT import *
from RooFitDecorators import *
gSystem.Load("libp2vv")
from math import pi

#######################
### Plot ICHEP Like ###
#######################

sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))
stash = []

def plotPdfComponents( pdf, obs, *args  ) :
    pdf.plotOn(obs,RooFit.Components("psibkg"), bkgcolor,dashed,lw)
    pdf.plotOn(obs,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
    pdf.plotOn(obs,RooFit.Components("sig"), sigcolor,dashed,lw)
    pdf.plotOn(obs,lw)


def plotAngles(data,pdf,angles ) :
    c2 = TCanvas('Angles','Angles',900,700)
    c2.Divide(3,4)
    for (i,a) in enumerate( angles ) :
        #==========================================================================================================
        c2.cd(1+i)
        _f = a.frame(RooFit.Bins(15),RooFit.Title('%s full mass region' % a.GetTitle() ))
        data.plotOn(_f)
        pdf.plotOn(_f,RooFit.Components('sig'),sigcolor,dashed,lw)
        pdf.plotOn(_f,RooFit.Components('bkg'),bkgcolor,dashed,lw)
        pdf.plotOn(_f,lw)
        _f.Draw()
        c2.Update()
        #==========================================================================================================
        c2.cd(4+i)
        _f = a.frame(RooFit.Bins(15),RooFit.Title('%s signal region' % a.GetTitle()))
        data.plotOn(_f,RooFit.CutRange('sigRegion'))
        pdf.plotOn(_f,RooFit.ProjectionRange('sigRegion'),RooFit.Components('sig'),sigcolor,dashed,lw)
        pdf.plotOn(_f,RooFit.ProjectionRange('sigRegion'),RooFit.Components('bkg'),bkgcolor,dashed,lw)
        pdf.plotOn(_f,RooFit.ProjectionRange('sigRegion'),lw)
        _f.Draw()
        c2.Update()
        #==========================================================================================================
        c2.cd(7+i)
        _f = a.frame(RooFit.Bins(15),RooFit.Title('%s sidebands'% a.GetTitle()))
        data.plotOn(_f,RooFit.CutRange('leftSideband,rightSideband'))
        pdf.plotOn(_f,RooFit.ProjectionRange('leftSideband,rightSideband'),lw)
        _f.Draw()
        c2.Update()
        #==========================================================================================================
        c2.cd(10+i)
        _f = a.frame(RooFit.Bins(15),RooFit.Title('%s right sideband'% a.GetTitle() ))
        data.plotOn(_f,RooFit.CutRange('rightSideband'))
        pdf.plotOn(_f,RooFit.ProjectionRange('rightSideband'),lw)
        _f.Draw()
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

sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))

declareObservables(ws)
file = TFile('Bs2JpsiPhiTuple.root')
obs = ws.set('observables')
angles = ws.set('helicityangles')
obs.add( angles )
data = RooDataSet('data','data',file.Get('dataset'),obs,'t==t && m==m')


#data
#sigmam_sigdata  = RooDataSet("sigmam_sig","sigmam_sig",data,data.get(),"(t>0.3)& (5365-15)<m&m<(5365+15)")
#sigmam_bkgdata  = RooDataSet("sigmam_bkg","sigmam_bkg",data,data.get(),"(t<0.3)& ( m<(5365-20) | (5365+20<m) )
#sigmam_data = RooDataHist("sigmam_%s_data"%sample,"hist m Err Per Ev",RooArgSet(sigmam),data)
#sigmat_pdf = ws.put(RooHistPdf("sigmam_%s"%sample,"sigmat_%s"%sample,RooArgSet(sigmat),stdata[sample])) # ,2)

# todo: make this a J/psi phi builder, so that we can also have a J/psi K* one ;-]
#       and we can tweak the initial values of the masses, set the right mass windows...
mb = MassPdfBuilder(ws,ws['m'],ws['mdau1'],ws['mdau2']) 
                                            

####
c = TCanvas()
c.Divide(2,3)

### TODO: move these plots into the mass PDF builder...
if True :
    c.cd(1)
    mdau1pdf = mb.dau1Pdf()
    mdau1pdf.fitTo(data,RooFit.NumCPU(7))
    for i in mdau1pdf.getParameters(data): i.setConstant(True)
    f_dau1 = mb.dau1Obs().frame()
    data.plotOn(f_dau1)
    mdau1pdf.plotOn(f_dau1,sigcolor)
    mdau1pdf.plotOn(f_dau1,nonpsicolor,RooFit.Components("*bkg*"))
    f_dau1.Draw()

if True:
    c.cd(2)
    mdau2pdf = mb.dau2Pdf()
    mdau2pdf.fitTo(data,RooFit.NumCPU(7))
    for i in mdau2pdf.getParameters(data) : i.setConstant(True)
    f_dau2 = mb.dau2Obs().frame()
    data.plotOn(f_dau2)
    mdau2pdf.plotOn(f_dau2,sigcolor)
    mdau2pdf.plotOn(f_dau2,nonpsicolor,RooFit.Components("*bkg*"))
    f_dau2.Draw()

if True :
    mpdf = mb.Pdf()
    ws['N_sig'].setVal(1000)
    ws['N_psibkg'].setVal( (data.numEntries()- ws['N_sig'].getVal())/2 )
    ws['N_nonpsibkg'].setVal( (data.numEntries()- ws['N_sig'].getVal())/2 )
    mpdf.fitTo(data,RooFit.NumCPU(7))
    for i in mpdf.getParameters(data) : i.setConstant(True)
    ## TODO: at some point, we need to add the KK mass to the story...
    for i,obs in enumerate( [ mb.Obs(), mb.dau1Obs()])  : #  , mb.dau2Obs() ] ):
        # TODO: plot current in signal of other, sideband of other....
        c.cd(3+i)
        frame = obs.frame()
        data.plotOn(frame)
        mpdf.plotOn(frame,RooFit.Components("m_sig"), sigcolor,dashed,lw)
        mpdf.plotOn(frame,RooFit.Components("m_psibkg"), bkgcolor,dashed,lw)
        mpdf.plotOn(frame,RooFit.Components("m_nonpsibkg"), nonpsicolor,dashed,lw)
        mpdf.plotOn(frame,lw)
        frame.Draw()
#import sys
#sys.exit(0)

c.Print("InitialMassFits.eps")

### Now that we have the masses fitted, let's make angle SPlots...
## TODO: move this into the mass PDF builder...
for i in mpdf.getParameters(data) : i.setConstant( i not in mb.yields() )
mpdf.Print("T")
mpdf.getParameters(data).Print("V")
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
      #stdata[sample].plotOn(p) # python has trouble finding the right plotOn...
      stpdf[sample].plotOn(p,RooFit.LineColor( [kRed,kBlue,kBlack][f]))
p.Draw()
      
x = TimeResolutionBuilder(ws)
y = BkgTimePdfBuilder(ws, x, stpdf)

if True :
    # create a dummy signal PDF...
    ws.factory("PROD:t_sig( Decay(t,t_sig_tau[1.5,1.2,1.8],tres_sig,SingleSided)|sigmat,sigmat_sig)")
    ### and now multiply and sum together...
    ws.factory("PROD:sig(       m_sig,       t_sig        )")
    ws.factory("PROD:psibkg(    m_psibkg,    t_psibkg     )")
    ws.factory("PROD:nonpsibkg( m_nonpsibkg, t_nonpsibkg  )")

    # TODO: create an intermediate bkg which is SUM(fpsi*psibkg,nonpsibkg)
    # as we can then make nicer plots -- alternative, we just add psi and nonpsi
    # by hand... which may be faster anyway ;-)
else :
    # and now, build the angular distributions for psi and non-psi background, using the SWeights we
    # just got...
    ab = abasis(ws,angles)
    x = BkgAnglePdfBuilder( ws, ab, data
                          , { 'psibkg'    : { 'ranges' : (range(4),range(8),range(-4,5))  , 'weight' : 'N_psibkg_sw'    }
                            , 'nonpsibkg' : { 'ranges' : (range(4),range(12),range(-4,5)) , 'weight' : 'N_nonpsibkg_sw' } 
                            } )
    (c1,c2) = x.makeplots()
    ws.factory("wtag[0.5]")
    ws.factory("{S[0],D[1],C[0]}")
    ws.factory("{t_sig_tau[1.4,1.2,1.8],t_sig_dm[17.7],t_sig_dG[0.05,-0.3,0.3]}")
    definePolarAngularAmplitudes(ws)
    print 'about to build J/psi phi signal'
    sigjpsiphi = buildJpsiphi(ws,'jpsiphisignal',False)
    print 'about to createProjection of J/psi phi signal over tagdecision'
    t_sig = sigjpsiphi # ws.put(sigjpsiphi.createProjection( ws.argSet('tagdecision') ))
    print 'about to multiply signal'
    ws.factory("PROD:ta_sig( %s |sigmat,sigmat_sig)" % t_sig.GetName() )
    ### and now multiply and sum together...
    print 'about to multiply background'
    ws.factory("PROD:sig(       m_sig,       ta_sig         )")
    ws.factory("PROD:psibkg(    m_psibkg,    t_psibkg,    %s )" % x.psibkgPdf().GetName())
    ws.factory("PROD:nonpsibkg( m_nonpsibkg, t_nonpsibkg, %s )" % x.nonpsibkgPdf().GetName())
    ws['N_sig'].setVal(1000)
    ws['N_psibkg'].setVal( (data.numEntries()- ws['N_sig'].getVal())/2 )
    ws['N_nonpsibkg'].setVal( (data.numEntries()- ws['N_sig'].getVal())/2 )

ws.factory("SUM::pdf(N_sig*sig,N_psibkg*psibkg,N_nonpsibkg*nonpsibkg)")
pdf = ws['pdf']
pdf.getParameters(data).readFromFile('initialvalues.txt')
### Maybe we should first fit the three time splots with the three components to get them 
### to a reasonable starting point??
result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Save(True),RooFit.Minos(false))
pdf.getParameters(data).writeToFile('fitresult.txt')





t = ws['t']
m = ws['m']
t.setRange('largeTime',0.3,t.getMax())

c = TCanvas()
c.Divide(4,2)

#===========================================================================================================
c.cd(1)
mpsiplot = ws.var('mdau1').frame(RooFit.Bins(50))
data.plotOn(mpsiplot,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mpsiplot,RooFit.Components("sig"),sigcolor,dashed,lw)
pdf.plotOn(mpsiplot,RooFit.Components("psibkg"),bkgcolor,dashed,lw)
pdf.plotOn(mpsiplot,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
pdf.plotOn(mpsiplot,lw)
mpsiplot.Draw()
c.Update()

#===========================================================================================================
c.cd(2)
_m = m.frame(RooFit.Bins(50),RooFit.Title('m'))
data.plotOn(_m,RooFit.MarkerSize(0.7),xes)
pdf.plotOn(_m,RooFit.Components("psibkg"),bkgcolor,dashed,lw)
pdf.plotOn(_m,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
pdf.plotOn(_m,RooFit.Components("sig"),sigcolor,dashed,lw)
pdf.plotOn(_m,lw)
_m.Draw() 
c.Update()

#===========================================================================================================
c.cd(3)
_m = m.frame(RooFit.Bins(50),RooFit.Title('m (t>0.3)'))
data.plotOn(_m,RooFit.MarkerSize(0.7),xes,RooFit.CutRange("largeTime"))
pdf.plotOn(_m,RooFit.Components("psibkg"),RooFit.ProjectionRange("largeTime"), bkgcolor,dashed,lw)
pdf.plotOn(_m,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange("largeTime"), nonpsicolor,dashed,lw)
pdf.plotOn(_m,RooFit.Components("sig"),RooFit.ProjectionRange("largeTime"), sigcolor,dashed,lw)
pdf.plotOn(_m,lw,RooFit.ProjectionRange("largeTime"))
_m.Draw() 
c.Update()

#===========================================================================================================
c.cd(5)
_tb = t.frame(-0.4,0.4,100)
data.plotOn(_tb,RooFit.MarkerSize(0.5),xes,err)
pdf.plotOn(_tb,RooFit.Components("sig"),sigcolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("psibkg"),bkgcolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("nonpsibkg"),nonpsicolor,dashed,lw)
pdf.plotOn(_tb,lw)
_tb.SetMinimum(0.1) 
_tb.SetTitle("proper time full mass range")
_tb.Draw()
c.Update()

#===========================================================================================================
c.cd(6)
gPad.SetLogy()

_tb = t.frame(-4,12.,50)
data.plotOn(_tb,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,RooFit.Components("sig"),RooFit.ProjectionRange('sigRegion'), sigcolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("psibkg"),RooFit.ProjectionRange('sigRegion'), bkgcolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange('sigRegion'), nonpsicolor,dashed,lw)
pdf.plotOn(_tb,lw,RooFit.ProjectionRange('sigRegion'))

_tb.SetMinimum(0.1) 
_tb.SetTitle("proper time signal region")
_tb.Draw()
c.Update()


#===========================================================================================================
c.cd(7)
gPad.SetLogy()

_tb = t.frame(-4,12,50)
data.plotOn(_tb,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange('leftSideband'), nonpsicolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("psibkg"),RooFit.ProjectionRange('leftSideband'), bkgcolor,dashed,lw)
pdf.plotOn(_tb,lw,RooFit.ProjectionRange('leftSideband'))

_tb.SetMinimum(0.1) 
_tb.SetTitle("proper time left sideband")
_tb.Draw()
c.Update()
#===========================================================================================================
c.cd(8)
gPad.SetLogy()

_tb = t.frame(-4,12,50)
data.plotOn(_tb,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,RooFit.Components("nonpsibkg"),RooFit.ProjectionRange('rightSideband'), nonpsicolor,dashed,lw)
pdf.plotOn(_tb,RooFit.Components("psibkg"),RooFit.ProjectionRange('rightSideband'), bkgcolor,dashed,lw)
pdf.plotOn(_tb,lw,RooFit.ProjectionRange('rightSideband'))

_tb.SetMinimum(0.1) 
_tb.SetTitle("proper time right sideband")
_tb.Draw()
c.Update()




if False :

    # pick up the psi and nonpsi mass PDFs and multiply them by their angles...

    ws.factory("PROD::psibkg(%s,%s)"%(mb.psibkgPdf().GetName(),abf.psibkgPdf().GetName() ) )
    ws.factory("PROD::nonpsibkg(%s,%s)"%(mb.nonpsibkgPdf().GetName(),abf.nonpsibkgPdf().GetName() ))
    ws.factory("PROD::sig(%s,%s)"%(mb.sigPdf().GetName(),abf.dummysigPdf().GetName() ))
    ws.factory("SUM::bkg(fpsi[0.5,0.01,0.99]*psibkg,nonpsibkg)")
    ws.factory("SUM::pdf(Nsig[12000,0,1000]*sig,Nbkg[81200,0,1000000]*bkg)")

    pdf = ws['pdf']
    pdf.fitTo(data,RooFit.NumCPU(8))
    plotAngles( data,pdf, angles )

