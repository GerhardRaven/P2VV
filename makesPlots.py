#######################################
### Author: Daan van Eijk
### Updated on: Jun 8 11
### Description: This script shows how to get sPlots weighted data from the original data.
###              I implemented a new (simpler) mass model, the one that is also used in Note 6 (single gaussian, exp. bkg)
###              This is now in ModelBuilders as mpdf_fit()
########################################
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

#angles
ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f]}"%(-pi,pi))
#time
ws.factory("{t[-2,20.],sigmat[0.005,0.1}")

#tagging
ws.factory("{etatag[0.,0.,0.501], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}")

#masses
ws.factory("{ m[%f,%f], mdau1[%f,%f], mdau2[%f,%f]}"%(5200,5550,3097-60,3097+60,1019.455-12,1019.455+12))

#read data (all observables you want sPlots for)

datafile = TFile('/data/bfys/dveijk/DataJpsiPhi/Bs2JpsiPhiForTaggedFit.root')
NTupletree = datafile.Get('MyTree')

ws.defineSet("observables","trcospsi,trcostheta,trphi,t,sigmat,tagdecision,etatag,m,mdau1,mdau2")

#Make the time cut to supress bkg
t = ws.var('t')
tmin = 0.3
tmax = 14.
t.setRange(tmin,tmax)

data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m==m && trcospsi==trcospsi && trcostheta==trcostheta && trphi==trphi && tagdecision==tagdecision && etatag==etatag')
ws.put(data)

# Build mass pdf
mb = MassPdfBuilder(ws,ws['m'],ws['mdau1'],ws['mdau2'],'Bs2Jpsiphi')
m_pdf_fit = mb.mpdf_fit()

#Fit yields
print 'These are the initiated yields:', mb.yields_fit().Print('v')
m_pdf_fit.fitTo(data)
print 'These are the fitted yields:', mb.yields_fit().Print('v')

masscanvas = TCanvas('masscanvas','masscanvas')
mframe = ws['m'].frame()
data.plotOn(mframe)
m_pdf_fit.plotOn(mframe)
mframe.Draw()

# Make sPlots
for i in m_pdf_fit.getParameters(data) : i.setConstant( i not in mb.yields_fit() )
#m_pdf_fit.Print("T")
#m_pdf_fit.getParameters(data).Print("V")
splot = RooStats.SPlot("splotdata","splotdata",data,m_pdf_fit,mb.yields_fit())
#wdata = splot.GetSDataSet()

if True :
    c = TCanvas('all vars from sPlots','all vars from sPlots')
    c2 = TCanvas('zoom t from sPlots','zoom t from sPlots')
    stash.append(c)
    stash.append(c2)
    observables = RooArgSet(ws['trcospsi'],ws['trcostheta'],ws['trphi'])
    observables.add(ws['t'])
    observables.add(ws['sigmat'])
    observables.add(ws['etatag'])
    #observables.add(ws['tagdecision'])
    observables.add(ws['mdau1'])
    observables.add(ws['mdau2'])
    c.Divide(len(observables),2)
    c2.Divide(1,2)
    for (i,sample) in enumerate([ j.GetName() for j in mb.yields_fit() ]):
        dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","%s_sw"%sample) # need a dummy cut, as passing a (const char*)0 is kind of difficult...
        # should only make splots for observables not used in fit...
        for (j,observable) in enumerate( observables ) :
            f = observable.frame(RooFit.Bins(10))
            dataw.plotOn(f, RooFit.DataError(RooAbsData.SumW2) )
            c.cd(i*len(observables)+j+1)
            f.Draw()
        c2.cd(1+i)
        f = ws['t'].frame(-2.,2.,100)
        dataw.plotOn(f)
        f.Draw()

#make sigmat RooHistPdf's
sigmat = ws['sigmat']
sigmat.setBins(40)
stdata = {}
stpdf = {}
c = TCanvas('sigmat from sPlots','sigmat from sPlots')
c.Divide(1,2)
for (f,sample) in enumerate([ j.GetName() for j in mb.yields_fit() ]):
      dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","%s_sw"%sample) 
      stdata[sample] = RooDataHist("sigmat_%s_data"%sample,"hist Err Per Ev",RooArgSet(sigmat),dataw)
      stpdf[sample] = ws.put(RooHistPdf("sigmat_%s"%sample,"sigmat_%s"%sample,RooArgSet(sigmat),stdata[sample])) # ,2)
      c.cd(f+1)
      p= sigmat.frame()
      dataw.plotOn(p) # python has trouble finding the right plotOn...
      stpdf[sample].plotOn(p,RooFit.LineColor( [kBlue,kRed,kBlack][f]))
      p.Draw()

ws['tagdecision'].setRange('B','Bs_Jpsiphi')
ws['tagdecision'].setRange('Bbar','Bsbar_Jpsiphi')
ws['tagdecision'].setRange('untagged','untagged')

#make etatag plots from data
etatag = ws['etatag']
etatag.setBins(40)

etatagcanvas = TCanvas('etatag from data','etatag from data')
etatagcanvas.Divide(4,1)

etatagcanvas.cd(1)
p = etatag.frame(RooFit.Title('etatag'))
data.plotOn(p)
etatagdata = RooDataHist("etatag_data","hist etatag Per Ev",RooArgSet(etatag),data)
etatagpdf = RooHistPdf("etatag_pdf","etatag_pdf",RooArgSet(etatag),etatagdata)
ws.put(etatagpdf)
etatagpdf.plotOn(p)
p.Draw()
etatagcanvas.cd(2)
p = etatag.frame(RooFit.Title('tagged as B'))
data.plotOn(p,RooFit.CutRange('B'))
p.Draw()
etatagcanvas.cd(3)
p = etatag.frame(RooFit.Title('tagged as Bbar'))
data.plotOn(p,RooFit.CutRange('Bbar'))
p.Draw()
etatagcanvas.cd(4)
p = etatag.frame(RooFit.Title('untagged'))
data.plotOn(p,RooFit.CutRange('untagged'))
p.Draw()

writefile = TFile('ws_etatagPDF.root','RECREATE')
ws.Write('ws')
writefile.Close()

#make etatag plots from sdata
datawsig = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","N_sig_sw") 
#etatagdatasig = RooDataHist("etatag_N_sig_data","hist etatag Per Ev",RooArgSet(etatag),datawsig)
#etatagpdfsig = RooHistPdf("etatag_N_sig","etatag_N_sig",RooArgSet(etatag),etatagdatasig)

datawbkg = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","N_bkg_sw") 
#etatagdatabkg = RooDataHist("etatag_N_bkg_data","hist etatag Per Ev",RooArgSet(etatag),datawbkg)
#etatagpdfbkg = RooHistPdf("etatag_N_bkg","etatag_N_bkg",RooArgSet(etatag),etatagdatabkg)

summarycanvas = TCanvas('etatag from sPlots','etatag from sPlots')
summarycanvas.Divide(4,2)

summarycanvas.cd(1)
p= etatag.frame(RooFit.Title('signal'))
datawsig.plotOn(p)
p.Draw()
summarycanvas.cd(2)
p= etatag.frame(RooFit.Title('signal, tagged as B'))
datawsig.plotOn(p,RooFit.CutRange("B"))
p.Draw()
summarycanvas.cd(3)
p= etatag.frame(RooFit.Title('signal, tagged as Bbar'))
datawsig.plotOn(p,RooFit.CutRange('Bbar'))
p.Draw()
summarycanvas.cd(4)
p= etatag.frame(RooFit.Title('signal, untagged'))
datawsig.plotOn(p,RooFit.CutRange('untagged'))
p.Draw()

summarycanvas.cd(5)
p= etatag.frame(RooFit.Title('bkg'))
datawbkg.plotOn(p)
p.Draw()
summarycanvas.cd(6)
p= etatag.frame(RooFit.Title('bkg, tagged as B'))
datawbkg.plotOn(p,RooFit.CutRange("B"))
p.Draw()
summarycanvas.cd(7)
p= etatag.frame(RooFit.Title('bkg, tagged as Bbar'))
datawbkg.plotOn(p,RooFit.CutRange("Bbar"))
p.Draw()
summarycanvas.cd(8)
p= etatag.frame(RooFit.Title('bkg, untagged'))
datawbkg.plotOn(p,RooFit.CutRange("untagged"))
p.Draw()
