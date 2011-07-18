########################################
### Author: Daan van Eijk
### Updated on: Jun 5 11
### Description: This script reads the root file with the workspace for the tagged fit and fits
###              The MakeProfile function is implemented in RooFitDecorators.py, but it makes a DLL profile in one go.
###              For big grids this becomes too time-consuming. So there are now scripts to make ganga jobs for profile production:
###              These are TaggedProfiles.py and SubmitTaggedProfiles.py
########################################

from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi
from array import array

from RooFitDecorators import *

#import rootStyle
#from ROOT import (gROOT,gStyle,TStyle)
#myStyle = rootStyle.plainstyle()
#gROOT.SetStyle(myStyle.GetName())
#gROOT.ForceStyle()
#gStyle.UseCurrentStyle()

angcorr = False
taggingsyst = True
blinded = True
phisparam = True

name = 'BsJpsiPhi2011_fit_phis'

wsfile = TFile('TaggedWS.root')
ws = wsfile.Get('ws')

if taggingsyst:
    if not angcorr:
        pdfbeforeblinding = ws.pdf('pdf_ext_inc_tag_syst')
    if angcorr:
        pdfbeforeblinding = ws.pdf('pdf_ext_angcorrpdf_inc_tag_syst')
    
else:
    pdfbeforeblinding = ws.pdf('pdf_ext_angcorrpdf')

data = ws.data('data')

################
### Blinding ###
################

if blinded:
    #Building blinded parameters
    if phisparam:
        ws.factory("RooUnblindUniform::t_sig_dG_blind('BsRooBarb',0.2,t_sig_dG)")
        ws.factory("RooUnblindUniform::phis_blind('BsCustard',3.,phis)")    

        #calling RooCustomizer
        customizer = RooCustomizer(pdfbeforeblinding,'blinded')
        customizer.replaceArg( ws['t_sig_dG'], ws['t_sig_dG_blind'] )
        customizer.replaceArg( ws['phis'], ws['phis_blind'] )

    else:
        ws.factory("RooUnblindUniform::t_sig_dG_blind('BsRooBarb',0.2,t_sig_dG)")
        ws.factory("RooUnblindUniform::S_blind('BsMickeyMouseLP2011',2.,S)")
        #ws.factory("RooUnblindUniform::D_blind('BsMiniMouseLP2011',2.,D)")
        #ws.factory("RooUnblindUniform::C_blind('BsDonaldDuckLP2011',2.,C)")

        #calling RooCustomizer
        customizer = RooCustomizer(pdfbeforeblinding,'blinded')
        customizer.replaceArg( ws['t_sig_dG'], ws['t_sig_dG_blind'] )
        customizer.replaceArg( ws['S'], ws['S_blind'] )
        #customizer.replaceArg( ws['D'], ws['D_blind'] )
        #customizer.replaceArg( ws['C'], ws['C_blind'] )

    blindedpdf = customizer.build()
    ws.put(blindedpdf)

if blinded:
    pdf = blindedpdf
else:
    pdf = pdfbeforeblinding


#Set boundaries of some variables wider for blinding
if phisparam:
    ws.var('t_sig_dG').setVal(0.060)
    ws.var('t_sig_dG').setMin(-4)
    ws.var('t_sig_dG').setMax(4)
    
    ws.var('phis').setVal(-1.50)
    ws.var('phis').setMin(-2*pi)
    ws.var('phis').setMax(2*pi)
else:
    ws.var('t_sig_dG').setVal(0.060)
    ws.var('t_sig_dG').setMin(-4)
    ws.var('t_sig_dG').setMax(4)
    
    #ws.var('S').setVal(-0.04)
    #ws.var('S').setMin(-4)
    #ws.var('S').setMax(4)

    #ws.var('D').setVal(1)
    #ws.var('D').setMin(-4)
    #ws.var('D').setMax(4)


#ws['rs2'].setVal(0.)
#ws['rs2'].setConstant()
#ws['deltas'].setConstant()

#Constrain deltams
ws.factory("Gaussian::dmsconstraint(t_sig_dm,t_sig_dm_mean[17.77],t_sig_dm_sigma[0.12])")

if taggingsyst:
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))
else:
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true),RooFit.ExternalConstraints(RooArgSet(ws.pdf('dmsconstraint'))))

######PLOT######
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

signame = 'sig_pdf_angcorrpdf_inc_tag_syst'
bkgname = 'bkg_pdf_angcorrpdf'

sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))

msigmin = 5345.
msigmax = 5387.
mmin = 5200.
mmax = 5550.
ws['m'].setRange('sigRegion',msigmin,msigmax)
ws['m'].setRange('leftSideband',mmin,msigmin)
ws['m'].setRange('rightSideband',msigmax,mmax)

#Make Sanity plots
C2 = TCanvas('C2','C2')
C2.Divide(3,2)

C2.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw)
pdf.plotOn(tframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(tframe,RooFit.Components(signame),sigcolor,dashed,lw)
tframe.Draw()

C2.cd(2)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw)
pdf.plotOn(mframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(mframe,RooFit.Components(signame),sigcolor,dashed,lw)
mframe.Draw()

C2.cd(4)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(30))
data.plotOn(trcosthetaframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcosthetaframe,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(signame),sigcolor,dashed,lw)
trcosthetaframe.Draw()

C2.cd(5)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(30))
data.plotOn(trcospsiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcospsiframe,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(signame),sigcolor,dashed,lw)
trcospsiframe.Draw()

C2.cd(6)
trphiframe = ws['trphi'].frame(RooFit.Bins(30))
data.plotOn(trphiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trphiframe,lw)
pdf.plotOn(trphiframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(trphiframe,RooFit.Components(signame),sigcolor,dashed,lw)
trphiframe.Draw()

C3 = TCanvas('C3','C3')
C3.Divide(3,2)

C3.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(tframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(tframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
tframe.Draw()

C3.cd(2)
mframe = ws['m'].frame(RooFit.Range(msigmin,msigmax),RooFit.Bins(30))
data.plotOn(mframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(mframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(mframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
mframe.Draw()

C3.cd(4)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(30))
data.plotOn(trcosthetaframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcosthetaframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(trcosthetaframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
trcosthetaframe.Draw()

C3.cd(5)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(30))
data.plotOn(trcospsiframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcospsiframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(trcospsiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
trcospsiframe.Draw()

C3.cd(6)
trphiframe = ws['trphi'].frame(RooFit.Bins(30))
data.plotOn(trphiframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trphiframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(trphiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(trphiframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
trphiframe.Draw()

writeFitParamsLatex(result,name,False)
dict = writeCorrMatrixLatex(result,name)

assert False
################
### Profiles ###
################

#setting back values
#ws.var('#Gamma').setVal(0.68)
#ws.var('t_sig_dG').setVal(0.060)
#ws.var('phis').setVal(0.0)

phis = ws.var('phis')
deltaGamma = ws.var('t_sig_dG')

#MakeProfile('ProfiledGamma_phis_tagged',data,pdf,15,phis,0,2*pi,deltaGamma,-0.7,0.7)

etataghist = RooDataHist('etataghist','etataghist',RooArgSet(ws.var('etatag')),data)
etatagpdf = RooHistPdf('etatagpdf','etatagpdf',RooArgSet(ws.var('etatag')),etataghist)
etatagframe = ws.var('etatag').frame(RooFit.Bins(20))
data.plotOn(etatagframe)
etatagpdf.plotOn(etatagframe)
etatagframe.Draw()

data.table(ws.cat('tagdecision')).Print('v')
