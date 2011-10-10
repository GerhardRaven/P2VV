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

name = 'BsJpsiPhi2011_fit_phis_biased'
phisparam = True
angcorr = True

wsfile = TFile('BiasedWS.root')
ws = wsfile.Get('ws')

if angcorr:
    accpdfbeforeblinding = ws['accpdf_angcorr_inc_tag_syst']
else:
    accpdfbeforeblinding = ws['accpdf_inc_tag_syst']

################
### Blinding ###
################

#Building blinded parameters
if phisparam:
    ws.factory("RooUnblindUniform::t_sig_dG_blind('BsRooBarb',0.2,t_sig_dG)")
    ws.factory("RooUnblindUniform::phis_blind('BsCustard',3.,phis)")    
    #calling RooCustomizer
    customizer = RooCustomizer(accpdfbeforeblinding,'blinded')
    customizer.replaceArg( ws['t_sig_dG'], ws['t_sig_dG_blind'] )
    customizer.replaceArg( ws['phis'], ws['phis_blind'] )
    accblindedpdf = customizer.build()

#Constrain deltams
ws.factory("Gaussian::dmsconstraint(t_sig_dm,t_sig_dm_mean[17.77],t_sig_dm_sigma[0.12])")

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

pdf = accblindedpdf
data = ws['fullybiaseddata']

#ws['t_bkg_tau'].setMin(0.)
#ws['t_bkg_tau'].setMax(1.5)
#ws['t_bkg_tau'].setVal(0.5)

## #Check
## ws['m_pdf'].fitTo(data,RooFit.NumCPU(8))

######PLOT######
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

ws['rs2'].setVal(0.1)
#ws['rs2'].setConstant()
#ws['deltas'].setConstant()

## ws['m_bkg_exp'].setConstant()
## ws['m_sig_sigma_1'].setConstant()
## ws['m_sig_mean'].setConstant()
ws['Nbkg'].setMax(5000)
ws['Nsig'].setMax(5000)
sw = TStopwatch()

sw.Start()
#result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))
result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))
sw.Stop() 

signame = 'sigpdf_angcorr_inc_tag_syst_blinded'
bkgname = 'bkgpdf_angcorr'

#signame = 'sig_pdf__inc_tag_syst_blinded'
#bkgname = 'bkg_pdf'

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
C2.Divide(2,2)

C2.cd(1)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
ws['effpdf'].plotOn(tframe)
tframe.Draw()

C2.cd(2)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

C2.cd(3)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
pdf.plotOn(tframe,lw)
tframe.Draw()

C2.cd(4)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw)
#pdf.plotOn(tframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
#pdf.plotOn(tframe,RooFit.Components(signame),sigcolor,dashed,lw)
tframe.Draw()

#Make Sanity plots
C3 = TCanvas('C3','C3')
C3.Divide(3,1)

C3.cd(1)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
mframe.Draw()

C3.cd(2)
mframe = ws['m'].frame(RooFit.Bins(30))
pdf.plotOn(mframe,lw)
mframe.Draw()

C3.cd(3)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw)
#pdf.plotOn(mframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
#pdf.plotOn(mframe,RooFit.Components(signame),sigcolor,dashed,lw)
mframe.Draw()

assert False

writeFitParamsLatex(result,name,False)
dict = writeCorrMatrixLatex(result,name)
