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

blinded = False
angcorr = False

name = 'SimFit'

wsfile = TFile('SimWS.root')
ws = wsfile.Get('ws')

if angcorr:
    pdfbeforeblinding = ws['simpdf_angcorr_inc_tag_syst']
else:
    pdfbeforeblinding = ws['simpdf_inc_tag_syst']

################
### Blinding ###
################

#Building blinded parameters
if blinded:
    ws.factory("RooUnblindUniform::t_sig_dG_blind('BsRooBarb',0.2,t_sig_dG)")
    ws.factory("RooUnblindUniform::phis_blind('BsCustard',3.,phis)")    

    #calling RooCustomizer
    customizer = RooCustomizer(pdfbeforeblinding,'blinded')
    customizer.replaceArg( ws['t_sig_dG'], ws['t_sig_dG_blind'] )
    customizer.replaceArg( ws['phis'], ws['phis_blind'] )
    blindedpdf = customizer.build()

    ws.put(blindedpdf)

    if angcorr:
        pdf = ws['simpdf_angcorr_inc_tag_syst_blinded']
    else:
        pdf = ws['simpdf_inc_tag_syst_blinded']

else:
    pdf = pdfbeforeblinding

ws['#Gamma'].setVal(0.6761)
ws['t_sig_dG'].setVal(0.106)
ws['rperp2'].setVal(0.260)
ws['rz2'].setVal(0.484)
ws['rs2'].setVal(0.045)
ws['deltapar'].setVal(3.26)
ws['deltaperp'].setVal(2.61)
ws['deltas'].setVal(2.93)
ws['phis'].setVal(0.151)
ws['Nsig_UB'].setVal(7280)
ws['Nbkg_UB'].setVal(3746)

data = ws['jointdata']

#ws['rs2'].setVal(0.)
#ws['rs2'].setConstant()
#ws['deltas'].setVal(0.)
#ws['deltas'].setConstant()

sw = TStopwatch()
sw.Start()
result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'),ws.pdf('tres_SFconstraint'))))
sw.Stop()

writeFitParamsLatex(result,name,False)
dict = writeCorrMatrixLatex(result,name)

#mcanvas = TCanvas('mcanvas','mcanvas')
#mframe = ws['m'].frame()
#data.plotOn(mframe)
#pdf.plotOn(mframe,RooFit.ProjWData(data))
#mframe.Draw()

assert False
######PLOT######
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

signame = 'accsig_inc_tag_syst_blinded'
bkgname = 'accbkg'

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
