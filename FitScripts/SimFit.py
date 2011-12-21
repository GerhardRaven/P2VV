########################################
### Author: Daan van Eijk
### Updated on: Jun 5 11
### Description: This script reads the root file with the workspace for the tagged fit and fits
###              The MakeProfile function is implemented in RooFitDecorators.py, but it makes a DLL profile in one go.
###              For big grids this becomes too time-consuming. So there are now scripts to make ganga jobs for profile production:
###              These are TaggedProfiles.py and SubmitTaggedProfiles.py
########################################

from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi
from array import array

from RooFitDecorators import *

nbins_timeacc = 100

name = 'SimFit_%s_bins'%(nbins_timeacc)

wsfile = TFile('SimWS_%s_bins.root'%(nbins_timeacc))

ws = wsfile.Get('ws')

#resultfile = TFile('result_%s.root'%(name))
#prelimresult =resultfile.Get('fitresult_simpdf_inc_tag_syst_jointdata')

blinded = False
angcorr = False

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

data = ws['jointdata']

#For non-extended fit
ws['#Gamma'].setVal(0.6761)
ws['t_sig_dG'].setVal(0.106)
ws['rperp2'].setVal(0.260)
ws['rz2'].setVal(0.484)
ws['rs2'].setVal(0.045)
ws['deltapar'].setVal(3.26)
ws['deltaperp'].setVal(2.61)
ws['deltas'].setVal(2.93)
ws['phis'].setVal(0.151)

#set = prelimresult.floatParsFinal()
#for i in set:
#    parname = i.GetName()
#    value = i.getVal()
#    error = i.getError()
#    ws[parname].setVal(value)
#    ws[parname].setError(error)
    
fit = True

if fit:
    sw = TStopwatch()
    sw.Start()
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'),ws.pdf('tres_SFconstraint'))))
    #result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'),ws.pdf('tres_SFconstraint'))))
    sw.Stop()

    result.SaveAs('result_%s.root'%(name))
    result.writepars(name,False)
    result.writecorr(name)

assert False

######PLOT######
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

signame_UB = 'sig_pdf_inc_tag_syst'
bkgname_UB = 'bkg_pdf_UB'
signame_B  = 'acc_sig_pdf_inc_tag_syst'
bkgname_B  = 'acc_bkg_pdf'

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

projectdata = RooDataSet('projectdata','projectdata',data,RooArgSet(ws['fitcat']))

#Make Sanity plots
masscanvas = TCanvas('masscanvas','masscanvas')
masscanvas.Divide(1,2)

masscanvas.cd(1)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.Cut("fitcat==fitcat::AllUnbiased"))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(signame_UB),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(bkgname_UB),RooFit.ProjWData(projectdata),dashed,bkgcolor) 
mframe.Draw()

masscanvas.cd(2)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.Cut("fitcat==fitcat::FullyBiased"))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(signame_B),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(bkgname_B),RooFit.ProjWData(projectdata),dashed,bkgcolor)
mframe.Draw()

masscanvas.SaveAs('plaatjemass.eps')

timecanvas = TCanvas('timecanvas','timecanvas')
timecanvas.Divide(1,2)

timecanvas.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.Cut("fitcat==fitcat::AllUnbiased"))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(signame_UB),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(bkgname_UB),RooFit.ProjWData(projectdata),dashed,bkgcolor) 
tframe.Draw()

timecanvas.cd(2)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.Cut("fitcat==fitcat::FullyBiased"))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(signame_B),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(bkgname_B),RooFit.ProjWData(projectdata),dashed,bkgcolor) 
tframe.Draw()

timecanvas.SaveAs('plaatjetime.eps')
