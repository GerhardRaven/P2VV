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

import rootStyle
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

#from lhcbStyle import *
#gStyle.SetOptTitle(0)
#lhcbname = printLHCb(gStyle,'L','Prelim','#sqrt{s}=7 TeV Data, L = 370 pb^{-1}')
#style = lhcbStyle()

angcorr = True
taggingsyst = True
blinded = False
phisparam = True
splitmKK = True

name = 'LargePhiWindow'
wsfile = TFile('UnbiasedWS_LargePhiWindow.root')

#name = 'FitXC_angcorr_withbkg'
#wsfile = TFile('UnbiasedWS_withbkg.root')

ws = wsfile.Get('ws')

if taggingsyst:
    if not angcorr:
        pdfbeforeblinding = ws['pdf_ext_inc_tag_syst']
    if angcorr:
        pdfbeforeblinding = ws['pdf_ext_angcorr_inc_tag_syst']
    
else:
    pdfbeforeblinding = ws['pdf_ext_angcorr']

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
if blinded:
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

ws['rz2'].setVal(0.485)
ws['rz2'].setMin(0.)
ws['rz2'].setMax(1.)

ws['rs2'].setVal(0.05)
ws['rs2'].setMin(0.)
ws['rs2'].setMax(1.)
#ws['rs2'].setConstant()
#ws['deltas'].setConstant()

ws['rperp2'].setVal(0.259)
ws['rperp2'].setMin(0.)
ws['rperp2'].setMax(1.)

ws['deltas'].setVal(0.5)
ws['deltas'].setMin(-2*pi)
ws['deltas'].setMax(2*pi)

ws['deltaperp'].setVal(1.)
ws['deltaperp'].setMin(0.)
ws['deltaperp'].setMax(2.*pi)

ws['deltapar'].setVal(3.)
ws['deltapar'].setMin(-2.*pi)
ws['deltapar'].setMax(2.*pi)

ws['Nbkg'].setMax(50000)
ws['Nsig'].setMax(10000)

ws['phis'].setVal(0.1)
ws['phis'].setMin(-pi)
ws['phis'].setMax(pi)

#ws.var('t_sig_dG').setVal(0.1)

if splitmKK:
    #################################
    ### TEST OF SPLIT IN MKK BINS ###
    #################################

    mdau2bound1 = ws['mdau2'].getMin()
    mdau2bound2 = 1008
    mdau2bound3 = 1020
    mdau2bound4 = 1032
    mdau2bound5 = ws['mdau2'].getMax()
    mKKcat = RooThresholdCategory('mKKcat','mKKcat',ws['mdau2'],'Bin1') #Can't make it using factory: something goes wrong when trying to make the default value (LeftSideBand in this case)
    ws.put(mKKcat)
    ws['mKKcat'].addThreshold(mdau2bound2,"Bin1")
    ws['mKKcat'].addThreshold(mdau2bound3,"Bin2")
    ws['mKKcat'].addThreshold(mdau2bound4,"Bin3")
    ws['mKKcat'].addThreshold(mdau2bound5,"Bin4")

    ws['unbiaseddata'].table(ws['mKKcat']).Print('v')

    mdau2canvas = TCanvas('mdau2canvas','mdau2canvas')
    mdau2frame = ws['mdau2'].frame(RooFit.Bins(40))
    ws['unbiaseddata'].plotOn(mdau2frame)
    
    
    mdau2frame.Draw()

    print "That's great: I have to print mdau2canvas.GetUymax() first, otherwise the boundaries aren't drawn nicely...."
    print mdau2canvas.GetUymax()

    bin1rightline = TLine(mdau2bound2,0,mdau2bound2,mdau2canvas.GetUymax())
    bin1rightline.SetLineColor(2)
    bin2rightline = TLine(mdau2bound3,0,mdau2bound3,mdau2canvas.GetUymax())
    bin2rightline.SetLineColor(2)
    bin3rightline = TLine(mdau2bound4,0,mdau2bound4,mdau2canvas.GetUymax())
    bin3rightline.SetLineColor(2)

    bin1rightline.Draw()
    bin2rightline.Draw()
    bin3rightline.Draw()

    mdau2canvas.SaveAs('mdau2plot.eps')
    
    addcolumn = ws['unbiaseddata'].addColumn(ws['mKKcat'])
    addcolumn.SetName('thismKKcat')
    ws.put(addcolumn)

    #ws.factory("SIMCLONE::simpdf(pdf_ext_inc_tag_syst_blinded,$SplitParam({rs2,deltas,Nbkg,Nsig},thismKKcat))")
    ws.factory("SIMCLONE::simpdf(pdf_ext_inc_tag_syst,$SplitParam({rs2,deltas,Nbkg,Nsig},thismKKcat))")

    print 'Going to fit the Simultaneous PDF....'
    result = ws['simpdf'].fitTo(ws['unbiaseddata'],RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))

    ### Plot roocategory against parameters ###
    NmKKBins = 4
    xlist = [0,1,2,3]
    xerrlist = [0.,0.,0.,0.]
    x = array('f', xlist)
    xerr = array('f', xerrlist)

    rs2namelist = ['rs2_Bin1','rs2_Bin2','rs2_Bin3','rs2_Bin4']
    deltasnamelist = ['deltas_Bin1','deltas_Bin2','deltas_Bin3','deltas_Bin4']
    rs2list = []
    rs2errlist = []
    deltaslist = []
    deltaserrlist = []

    for i in range(0,NmKKBins):
        rs2list.append(ws[rs2namelist[i]].getVal())
        rs2errlist.append(ws[rs2namelist[i]].getError())
        deltaslist.append(ws[deltasnamelist[i]].getVal())
        deltaserrlist.append(ws[deltasnamelist[i]].getError())

    rs2y = array('f',rs2list)
    rs2yerr = array('f',rs2errlist)
    rs2_v_cat = TGraphErrors(NmKKBins, x, rs2y, xerr, rs2yerr);
    rs2_v_cat.SetTitle("|A_{s}|^{2} vs. mKK category")
    rs2_v_cat.SetMarkerStyle(20)
    rs2_v_cat.SetMarkerColor(4)
    rs2_v_cat.GetXaxis().SetTitle("mKK category")
    rs2_v_cat.GetYaxis().SetTitle("|A_{s}|^{2}")
    #rs2_v_cat.SetMaximum(1.3)

    deltasy = array('f',deltaslist)
    deltasyerr = array('f',deltaserrlist)
    deltas_v_cat = TGraphErrors(NmKKBins, x, deltasy, xerr, deltasyerr);
    deltas_v_cat.SetTitle("\delta_{s} vs. mKK category")
    deltas_v_cat.SetMarkerStyle(20)
    deltas_v_cat.SetMarkerColor(4)
    deltas_v_cat.GetXaxis().SetTitle("mKK category")
    deltas_v_cat.GetYaxis().SetTitle("\delta_{s}")
    #deltas_v_cat.SetMaximum(1.3)

    ParamCanvas = TCanvas('ParamCanvas','ParamCanvas')
    ParamCanvas.Divide(2)

    ParamCanvas.cd(1)
    rs2_v_cat.GetXaxis().SetLimits(-1,NmKKBins)
    rs2_v_cat.GetXaxis().SetBinLabel(20,'Bin1')
    rs2_v_cat.GetXaxis().SetBinLabel(40,'Bin2')
    rs2_v_cat.GetXaxis().SetBinLabel(60,'Bin3')
    rs2_v_cat.GetXaxis().SetBinLabel(80,'Bin4')
    rs2_v_cat.GetXaxis().LabelsOption('h')
    rs2_v_cat.Draw('AP')

    ParamCanvas.cd(2)
    deltas_v_cat.GetXaxis().SetLimits(-1,NmKKBins)
    deltas_v_cat.GetXaxis().SetBinLabel(20,'Bin1')
    deltas_v_cat.GetXaxis().SetBinLabel(40,'Bin2')
    deltas_v_cat.GetXaxis().SetBinLabel(60,'Bin3')
    deltas_v_cat.GetXaxis().SetBinLabel(80,'Bin4')
    deltas_v_cat.GetXaxis().LabelsOption('h')
    deltas_v_cat.Draw('AP')

    #Show sig/bkg per category:
    masscanvas = TCanvas('masscanvas','masscanvas')
    masscanvas.Divide(2,2)

    ws['mdau2'].setRange('Bin1',mdau2bound1,mdau2bound2)
    ws['mdau2'].setRange('Bin2',mdau2bound2,mdau2bound3)
    ws['mdau2'].setRange('Bin3',mdau2bound3,mdau2bound4)
    ws['mdau2'].setRange('Bin4',mdau2bound4,mdau2bound5)

    argset = RooArgSet(ws['thismKKcat'])
    #bin1 = 'pdf_ext_inc_tag_syst_blinded_Bin1'
    #bin2 = 'pdf_ext_inc_tag_syst_blinded_Bin2'
    #bin3 = 'pdf_ext_inc_tag_syst_blinded_Bin3'
    #bin4 = 'pdf_ext_inc_tag_syst_blinded_Bin4'

    #catdataset = RooDataSet('catdataset','catdataset',ws['unbiaseddata'],RooArgSet(ws['thismKKcat']))
    #ws.put(catdataset)

    masscanvas.cd(1)
    mframe1 = ws['m'].frame(RooFit.Bins(30))
    ws['unbiaseddata'].plotOn(mframe1,RooFit.CutRange('Bin1'))
    ws['simpdf'].plotOn(mframe1,RooFit.Slice(ws['thismKKcat'],'Bin1'),RooFit.ProjWData(argset,ws['unbiaseddata']))
    mframe1.Draw()

    masscanvas.cd(2)
    mframe2 = ws['m'].frame(RooFit.Bins(30))
    ws['unbiaseddata'].plotOn(mframe2,RooFit.CutRange('Bin2'))
    ws['simpdf'].plotOn(mframe2,RooFit.Slice(ws['thismKKcat'],'Bin2'),RooFit.ProjWData(argset,ws['unbiaseddata']))
    mframe2.Draw()

    masscanvas.cd(3)
    mframe3 = ws['m'].frame(RooFit.Bins(30))
    ws['unbiaseddata'].plotOn(mframe3,RooFit.CutRange('Bin3'))
    ws['simpdf'].plotOn(mframe3,RooFit.Slice(ws['thismKKcat'],'Bin3'),RooFit.ProjWData(argset,ws['unbiaseddata']))
    mframe3.Draw()

    masscanvas.cd(4)
    mframe4 = ws['m'].frame(RooFit.Bins(30))
    ws['unbiaseddata'].plotOn(mframe4,RooFit.CutRange('Bin4'))
    ws['simpdf'].plotOn(mframe4,RooFit.Slice(ws['thismKKcat'],'Bin4'),RooFit.ProjWData(argset,ws['unbiaseddata']))
    mframe4.Draw()
    
    ##################################################################################

pdf = pdf
data = ws['unbiaseddata']

if taggingsyst:
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws['p0'],ws['p1'],ws['dmsconstraint'],ws['tres_SFconstraint'])))
else:
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws['dmsconstraint'],ws['tres_SFconstraint'])))

writeFitParamsLatex(result,name,False)
dict = writeCorrMatrixLatex(result,name)


######PLOT######
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

signame = 'sig_pdf_inc_tag_syst_blinded'
bkgname = 'bkg_pdf'

sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))

msigmin = 5310.
msigmax = 5420.
mmin = 5200.
mmax = 5550.
ws['m'].setRange('sigRegion',msigmin,msigmax)
ws['m'].setRange('leftSideband',mmin,msigmin)
ws['m'].setRange('rightSideband',msigmax,mmax)

#Make Sanity plots
C2 = TCanvas('fullmasswindow','fullmasswindow')
C2.Divide(3,2)

C2.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw)
pdf.plotOn(tframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(tframe,RooFit.Components(signame),sigcolor,dashed,lw)
tframe.Draw()
lhcbname.Draw('same')

C2.cd(2)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw)
pdf.plotOn(mframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(mframe,RooFit.Components(signame),sigcolor,dashed,lw)
mframe.Draw()
lhcbname.Draw('same')

C2.cd(4)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(30))
data.plotOn(trcosthetaframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcosthetaframe,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(signame),sigcolor,dashed,lw)
trcosthetaframe.Draw()
lhcbname.Draw('same')

C2.cd(5)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(30))
data.plotOn(trcospsiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcospsiframe,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(signame),sigcolor,dashed,lw)
trcospsiframe.Draw()
lhcbname.Draw('same')

C2.cd(6)
trphiframe = ws['trphi'].frame(RooFit.Bins(30))
data.plotOn(trphiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trphiframe,lw)
pdf.plotOn(trphiframe,RooFit.Components(bkgname),bkgcolor,dashed,lw)
pdf.plotOn(trphiframe,RooFit.Components(signame),sigcolor,dashed,lw)
trphiframe.Draw()
lhcbname.Draw('same')

C2.SaveAs(name+'_fullmasswindow.eps')

C3 = TCanvas('signalregion','signalregion')
C3.Divide(3,2)

C3.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(tframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(tframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
tframe.Draw()
lhcbname.Draw('same')

C3.cd(2)
mframe = ws['m'].frame(RooFit.Range(msigmin,msigmax),RooFit.Bins(30))
data.plotOn(mframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(mframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(mframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
mframe.Draw()
lhcbname.Draw('same')

C3.cd(4)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(30))
data.plotOn(trcosthetaframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcosthetaframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(trcosthetaframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
trcosthetaframe.Draw()
lhcbname.Draw('same')

C3.cd(5)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(30))
data.plotOn(trcospsiframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcospsiframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(trcospsiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
trcospsiframe.Draw()
lhcbname.Draw('same')

C3.cd(6)
trphiframe = ws['trphi'].frame(RooFit.Bins(30))
data.plotOn(trphiframe,RooFit.CutRange('sigRegion'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trphiframe,lw,RooFit.ProjectionRange('sigRegion'))
pdf.plotOn(trphiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('sigRegion'),bkgcolor,dashed,lw)
pdf.plotOn(trphiframe,RooFit.Components(signame),RooFit.ProjectionRange('sigRegion'),sigcolor,dashed,lw)
trphiframe.Draw()
lhcbname.Draw('same')

C3.SaveAs(name+'_signalregion.eps')

C4 = TCanvas('leftsideband','leftsideband')
C4.Divide(3,2)

C4.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw,RooFit.ProjectionRange('leftSideband'))
pdf.plotOn(tframe,RooFit.Components(bkgname),RooFit.ProjectionRange('leftSideband'),bkgcolor,dashed,lw)
pdf.plotOn(tframe,RooFit.Components(signame),RooFit.ProjectionRange('leftSideband'),sigcolor,dashed,lw)
tframe.Draw()
lhcbname.Draw('same')

C4.cd(2)
mframe = ws['m'].frame(RooFit.Range(mmin,msigmin),RooFit.Bins(30))
data.plotOn(mframe,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw,RooFit.ProjectionRange('leftSideband'))
pdf.plotOn(mframe,RooFit.Components(bkgname),RooFit.ProjectionRange('leftSideband'),bkgcolor,dashed,lw)
pdf.plotOn(mframe,RooFit.Components(signame),RooFit.ProjectionRange('leftSideband'),sigcolor,dashed,lw)
mframe.Draw()
lhcbname.Draw('same')

C4.cd(4)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(30))
data.plotOn(trcosthetaframe,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcosthetaframe,lw,RooFit.ProjectionRange('leftSideband'))
pdf.plotOn(trcosthetaframe,RooFit.Components(bkgname),RooFit.ProjectionRange('leftSideband'),bkgcolor,dashed,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(signame),RooFit.ProjectionRange('leftSideband'),sigcolor,dashed,lw)
trcosthetaframe.Draw()
lhcbname.Draw('same')

C4.cd(5)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(30))
data.plotOn(trcospsiframe,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcospsiframe,lw,RooFit.ProjectionRange('leftSideband'))
pdf.plotOn(trcospsiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('leftSideband'),bkgcolor,dashed,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(signame),RooFit.ProjectionRange('leftSideband'),sigcolor,dashed,lw)
trcospsiframe.Draw()
lhcbname.Draw('same')

C4.cd(6)
trphiframe = ws['trphi'].frame(RooFit.Bins(30))
data.plotOn(trphiframe,RooFit.CutRange('leftSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trphiframe,lw,RooFit.ProjectionRange('leftSideband'))
pdf.plotOn(trphiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('leftSideband'),bkgcolor,dashed,lw)
pdf.plotOn(trphiframe,RooFit.Components(signame),RooFit.ProjectionRange('leftSideband'),sigcolor,dashed,lw)
trphiframe.Draw()
lhcbname.Draw('same')

C4.SaveAs(name+'_leftSideband.eps')

C5 = TCanvas('rightsideband','rightsideband')
C5.Divide(3,2)

C5.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw,RooFit.ProjectionRange('rightSideband'))
pdf.plotOn(tframe,RooFit.Components(bkgname),RooFit.ProjectionRange('rightSideband'),bkgcolor,dashed,lw)
pdf.plotOn(tframe,RooFit.Components(signame),RooFit.ProjectionRange('rightSideband'),sigcolor,dashed,lw)
tframe.Draw()
lhcbname.Draw('same')

C5.cd(2)
mframe = ws['m'].frame(RooFit.Range(msigmax,mmax),RooFit.Bins(30))
data.plotOn(mframe,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw,RooFit.ProjectionRange('rightSideband'))
pdf.plotOn(mframe,RooFit.Components(bkgname),RooFit.ProjectionRange('rightSideband'),bkgcolor,dashed,lw)
pdf.plotOn(mframe,RooFit.Components(signame),RooFit.ProjectionRange('rightSideband'),sigcolor,dashed,lw)
mframe.Draw()
lhcbname.Draw('same')

C5.cd(4)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(30))
data.plotOn(trcosthetaframe,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcosthetaframe,lw,RooFit.ProjectionRange('rightSideband'))
pdf.plotOn(trcosthetaframe,RooFit.Components(bkgname),RooFit.ProjectionRange('rightSideband'),bkgcolor,dashed,lw)
pdf.plotOn(trcosthetaframe,RooFit.Components(signame),RooFit.ProjectionRange('rightSideband'),sigcolor,dashed,lw)
trcosthetaframe.Draw()
lhcbname.Draw('same')

C5.cd(5)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(30))
data.plotOn(trcospsiframe,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcospsiframe,lw,RooFit.ProjectionRange('rightSideband'))
pdf.plotOn(trcospsiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('rightSideband'),bkgcolor,dashed,lw)
pdf.plotOn(trcospsiframe,RooFit.Components(signame),RooFit.ProjectionRange('rightSideband'),sigcolor,dashed,lw)
trcospsiframe.Draw()
lhcbname.Draw('same')

C5.cd(6)
trphiframe = ws['trphi'].frame(RooFit.Bins(30))
data.plotOn(trphiframe,RooFit.CutRange('rightSideband'),RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trphiframe,lw,RooFit.ProjectionRange('rightSideband'))
pdf.plotOn(trphiframe,RooFit.Components(bkgname),RooFit.ProjectionRange('rightSideband'),bkgcolor,dashed,lw)
pdf.plotOn(trphiframe,RooFit.Components(signame),RooFit.ProjectionRange('rightSideband'),sigcolor,dashed,lw)
trphiframe.Draw()
lhcbname.Draw('same')

C5.SaveAs(name+'_rightSideband.eps')
