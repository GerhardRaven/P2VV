#Do 2011 Unbiased Fit with DSC parameterization

from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi
from array import array

from RooFitDecorators import *

from lhcbStyle import *
gStyle.SetOptTitle(0)
lhcbname = printLHCb(gStyle,'L','Prelim','#sqrt{s}=7 TeV Data, L = 370 pb^{-1}')
#style = lhcbStyle()

angcorr = True
taggingsyst = True
blinded = True
phisparam = False
splitmKK = False

#name = 'FitXC_angcorr_flatbkg'
#wsfile = TFile('UnbiasedWS_flatbkg.root')

name = 'angcorr_DSC'
wsfile = TFile('UnbiasedWS_DSC.root')

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
        ws.factory("RooUnblindUniform::D_blind('BsMiniMouseLP2011',2.,D)")
        #ws.factory("RooUnblindUniform::C_blind('BsDonaldDuckLP2011',2.,C)")

        #calling RooCustomizer
        customizer = RooCustomizer(pdfbeforeblinding,'blinded')
        customizer.replaceArg( ws['t_sig_dG'], ws['t_sig_dG_blind'] )
        customizer.replaceArg( ws['S'], ws['S_blind'] )
        customizer.replaceArg( ws['D'], ws['D_blind'] )
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

ws['rs2'].setVal(0.03)
ws['rs2'].setMin(0.)
ws['rs2'].setMax(1.)
#ws['rs2'].setConstant()
#ws['deltas'].setConstant()

ws['rz2'].setVal(0.485)
ws['rz2'].setMin(0.)
ws['rz2'].setMax(1.)

ws['rs2'].setVal(0.05)
ws['rs2'].setMin(0.)
ws['rs2'].setMax(1.)

ws['rperp2'].setVal(0.259)
ws['rperp2'].setMin(0.)
ws['rperp2'].setMax(1.)

ws['deltas'].setVal(3.0)
ws['deltas'].setMin(-2.*pi)
ws['deltas'].setMax(2.*pi)

ws['deltaperp'].setVal(2.79)
ws['deltaperp'].setMin(-2.*pi)
ws['deltaperp'].setMax(2.*pi)

ws['deltapar'].setVal(3.24)
ws['deltapar'].setMin(-2.*pi)
ws['deltapar'].setMax(2.*pi)

ws['S'].setMin(-0.5)
ws['S'].setMax(2.5)

ws['D'].setMin(-0.5)
ws['D'].setMax(2.5)

if splitmKK:
    #################################
    ### TEST OF SPLIT IN MKK BINS ###
    #################################
    mdau2canvas = TCanvas('mdau2canvas','mdau2canvas')
    mdau2frame = ws['mdau2'].frame()
    ws['unbiaseddata'].plotOn(mdau2frame)
    mdau2frame.Draw()

    mdau2min = ws['mdau2'].getMin()
    mdau2sigmin = 1014
    mdau2sigmax = 1026
    mdau2max = ws['mdau2'].getMax()
    mKKcat = RooThresholdCategory('mKKcat','mKKcat',ws['mdau2'],'LeftSideBand') #Can't make it using factory: something goes wrong when trying to make the default value (LeftSideBand in this case)
    ws.put(mKKcat)
    ws['mKKcat'].addThreshold(mdau2sigmin,"LeftSideBand")
    ws['mKKcat'].addThreshold(mdau2sigmax,"Signal")
    ws['mKKcat'].addThreshold(mdau2max,"RightSideBand")
    ws['unbiaseddata'].table(ws['mKKcat']).Print('v')

    addcolumn = ws['unbiaseddata'].addColumn(ws['mKKcat'])
    addcolumn.SetName('thismKKcat')
    ws.put(addcolumn)

    #ws.factory("SIMCLONE::simpdf(pdf_ext_inc_tag_syst_blinded,$SplitParam({rs2,deltas,Nbkg,Nsig},thismKKcat))")
    ws.factory("SIMCLONE::simpdf(pdf_ext_inc_tag_syst_blinded,$SplitParam({rs2,deltas,Nbkg,Nsig},thismKKcat))")

    print 'Going to fit the Simultaneous PDF....'
    result = ws['simpdf'].fitTo(ws['unbiaseddata'],RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))

    ### Plot roocategory against parameters ###
    NmKKBins = 3
    xlist = [0,1,2]
    xerrlist = [0.,0.,0.]
    x = array('f', xlist)
    xerr = array('f', xerrlist)

    rs2namelist = ['rs2_LeftSideBand','rs2_Signal','rs2_RightSideBand']
    deltasnamelist = ['deltas_LeftSideBand','deltas_Signal','deltas_RightSideBand']
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
    rs2_v_cat.SetTitle("|A_{s}|^{2} vs. category")
    rs2_v_cat.SetMarkerStyle(20)
    rs2_v_cat.SetMarkerColor(4)
    rs2_v_cat.GetXaxis().SetTitle("Category")
    rs2_v_cat.GetYaxis().SetTitle("|A_{s}|^{2}")
    #rs2_v_cat.SetMaximum(1.3)

    deltasy = array('f',deltaslist)
    deltasyerr = array('f',deltaserrlist)
    deltas_v_cat = TGraphErrors(NmKKBins, x, deltasy, xerr, deltasyerr);
    deltas_v_cat.SetTitle("\delta_{s} vs. category")
    deltas_v_cat.SetMarkerStyle(20)
    deltas_v_cat.SetMarkerColor(4)
    deltas_v_cat.GetXaxis().SetTitle("Category")
    deltas_v_cat.GetYaxis().SetTitle("\delta_{s}")
    #deltas_v_cat.SetMaximum(1.3)

    ParamCanvas = TCanvas('ParamCanvas','ParamCanvas')
    ParamCanvas.Divide(2)

    ParamCanvas.cd(1)
    rs2_v_cat.GetXaxis().SetLimits(-1,3)
    rs2_v_cat.GetXaxis().SetBinLabel(25,'LeftSideBand')
    rs2_v_cat.GetXaxis().SetBinLabel(50,'Signal')
    rs2_v_cat.GetXaxis().SetBinLabel(75,'RightSideBand')
    rs2_v_cat.GetXaxis().LabelsOption('h')
    rs2_v_cat.Draw('AP')

    ParamCanvas.cd(2)
    deltas_v_cat.GetXaxis().SetLimits(-1,3)
    deltas_v_cat.GetXaxis().SetBinLabel(25,'LeftSideBand')
    deltas_v_cat.GetXaxis().SetBinLabel(50,'Signal')
    deltas_v_cat.GetXaxis().SetBinLabel(75,'RightSideBand')
    deltas_v_cat.GetXaxis().LabelsOption('h')
    deltas_v_cat.Draw('AP')

    #Show sig/bkg per category:
    masscanvas = TCanvas('masscanvas','masscanvas')
    masscanvas.Divide(3,1)

    ws['mdau2'].setRange('LeftSideBand',mdau2min,mdau2sigmin)
    ws['mdau2'].setRange('Signal',mdau2sigmin,mdau2sigmax)
    ws['mdau2'].setRange('RightSideBand',mdau2sigmax,mdau2max)

    argset = RooArgSet(ws['thismKKcat'])
    leftsideband = 'pdf_ext_inc_tag_syst_blinded_LeftSideBand'
    signal = 'pdf_ext_inc_tag_syst_blinded_Signal'
    rightsideband = 'pdf_ext_inc_tag_syst_blinded_RightSideBand'
    #catdataset = RooDataSet('catdataset','catdataset',ws['unbiaseddata'],RooArgSet(ws['thismKKcat']))
    #ws.put(catdataset)

    masscanvas.cd(1)
    mframe1 = ws['m'].frame()
    ws['unbiaseddata'].plotOn(mframe1,RooFit.CutRange('LeftSideBand'))
    ws['simpdf'].plotOn(mframe1,RooFit.Slice(ws['thismKKcat'],'LeftSideBand'),RooFit.ProjWData(argset,ws['unbiaseddata']))
    mframe1.Draw()

    masscanvas.cd(2)
    mframe2 = ws['m'].frame()
    ws['unbiaseddata'].plotOn(mframe2,RooFit.CutRange('Signal'))
    ws['simpdf'].plotOn(mframe2,RooFit.Slice(ws['thismKKcat'],'Signal'),RooFit.ProjWData(argset,ws['unbiaseddata']))
    mframe2.Draw()

    masscanvas.cd(3)
    mframe3 = ws['m'].frame()
    ws['unbiaseddata'].plotOn(mframe3,RooFit.CutRange('RightSideBand'))
    ws['simpdf'].plotOn(mframe3,RooFit.Slice(ws['thismKKcat'],'RightSideBand'),RooFit.ProjWData(argset,ws['unbiaseddata']))
    mframe3.Draw()
    
    ##################################################################################
pdf = pdf
data = ws['unbiaseddata']

#Construct unbinned likelihood
nll = pdf.createNLL(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.ExternalConstraints(RooArgSet(ws['p0'],ws['p1'],ws['dmsconstraint'],ws['tres_SFconstraint'])))

print 'NLL value =',nll.getVal()

if taggingsyst:
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws['p0'],ws['p1'],ws['dmsconstraint'],ws['tres_SFconstraint'])))
else:
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws['dmsconstraint'],ws['tres_SFconstraint'])))

print 'NLL value =',nll.getVal()

writeFitParamsLatex(result,name,False)
dict = writeCorrMatrixLatex(result,name)

#################
### Single LL ###
#################
CanL = TCanvas('CanL','CanL')
CanL.Divide(3,1)

#Minimize likelihood w.r.t all parameters before making plots
#minuit = RooMinuit(nll)
#minuit.migrad()

CanL.cd(1)
gammaframe = ws['t_sig_dG'].frame(RooFit.Bins(10),RooFit.Title("LL in #Delta#Gamma"))#Range(0.01,0.95)
gammaframe.SetXTitle('#Delta#Gamma')
nll.plotOn(gammaframe,RooFit.ShiftToZero())
gammaframe.Draw()

CanL.cd(2)
Sframe = ws['S'].frame(RooFit.Bins(10),RooFit.Title("LL in S"))#Range(3.3,5.0),
Sframe.SetXTitle('S')
nll.plotOn(Sframe,RooFit.ShiftToZero())
Sframe.Draw()

CanL.cd(3)
Dframe = ws['D'].frame(RooFit.Bins(10),RooFit.Title("LL in D"))#Range(3.3,5.0),
Dframe.SetXTitle('D')
nll.plotOn(Dframe,RooFit.ShiftToZero())
Dframe.Draw()

######PLOT######
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

signame = 'sig_pdf_angcorr_inc_tag_syst_blinded'
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
