samples = [ 'sigRange', 'signal', 'peakBkg' ]
vars = [ 't', 'ctk', 'ctl', 'phi', 'd', 'mB', 'mBErr', 'p', 'pT', 'st', 'wOS', 'wSS', 'qOS', 'qSS', 'cOS', 'cSS', 'KKCat' ]
selStr = ''
# 26468.096538697835, 44234.582149922244, 25140.64489205761
histType = dict( t = 'd', ctk = 'd', ctl = 'd', phi = 'd', mB = 'd', mBErr = 'd', d = 'd', p = 'd', pT = 'd', st = 'd', wOS = 'd', wSS = 'd'
                , qOS = 'i', qSS = 'i', cOS = 'i', cSS = 'i', KKCat = 'i' )
scaleFacs = dict( t = 1., ctk = 1., ctl = 1., phi = 1., d = 2.99792458e8 / 5366.77 * 1.e-9, mB = 1., mBErr = 1., p = 1.e-3, pT = 1.e-3, st = 1.
                 , wOS = 1., wSS = 1., qOS = 1., qSS = 1., cOS = 1., cSS = 1., KKCat = 1. )
yScaleMin = dict( sigRange = { }, signal = dict( mB = -500. ), peakBkg = { } )
yScaleMax = dict(  sigRange = dict( t = 16000., ctk = 4500., ctl = 4500., phi = 4500., d = 9000., mB = 16000., mBErr = 7000., p = 7000., pT = 5000.
                                   , st = 14000., wOS = 2500., wSS = 12000., qOS = 65000., qSS = 65000., cOS = 85000., cSS = 65000.
                                   , KKCat = 40000. )
                 , signal   = dict( t =  6000., ctk = 3500., ctl = 3500., phi = 3500., d = 5200., mB = 18000., mBErr = 2000. *
                                    25140.6/26468.096538697835, p = 5200., pT = 4000.
                                   , st = 11000., wOS = 1800., wSS = 10000., qOS = 60000., qSS = 60000., cOS = 80000., cSS = 60000.
                                   , KKCat = 40000. )
                 , peakBkg  = dict( t = 450., ctk = 350., ctl = 350., phi = 350., d = 450., mB = 500., mBErr = 500., p = 450., pT = 300.
                                   , st = 900., wOS = 120., wSS = 1000., qOS = 4000., qSS = 4000., cOS = 5000., cSS = 4000.
                                   , KKCat = 3500. )
                )
varUnits = dict( t = 'ps', ctk = '', ctl = '', phi = 'rad', d = 'mm', mB = 'MeV/c^{2}', mBErr = 'MeV/c^{2}', p = 'GeV/c', pT = 'GeV/c', st = 'ps'
                , wOS = '', wSS = '', KKCat = '' )
binPrec = dict( t = 2, ctk = 2, ctl = 2, phi = 2, d = 2, mB = 2, mBErr = 2, p = 2, pT = 2, st = 2, wOS = 2, wSS = 2, qOS = 2, qSS = 2
               , cOS = 2, cSS = 2, KKCat = 2 )
varTitles = dict(  t     = 'Decay time [%s]' % varUnits['t']
                 , ctk   = 'cos#kern[0.1]{#theta_{K}}'
                 , ctl   = 'cos#kern[0.1]{#theta_{#mu}}'
                 , phi   = '#varphi_{h} [rad]'
                 , d     = 'Flight distance [%s]' % varUnits['d']
                 , mB    = 'm(J/#psi K^{+}K^{-}) [%s]' % varUnits['mB']
                 , mBErr = '#sigma(J/#psi K^{+}K^{-}) [%s]' % varUnits['mBErr']
                 , p     = 'Momentum [%s]' % varUnits['p']
                 , pT    = 'Transverse momentum [%s]' % varUnits['pT']
                 , st    = 'Estimated decay-time resolution [%s]' % varUnits['st']
                 , wOS   = '#eta_{OS}'
                 , wSS   = '#eta_{SS}'
                 , qOS   = 'OS flavour tag'
                 , qSS   = 'SS flavour tag'
                 , cOS   = 'OS flavour tagging category'
                 , cSS   = 'SS flavour tagging category'
                 , KKCat = 'KK-mass category'
                )
plotsFilePath = 'generalPlots.pdf'
plotObjFilePath = 'generalPlots.root'
nTupleFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/fitNTuple_peakBkg_2011_2012_Reco14_TOS_HLT2B_20140309.root'
#nTupleFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/nTupleC_merged.root'

from ROOT import TFile
nTupleFile = TFile.Open(nTupleFilePath)
nTuple = nTupleFile.Get('DecayTree')
nTupleVars = dict( t = 'time', ctk = 'helcosthetaK', ctl = 'helcosthetaL', phi = 'helphi', d = 'B_P*time', mB = 'mass', mBErr = 'sigmam', p = 'B_P'
                  , pT = 'B_Pt', st = 'sigmat', wOS = 'tagomega_os_cb', wSS = 'tagomega_ss_nn', qOS = 'iTagOS', qSS = 'iTagSS'
                  , cOS = 'tagCatP2VVOS', cSS = 'tagCatP2VVSS', KKCat = 'KKMassCat' )

from math import pi
from P2VV.Load import LHCbStyle
plotFile = TFile.Open( plotObjFilePath, 'RECREATE' )
from ROOT import TH1, TH1D
TH1.SetDefaultSumw2(True)
hists = dict( [ ( sample, dict(  t     = TH1D( 'tHist_%s'     % sample, 'tHist',     47,  0.3,   5.   )
                               , ctk   = TH1D( 'ctkHist_%s'   % sample, 'ctkHist',   40, -1,     1.   )
                               , ctl   = TH1D( 'ctlHist_%s'   % sample, 'ctlHist',   40, -1,     1.   )
                               , phi   = TH1D( 'phiHist_%s'   % sample, 'phiHist',   40, -pi,    pi   )
                               , d     = TH1D( 'dHist_%s'     % sample, 'dHist',     50,  0.,    25.  )
                               , mBErr = TH1D( 'mBErrHist_%s' % sample, 'mBErrHist', 50,  2.,     8.  )
                               , p     = TH1D( 'pHist_%s'     % sample, 'pHist',     50,  0.,    250. )
                               , pT    = TH1D( 'pTHist_%s'    % sample, 'pTHist',    53,  0.,    15.9 )
                               , st    = TH1D( 'stHist_%s'    % sample, 'stHist',    50,  0.,    0.12 )
                               , wOS   = TH1D( 'wOSHist_%s'   % sample, 'wOSHist',   50,  0.,    0.5  )
                               , wSS   = TH1D( 'wSSHist_%s'   % sample, 'wSSHist',   50,  0.,    0.5  )
                               , qOS   = TH1D( 'qOSHist_%s'   % sample, 'qOSHist',    3, -1.5,   1.5  )
                               , qSS   = TH1D( 'qSSHist_%s'   % sample, 'qSSHist',    3, -1.5,   1.5  )
                               , cOS   = TH1D( 'cOSHist_%s'   % sample, 'cOSHist',    2, -0.5,   1.5  )
                               , cSS   = TH1D( 'cSSHist_%s'   % sample, 'cSSHist',    2, -0.5,   1.5  )
                               , KKCat = TH1D( 'KKCatHist_%s' % sample, 'KKCat',      6, -0.5,   5.5  )
                              )
              ) for sample in samples ]
            )
hists['sigRange']['mB'] = TH1D( 'mBHist_sigRange', 'mBHist', 24, 5338., 5398. )
hists['signal']['mB']   = TH1D( 'mBHist_signal',   'mBHist', 56, 5300., 5440. )
hists['peakBkg']['mB']  = TH1D( 'mBHist_peakBkg',  'mBHist', 35, 5200., 5550. )

labelText = ''
if labelText :
    from ROOT import TLatex
    label = TLatex()
    label.SetTextAlign(32)
    label.SetTextSize(0.072)

import P2VV.RooFitWrappers
from ROOT import gStyle
gStyle.SetColorModelPS(1)

from ROOT import TCanvas, kFullDotLarge, kBlack as markCol
for sample in samples :
    plotsFile = '%s_%s.%s' % ( plotsFilePath.split('.')[0], sample, plotsFilePath.split('.')[1] )
    canvs = dict( [ ( var, TCanvas( var + 'Canv_%s' % sample ) ) for var in vars ] )
    canvs[ vars[0] ].Print( plotsFile + '[' )
    for var in vars :
        hist = hists[sample][var]
        hist.SetMinimum( yScaleMin[sample][var] if var in yScaleMin[sample] else 0. )
        if var in yScaleMax[sample] : hist.SetMaximum( yScaleMax[sample][var] )
        hist.SetMarkerStyle(kFullDotLarge)
        hist.SetMarkerSize(0.6)
        hist.SetMarkerColor(markCol)
        hist.SetLineColor(markCol)
        hist.SetLineWidth(2)
        hist.SetXTitle( varTitles[var] )
        yTitle = 'Decays' if sample == 'signal' else 'Candidates'
        if var == 'phi' :
          yTitle += ( ' / (%%.%dg#kern[0.3]{#pi} %s)' % ( binPrec[var], varUnits[var] ) ) % ( hist.GetBinWidth(1) / pi )
        elif histType[var] == 'd' and varUnits[var] :
          yTitle += ( ' / (%%.%dg %s)' % ( binPrec[var], varUnits[var] ) ) % hist.GetBinWidth(1)
        elif histType[var] == 'd' :
          yTitle += ( ' / %%.%dg' % binPrec[var] ) % hist.GetBinWidth(1)
        hist.SetYTitle( yTitle )
        hist.SetTitleOffset( 1.15 if var == 'mB' else 1.10, 'x' )
        hist.SetTitleOffset( 1.00 if sample == 'peakBkg' else 1.15, 'y' )
        if var == 'mB' : hist.SetLabelOffset( 0.03, 'x' )
        if histType[var] == 'i' : hist.GetXaxis().SetNdivisions( hist.GetNbinsX(), 0, 0 )
        hist.GetYaxis().SetNdivisions( 5, 5, 0 )

        canv = canvs[var].cd()
        canv.SetLeftMargin(0.18)
        canv.SetRightMargin(0.05)
        canv.SetBottomMargin(0.18)
        canv.SetTopMargin(0.05)

        if sample == 'peakBkg' :
            sel = ( '-wMC*(wMC < 0. && %s)' % selStr ) if selStr else '-wMC*(wMC < 0.)'
        elif sample == 'signal' :
            sel = ( 'sWeights_ipatia*(%s)' % selStr ) if selStr else 'sWeights_ipatia'
        else :
            sel = ( 'fabs(mass - 5368.) < 30. && %s' % selStr ) if selStr else 'fabs(mass - 5368.) < 30.'
        nTuple.Draw( '%s*%f >> %s' % ( nTupleVars[var], scaleFacs[var], hist.GetName() ), sel, 'E1' )
        if var not in yScaleMin[sample] or yScaleMin[sample][var] >= 0. :
            for bin in range( hist.GetNbinsX() ) :
                if hist.GetBinContent( bin + 1 ) < 0. : hist.SetBinContent( bin + 1, 0. )
        canv.Update()

        xMin = hist.GetXaxis().GetXmin()
        xMax = hist.GetXaxis().GetXmax()
        yMax = hist.GetMaximum()
        if labelText : label.DrawLatex( xMin + 0.90 * ( xMax - xMin ), 0.85 * yMax, labelText )
        canv.Print(plotsFile)
    canvs[ vars[0] ].Print( plotsFile + ']' )

plotsFile = '%s_bkgSub.%s' % ( plotsFilePath.split('.')[0], plotsFilePath.split('.')[1] )
canvs = dict( [ ( var, TCanvas( var + 'Canv_bkgSub' ) ) for var in vars ] )
canvs[ vars[0] ].Print( plotsFile + '[' )
from ROOT import kRed, kBlue
for var in vars :
    if var == 'mB' : continue
    hist1 = hists['sigRange'][var]
    hist2 = hists['signal'][var]
    hist1.SetMarkerColor(kRed)
    hist1.SetLineColor(kRed)
    hist2.SetMarkerColor(kBlue)
    hist2.SetLineColor(kBlue)

    canv = canvs[var].cd()
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)

    hist1.Draw('E1')
    hist2.Draw('E1 SAMES')
    canv.Print(plotsFile)

canvs[ vars[0] ].Print( plotsFile + ']' )

from ROOT import TObject
plotFile.Write( plotObjFilePath, TObject.kOverwrite )
plotFile.Close()
nTupleFile.Close()
