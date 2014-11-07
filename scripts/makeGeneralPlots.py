vars = [ 't', 'd', 'p', 'pT', 'st', 'wOS', 'wSS', 'qOS', 'qSS', 'cOS', 'cSS', 'KKCat' ]
histType = dict( t = 'd', d = 'd', p = 'd', pT = 'd', st = 'd', wOS = 'd', wSS = 'd', qOS = 'i', qSS = 'i', cOS = 'i', cSS = 'i'
                , KKCat = 'i' )
scaleFacs = dict( t = 1., d = 2.99792458e8 / 5366.77 * 1.e-9, p = 1.e-3, pT = 1.e-3, st = 1., wOS = 1., wSS = 1., qOS = 1., qSS = 1.
                 , cOS = 1., cSS = 1., KKCat = 1. )
yScale = dict( t = 6000., d = 5200., p = 5200., pT = 4000., st = 11000., wOS = 1800., wSS = 7000., qOS = 60000., qSS = 60000.
              , cOS = 80000., cSS = 60000., KKCat = 40000. )
varUnits = dict( t = 'ps', d = 'mm', p = 'GeV/c', pT = 'GeV/c', st = 'ps', wOS = '', wSS = '', KKCat = '' )
binPrec = dict( t = 2, d = 2, p = 2, pT = 2, st = 2, wOS = 2, wSS = 2, qOS = 2, qSS = 2, cOS = 2, cSS = 2, KKCat = 2 )
varTitles = dict(  t     = 'Decay time [%s]' % varUnits['t']
                 , d     = 'Flight distance [%s]' % varUnits['d']
                 , p     = 'Momentum [%s]' % varUnits['p']
                 , pT    = 'Transverse momentum [%s]' % varUnits['pT']
                 , st    = 'Estimated decay-time resolution [%s]' % varUnits['st']
                 , wOS   = 'Estimated OS wrong-tag probability'
                 , wSS   = 'Estimated SS wrong-tag probability'
                 , qOS   = 'OS flavour tag'
                 , qSS   = 'SS flavour tag'
                 , cOS   = 'OS flavour tagging category'
                 , cSS   = 'SS flavour tagging category'
                 , KKCat = 'KK-mass category'
                )
plotsFilePath = 'generalPlots.pdf'
plotObjFilePath = 'generalPlots.root'
nTupleFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/fitNTuple_peakBkg_2011_2012_Reco14_TOS_HLT2B_20140309.root'

from ROOT import TFile
nTupleFile = TFile.Open(nTupleFilePath)
nTuple = nTupleFile.Get('DecayTree')
nTupleVars = dict( t = 'time', d = 'B_P*time', p = 'B_P', pT = 'B_Pt', st = 'sigmat', wOS = 'tagomega_os_cb', wSS = 'tagomega_ss_nn'
                  , qOS = 'iTagOS', qSS = 'iTagSS', cOS = 'tagCatP2VVOS', cSS = 'tagCatP2VVSS', KKCat = 'KKMassCat' )

from P2VV.Load import LHCbStyle
plotFile = TFile.Open( plotObjFilePath, 'RECREATE' )
from ROOT import TH1, TH1D
TH1.SetDefaultSumw2(True)
hists = dict(  t     = TH1D( 'tHist',   'tHist',   50,  0.,  5.   )
             , d     = TH1D( 'dHist',   'dHist',   50,  0.,  25.  )
             , p     = TH1D( 'pHist',   'pHist',   50,  0.,  250. )
             , pT    = TH1D( 'pTHist',  'pTHist',  53,  0.,  15.9 )
             , st    = TH1D( 'stHist',  'stHist',  50,  0.,  0.12 )
             , wOS   = TH1D( 'wOSHist', 'wOSHist', 50,  0.,  0.5  )
             , wSS   = TH1D( 'wSSHist', 'wSSHist', 50,  0.,  0.5  )
             , qOS   = TH1D( 'qOSHist', 'qOSHist',  3, -1.5, 1.5  )
             , qSS   = TH1D( 'qSSHist', 'qSSHist',  3, -1.5, 1.5  )
             , cOS   = TH1D( 'cOSHist', 'cOSHist',  2, -0.5, 1.5  )
             , cSS   = TH1D( 'cSSHist', 'cSSHist',  2, -0.5, 1.5  )
             , KKCat = TH1D( 'KKCat',   'KKCat',    6, -0.5, 5.5  )
            )
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
canvs = dict( [ ( var, TCanvas( var + 'Canv' ) ) for var in vars ] )
canvs[ vars[0] ].Print( plotsFilePath + '[' )
for var in vars :
    hists[var].SetMinimum(0.)
    hists[var].SetMaximum( yScale[var] )
    hists[var].SetMarkerStyle(kFullDotLarge)
    hists[var].SetMarkerSize(0.6)
    hists[var].SetMarkerColor(markCol)
    hists[var].SetLineColor(markCol)
    hists[var].SetLineWidth(2)
    hists[var].SetXTitle( varTitles[var] )
    if histType[var] == 'd' :
      yTitle = ( 'Decays / (%%.%dg%s)' % ( binPrec[var], ( ' ' + varUnits[var] ) if varUnits[var] else '' ) ) % hists[var].GetBinWidth(1)
    else :
      yTitle = 'Decays'
    hists[var].SetYTitle( yTitle )
    hists[var].SetTitleOffset( 1.10, 'x' )
    hists[var].SetTitleOffset( 1.15, 'y' )
    if histType[var] == 'i' : hists[var].GetXaxis().SetNdivisions( hists[var].GetNbinsX(), 0, 0 )

    canv = canvs[var].cd()
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)

    nTuple.Draw( '%s*%f >> %s' % ( nTupleVars[var], scaleFacs[var], hists[var].GetName() ), 'sWeights_ipatia', 'E1' )
    for bin in range( hists[var].GetNbinsX() ) :
        if hists[var].GetBinContent( bin + 1 ) < 0. : hists[var].SetBinContent( bin + 1, 0. )
    canv.Update()

    xMin = hists[var].GetXaxis().GetXmin()
    xMax = hists[var].GetXaxis().GetXmax()
    yMax = hists[var].GetMaximum()
    if labelText : label.DrawLatex( xMin + 0.90 * ( xMax - xMin ), 0.85 * yMax, labelText )
    canv.Print(plotsFilePath)
canvs[ vars[0] ].Print( plotsFilePath + ']' )

from ROOT import TObject
plotFile.Write( plotObjFilePath, TObject.kOverwrite )
plotFile.Close()
nTupleFile.Close()
