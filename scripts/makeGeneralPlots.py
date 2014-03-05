vars = [ 't', 'd', 'p', 'pT' ]
scaleFacs = dict( t = 1., d = 2.99792458e8 / 5366.77 * 1.e-9, p = 1.e-3, pT = 1.e-3 )
yScale = dict( t = 6000., d = 5200., p = 5200., pT = 4000. )
varUnits = dict( t = 'ps', d = 'mm', p = 'GeV/c', pT = 'GeV/c' )
binPrec = dict( t = 2, d = 2, p = 2, pT = 2 )
varTitles = dict(  t  = 'Decay time (%s)' % varUnits['t']
                 , d  = 'Flight distance (%s)' % varUnits['d']
                 , p  = 'Momentum (%s)' % varUnits['p']
                 , pT = 'Transverse momentum (%s)' % varUnits['pT']
                )
plotsFilePath = 'generalPlots.ps'
plotObjFilePath = 'generalPlots.root'
nTupleFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/fitNTuple_peakBkg_2011_2012_Reco14_TOS_HLT2B_20140215.root'

from ROOT import TFile
nTupleFile = TFile.Open(nTupleFilePath)
nTuple = nTupleFile.Get('DecayTree')
nTupleVars = dict( t = 'time', d = 'B_P*time', p = 'B_P', pT = 'B_Pt' )

from P2VV.Load import LHCbStyle
plotFile = TFile.Open( plotObjFilePath, 'RECREATE' )
from ROOT import TH1, TH1D
TH1.SetDefaultSumw2(True)
hists = dict(  t  = TH1D( 'tHist',  'tHist',  50, 0., 5.   )
             , d  = TH1D( 'dHist',  'dHist',  50, 0., 25.  )
             , p  = TH1D( 'pHist',  'pHist',  50, 0., 250. )
             , pT = TH1D( 'pTHist', 'pTHist', 53, 0., 15.9 )
            )

from ROOT import TLatex
label = TLatex()
label.SetTextAlign(32)
label.SetTextSize(0.072)

from ROOT import TCanvas, kFullDotLarge, kBlack as markCol
canvs = dict( t = TCanvas('tCanv'), d = TCanvas('dCanv'), p = TCanvas('pCanv'), pT = TCanvas('pTCanv') )
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
    hists[var].SetYTitle( ( 'Decays / (%%.%dg %s)' % ( binPrec[var], varUnits[var] ) ) % hists[var].GetBinWidth(1) )
    hists[var].SetTitleOffset( 1.10, 'x' )
    hists[var].SetTitleOffset( 1.15, 'y' )

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
    label.DrawLatex( xMin + 0.90 * ( xMax - xMin ), 0.85 * yMax, 'LHCb unofficial' )
    canv.Print(plotsFilePath)
canvs[ vars[0] ].Print( plotsFilePath + ']' )

assert False
from ROOT import TObject
plotFile.Write( plotObjFilePath, TObject.kOverwrite )
plotFile.Close()
nTupleFile.Close()
