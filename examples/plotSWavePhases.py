from P2VVLoad import ROOTStyle

plotsFilePath   = 'plotSWavePhases.ps'
deltaSAxisRange = ( -1.2, 4.2 )
KKMassTitle     = 'm_{KK} (MeV)'
deltaSTitle     = '#delta_{S} - #delta_{#perp}    (rad)'
CFitTitle       = 'S-Wave Phases Classic Fit'
SFitTitle       = 'S-Wave Phases S-Fit'

from array import array
KKMass        = array( 'd', [ 1020. - 21., 1020. - 6., 1020. + 6., 1020. + 21. ] )
KKMassLowErr  = array( 'd', [ 9.,          6.,         6.,         9.          ] )
KKMassHighErr = KKMassLowErr

from ROOT import TGraphAsymmErrors
deltaSCFit1        = array( 'd', [ 9.5903e-01, 3.7027e-01, -5.7206e-01, -5.5373e-01 ] )
deltaSCFit1LowErr  = array( 'd', [ 4.03e-01,   1.39e-01,    2.33e-01,    1.71e-01   ] )
deltaSCFit1HighErr = deltaSCFit1LowErr
deltaSCFit1Graph = TGraphAsymmErrors(len(KKMass), KKMass, deltaSCFit1, KKMassLowErr, KKMassHighErr, deltaSCFit1LowErr, deltaSCFit1HighErr)

deltaSCFit2        = array( 'd', [ 2.1797e+00, 2.7708e+00,  3.7159e+00,  3.6958e+00 ] )
deltaSCFit2LowErr  = array( 'd', [ 4.02e-01,   1.34e-01,    2.70e-01,    1.72e-01   ] )
deltaSCFit2HighErr = deltaSCFit2LowErr
deltaSCFit2Graph = TGraphAsymmErrors(len(KKMass), KKMass, deltaSCFit2, KKMassLowErr, KKMassHighErr, deltaSCFit2LowErr, deltaSCFit2HighErr)

deltaSSFit1        = array( 'd', [ 1.3853e+00, 3.2835e-01, -4.9483e-01, -6.0078e-01 ] )
deltaSSFit1LowErr  = array( 'd', [ 3.94e-01,   1.15e-01,    1.82e-01,    1.69e-01   ] )
deltaSSFit1HighErr = deltaSSFit1LowErr
deltaSSFit1Graph = TGraphAsymmErrors(len(KKMass), KKMass, deltaSSFit1, KKMassLowErr, KKMassHighErr, deltaSSFit1LowErr, deltaSSFit1HighErr)

deltaSSFit2        = array( 'd', [ 1.7570e+00, 2.8132e+00,  3.6369e+00,  3.7419e+00 ] )
deltaSSFit2LowErr  = array( 'd', [ 3.94e-01,   1.15e-01,    1.84e-01,    1.69e-01   ] )
deltaSSFit2HighErr = deltaSSFit2LowErr
deltaSSFit2Graph = TGraphAsymmErrors(len(KKMass), KKMass, deltaSSFit2, KKMassLowErr, KKMassHighErr, deltaSSFit2LowErr, deltaSSFit2HighErr)

from ROOT import kBlack, kBlue, kFullCircle, kFullSquare
from operator import methodcaller
from itertools import product, chain, imap
for (obj,opt) in chain( product( [ deltaSCFit1Graph, deltaSSFit1Graph ] 
                               , [ ('SetLineColor', kBlue) ,               ('SetMarkerColor', kBlue)
                                 , ('SetLineWidth', 2) ,                   ('SetMarkerStyle', kFullCircle)
                                 , ('SetMinimum',   deltaSAxisRange[0] ) , ('SetMaximum',     deltaSAxisRange[1] ) ] 
                               )
                      , product( [ deltaSCFit2Graph, deltaSSFit2Graph ]  
                               , [ ('SetLineColor',kBlack) ,               ('SetMarkerColor',kBlack)
                                 , ('SetLineWidth',2) ,                    ('SetMarkerStyle',kFullSquare) ]
                               ) 
                      , product( imap( methodcaller('GetYaxis'), [ deltaSCFit1Graph, deltaSSFit1Graph ] )
                               , [ ('SetTitleOffset', 1 ), ( 'SetTitle', deltaSTitle ) ] 
                               )
                      , product( imap( methodcaller('GetXaxis'), [ deltaSCFit1Graph, deltaSSFit1Graph ] ) 
                               , [ ('SetTitle',KKMassTitle ) ]
                               )
                      ) :
    methodcaller( *opt )( obj )


deltaSCFit1Graph.SetTitle(CFitTitle)
deltaSSFit1Graph.SetTitle(SFitTitle)

from ROOT import TLegend
leg = TLegend( 0.52, 0.43, 0.90, 0.63 )
leg.AddEntry( deltaSCFit1Graph, 'solution I  (#Delta#Gamma_{s} > 0)', 'LPE' )
leg.AddEntry( deltaSCFit2Graph, 'solution II (#Delta#Gamma_{s} < 0)', 'LPE' )
leg.SetBorderSize(0)
leg.SetFillStyle(0)

from ROOT import TCanvas
CFitCanv = TCanvas( 'CFitCanv', 'S-Wave Phases C-Fit' )
CFitCanv.SetLeftMargin(0.1)
deltaSCFit1Graph.Draw('AP')
deltaSCFit2Graph.Draw('P SAMES')
leg.Draw()
CFitCanv.Print( plotsFilePath + '(' )

SFitCanv = TCanvas( 'SFitCanv', 'S-Wave Phases S-Fit' )
SFitCanv.SetLeftMargin(0.1)
deltaSSFit1Graph.Draw('AP')
deltaSSFit2Graph.Draw('P sames')
leg.Draw()
SFitCanv.Print( plotsFilePath + ')' )