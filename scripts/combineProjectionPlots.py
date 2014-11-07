plotVar = 'timeLin'
ROOTFilePaths = [ 'plots/Reco14/projPlot_%s_%s.root' % ( plotVar, comp ) for comp in [ 'total', 'even', 'odd', 'S' ] ]
plotFilePath = 'plots/Reco14/projPlot_%s_LHCb.pdf' % plotVar
minMax = dict(  timeLin = ( 0.,    6000. )
              , timeLog = ( 1.e-1, 2.e4  )
              , ctk     = ( 0.,    3500. )
              , ctl     = ( 0.,    3500. )
              , phi     = ( 0.,    3500. )
             )
drawLabel = True

from ROOT import gStyle
from P2VV.Load import LHCbStyle
gStyle.SetLineStyleString( 5, ' 40 20 10 20' )
gStyle.SetLineStyleString( 7, ' 40 20'       )
gStyle.SetLineStyleString( 9, ' 100 20'      )
gStyle.SetColorModelPS(0)

projPlot = None
LHCbLabel = None
from ROOT import TFile, RooPlot, TLatex
for it, path in enumerate(ROOTFilePaths) :
    plot = None
    ROOTFile = TFile.Open(path)
    print 'reading file %s' % ROOTFile.GetName()
    for key in ROOTFile.GetListOfKeys() :
        obj = key.ReadObj()
        if it == 0 and type(obj) == RooPlot :
            projPlot = obj
        elif type(obj) == RooPlot :
            obj.getCurve('pdf').SetName( 'pdf%d' % it )
            projPlot.addPlotable( obj.getCurve( 'pdf%d' % it ), obj.getDrawOptions( 'pdf%d' % it ).Data() )
            obj.remove( 'pdf%d' % it, False )
        elif it == 0 and type(obj) == TLatex :
            LHCbLabel = obj

projPlot.SetMinimum( minMax[plotVar][0] )
projPlot.SetMaximum( minMax[plotVar][1] )
if not drawLabel :
  projPlot.SetYTitle( projPlot.GetYaxis().GetTitle().replace( 'Candidates', 'Decays' ) )

# draw plot
from ROOT import TCanvas
canv0 = TCanvas('plotCanv')
canv0.SetLeftMargin(0.18)
canv0.SetRightMargin(0.05)
canv0.SetBottomMargin(0.18)
canv0.SetTopMargin(0.05)
if plotVar == 'timeLog' : canv0.SetLogy()

projPlot.Draw()
if drawLabel : LHCbLabel.Draw()
canv0.Print( plotFilePath + '(' )

# draw plot with residuals
residHist = projPlot.residHist( 'data', 'pdf', True )
residFrame = projPlot.emptyClone( projPlot.GetName() + '_resid' )
residFrame.SetXTitle( projPlot.GetXaxis().GetTitle() )
residFrame.SetBarWidth(1.02)
residFrame.addPlotable( residHist, 'BX' )
residFrame.SetMinimum(-5.)
residFrame.SetMaximum(+5.)

from ROOT import TLine, kRed
zeroLine = TLine( residFrame.GetXaxis().GetXmin(), 0., residFrame.GetXaxis().GetXmax(), 0. )
zeroLine.SetLineColor(kRed)
residFrame.addObject(zeroLine)

canv1 = TCanvas('plotCanvResid')
canv1.SetLeftMargin(0.18)
canv1.SetRightMargin(0.05)
canv1.SetBottomMargin(0.18)
canv1.SetTopMargin(0.05)

from ROOT import TPad
canv1.cd()
plotPad = TPad( 'plotPad', 'plotPad', 0., 0.32, 1., 1. )
if plotVar == 'timeLog' : plotPad.SetLogy()
plotPad.SetNumber(1)
plotPad.SetLeftMargin(0.18)
plotPad.SetRightMargin(0.12)
plotPad.SetBottomMargin(0.05)
plotPad.SetTopMargin(0.06)
plotPad.Draw()
canv1.cd(1)

projPlot.SetLabelSize( 0.92 / 0.68 * 0.060, 'y' )
projPlot.SetTitleSize( 0.92 / 0.68 * 0.072, 'y' )
projPlot.SetLabelOffset( 1.1,  'x' )
projPlot.SetLabelOffset( 0.01, 'y' )
projPlot.SetTitleOffset( 0.7 if plotVar == 'timeLog' else 0.8,  'y' )
projPlot.Draw()
if drawLabel :
    LHCbLabel.SetTextSize( 0.92 / 0.68 * 0.072 )
    LHCbLabel.SetX( 0.93 * LHCbLabel.GetX() )
    LHCbLabel.SetY( 0.97 * LHCbLabel.GetY() )
    LHCbLabel.Draw()

canv1.cd()
residPad = TPad( 'residPad', 'residPad', 0., 0., 1., 0.32 )
residPad.SetNumber(2)
residPad.SetLeftMargin(0.18)
residPad.SetRightMargin(0.12)
residPad.SetBottomMargin(0.54)
residPad.SetTopMargin(0.06)
residPad.Draw()
canv1.cd(2)

residFrame.GetYaxis().SetNdivisions( 4, 5, 0 )
residFrame.SetLabelSize( 0.92 / 0.32 * 0.060, 'x' )
residFrame.SetLabelSize( 0.92 / 0.32 * 0.060, 'y' )
residFrame.SetTitleSize( 0.92 / 0.32 * 0.072, 'x' )
residFrame.SetLabelOffset( 0.06, 'x' )
residFrame.SetLabelOffset( 0.01, 'y' )
residFrame.SetTitleOffset( 1.1, 'x' )
residFrame.Draw()

canv1.Print( plotFilePath + ')' )
