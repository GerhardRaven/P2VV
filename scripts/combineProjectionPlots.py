plotVar = 'timeLog'
ROOTFilePaths = [ 'plots/Reco14/projPlot_%s_%s_moreBins.root' % ( plotVar, comp ) for comp in [ 'total' ] ]#, 'even', 'odd', 'S' ] ]
plotFilePath = 'plots/Reco14/projPlot_%s_moreBins.ps' % plotVar
minMax = dict(  timeLin = ( 0.,    6000. )
              , timeLog = ( 1.e-1, 2.e4  )
              , ctk     = ( 0.,    3500. )
              , ctl     = ( 0.,    3500. )
              , phi     = ( 0.,    3500. )
             )

from ROOT import gStyle
from P2VV.Load import LHCbStyle
gStyle.SetLineStyleString( 5, ' 40 20 10 20' )
gStyle.SetLineStyleString( 7, ' 40 20'       )
gStyle.SetLineStyleString( 9, ' 100 20'      )

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

from ROOT import TCanvas
canv = TCanvas()
canv.SetLeftMargin(0.18)
canv.SetRightMargin(0.05)
canv.SetBottomMargin(0.18)
canv.SetTopMargin(0.05)
if plotVar == 'timeLog' : canv.SetLogy()

projPlot.Draw()
LHCbLabel.Draw()

canv.Print(plotFilePath)
