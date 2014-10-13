mode = ''
if mode == 'polarizations' :
    plotFilePath = 'asymPlot_polarDep_realBins_points_polar.pdf'
    ROOTFilePaths = [ 'asymPlot_polarDep%s_realBins_points.root' % model for model in [ '', '_even', '_odd', '_S' ] ]
    minMax = ( -0.07, +0.07 )
else :
    plotFilePath = 'asymPlot_realBins_points.pdf'
    ROOTFilePaths = [ 'asymPlot_%s_realBins_points.root' % model for model in [ 'polarDep', 'lamb_phi', 'phi' ] ]
    minMax = ( -0.07, +0.07 )
drawLabel = True

from ROOT import gStyle
from P2VV.Load import LHCbStyle
gStyle.SetLineStyleString( 5, ' 40 20 10 20' )
gStyle.SetLineStyleString( 7, ' 40 20'       )
gStyle.SetLineStyleString( 9, ' 100 20'      )
gStyle.SetColorModelPS(1)

# function to set graph attributes
def setGraphAtts( graph, colour, markerSize, lineStyle = None ) :
    from ROOT import kFullDotLarge
    graph.SetMarkerStyle(kFullDotLarge)
    graph.SetLineColor(colour)
    graph.SetMarkerColor(colour)
    graph.SetLineWidth(3)
    graph.SetMarkerSize(markerSize)
    if lineStyle != None : graph.SetLineStyle(lineStyle)
    graph.GetXaxis().SetTitle('Decay time (modulo 2#pi/#Deltam_{s}) [ps]')
    graph.GetYaxis().SetTitle('B/#bar{B}-tag asymmetry')
    graph.GetXaxis().SetTitleOffset(1.2)
    graph.GetYaxis().SetTitleOffset(1.2)
    graph.SetMinimum( minMax[0] )
    graph.SetMaximum( minMax[1] )
    graph.GetXaxis().SetNdivisions( 10, 5, 0, True )
    graph.GetYaxis().SetNdivisions(  5, 5, 0, True )

from ROOT import TGraphErrors as Graph
Graph.setAtts = setGraphAtts

dummyGraph = None
dataGraph = None
pdfGraphs = [ ]
LHCbLabel = None
from ROOT import TFile, TCanvas, TGraphErrors, TLine, TLatex, kBlack, kBlue, kRed, kGreen, kMagenta, kSolid
lineColors = [ kBlue, kRed, kGreen + 3, kMagenta + 3 ] if mode == 'polarizations' else [ kBlue, kGreen + 3, kRed ]
lineStyles = [ kSolid, 9, 7, 5 ] if mode == 'polarizations' else [ kSolid, 5, 7 ]
for it, path in enumerate(ROOTFilePaths) :
    plot = None
    ROOTFile = TFile.Open(path)
    print 'reading file %s' % ROOTFile.GetName()
    for key in ROOTFile.GetListOfKeys() :
        obj = key.ReadObj()
        if type(obj) == TGraphErrors :
            if obj.GetName() == 'data' and it == 0 :
                dataGraph = obj
                dataGraph.setAtts( kBlack, 1. )
            elif obj.GetName() == 'pdf' :
                pdfGraphs.append(obj)
                pdfGraphs[-1].setAtts( lineColors[it], 0.1, lineStyles[it] )
            elif obj.GetName() == 'dummy' and it == 0 :
                dummyGraph = obj
                dummyGraph.setAtts( kBlack, 0. )
        elif type(obj) == TLine and it == 0 :
            zeroLine = obj
assert dummyGraph and dataGraph and pdfGraphs

canv = TCanvas('canv')
canv.SetLeftMargin(0.18)
canv.SetRightMargin(0.05)
canv.SetBottomMargin(0.18)
canv.SetTopMargin(0.05)

dummyGraph.Draw('AP')
zeroLine.Draw()
dataGraph.Draw('P SAMES')
for graph in reversed(pdfGraphs) : graph.Draw('L SAMES')
canv.Print(plotFilePath)
