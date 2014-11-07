dataSetPath   = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140309.root'
dataSetName   = 'JpsiKK_sigSWeight'
plotsFilePath = 'sigWeights.pdf'
massName      = 'mass'
hlt2BName     = 'hlt2_biased'
wMCName       = 'wMC'
plotSamples   = [ 'all', 'realB', 'realNB', 'MCB', 'MCNB' ]
wMinMax       = [ ( -5.4, +2.0 ), ( -1.6, +2.0 ) ]
#wMinMax       = [ ( None, None ), ( None, None ) ]

# read data
from ROOT import TFile
dataSetFile = TFile.Open(dataSetPath)
dataSet = dataSetFile.Get(dataSetName)
dataSetFile.Close()

# get signal weights as a function of mass
import P2VV.RooFitWrappers
massVals = dict( [ ( sample, [ ] ) for sample in plotSamples ] )
weightVals  = dict( [ ( sample, [ ] ) for sample in plotSamples ] )
mass = dataSet.get().find(massName)
for obsSet in dataSet :
    if 'all' in plotSamples :
        massVals['all'].append( mass.getVal() )
        weightVals['all'].append( dataSet.weight() )
    if 'realB' in plotSamples and obsSet.getRealValue(wMCName) > 0. and obsSet.getCatIndex(hlt2BName) == 1 :
        massVals['realB'].append( mass.getVal() )
        weightVals['realB'].append( dataSet.weight() )
    if 'realNB' in plotSamples and obsSet.getRealValue(wMCName) > 0. and obsSet.getCatIndex(hlt2BName) == 0 :
        massVals['realNB'].append( mass.getVal() )
        weightVals['realNB'].append( dataSet.weight() )
    if 'MCB' in plotSamples and obsSet.getRealValue(wMCName) < 0. and obsSet.getCatIndex(hlt2BName) == 1:
        massVals['MCB'].append( mass.getVal() )
        weightVals['MCB'].append( dataSet.weight() )
    if 'MCNB' in plotSamples and obsSet.getRealValue(wMCName) < 0. and obsSet.getCatIndex(hlt2BName) == 0:
        massVals['MCNB'].append( mass.getVal() )
        weightVals['MCNB'].append( dataSet.weight() )

# create weights vs. mass graphs
from ROOT import gStyle
from P2VV.Load import LHCbStyle
gStyle.SetColorModelPS(1)
from array import array
massArrs   = dict( [ ( sample, array( 'd', massVals[sample]   ) ) for sample in plotSamples ] )
weightArrs = dict( [ ( sample, array( 'd', weightVals[sample] ) ) for sample in plotSamples ] )
from ROOT import TGraph, kFullDotLarge, kBlue, kRed, kMagenta, kCyan, kOrange
graphs = dict( [ ( sample, TGraph( len( massArrs[sample] ), massArrs[sample], weightArrs[sample] ) )\
                 for sample in plotSamples ] )
for sample, graph in graphs.iteritems() :
    graph.SetMarkerStyle(kFullDotLarge)
    graph.SetMarkerSize(0.2)
    graphColor = kRed + 3 if sample == 'MCNB' else kOrange + 2 if sample == 'MCB' else kCyan + 2 if sample == 'realNB' else kBlue
    graph.SetMarkerColor(graphColor)

    if wMinMax[0][0] != None : graph.SetMinimum( wMinMax[0][0] )
    if wMinMax[0][1] != None : graph.SetMaximum( wMinMax[0][1] )
    graph.GetXaxis().SetLimits( mass.getMin(), mass.getMax() )
    graph.GetXaxis().SetTitle('m(J/#psi K^{+}K^{-}) [MeV/c^{2}]')
    graph.GetYaxis().SetTitle('signal weight')
    graph.GetXaxis().SetTitleOffset(1.2)
    graph.GetYaxis().SetTitleOffset(1.0)
    graph.GetYaxis().SetNdivisions(6, 5, 0)

# draw graphs
from ROOT import TCanvas, TLine
zeroLine = TLine()
canvs = dict( [ ( sample, TCanvas(sample) ) for sample in plotSamples ] )
if any( sample in plotSamples for sample in [ 'realB', 'realNB', 'MCB', 'MCNB' ] ) :
    plotSamples.append('allSplit')
    canvs['allSplit'] = TCanvas('allSplit')
for sample, canv in canvs.iteritems() :
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)
    canv.cd()
    if sample == 'allSplit' :
        drawOpts = 'AP'
        for subSample in [ 'realB', 'realNB', 'MCB', 'MCNB' ] :
            if not subSample in plotSamples : continue
            graphs[subSample].Draw(drawOpts)
            drawOpts = 'P SAMES'
        sample = 'realB'
    else :
        graphs[sample].Draw('AP')
    zeroLine.DrawLine( graphs[sample].GetXaxis().GetXmin(), 0., graphs[sample].GetXaxis().GetXmax(), 0. )

for it, sample in enumerate(plotSamples) :
    canvs[sample].Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len(plotSamples) - 1 else '' ) )

# draw zoomed graphs
for graph in graphs.itervalues() :
    graph.GetXaxis().SetLimits( mass.getMin('Signal'), mass.getMax('Signal') )
    if wMinMax[1][0] != None : graph.SetMinimum( wMinMax[1][0] )
    if wMinMax[1][1] != None : graph.SetMaximum( wMinMax[1][1] )
canvs = dict( [ ( sample, TCanvas( sample + '_zoom' ) ) for sample in plotSamples ] )
for sample, canv in canvs.iteritems() :
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)
    canv.cd()
    if sample == 'allSplit' :
        drawOpts = 'AP'
        for subSample in [ 'realB', 'realNB', 'MCB', 'MCNB' ] :
            if not subSample in plotSamples : continue
            graphs[subSample].Draw(drawOpts)
            drawOpts = 'P SAMES'
        sample = 'realB'
    else :
        graphs[sample].Draw('AP')
    zeroLine.DrawLine( graphs[sample].GetXaxis().GetXmin(), 0., graphs[sample].GetXaxis().GetXmax(), 0. )

fPath = plotsFilePath.split('.')
for it, sample in enumerate(plotSamples) :
    canvs[sample].Print( fPath[0] + '_zoom.' + fPath[1] + ( '(' if it == 0 else ')' if it == len(plotSamples) - 1 else '' ) )
