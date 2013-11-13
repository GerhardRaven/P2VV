testParName = 'ASOddPhase_bin1'
testParTitle = '#delta_{S_{1}}'

nllFilePaths = [ 'nllPars_12345.par' ]
fitFilePaths = [ 'fitPars_1234%d_00%d.par' % ( id, it ) for id in range( 5, 7 ) for it in range(3) ]
plotFilePath = 'nllPlots.ps'
parRange = ( 0.5, 2.8 )
nllRange = ( 0., 8. )

nPointsPara = 10000
meanPara = 2.32
sigmaPara = 0.20

# get parabola values
parValsPara = array( 'd', [  ] )
nllValsPara = array( 'd', [ ] )


# get NLL values
parValsNLL = [ ]
nllVals = [ ]
for filePath in nllFilePaths :
    try :
        nllFile = open(filePath)
    except :
        print 'plotNLL: ERROR: unable to open file "%s"' % filePath

    valPos = -1
    while True :
        line = nllFile.readline()
        if not line : break

        line = line.strip().split()
        if valPos >= 0 :
            parValsNLL.append( float( line[0] ) )
            nllVals.append( float( line[2] ) )

        elif len(line) >= 3 and line[0] == 'test' and line[1] == 'parameters:' and testParName in line[ 2 : ] :
            valPos = line.index(testParName) - 2

    if valPos < 0 :
        print 'plotNLL: ERROR: file "%s" is not an NLL file for parameter "%s"' % ( filePath, testParName )

from array import array
if parValsNLL :
    minNLL = min(nllVals)
    parValsNLLArr = array( 'd', parValsNLL )
    nllValsArr = array( 'd', [ val - minNLL for val in nllVals ] )

# get fit values
parValsFit = [ ]
nllValsFit = [ ]

from P2VV.Parameterizations.FullPDFs import PdfConfiguration
for filePath in fitFilePaths :
    pdfConfig = PdfConfiguration()
    fitStatus = pdfConfig.readParametersFromFile( filePath, Verbose = False )
    assert fitStatus[0] == 0, 'plotNLL: ERROR: fit status is "%d" for fit "%s"' % ( fitStatus[0], filePath )
    assert testParName in pdfConfig.parameters(), 'plotNLL: ERROR: "%s" not in parameters of file "%s"' % ( testParName, filePath )

    parValsFit.append( pdfConfig.parameters()[testParName][0] )
    nllValsFit.append( fitStatus[1] )

if parValsFit :
    minNLL = min(nllValsFit)
    parValsFitArr = array( 'd', parValsFit )
    nllValsFitArr = array( 'd', [ val - minNLL for val in nllValsFit ] )

# draw graphs
from P2VV.Load import LHCbStyle

graphs = [ ]
from ROOT import TGraph
graphs.append( None )
graphs.append( TGraph( len(parValsNLL), parValsNLLArr, nllValsArr )    if parValsNLL else None )
graphs.append( TGraph( len(parValsFit), parValsFitArr, nllValsFitArr ) if parValsFit else None )

from ROOT import kBlack, kRed, kBlue
for graphIt, graph in enumerate(graphs) :
    if not graph : continue

    graph.SetLineWidth(3)
    graph.SetLineColor( kBlue if graphIt == 1 else kRed if graphIt == 2 else kBlack )
    graph.SetMinimum( nllRange[0] )
    graph.SetMaximum( nllRange[1] )
    graph.GetXaxis().SetLimits( parRange[0], parRange[1] )
    graph.GetXaxis().SetTitle(testParTitle)
    graph.GetYaxis().SetTitle('#Delta log(L)')
    graph.GetXaxis().SetTitleOffset(1.1)
    graph.GetYaxis().SetTitleOffset(0.9)

from ROOT import TCanvas
canv = TCanvas('canv')
canv.SetLeftMargin(0.18)
canv.SetRightMargin(0.05)
canv.SetBottomMargin(0.18)
canv.SetTopMargin(0.05)

firstGraph = True
for graph in graphs :
    if graph and firstGraph :
        graph.Draw('AC')
        firstGraph = False
    elif graph :
        graph.Draw('SAMES C')
canv.Print(plotFilePath)
