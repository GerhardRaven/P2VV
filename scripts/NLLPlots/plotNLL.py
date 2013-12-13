amp = 'AS'
testParName = '__phiCP_%s__' % amp
testParTitle = '#phi_{s}^{%s}' % ( '0' if amp == 'A0' else '#parallel' if amp == 'Apar' else '#perp' if amp == 'Aperp' else 'S' )

from glob import glob
nllFilePaths = [ 'nllVals/nllPars_phiCP_%s.par' % amp ]
fitFilePaths = glob( 'nllVals/fitPars_phiCP_%s_*.par' % amp )
plotFilePath = 'phiCP_%s.ps' % amp
parRange = ( -0.60, 1.0 )
nllRange = ( 0., 14. )

nPointsPara = 1000
meanPara =  0.1909084  # 0.1909084  # -0.23525682 # -0.28361833 # -0.26138806
sigmaPara = 0.14208314 # 0.14208314 # 0.066643948 # 0.063555006 # 0.054327098

# get parabola values
from math import sqrt
minXPara = max( parRange[0], meanPara - sigmaPara * sqrt( 2. * nllRange[1] ) )
maxXPara = min( parRange[1], meanPara + sigmaPara * sqrt( 2. * nllRange[1] ) )
from array import array
parValsParaArr = array( 'd', [ minXPara + float(it) / float(nPointsPara - 1) * ( maxXPara - minXPara ) for it in range(nPointsPara) ] )
nllValsParaArr = array( 'd', [ 0.5 * ( ( parVal - meanPara ) / sigmaPara )**2 for parVal in parValsParaArr ] )

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
    if fitStatus[0] != 0 :
        print 'plotNLL: WARNING: fit status is "%d" for fit "%s"' % ( fitStatus[0], filePath )
    assert testParName in pdfConfig.parameters(), 'plotNLL: ERROR: "%s" not in parameters of file "%s"' % ( testParName, filePath )

    parVal = pdfConfig.parameters()[testParName][0]
    pos = 0
    for it in range( len(parValsFit) ) :
        if parVal > parValsFit[it] :
            break
        pos += 1
    parValsFit.insert( pos, parVal )
    nllValsFit.insert( pos, fitStatus[1] )

print parValsFit
if parValsFit :
    minNLL = min(nllValsFit)
    parValsFitArr = array( 'd', parValsFit )
    nllValsFitArr = array( 'd', [ val - minNLL for val in nllValsFit ] )

# draw graphs
from P2VV.Load import LHCbStyle

graphs = [ ]
from ROOT import TGraph
graphs.append( TGraph( len(parValsParaArr), parValsParaArr, nllValsParaArr ) if parValsParaArr else None )
graphs.append( TGraph( len(parValsNLL),     parValsNLLArr,  nllValsArr     ) if parValsNLL     else None )
graphs.append( TGraph( len(parValsFit),     parValsFitArr,  nllValsFitArr  ) if parValsFit     else None )

from ROOT import kBlack, kRed, kBlue, kFullDotLarge
for graphIt, graph in enumerate(graphs) :
    if not graph : continue

    graph.SetLineWidth( 3 if graphIt > 0 else 2 )
    graph.SetMarkerStyle(kFullDotLarge)
    graph.SetMarkerSize(0.6)
    graph.SetLineColor( kBlue if graphIt == 1 else kRed if graphIt == 2 else kBlack )
    graph.SetMarkerColor( kBlue if graphIt == 1 else kRed if graphIt == 2 else kBlack )
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
for it, graph in enumerate(graphs) :
    if graph and firstGraph :
        graph.Draw('AL')
        firstGraph = False
    elif graph :
        graph.Draw( 'SAMES' + ( ' PC' if it == 2 else ' L' ) )
canv.Print(plotFilePath)
