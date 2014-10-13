model             = 'polarDep'
plotComp          = ''
createDataSets    = False
plotData          = True
plotPDFInt        = False
plotPDF           = True
blind             = False
applyDilWeights   = True
applyAngWeights   = '' # 'para_perp'
dilutionCut       = 0.01
dataPath          = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
plotsFilePath     = 'asymPlot_%s%s_realBins_points.pdf'  % ( model, '_' + plotComp if plotComp else '' ) # 'asymPlot_08bins_m1p590_f0050_b.ps' # 'asymPlot_08bins_f0050_b.ps' # 'asymPlot_04bins_f00255075_b.ps'
plotsROOTFilePath = 'asymPlot_%s%s_realBins_points.root' % ( model, '_' + plotComp if plotComp else '' )
dataSetName       = 'JpsiKK_sigSWeight'
dataSetFileIn     = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'
dataSetFileOut    = 'asymmetryData.root'
pdfFilePaths      = 'pdfVals/pdfVals_%s%s_08bins_6p760_b*_p*_f100p_b.par'
#pdfFilePaths      = 'pdfVals/pdfVals_%s_08bins_m3p234_b*_p*_f100p_b.par' % model
#pdfFilePaths      = 'pdfVals/pdfVals_08bins_m1p616_b*_p*_f0050_b.par'
timeFracs         = [ '%.2f' % ( float(it) / 100. ) for it in range(100) ] # [ '0.00', '0.50' ] # [ '0.00', '0.25', '0.50', '0.75' ]
labelText = '' # 'LHCb'
if blind :
    addText = 'blind: C #approx +0.04; S #approx -0.07'
else :
    addText = ''

from math import pi, sqrt, sin, cos, atan2
numTimeBins = 8
asymMax     = 0.07
#Deltam      = 17.723
Deltam      = 17.69654
oscPeriod   = 2. * pi / Deltam

#phis        = -0.057 #+0.073
#lamb        = 0.963
#AperpSq     = 0.251
#fS          = ( 2443. * 0.426 + 9885. * 0.0587 + 34190. * 0.0096 + 28012. * 0.0094 + 12428. * 0.0482 + 7676. * 0.1920 ) / 94634.
#binOffset   = -0.5 - atan2( 1. - lamb**2, 2. * ( (1. - fS) * (1. - 2. * AperpSq) - fS ) * lamb * sin(phis) ) / 2. / pi * float(numTimeBins)
#periodShift = 0

binWidth = 2. * pi / Deltam / float(numTimeBins)
binOffset = 0.3 / binWidth
periodShift = -1

# check plot component
assert plotComp in [ '', 'even', 'odd', 'S' ]

# check period shift
assert type(periodShift) == int
periodShift = float(periodShift)

# function to set graph attributes
def setGraphAtts( graph, colour, markerSize, minimum, maximum ) :
    from ROOT import kFullDotLarge
    graph.SetMarkerStyle(kFullDotLarge)
    graph.SetLineColor(colour)
    graph.SetMarkerColor(colour)
    graph.SetLineWidth(3)
    graph.SetMarkerSize(markerSize)
    graph.GetXaxis().SetTitle('Decay time (modulo 2#pi/#Deltam_{s}) [ps]')
    graph.GetYaxis().SetTitle('B/#bar{B}-tag asymmetry')
    graph.GetXaxis().SetTitleOffset(1.2)
    graph.GetYaxis().SetTitleOffset(1.2)
    graph.GetXaxis().SetLimits( ( periodShift + binOffset / float(numTimeBins) ) * oscPeriod
                               , ( periodShift + ( 1. + binOffset / float(numTimeBins) ) ) * oscPeriod )
    graph.SetMinimum(minimum)
    graph.SetMaximum(maximum)
    graph.GetXaxis().SetNdivisions( 10, 5, 0, True )
    graph.GetYaxis().SetNdivisions(  5, 5, 0, True )

from ROOT import TGraphErrors as Graph
Graph.setAtts = setGraphAtts

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

from ROOT import TFile, RooDataSet, TObject
from P2VV.RooFitWrappers import __dref__
if createDataSets :
    # read data set from file
    dataFile = TFile.Open(dataSetFileIn)
    dataSet = dataFile.Get(dataSetName)
    dataFile.Close()
    dataSet.Print()

    # get set of observables in data set
    from ROOT import RooArgSet
    from P2VV.RooFitWrappers import RealVar, Category
    obsSet = RooArgSet( dataSet.get() )
    sigWeight  = RealVar( Name = 'sigWeight',  Value = 1. )
    dilution = RealVar( Name = 'dilution', Value = 1. )
    asymCat = Category( Name = 'asymCat', States = dict( plus = +1, minus = -1 ) )
    obsSet.add( __dref__(sigWeight)  )
    obsSet.add( __dref__(dilution) )
    obsSet.add( __dref__(asymCat)    )
    if applyAngWeights :
        angWeight = RealVar( Name = 'angWeight_%s' % applyAngWeights, Value = 1. )
        obsSet.add( __dref__(angWeight)  )

    # create data set with events in two asymmetry categories
    print 'plotAsymmetry: creating dataset with events in two asymmetry categories'
    dataSetAsym = RooDataSet( 'asymData', 'asymData', obsSet )
    from math import exp
    for evSet in dataSet :
        # get dilution from resolution
        sigmat = 1.45 * evSet.getRealValue('sigmat')
        resDil = exp( -0.5 * Deltam**2 * sigmat**2 )

        # get dilution from tagging
        tagCat = ( int( evSet.getCatIndex('tagCatP2VVOS') > 0 ), int( evSet.getCatIndex('tagCatP2VVSS') > 0 ) )
        tags   = ( evSet.getCatIndex('iTagOS'),       evSet.getCatIndex('iTagSS')       )
        assert abs(tags[0]) == 1 and abs(tags[1]) == 1

        tag = tags[0]
        tagDil = 0.
        if tagCat[0] == 1 :
            tagDil = 1. - 2. * evSet.getRealValue('tagomega_os_cb')
            if tagCat[1] == 1 :
                tagDilAlt = 1. - 2. * evSet.getRealValue('tagomega_ss_nn')
                if tags[0] == tags[1] :
                    tagDil = ( tagDil + tagDilAlt ) / ( 1 + tagDil * tagDilAlt )
                else :
                    tagDil = ( tagDil - tagDilAlt ) / ( 1 - tagDil * tagDilAlt )
                    if tagDil < 0. :
                        tag = -tag
                        tagDil = -tagDil
        elif tagCat[1] == 1 :
            tag = tags[1]
            tagDil = 1. - 2. * evSet.getRealValue('tagomega_ss_nn')

        if resDil * tagDil < dilutionCut : continue

        # set signal weight
        sigWeight.setVal( dataSet.weight() )

        # set tagging observables
        dilution.setVal( resDil * tagDil )
        asymCat.setIndex(tag)

        # calculate angular weight
        if applyAngWeights :
            ctk = evSet.getRealValue('helcosthetaK')
            ctl = evSet.getRealValue('helcosthetaL')
            phi = evSet.getRealValue('helphi')
        if applyAngWeights == 'ang' :
            angWeight.setVal( 2. - 5. * ( 1. - ctl**2 ) * sin(phi)**2 )
        elif applyAngWeights == 'para_perp' :
            angWeight.setVal( ( 9. - 15. * ctk**2 ) * ( 1. - ctl**2 ) * sin(2. * phi) )
            #angWeight.setVal( ( 1. - ctl**2 ) * sin(2. * phi) )

        # add event to dataset
        dataSetAsym.add(obsSet)

    dataFile = TFile.Open( dataSetFileOut, 'RECREATE' )
    dataFile.Add(dataSetAsym)
    dataFile.Write( dataSetFileOut, TObject.kOverwrite )
    dataFile.Close()

else :
    # read data set with events in two asymmetry categories
    print 'plotAsymmetry: reading dataset with events in two asymmetry categories'
    dataFile = TFile.Open(dataSetFileOut)
    dataSetAsym = dataFile.Get('asymData')
    dataFile.Close()

dataSetAsym.Print()

# create weighted data set
from ROOT import RooProduct, RooArgSet, RooArgList
obsSet = RooArgSet( dataSetAsym.get() )
prodList = RooArgList( obsSet.find('sigWeight') )
if applyDilWeights : prodList.add( obsSet.find('dilution') )
if applyAngWeights : prodList.add( obsSet.find( 'angWeight_%s' % applyAngWeights ) )
weightVar = RooProduct( 'weightVar', 'weightVar', prodList )
weightVar = dataSetAsym.addColumn(weightVar)
obsSet.add(weightVar)

dataSetAsymW = RooDataSet( 'asymDataW', 'asymDataW', obsSet, Import = dataSetAsym, WeightVar = ( 'weightVar', True ) )
del dataSetAsym
ws.put(dataSetAsymW)
del dataSetAsymW
dataSetAsymW = ws['asymDataW']
obsSet = RooArgSet( dataSetAsymW.get() )
dataSetAsymW.Print()

# get sums of weights
sumW = dict( plus = 0., minus = 0. )
for evSet in dataSetAsymW :
    if evSet.getCatIndex('asymCat') == 1 :
        sumW['plus'] += dataSetAsymW.weight()
    else :
        sumW['minus'] += dataSetAsymW.weight()
assert abs( dataSetAsymW.sumEntries() - sumW['plus'] - sumW['minus'] ) < 1.e-5
if applyAngWeights :
    ASumW = 2. * sumW['plus'] - ( sumW['plus'] + sumW['minus'] )
else :
    ASumW = 2. * sumW['plus'] / ( sumW['plus'] + sumW['minus'] ) - 1.

# create arrays of time bins
from array import array
timeArr = array( 'd', [ ( periodShift + ( float(it) + binOffset + 0.5 ) / float(numTimeBins) ) * oscPeriod\
                       for it in range( numTimeBins + 1 ) ] )
timeErrArr = array( 'd', [ 0.5 / float(numTimeBins) * oscPeriod ] * ( numTimeBins + 1 ) )
timeArrPdf = array( 'd', [ ( periodShift + ( float(it) + binOffset + float(frac) ) / float(numTimeBins) ) * oscPeriod\
                           for it in range( numTimeBins + 1 ) for frac in timeFracs ]\
                         + [ ( periodShift + 1. + ( 1. + binOffset + float( timeFracs[0] ) ) / float(numTimeBins) ) * oscPeriod ] )
timeErrArrPdf = array( 'd', [ 0. ] * len(timeArrPdf) )

pdfBlindVals    = { }
pdfIntBlindVals = { }
if blind or plotPDFInt or plotPDF :
    # function to read files with PDF values
    def readPDFValues(filePathsString) :
        from glob import glob
        filePaths = glob(filePathsString)
        _pdfIntVals = dict( plus = [ 0. ] * numTimeBins, minus = [ 0. ] * numTimeBins )
        _pdfIntBlindVals = dict( plus = [ 0. ] * numTimeBins, minus = [ 0. ] * numTimeBins )
        _pdfVals = dict(  plus  = [ dict( ( frac, 0. ) for frac in timeFracs ) for bin in range(numTimeBins) ]
                        , minus = [ dict( ( frac, 0. ) for frac in timeFracs ) for bin in range(numTimeBins) ]
                       )
        _pdfBlindVals = dict(  plus  = [ dict( ( frac, 0. ) for frac in timeFracs ) for bin in range(numTimeBins) ]
                             , minus = [ dict( ( frac, 0. ) for frac in timeFracs ) for bin in range(numTimeBins) ]
                            )
        fileCount = 0
        for path in filePaths :
            valFile = open(path)
            bin = -1
            plusIntBlindVal = None
            minIntBlindVal  = None
            readBlindPdf    = False
            while True :
                line = valFile.readline().split()
                if not line : break
                if bin < 0 :
                    assert line[0] == 'bin'
                    bin = int( line[1] )
                elif line[0] == 'integrals' :
                    plusIntVal = float( line[2] )
                    minIntVal  = float( line[3] )
                    if len(line) > 4 :
                        plusIntBlindVal = float( line[4] )
                        minIntBlindVal  = float( line[5] )
                else :
                    _pdfVals['plus'][bin][ line[0] ]  += float( line[3] )
                    _pdfVals['minus'][bin][ line[0] ] += float( line[4] )
                    if len(line) > 5 :
                        readBlindPdf = True
                        _pdfBlindVals['plus'][bin][ line[0] ]  += float( line[5] )
                        _pdfBlindVals['minus'][bin][ line[0] ] += float( line[6] )

            _pdfIntVals['plus'][bin]  += plusIntVal
            _pdfIntVals['minus'][bin] += minIntVal
            if blind or readBlindPdf :
                assert plusIntBlindVal != None and minIntBlindVal != None and readBlindPdf
                _pdfIntBlindVals['plus'][bin]  += plusIntBlindVal
                _pdfIntBlindVals['minus'][bin] += minIntBlindVal

            valFile.close()
            fileCount += 1
        print 'plotAsymmetry: read data from %d PDF values files (%s)' % ( fileCount, filePathsString )
        return ( _pdfIntVals, _pdfIntBlindVals, _pdfVals, _pdfBlindVals )

    ( pdfIntVals, pdfIntBlindVals, pdfVals, pdfBlindVals ) = readPDFValues( pdfFilePaths % ( model, '_' + plotComp if plotComp else '' ) )
    if plotComp :
        ( pdfIntValsTot, pdfIntBlindValsTot, pdfValsTot, pdfBlindValsTot ) = readPDFValues( pdfFilePaths % ( model, '' ) )

from ROOT import gStyle
from P2VV.Load import LHCbStyle
gStyle.SetColorModelPS(1)
pdfGraph = None
if plotPDF :
    # plot PDF asymmetry in time bins
    if blind or not pdfBlindVals :
        vals = pdfVals
        valsTot = pdfValsTot if plotComp else pdfVals
    else :
        vals = pdfBlindVals
        valsTot = pdfBlindValsTot if plotComp else pdfBlindVals
        
    pdfArr = array( 'd', [ ( pVals[frac] - mVals[frac] ) / ( pValsTot[frac] + mValsTot[frac] )\
                           for pVals, mVals, pValsTot, mValsTot in zip( vals['plus'], vals['minus'], valsTot['plus'], valsTot['minus'] )\
                           for frac in timeFracs ] )
    for it in range( len(timeFracs) + 1 ) : pdfArr.append( pdfArr[it] )
    pdfErrArr = array( 'd', [ 0. ] * len(pdfArr) )
    pdfGraph = Graph( len(timeArrPdf), timeArrPdf, pdfArr, timeErrArrPdf, pdfErrArr )
    pdfGraph.SetName('pdf')

pdfIntGraph = None
if plotPDFInt :
    # plot integrated-PDF asymmetry in time bins
    if blind or not pdfIntBlindVals :
        intVals = pdfIntVals
    else :
        intVals = pdfIntBlindVals

    pdfIntArr = array( 'd', [ 2. * pVal / ( pVal + mVal ) - 1. for pVal, mVal in zip( intVals['plus'], intVals['minus'] ) ] )
    pdfIntArr.append( pdfIntArr[0] )
    pdfIntErrArr = array( 'd', [ 0. ] * len(pdfIntArr) )
    pdfIntGraph = Graph( len(timeArr), timeArr, pdfIntArr, timeErrArr, pdfIntErrArr )
    pdfIntGraph.SetName('pdfInt')

dataGraph = None
if plotData :
    # get data asymmetries in time bins
    eventSums = dict( [ ( var, [ 0. ] * numTimeBins ) for var in ( 'm0', 'm1', 'n0', 'n1' ) ] )
    timeBins = [ ( 1. + float(it) + binOffset ) / float(numTimeBins) * oscPeriod for it in range(numTimeBins) ]
    for evSet in dataSetAsymW :
        time = evSet.getRealValue('time')
        bin = 0
        iter = 0
        while time >= float(iter) * oscPeriod + timeBins[bin] :
            if bin < len(timeBins) - 1 :
                bin += 1
            else :
                bin = 0
                iter += 1

        weight = dataSetAsymW.weight()
        eventSums['m0'][bin] += weight
        eventSums['n0'][bin] += weight**2
        if evSet.getCatIndex('asymCat') == 1 :
            eventSums['m1'][bin] += weight
            eventSums['n1'][bin] += weight**2
    assert abs( sumW['plus'] - sum( eventSums['m1'] ) ) < 1.e-5 and abs( sumW['plus'] + sumW['minus'] - sum( eventSums['m0'] ) ) < 1.e-5

    # plot data asymmetry as a function of time
    if applyAngWeights :
        dataArr = array( 'd', [ 2. * m1 - m0 for m0, m1 in zip( eventSums['m0'], eventSums['m1'] ) ] )
        dataErrArr = array( 'd', [ sqrt(n0) for n0 in eventSums['n0'] ] )
    else :
        dataArr = array( 'd', [ 2. * m1 / m0 - 1. for m0, m1 in zip( eventSums['m0'], eventSums['m1'] ) ] )
        if blind :
            for bin, ( plusVal, minVal, plusValB, minValB )\
                    in enumerate( zip( pdfIntVals['plus'], pdfIntVals['minus'], pdfIntBlindVals['plus'], pdfIntBlindVals['minus'] ) ) :
                dataArr[bin] += 2. * ( plusVal / ( plusVal + minVal ) - plusValB / ( plusValB + minValB ) )
        dataErrArr = array( 'd', [ 2. / m0 * sqrt( n0 * ( 1. + nu )**2 / 4. - n1 * nu )\
                                  for nu, m0, n0, n1 in zip( dataArr, eventSums['m0'], eventSums['n0'], eventSums['n1'] ) ] )
    dataArr.append( dataArr[0] )
    dataErrArr.append( dataErrArr[0] )
    dataGraph = Graph( len(timeArr), timeArr, dataArr, timeErrArr, dataErrArr )
    dataGraph.SetName('data')
    if not asymMax :
        asymMax = 1.1 * max( max( abs( val + err ), abs( val - err ) ) for val, err in zip( dataArr, dataErrArr ) )

# set graph attributes
if not asymMax : asymMax = 1.
from ROOT import kBlack, kBlue, kRed
dummyGraph = Graph( numTimeBins + 1 )
dummyGraph.SetName('dummy')
dummyGraph.setAtts( kBlack, 0., -asymMax, +asymMax )
if dataGraph : dataGraph.setAtts( kBlack, 1., -asymMax, +asymMax )
if pdfIntGraph : pdfIntGraph.setAtts( kBlue, 1., -asymMax, +asymMax )
if pdfGraph : pdfGraph.setAtts( kRed if pdfIntGraph else kBlue, 0.1, -asymMax, +asymMax )

# create line
from ROOT import TLine, kDotted
dotLine = TLine()
dotLine.SetLineStyle(kDotted)

# create label
from ROOT import TLatex
label = TLatex()
label.SetTextAlign(12)
label.SetTextSize(0.072)

if addText :
    # create additional text
    text = TLatex()
    text.SetTextAlign(22)
    text.SetTextSize(0.060)

# draw graphs
tMin = dummyGraph.GetXaxis().GetXmin()
tMax = dummyGraph.GetXaxis().GetXmax()
from ROOT import TCanvas
canv = TCanvas('canv')
canv.SetLeftMargin(0.18)
canv.SetRightMargin(0.05)
canv.SetBottomMargin(0.18)
canv.SetTopMargin(0.05)
dummyGraph.Draw('AP')
horLine = dotLine.DrawLine( tMin, ASumW, tMax, ASumW )
#vertLine0 = dotLine.DrawLine( 0., dummyGraph.GetMinimum(), 0., dummyGraph.GetMaximum() )
#vertLine1 = dotLine.DrawLine( oscPeriod, dummyGraph.GetMinimum(), oscPeriod, dummyGraph.GetMaximum() )
if dataGraph : dataGraph.Draw('P SAMES')
if pdfIntGraph : pdfIntGraph.Draw('P SAMES')
if pdfGraph : pdfGraph.Draw('L SAMES')
if labelText : label.DrawLatex( tMin + 0.77 * ( tMax - tMin ), 0.76 * asymMax, labelText )
if addText : text.DrawLatex( 0.5 * ( tMax + tMin ), -0.7 * asymMax, addText )
canv.Print(plotsFilePath)

# save graphs and canvas
plotsROOTFile = TFile.Open( plotsROOTFilePath, 'RECREATE' )
plotsROOTFile.Add(dummyGraph)
plotsROOTFile.Add(horLine)
#plotsROOTFile.Add(vertLine0)
#plotsROOTFile.Add(vertLine1)
if dataGraph : plotsROOTFile.Add(dataGraph)
if pdfIntGraph : plotsROOTFile.Add(pdfIntGraph)
if pdfGraph : plotsROOTFile.Add(pdfGraph)
plotsROOTFile.Add(canv)
plotsROOTFile.Write( plotsROOTFilePath, TObject.kOverwrite)
plotsROOTFile.Close()
