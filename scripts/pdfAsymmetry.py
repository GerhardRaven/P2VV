dataPath         = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetFilePath  = 'asymmetryData.root'
parFilePath      = '20112012Reco14DataPlotValues_6KKMassBins.par'
applyPlotWeights = True
blindVars        = [ 'phiCP', 'dGamma' ]

from math import pi, sqrt
Deltam         = 17.762
numTimeBinsTot = 8
timeBins       = [ 0 ]
periods        = [ 0 ]
binOffset      = -0.5
timeFracs      = [ 0., 0.5 ]

import sys
if len(sys.argv) > 1 :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( 'numTimeBinsTot', type = int )
    parser.add_argument( 'startBin',       type = int )
    parser.add_argument( 'numTimeBins',    type = int )
    parser.add_argument( 'startPeriod',    type = int )
    parser.add_argument( 'numPeriods',     type = int )
    parser.add_argument( 'binOffset',      type = float )
    parser.add_argument( 'Deltam',         type = float )
    parser.add_argument( 'timeFracs',      type = float, nargs = '+' )

    args = parser.parse_args()
    numTimeBinsTot = args.numTimeBinsTot
    timeBins = [ args.startBin + it for it in range(args.numTimeBins) ]
    periods = [ args.startPeriod + it for it in range(args.numPeriods) ]
    binOffset = args.binOffset
    Deltam = args.Deltam
    timeFracs = [ val for val in args.timeFracs ]

jobID = '%02dbins_%s_b%s_p%s_f%s%s' % ( numTimeBinsTot, ( '%.3f' % binOffset ).replace( '-', 'm' ).replace( '.', 'p' )
                                       , ''.join( '%02d' % bin for bin in timeBins ), ''.join( '%02d' % per for per in periods )
                                       , ''.join( '%02.0f' % ( 100. * frac ) for frac in timeFracs ), '_b' if blindVars else '' )
pdfValsFilePath = 'pdfVals/pdfVals_%s.par' % jobID
oscPeriod = 2. * pi / Deltam

print 'pdfAsymmetry: job %s parameters:' % jobID
print '  dataset file: ' + dataSetFilePath
print '  parameter file: ' + parFilePath
print '  apply plot weights: %s' % ( 'yes' if applyPlotWeights else 'no' )
print '  blinded parameters: %s' % ', '.join( par for par in blindVars )
print '  Deltam = %.5f / ps' % Deltam
print '  oscillation period = %.3f ps' % oscPeriod
print '  total number of time bins: %d' % numTimeBinsTot
print '  time bins: %s' % ', '.join( '%d' % bin for bin in timeBins )
print '  periods: %s' % ', '.join( '%d' % per for per in periods )
print '  time offset (bins) = %.3f' % binOffset
print '  PDF time points: %s' % ', '.join( '%.2f' % frac for frac in timeFracs )
print '  output file: ' + pdfValsFilePath

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set with events in two asymmetry categories
print 'pdfAsymmetry: reading dataset with events in two asymmetry categories'
from ROOT import TFile
dataFile = TFile.Open(dataSetFilePath)
dataSetAsym = dataFile.Get('asymData')
dataFile.Close()
dataSetAsym.Print()

# create weighted data set
from ROOT import RooProduct, RooArgSet, RooArgList
obsSet = RooArgSet( dataSetAsym.get() )
prodList = RooArgList( obsSet.find('sigWeight') )
if applyPlotWeights : prodList.add( obsSet.find('dilution') )
weightVar = RooProduct( 'weightVar', 'weightVar', prodList )
weightVar = dataSetAsym.addColumn(weightVar)
obsSet.add(weightVar)

from ROOT import RooDataSet
dataSetAsymW = RooDataSet( 'asymDataW', 'asymDataW', obsSet, Import = dataSetAsym, WeightVar = ( 'weightVar', True ) )
del dataSetAsym
ws.put(dataSetAsymW)
del dataSetAsymW
dataSetAsymW = ws['asymDataW']
obsSet = RooArgSet( dataSetAsymW.get() )
dataSetAsymW.Print()

# build PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig( RunPeriods = '3fb' )

timeEffFile2011 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
timeEffFile2012 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeEffFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeEffFile2012
pdfConfig['anglesEffType'] = 'basisSig4'
pdfConfig['angEffMomsFiles'] = dataPath + 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Basis'
pdfConfig['signalData'] = dataSetAsymW
pdfConfig['readFromWS'] = True

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()
pdfConfig.readParametersFromFile( filePath = parFilePath )
pdfConfig.setParametersInPdf(pdf)

# create projection data sets
projSet = [ obs for obs in pdf.ConditionalObservables() ] + [ cat for cat in pdf.indexCat().inputCatList() ]
projSet = RooArgSet( obsSet.find( var.GetName() ) for var in projSet )
projDataSets = dict(  plus  = dataSetAsymW.reduce( Name = 'asymDataWPlus',  SelectVars = projSet, Cut = 'asymCat==+1' )
                    , minus = dataSetAsymW.reduce( Name = 'asymDataWMinus', SelectVars = projSet, Cut = 'asymCat==-1' )
                   )
print 'pdfAsymmetry: projection data sets:'
projDataSets['plus'].Print()
projDataSets['minus'].Print()

# create time ranges for time bins
rangeNames = dict( )
timePoints = dict( )
lastBin = False
if not timeBins : timeBins = range(8)
if not periods : periods = range( int( ws['time'].getMax() / oscPeriod - binOffset / float(numTimeBinsTot) ) + 1 )
for per in periods :
    for bin in timeBins :
        timeLo = ( float(per) + ( float(bin) + binOffset ) / float(numTimeBinsTot) ) * oscPeriod
        timeHi = timeLo + oscPeriod / float(numTimeBinsTot)
        if timeLo < ws['time'].getMin() :
            timeLo = ws['time'].getMin()
        if timeHi <= ws['time'].getMin() :
            continue
        if timeHi > ws['time'].getMax() :
            timeHi = ws['time'].getMax()
            lastBin = True
        ws['time'].setRange( 'bin%02d_%02d' % ( bin, per ), timeLo, timeHi )
        if bin not in rangeNames : rangeNames[bin] = [ ]
        rangeNames[bin].append( 'bin%02d_%02d' % ( bin, per ) )

        if bin not in timePoints : timePoints[bin] = [ [ ] for it in range( len(timeFracs) ) ]
        for it, frac in enumerate(timeFracs) :
            timeVal = ( float(per) + ( float(bin) + binOffset + frac ) / float(numTimeBinsTot) ) * oscPeriod
            if timeVal < ws['time'].getMin() * ( 1. - 1.e-7 ) or timeVal > ws['time'].getMax() * ( 1. + 1.e-7 ) : continue
            timePoints[bin][it].append(timeVal)

        if lastBin : break

    if lastBin : break

#t = ws['time']
#for bin, ( rVals, tVals ) in enumerate( zip( rangeNames.itervalues(), timePoints.itervalues() ) ) :
#  print bin,
#  for val in rVals :
#      print '%s: %.3f - %.3f ' % ( val, t.getMin(val), t.getMax(val) ),
#  for vals in tVals :
#      print '[%s]' % ', '.join( '%.3f' % val for val in vals ),
#  print

# create PDF integrals in bins for plus and minus categories
emptySet = RooArgSet()
timeSet  = RooArgSet( ws['time'] )
angleSet = RooArgSet( ws[varName] for varName in [ 'helcosthetaK', 'helcosthetaL', 'helphi' ] )
intSet   = RooArgSet(timeSet)
intSet.add(angleSet)
from ROOT import RooExplicitNormPdf
asymPdfs = dict(  plus  = RooExplicitNormPdf( 'pdfPlus', 'pdfPlus', timeSet, angleSet, pdf._var, pdf._var
                                             , projDataSets['plus'].sumEntries() / dataSetAsymW.sumEntries(), projDataSets['plus'] )
                , minus = RooExplicitNormPdf( 'pdfMinus', 'pdfMinus', timeSet, angleSet, pdf._var, pdf._var
                                             , projDataSets['minus'].sumEntries() / dataSetAsymW.sumEntries(), projDataSets['minus'] )
               )
asymPdfInts = dict(  plus  = dict( [ ( bin, RooExplicitNormPdf( 'pdfPlus%02d' % bin, 'pdfPlus%02d' % bin, emptySet, intSet, pdf._var
                                                               , pdf._var, projDataSets['plus'].sumEntries() / dataSetAsymW.sumEntries()
                                                               , projDataSets['plus'], ','.join( name for name in rangeNames[bin] )
                                                              ) ) for bin in sorted( rangeNames.keys() )
                                 ] )
                   , minus = dict( [ ( bin, RooExplicitNormPdf( 'pdfMinus%02d' % bin, 'pdfMinus%02d' % bin, emptySet, intSet, pdf._var
                                                               , pdf._var, projDataSets['minus'].sumEntries() / dataSetAsymW.sumEntries()
                                                               , projDataSets['minus'], ','.join( name for name in rangeNames[bin] )
                                                              ) ) for bin in sorted( rangeNames.keys() )
                                 ] )
                  )

print 'pdfAsymmetry: PDF parameters:'
asymPdfs['plus'].parSet().Print('v')

def getPdfVals() :
    pdfVals = dict(  plus  = dict( [ ( bin, [ [ ] for it in range( len(timeFracs) ) ] ) for bin in rangeNames.iterkeys() ] )
                   , minus = dict( [ ( bin, [ [ ] for it in range( len(timeFracs) ) ] ) for bin in rangeNames.iterkeys() ] )
                  )
    for bin, tPoints in timePoints.iteritems() :
        for pointIt, tVals in enumerate(tPoints) :
            for tVal in tVals :
                ws['time'].setVal(tVal)
                pdfVals['plus'][bin][pointIt].append( asymPdfs['plus'].getVal() )
                pdfVals['minus'][bin][pointIt].append( asymPdfs['minus'].getVal() )
    return pdfVals

asymPdfBlindIntVals = { }
asymPdfBlindVals = { }
if blindVars :
    bName = lambda name : '__%s__' % name
    blindVals = [ 2. * ws[ bName(name) ].getVal() - ws[name].getVal() for name in blindVars ]
    asymPdfBlindIntVals = dict(  plus  = dict( [ ( bin, p.getVal() ) for bin, p in asymPdfInts['plus'].iteritems()  ] )
                               , minus = dict( [ ( bin, p.getVal() ) for bin, p in asymPdfInts['minus'].iteritems() ] )
                              )
    asymPdfBlindVals = getPdfVals()
    for name, val in zip( blindVars, blindVals ) : ws[ bName(name) ].setVal(val)

asymPdfIntVals = dict(  plus  = dict( [ ( bin, p.getVal() ) for bin, p in asymPdfInts['plus'].iteritems()  ] )
                      , minus = dict( [ ( bin, p.getVal() ) for bin, p in asymPdfInts['minus'].iteritems() ] )
                     )
asymPdfVals = getPdfVals()

valsFile = open( pdfValsFilePath, 'w' )
for bin in sorted( rangeNames.keys() ) :
    periodStr = ', '.join( '%.5f - %.5f' % ( ws['time'].getMin(name), ws['time'].getMax(name) ) for name in rangeNames[bin] )
    valsFile.write( 'bin %d : %s\n' % ( bin, periodStr ) )
    valsFile.write( '    integrals : %.5g  %.5g' % ( asymPdfIntVals['plus'][bin], asymPdfIntVals['minus'][bin] ) )
    if asymPdfBlindIntVals :
        valsFile.write( '  %.5g  %.5g' % ( asymPdfBlindIntVals['plus'][bin], asymPdfBlindIntVals['minus'][bin] ) )
    valsFile.write('\n')
    for it, frac in enumerate(timeFracs) :
        valsFile.write( '    %.2f (%.5f) : %.5g  %.5g'\
                        % ( frac, ( float(bin) + binOffset + frac ) / float(numTimeBinsTot) * oscPeriod
                           , sum( asymPdfVals['plus'][bin][it] ), sum( asymPdfVals['minus'][bin][it] ) ) )
        if asymPdfBlindVals :
            valsFile.write( '  %.5g  %.5g' % ( sum( asymPdfBlindVals['plus'][bin][it] ), sum( asymPdfBlindVals['minus'][bin][it] ) ) )
        valsFile.write('\n')
valsFile.close()
