###########################################################################################################################################
## set script parameters ##
###########################

# NLL parameters
testParNames = [ '__phiCP_A0__' ]
nllParVals = [ ] #[ [ -0.55 + 0.55 / 1000. * float(it) ] for it in range(1001) ]
fitParVals = [ ] # [ [ 2.8 ], [ 2.9 ], [ 3.1 ], [ 3.2 ] ]
outDirPath = '/project/bfys/jleerdam/softDevel/P2VV2/test/phis_NLL/nllVals/'
jobID = 'phiCP_A0'
startIndex = 0
indexWidth = 4

import sys
if len(sys.argv) > 1 :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('jobID')
    parser.add_argument( 'startIndex', type = int )
    parser.add_argument( 'indexWidth', type = int )
    parser.add_argument('parName')
    parser.add_argument( 'parVals', type = float, nargs = '+' )

    args = parser.parse_args()
    jobID = args.jobID
    startIndex = args.startIndex
    testParNames = [ args.parName ]
    fitParVals = [ [ val ] for val in args.parVals ]

print 'job parameters:'
print '  ID:', jobID
print '  start index:', startIndex
print '  index width:', indexWidth
print '  parameter names:', testParNames
print '  parameter values:', fitParVals

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
parFileIn   = '20112012Reco14DataPlotValues_6KKMassBins_CPVDecay.par'
dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetName = 'JpsiKK_sigSWeight'
dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2DiegoMass_6KKMassBins_2TagCats.root'

# fit options
fitOpts = dict(  NumCPU    = 1
               , Optimize  = 2
               , Minimizer = 'Minuit2'
              )
pdfConfig['fitOptions'] = fitOpts
corrSFitErrCats         = [ 'runPeriod', 'KKMassCat' ]
randomParVals           = ( ) # ( 1., 12345 )

# PDF options
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file']\
        = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file']\
        = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
pdfConfig['angEffMomsFiles'] = dataPath + 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'

pdfConfig['lambdaCPParam'] = 'lambPhi_CPVDecay'


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
pdfConfig['signalData'] = dataSet
pdfConfig['readFromWS'] = True

# build the PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if not 'Optimize' in fitOpts or fitOpts['Optimize'] < 2 :
    # unset cache-and-track
    for par in pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

# data set with weights corrected for background dilution: for phi_s fit only!
from P2VV.Utilities.DataHandling import correctWeights
fitData = correctWeights( dataSet, corrSFitErrCats )

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)

# get test parameters
testPars = [ pdfPars.find(name) for name in testParNames ]


###########################################################################################################################################
## get NLL values ##
####################

# fix values of some parameters
for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
    if not pdfConfig['SSTagging'] :
        CEvenOdds.setConstant('avgCEven.*')
        CEvenOdds.setConstant( 'avgCOdd.*', True )
    else :
        for CEvenOdd in CEvenOdds :
            CEvenOdd.setConstant('avgCEven.*')
            CEvenOdd.setConstant( 'avgCOdd.*', True )

pdfBuild['amplitudes'].setConstant('C_SP')

for par in testPars : par.setConstant(True)

# print parameters
print 120 * '='
print 'Bs2JpsiKKFit: fit data:'
fitData.Print()
print 'Bs2JpsiKKFit: observables in PDF:'
pdfObs.Print('v')
print 'Bs2JpsiKKFit: parameters in PDF:'
pdfPars.Print('v')
print 'Bs2JpsiKKFit: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()

if nllParVals :
    print 120 * '='
    print 'Bs2JpsiKKFit: computing NLL values'

    # create NLL variable
    nll = pdf.createNLL( fitData, **fitOpts )
    #from ROOT import RooAbsReal
    #RooAbsReal.setHideOffset(False)
    #nll.enableOffsetting(True)

    # open output file for NLL values
    try :
        nllFile = open( outDirPath + 'nllPars_' + jobID + '.par', 'w' )
    except :
        raise RuntimeError( 'Bs2JpsiKKFit: ERROR: unable to open file \"%s\"' % ( outDirPath + 'nllPars_' + jobID + '.par' ) )

    # write test-parameter names to file
    nllFile.write('test parameters: ')
    for parIt, par in enumerate(testPars) : nllFile.write( par.GetName() + ( ', ' if parIt < len(testPars) - 1 else '\n' ) )

    for valIt, parVals in enumerate(nllParVals) :
        if valIt % 100 == 0 :
            print 'JpsiKKFit: NLL iteration %d' % valIt
            import sys
            sys.stdout.flush()

        # set test-parameter values
        for parIt, par in enumerate(testPars) :
            par.setVal( parVals[parIt] )
            nllFile.write( '{0:+14.8g} '.format( par.getVal() ) )

        # write NLL value to file
        nllFile.write( ': {0:+.16g}\n'.format( nll.getVal() ) )

    nllFile.close()

# fit data in loop
maxLenFitParName = max( len( par.GetName() ) for par in testPars )
for valIt, parVals in enumerate(fitParVals) :
    print 120 * '='

    if randomParVals :
        # give parameters random offsets before fitting
        import random
        print 'Bs2JpsiKKFit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
              % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
        random.seed( randomParVals[1] if randomParVals[1] else None )
        for par in pdfPars :
            if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

    # set test-parameter values
    print 'Bs2JpsiKKFit: fit iteration %d' % valIt
    print 'Bs2JpsiKKFit: setting test-parameter values'
    for parIt, par in enumerate(testPars) :
        par.setVal( parVals[parIt] )
        print ( '  {0:%ds} = {1:+10.4g}' % maxLenFitParName ).format( par.GetName(), par.getVal() )

    # fit data
    print 'Bs2JpsiKKFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )
    fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, Hesse = False, Offset = False, **fitOpts )

    # print parameter values
    from P2VV.Imports import parNames, parValues
    print 'Bs2JpsiKKFit: parameters:'
    fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )

    # write parameters to file
    from math import log10
    filePF = ( '_{0:0%dd}' % indexWidth ).format( startIndex + valIt )
    pdfConfig.getParametersFromPdf( pdf, fitData )
    pdfConfig.writeParametersToFile( filePath = outDirPath + 'fitPars_' + jobID + filePF + '.par'
                                    , FitStatus = ( fitResult.status(), fitResult.minNll(), fitResult.edm() ) )

