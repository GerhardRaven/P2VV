###########################################################################################################################################
## set script parameters ##
###########################

testParNames = [ 'delSDelta' ]
nllParVals = [ ] # [ [ -3.8 + 4.3 / 1000. * float(it) ] for it in range(1001) ]
fitParVals = [ [ -0.7 ], [ -0.6 ], [ -0.5 ], [ -0.4 ], [ -0.3 ] ]
outDirPath = './nllVals/'
jobID = 'delSDelta'

indexWidth = 4
startIndex = 0

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
parFileIn  = '20112012Reco14DataFitValues_4KKMassBins.par'

dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetName = 'JpsiKK_sigSWeight'
dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_4KKMassBins_2TagCats.root'
deltaSDiff  = ( 'ASOddPhase_bin0', 'ASOddPhase_bin3' )

# fit options
fitOpts = dict(  NumCPU    = 8
               , Optimize  = 2
               , Minimizer = 'Minuit2'
              )
pdfConfig['fitOptions'] = fitOpts
corrSFitErrCats         = [ 'runPeriod', 'KKMassCat' ]
randomParVals           = ( ) # ( 1., 12345 )

# PDF options
pdfConfig['timeResType'] = 'eventNoMean'
pdfConfig['externalConstr']['timeResSigmaSF'] = ( 1.45, 0. )

pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file']\
        = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file']\
        = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
pdfConfig['angEffMomsFiles'] = dataPath + 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'

pdfConfig['KKMassBinBounds'] = [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['CSPValues']       = [ 0.966, 0.797, 0.797, 0.966 ]


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

if deltaSDiff :
    # replace last S-wave phase by difference with first S-wave phase
    from P2VV.RooFitWrappers import RealVar, Addition, EditPdf
    from ROOT import RooNumber
    RooInf = RooNumber.infinity()
    delSDelta = RealVar( Name = 'delSDelta', Value = -1.5, Error = 0.3, MinMax = ( -RooInf, +RooInf ) )
    delSm1New = Addition( Name = 'ASOddPhase_m1', Arguments = [ ws[ deltaSDiff[0] ], ws['delSDelta'] ] )
    pdf = EditPdf( Name = pdf.GetName() + '_delSDelta', Original = pdf, Rules = { ws[ deltaSDiff[1] ] : delSm1New } )

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
## fit data ##
##############

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
        if valIt % 1000 == 0 :
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
    nll.IsA().Destructor(nll)

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

    fitResult.IsA().Destructor(fitResult)
