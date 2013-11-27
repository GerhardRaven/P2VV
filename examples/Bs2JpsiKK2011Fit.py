###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
doFit = True

parFileIn  = ''
parFileOut = ''

dataSetName = 'JpsiKK_sigSWeight'
dataSetFile = 'data/P2VVDataSets2011Reco12_6KKMassBins_2TagCats.root'

# fit options
fitOpts = dict(  NumCPU    = 2
               , Optimize  = 2
               , Minimizer = 'Minuit2'
               , Offset    = True
#               , Hesse     = False
               , Timer     = True
#               , Verbose   = True
              )
pdfConfig['fitOptions'] = fitOpts

fitRange      = ''
randomParVals = ( ) # ( 1., 12345 )
MinosPars     = [#  'AparPhase'
                 #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
                 #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
                ]

# PDF options
pdfConfig['numTimeResBins'] = 40
pdfConfig['timeResType']    = 'eventNoMean'
pdfConfig['timeEffType']    = 'paper2012'

pdfConfig['timeEffHistFiles'] = dict(  file      = 'data/Bs_HltPropertimeAcceptance_Data-20120816.root'
                                     , hlt1UB    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                     , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                    )

pdfConfig['anglesEffType']   = 'weights'
pdfConfig['angEffMomsFiles'] = 'data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

pdfConfig['SSTagging']   = True
pdfConfig['condTagging'] = True
pdfConfig['contEstWTag'] = True

pdfConfig['paramKKMass']     = 'simultaneous'
pdfConfig['KKMassBinBounds'] = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ]
pdfConfig['CSPValues']       = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

KKMassPars = pdfConfig['obsDict']['KKMass']
pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                  , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )

pdfConfig['lambdaCPParam'] = 'lambPhi'

pdfConfig['externalConstr']['dM']             = ( 17.63, 0.11 )
pdfConfig['externalConstr']['timeResSigmaSF'] = (  1.45, 0.06 )


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
worksp = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

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
fitData = correctWeights( dataSet, 'KKMassCat' if pdfConfig['paramKKMass'] == 'simultaneous' else '' )

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)


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

if randomParVals :
    # give parameters random offsets
    import random
    print 'Bs2JpsiKK2011Fit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
          % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
    random.seed( randomParVals[1] if randomParVals[1] else None )
    for par in pdfPars :
        if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

# print parameters
print 120 * '='
print 'Bs2JpsiKK2011Fit: fit data:'
fitData.Print()
print 'Bs2JpsiKK2011Fit: observables in PDF:'
pdfObs.Print('v')
print 'Bs2JpsiKK2011Fit: parameters in PDF:'
pdfPars.Print('v')
print 'Bs2JpsiKK2011Fit: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()

if doFit :
    # fit data
    print 120 * '='
    print 'Bs2JpsiKK2011Fit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    RooMinPars = [ ]
    if MinosPars :
        print 'Bs2JpsiKK2011Fit: running Minos for parameters',
        for parName in MinosPars :
            RooMinPars.append( pdfPars.find(parName) )
            print '"%s"' % RooMinPars[-1],
        print

    fitResult = pdf.fitTo( fitData, SumW2Error = False, Minos = RooMinPars, Save = True, Range = fitRange, **fitOpts )

    # print parameter values
    from P2VV.Imports import parNames, parValues2011 as parValues
    print 'Bs2JpsiKK2011Fit: parameters:'
    fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    fitResult.covarianceMatrix().Print()
    fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

    if parFileOut :
        # write parameters to file
        pdfConfig.getParametersFromPdf( pdf, fitData )
        pdfConfig.writeParametersToFile( filePath = parFileOut, FitStatus = ( fitResult.status(), fitResult.minNll(), fitResult.edm() ) )
