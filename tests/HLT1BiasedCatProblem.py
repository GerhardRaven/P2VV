###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
from P2VV.Parameterizations.FullPDFs import SimulCatSettings
doFit                  = True #'NLL'
HLT1Biased             = True
corrSFitErr            = 'sumWeight' # '' / 'sumWeight' / ( 0.887, [ 0.566, 0.863, 0.956, 0.948, 0.855, 0.662 ] ) / 'matrix'
pdfConfig['blind']     = {#  'phiCP'  : ( 'UnblindUniform', 'BsPhis2013EPS',  0.2  )
                          #, 'dGamma' : ( 'UnblindUniform', 'BsDGs2013EPS',   0.02 )
                         }
pdfConfig['numEvents'] = 54755
pdfConfig['sigFrac']   = 0.504
dataPath               = '/project/bfys/jleerdam/data/Bs2Jpsiphi/'

parFileIn     = '' #'2011DataFitValues.par'
parFileOut    = ''

dataSetName = 'JpsiKK_sigSWeight'
#dataSetFile = dataPath + 'P2VVDataSets2011Reco12_4KKMassBins_2TagCats.root'
#dataSetFile = dataPath + 'P2VVDataSets2011Reco12_unbiased_4KKMassBins_2TagCats.root'
dataSetFile = dataPath + 'P2VVDataSets2011Reco12_exclBiased_4KKMassBins_2TagCats.root'

MinosPars = [#  'phiCP', 'lambdaCP'
             #, 'AparPhase', 'AperpPhase'
             #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3'#, 'ASOddPhase_bin4', 'ASOddPhase_bin5'
             #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3'#,        'f_S_bin4',        'f_S_bin5'
             #, 'timeResSigmaSF'
            ]
dllPars = [ ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

# fit options
fitRange = ''
fitOpts = dict(  NumCPU    = 8
               , Optimize  = 2
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
               , Strategy  = 1
               , Offset    = True
              )
pdfConfig['fitOptions'] = fitOpts

fastFit           = False
manualTagCatBins  = False
constTagCatCoefs  = True
constAvgCEvenOdd  = True
constWTagAsyms    = 'P1'
constCSP          = True
constAmplitudes   = False
constLambdaCP     = ''  # 'lamb'
equalAbsLambdaCPs = False

# PDF options
pdfConfig['numTimeBins'] = 30
pdfConfig['numTimeResBins']     = 40
pdfConfig['timeResType']        = 'eventNoMean' # 'eventDoubleGaussConstantFixedMean' # 'event' # 'eventNoMean'
pdfConfig['constrainTResScale'] = '' # 'fixed'  # 'constrain'
pdfConfig['timeEffType']        = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['constrainDeltaM']    = ''  # 'constrain' # fixed
pdfConfig['constrainBeta']      = 'noBeta'      # '' / 'constrain' / 'fixed' / 'noBeta'

pdfConfig['timeEffHistFiles'] = dict(  file      = dataPath + 'Bs_HltPropertimeAcceptance_Data-20120816.root'
                                     , hlt1UB    ='Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                     , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                    )
#timeEffHistFiles = SimulCatSettings('timeEffHistFiles')
#timeEffHistFiles.addSettings( 'default', [ ]
#                             , dict(  file      = dataPath + 'Bs_HltPropertimeAcceptance_Data-20120816.root'
#                                    , hlt1UB    ='Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#                                    , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
#                                   )
#                            )
#timeEffHistFiles.addSettings( [ 'KKMassCat' ], [ [ 'bin0', 'bin3' ] ]
#                             , dict(  file      = dataPath + 'Bs_HltPropertimeAcceptance_Data-20120816.root'
#                                    , hlt1UB    ='Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#                                    , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
#                                   )
#                            )
#timeEffHistFiles.addSettings( [ 'KKMassCat' ], [ [ 'bin1', 'bin2' ] ]
#                             , dict(  file      = dataPath + 'Bs_HltPropertimeAcceptance_Data-20120816.root'
#                                    , hlt1UB    ='Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#                                    , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
#                                   )
#                            )
#pdfConfig['timeEffHistFiles'] = timeEffHistFiles

pdfConfig['transAngles']   = False
pdfConfig['anglesEffType'] = 'weights'  # 'weights' # 'basis012' # 'basisSig4'
#pdfConfig['angEffMomsFiles'] = dataPath + 'hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
angEffMomsFiles = SimulCatSettings('angEffMomsFiles')
angEffMomsFiles.addSettings( 'default', [ ], dataPath + 'hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights' )
angEffMomsFiles.addSettings( [ 'hlt1_excl_biased_dec' ], [ [ 'exclBiased' ] ]
                            , dataPath + 'hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights' )
angEffMomsFiles.addSettings( [ 'hlt1_excl_biased_dec' ], [ [ 'notExclBiased' ] ]
                            , dataPath + 'hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights' )
pdfConfig['angEffMomsFiles'] = angEffMomsFiles

pdfConfig['angularRanges'] = dict( ctheta = [  ( 'ctlBin0', -1.,   -0.70 )
                                             , ( 'ctlBin1', -0.70, -0.25 )
                                             , ( 'ctlBin2', -0.25, +0.25 )
                                             , ( 'ctlBin3', +0.25, +0.70 )
                                             , ( 'ctlBin4', +0.70, +1.   )
                                            ]
                                 )
if not pdfConfig['SFit'] : pdfConfig['bkgAnglePdfType'] = 'hybrid'

pdfConfig['sigTaggingType']   = 'tagUntag'
pdfConfig['SSTagging']        = True
pdfConfig['condTagging']      = True
pdfConfig['contEstWTag']      = True
pdfConfig['constrainTagging'] = ''  # 'constrain'

pdfConfig['amplitudeParam']  = 'phasesSWaveFrac' # 'bank' # 'phasesSWaveFrac'
pdfConfig['ASParam']         = 'deltaPerp'  # 'deltaPerp'
pdfConfig['AparParam']       = 'phase' # 'Mag2ReIm' # 'phase'
pdfConfig['ambiguityPars']   = False

pdfConfig['paramKKMass']     = 'simultaneous'  # 'simultaneous'
pdfConfig['KKMassBinBounds'] = [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['CSPValues']       = [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.959, 0.498, 0.968 ]
KKMassPars = pdfConfig['obsDict']['KKMass']
pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                  , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )

pdfConfig['lambdaCPParam'] = 'lambPhi' # 'lambPhi_CPVDecay_PSWaves'  # 'lambPhi'

pdfConfig['splitParams'] = dict( KKMassCat = [ 'f_S', 'ASOddPhase' ] )

pdfConfig['externalConstr']['dM'] = ( 17.63, 0.11 )

from P2VV.Imports import extConstraintValues
extConstraintValues.setVal( 'DM', ( 17.63, 0.11 ) )

dGammaVal = 0.108
dMVal     = 17.647

A0Mag2Val     =  0.5214
APerpMag2Val  =  0.2532
f_SVal        =  0.0266
AparPhaseVal  =  3.333
AperpPhaseVal =  2.998
ASOddPhaseVal =  0.0291

lambCPSqVal = 1. # 0.959**2
phiCPVal    = 0.009


###########################################################################################################################################
## build PDF ##
###############

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
pdfConfig['signalData'] = dataSet
pdfConfig['readFromWS'] = True

# build PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

# get variables
obsSetP2VV = [ pdfBuild['observables'][obs] for obs in [ 'time', 'cpsi', 'ctheta', 'phi', 'iTagOS' ] ]
time       = obsSetP2VV[0]
angles     = obsSetP2VV[ 1 : 4 ]
iTagOS     = obsSetP2VV[4]
iTagSS     = pdfBuild['observables']['iTagSS']
BMass      = pdfBuild['observables']['mass']
mumuMass   = pdfBuild['observables']['mumuMass']
KKMass     = pdfBuild['observables']['KKMass']
estWTagOS  = pdfBuild['observables']['wTagOS']
timeRes    = pdfBuild['observables']['timeRes']

if not pdfConfig['SFit'] : obsSetP2VV.append(BMass)

tagCatP2VVOS = pdfBuild['observables']['tagCatOS']
tagCatP2VVSS = pdfBuild['observables']['tagCatSS']
obsSetP2VV.append(tagCatP2VVOS)

if not 'Optimize' in fitOpts or fitOpts['Optimize'] < 2 :
    # unset cache-and-track
    for par in pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

from P2VV.Utilities.DataHandling import correctSWeights
fitData = correctSWeights( dataSet, 'N_cbkgMass_sw'
                          , 'KKMassCat' if pdfConfig['paramKKMass'] == 'simultaneous' else ''
                          , CorrectionFactors = None if corrSFitErr == 'sumWeight' else corrSFitErr )

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)


###########################################################################################################################################
## fit data ##
##############

# float/fix values of some parameters
from math import sqrt
if 'lamb' in constLambdaCP.lower() :
    pdfBuild['lambdaCP'].setConstant( 'lambdaCPSq' if pdfConfig['lambdaCPParam'] == 'lambSqPhi' else 'lambdaCP' )
    pdfBuild['lambdaCP'].parameter('lambdaCPSq').setVal(lambCPSqVal) if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
        else pdfBuild['lambdaCP'].parameter('lambdaCP').setVal( sqrt(lambCPSqVal) )
if 'phi' in constLambdaCP.lower() :
    pdfBuild['lambdaCP'].setConstant('phiCP')
    pdfBuild['lambdaCP'].parameter('phiCP').setVal(phiCPVal)
for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
    if not pdfConfig['SSTagging'] :
        CEvenOdds.setConstant('avgCEven.*')
        if constAvgCEvenOdd : CEvenOdds.setConstant( 'avgCOdd.*', True )
    else :
        for CEvenOdd in CEvenOdds :
            CEvenOdd.setConstant('avgCEven.*')
            if constAvgCEvenOdd : CEvenOdd.setConstant( 'avgCOdd.*', True )

if not constTagCatCoefs : pdfBuild['taggingParams'].setConstant( 'tagCatCoef.*', False )

if constWTagAsyms and constWTagAsyms != 'P1' :
    pdfBuild['tagCatsOS'].parameter('wTagDelP0OS').setVal(0.)
    pdfBuild['tagCatsSS'].parameter('wTagDelP0SS').setVal(0.)
    pdfBuild['tagCatsOS'].setConstant('wTagDelP0')
    pdfBuild['tagCatsSS'].setConstant('wTagDelP0')
if constWTagAsyms and constWTagAsyms != 'P0' :
    pdfBuild['tagCatsOS'].parameter('wTagDelP1OS').setVal(0.)
    pdfBuild['tagCatsSS'].parameter('wTagDelP1SS').setVal(0.)
    pdfBuild['tagCatsOS'].setConstant('wTagDelP1')
    pdfBuild['tagCatsSS'].setConstant('wTagDelP1')

if pdfConfig['paramKKMass'] == 'amplitudes' :
    for par in pdfBuild['signalKKMass'].pdf().getParameters(fitData) : par.setConstant(True)
    if not pdfConfig['SFit'] :
        for par in pdfBuild['backgroundKKMass'].pdf().getParameters(fitData) : par.setConstant(True)

if constCSP : pdfBuild['amplitudes'].setConstant('C_SP')

if fastFit or constAmplitudes :
    pdfBuild['amplitudes'].setConstant('A0Mag2')
    pdfBuild['amplitudes'].setConstant('AperpMag2')
    pdfBuild['amplitudes'].setConstant('f_S')
    pdfBuild['amplitudes'].setConstant('AparPhase')
    pdfBuild['amplitudes'].setConstant('AperpPhase')
    pdfBuild['amplitudes'].setConstant('ASOddPhase')
    pdfBuild['amplitudes'].parameter('A0Mag2').setVal(A0Mag2Val)
    pdfBuild['amplitudes'].parameter('AperpMag2').setVal(APerpMag2Val)
    pdfBuild['amplitudes'].parameter('f_S').setVal(f_SVal)
    pdfBuild['amplitudes'].parameter('AparPhase').setVal(AparPhaseVal)
    pdfBuild['amplitudes'].parameter('AperpPhase').setVal(AperpPhaseVal)
    pdfBuild['amplitudes'].parameter('ASOddPhase').setVal(ASOddPhaseVal)

if fastFit :
    pdfBuild['lambdaCP'].setConstant('lambdaCPSq') if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
        else pdfBuild['lambdaCP'].setConstant('lambdaCP')
    pdfBuild['lambdaCP'].parameter('lambdaCPSq').setVal(lambCPSqVal) if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
        else pdfBuild['lambdaCP'].parameter('lambdaCP').setVal( sqrt(lambCPSqVal) )
    pdfBuild['lambdaCP'].setConstant('phiCP')
    pdfBuild['lambdaCP'].parameter('phiCP').setVal(phiCPVal)
    for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
        if not pdfConfig['SSTagging'] :
            CEvenOdds.setConstant('avgCEven.*|avgCOdd.*')
        else :
            for CEvenOdd in CEvenOdds : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')
    pdfBuild['tagCatsOS'].setConstant('.*')
    pdfBuild['tagCatsSS'].setConstant('.*')
    pdfBuild['lifetimeParams'].setConstant('dM')
    pdfBuild['lifetimeParams'].setConstant('dGamma')
    pdfBuild['lifetimeParams'].parameter('dM').setVal(dMVal)
    pdfBuild['lifetimeParams'].parameter('dGamma').setVal(dGammaVal)
    pdfBuild['timeResModel'].setConstant('.*')
    pdfBuild['signalBMass'].setConstant('.*')
    if not pdfConfig['SFit'] :
        pdfBuild['backgroundBMass'].setConstant('.*')
        pdfBuild['backgroundTime'].setConstant('.*')
        if hasattr( pdfBuild, '_bkgTaggingPdf' ) : pdfBuild['bkgTaggingPdf'].setConstant('.*')
    pdfBuild['amplitudes'].setConstant('C_SP')

if pdfConfig['lambdaCPParam'].startswith('lambPhi_CPVDecay') :
    if equalAbsLambdaCPs : pdfBuild['lambdaCP'].setConstant('rhoCP_A.*')
    else :                 pdfBuild['lambdaCP'].setConstant('rhoCP_m')

# set value of HLT1 biased category
ws['hlt1_excl_biased_dec'].setIndex( 1 if HLT1Biased else 0 )

# print parameters
print 120 * '='
print 'HLT1BiasedCatProblem: fit data:'
fitData.Print()
print 'HLT1BiasedCatProblem: observables in PDF:'
pdfObs.Print('v')
print 'HLT1BiasedCatProblem: parameters in PDF:'
pdfPars.Print('v')
print 'HLT1BiasedCatProblem: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()

if doFit :
    # fit data
    print 120 * '='
    print 'HLT1BiasedCatProblem: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    RooMinPars = [ ]
    if MinosPars :
        print 'HLT1BiasedCatProblem: running Minos for parameters',
        for parName in MinosPars :
            RooMinPars.append( pdfPars.find(parName) )
            print '"%s"' % RooMinPars[-1],
        print

    fitResult = pdf.fitTo( fitData, SumW2Error = True if corrSFitErr == 'matrix' else False
                          , Minos = RooMinPars, Save = True, Range = fitRange
                          , **fitOpts
                         )

    from P2VV.Imports import parNames, parValues4KKBins as parValues
    print 'HLT1BiasedCatProblem: parameters:'
    fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    fitResult.covarianceMatrix().Print()
    fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

if parFileOut :
    # write parameters to file
    pdfConfig.getParametersFromPdf( pdf, fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut )
