###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
from P2VV.Parameterizations.FullPDFs import SimulCatSettings
doFit      = True
HLT1Biased = False
dataPath   = '/project/bfys/jleerdam/data/Bs2Jpsiphi/'

pdfConfig['runPeriods'] = [ ] #[ 2011, 2012 ]

parFileIn     = ''
parFileOut    = ''

dataSetName = 'JpsiKK_sigSWeight'
dataSetFile = dataPath + 'Reco14/P2VVDataSets2011Reco14_I2Mass_exclBiased_narrowKKMass_noKKMassBins_2TagCats.root'

MinosPars = [#  'phiCP', 'lambdaCP'
             #, 'AparPhase', 'AperpPhase'
             #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3'#, 'ASOddPhase_bin4', 'ASOddPhase_bin5'
             #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3'#,        'f_S_bin4',        'f_S_bin5'
             #, 'timeResSigmaSF'
            ]

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
constTagCatCoefs  = True
constAvgCEvenOdd  = True
constWTagAsyms    = ''
constCSP          = True
constAmplitudes   = False
constLambdaCP     = ''  # 'lamb'

# PDF options
pdfConfig['numTimeBins'] = 30
pdfConfig['numTimeResBins']     = 40
pdfConfig['timeResType']        = 'eventNoMean' # 'eventDoubleGaussConstantFixedMean' # 'event' # 'eventNoMean'
pdfConfig['constrainTResScale'] = 'fixed' # 'fixed'  # 'constrain'
pdfConfig['timeEffType']        = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['constrainDeltaM']    = 'fixed'  # 'constrain' # fixed
pdfConfig['constrainBeta']      = 'noBeta'      # '' / 'constrain' / 'fixed' / 'noBeta'

pdfConfig['timeEffHistFiles'] = dict(  file      = dataPath + 'Reco14/Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
                                     , hlt1UB    = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                     , hlt1ExclB = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                    )
#timeEffHistFiles = SimulCatSettings('timeEffHistFiles')
#timeEffHistFiles.addSettings( [ 'runPeriod' ], [ [ 'p2011' ] ]
#                             , dict(  file      = dataPath + 'Reco14/Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
#                                    , hlt1UB    = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#                                    , hlt1ExclB = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
#                                   )
#                            )
#timeEffHistFiles.addSettings( [ 'runPeriod' ], [ [ 'p2012' ] ]
#                             , dict(  file      = dataPath + 'Reco14/Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
#                                    , hlt1UB    = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'
#                                    , hlt1ExclB = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
#                                   )
#                            )
#pdfConfig['timeEffHistFiles'] = timeEffHistFiles

pdfConfig['anglesEffType'] = 'weights'
pdfConfig['angEffMomsFiles'] = dataPath + 'Reco14/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
#angEffMomsFiles = SimulCatSettings('angEffMomsFiles')
#angEffMomsFiles.addSettings( [ 'hlt1_excl_biased_dec' ], [ [ 'notExclB' ] ]
#                            , dataPath + 'Reco14/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
#                           )
#angEffMomsFiles.addSettings( [ 'hlt1_excl_biased_dec' ], [ [ 'exclB' ] ]
#                            , dataPath + 'Reco14/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
#                           )
#pdfConfig['angEffMomsFiles'] = angEffMomsFiles

pdfConfig['sigTaggingType']   = 'tagUntag'
pdfConfig['SSTagging']        = True
pdfConfig['condTagging']      = True
pdfConfig['contEstWTag']      = True
pdfConfig['constrainTagging'] = 'fixed'  # 'constrain'

pdfConfig['amplitudeParam']  = 'phasesSWaveFrac' # 'bank' # 'phasesSWaveFrac'
pdfConfig['ASParam']         = 'deltaPerp'  # 'deltaPerp'
pdfConfig['AparParam']       = 'phase' # 'Mag2ReIm' # 'phase'
pdfConfig['ambiguityPars']   = False

pdfConfig['paramKKMass']     = ''
pdfConfig['KKMassBinBounds'] = [ 1008., 1050. ] #[ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['CSPValues']       = [ 0.498 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.959, 0.498, 0.968 ]
KKMassPars = pdfConfig['obsDict']['KKMass']
pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                  , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )

pdfConfig['lambdaCPParam'] = 'lambPhi'

pdfConfig['splitParams'] = dict(#  KKMassCat = [ 'f_S', 'ASOddPhase' ]
                                #, runPeriod = [ 'Gamma' ] # [ 'betaTimeEff' ] # [ 'tagCatCoef0_1' ]
                                hlt1_excl_biased_dec = [ 'tagCatCoef0_1' ]
                               )

pdfConfig['externalConstr'] = { }

from P2VV.Imports import extConstraintValues
extConstraintValues.setVal( 'DM',      ( 17.768, 0.024  ) )
extConstraintValues.setVal( 'P0OS',    (  0.380,  0.011, 0.380 ) )
extConstraintValues.setVal( 'P1OS',    (  1.000,  0.062 ) )
extConstraintValues.setVal( 'DelP0OS', (  0.010,  0.001 ) )
extConstraintValues.setVal( 'DelP1OS', (  0.070,  0.010 ) )
extConstraintValues.setVal( 'P0SS',    (  0.437,  0.008, 0.437 ) )
extConstraintValues.setVal( 'P1SS',    (  1.00,   0.12  ) )
extConstraintValues.setVal( 'DelP0SS', ( -0.016,  0.002 ) )
extConstraintValues.setVal( 'DelP1SS', (  0.015,  0.019 ) )

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

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(dataSet)
pdfPars = pdf.getParameters(dataSet)


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
    for par in pdfBuild['signalKKMass'].pdf().getParameters(dataSet) : par.setConstant(True)

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
    pdfBuild['amplitudes'].setConstant('C_SP')

ws['hlt1_excl_biased_dec'].setIndex( 1 if HLT1Biased else 0 )

# print parameters
print 120 * '='
print 'JvLFit: fit data:'
dataSet.Print()
print 'JvLFit: observables in PDF:'
pdfObs.Print('v')
print 'JvLFit: parameters in PDF:'
pdfPars.Print('v')
print 'JvLFit: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()

if doFit :
    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( dataSet.numEntries(), 'weighted' if dataSet.isWeighted() else 'not weighted' )

    RooMinPars = [ ]
    if MinosPars :
        print 'JvLFit: running Minos for parameters',
        for parName in MinosPars :
            RooMinPars.append( pdfPars.find(parName) )
            print '"%s"' % RooMinPars[-1],
        print

    fitResult = pdf.fitTo( dataSet, SumW2Error = False
                          , Minos = RooMinPars, Save = True, Range = fitRange
                          , **fitOpts
                         )

    from P2VV.Imports import parNames, parValues
    print 'JvLFit: parameters:'
    fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    #fitResult.covarianceMatrix().Print()
    #fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

if parFileOut :
    # write parameters to file
    pdfConfig.getParametersFromPdf( pdf, dataSet )
    pdfConfig.writeParametersToFile( filePath = parFileOut )
