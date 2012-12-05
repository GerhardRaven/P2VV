###########################################################################################################################################
## set script parameters ##
###########################

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
pdfConfig['dataSample'] = ''   # '' / 'Summer2011'
pdfConfig['selection']  = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
generateData            = False
doFit                   = True
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
pdfConfig['blind']      = False
pdfConfig['nominalPdf'] = False  # nominal PDF option does not work at the moment
corrSFitErr             = 'sumWeight' # [ 1., 0.700, 0.952, 0.938, 0.764 ] # '' / 'matrix' / 'sumWeight'
randomParVals           = ( ) # ( 1., 12346 ) # ( 2., 12345 )

#plotsFile = 'plots/JvLSFit_SWavePhases.ps'
plotsFile = 'plots/JvLSFit.ps' if pdfConfig['SFit']\
       else 'plots/JvLCFit.ps'
parameterFile = None # 'JvLSFit.par' if pdfConfig['SFit'] else 'JvLCFit.par'

if readData :
    pdfConfig['nTupleName'] = 'DecayTree'
    pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Pass3-version2_Bs_050711_nocut_Phi_P2VV.root'
    pdfConfig['nominalDataSet'] = False
else :
    pdfConfig['nTupleName'] = None
    pdfConfig['nTupleFile'] = None

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLSFit.root' if pdfConfig['SFit'] else 'JvLCFit.root'

MinosPars = [  ] # [ 'AparPhase', 'ASOddPhase' ] # [ 'AparPhase', 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5' ]
dllPars = [ ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

# fit options
fitOpts = dict(  NumCPU    = 6
               , Optimize  = 2
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
              )
pdfConfig['fitOptions'] = fitOpts

# plot options
from ROOT import kBlack, kBlue, kRed, kGreen
lineWidth     = 3
lineColor     = kBlue
markStyle     = 8
markSize      = 0.6
markColor     = kBlack
markLineWidth = 2

# PDF options
pdfConfig['transversityAngles'] = False  # default: False | nominal: True

pdfConfig['bkgAnglePdf']          = 'hybrid'  # default/nominal: ''
pdfConfig['sigTaggingPdf']        = 'tagUntag'  # default: 'tagUntag' | nominal: 'tagCats'
pdfConfig['bkgTaggingPdf']        = 'tagUntagRelative'  # default: 'tagUntagRelative' | 'tagCatsRelative'
pdfConfig['multiplyByTagPdf']     = False
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['timeEffType']          = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['multiplyByAngEff']     = 'basis012'  # default: 'basis012'
pdfConfig['parameterizeKKMass']   = ''  # default/nominal: 'simultaneous'
pdfConfig['ambiguityParameters']  = False
pdfConfig['lifetimeRange']        = ( 0.3, 14. )
pdfConfig['SWeightsType']         = ''  # default/nominal: 'simultaneousFreeBkg'
pdfConfig['KKMassBinBounds']      = [ 1008., 1032. ] # [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.33, 0.09), (0.073, 0.030), (0.009, 0.012), (0.012, 0.010), (0.061, 0.027), (0.18, 0.04) ]
#                                     , [ (1.1,  0.5 ), (0.7,   0.2  ), (0.4,   0.4  ), (-0.6,  0.3  ), (-0.4, 0.2   ), (-0.7, 0.2 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.28, 0.11), (0.06, 0.02), (0.04, 0.02), (0.27, 0.07) ]
#                                     , [ (2.7,  0.4 ), (0.22,   0.14), (-0.11, 0.17 ), (-0.97, 0.3 ) ] )
pdfConfig['CSPValues']            = [ 0.498 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ]

pdfConfig['sameSideTagging']    = False  # nominal: True
pdfConfig['conditionalTagging'] = True  # nominal: True
pdfConfig['continuousEstWTag']  = True  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 20
pdfConfig['constrainTagging']   = 'fixed'  # nominal: 'constrain'

pdfConfig['timeResType']           = 'eventNoMean' # 'event' # 'eventNoMean'
pdfConfig['numTimeResBins']        = 25
pdfConfig['constrainTimeResScale'] = 'constrain'  # nominal: 'constrain'

pdfConfig['numEvents'] = 10000
pdfConfig['signalFraction'] = 0.45
pdfConfig['massRangeBackground'] = True

pdfConfig['amplitudeParam'] = 'phasesSWaveFrac' # default: 'bank' | nominal: 'phasesSWaveFrac'
pdfConfig['ASParam']        = 'deltaPerp'  # default/nominal: 'deltaPerp'
pdfConfig['AparParam']      = 'phase' # default: 'Mag2ReIm' | nominal: 'phase'

pdfConfig['constrainDeltaM'] = 'fixed'  # nominal: 'constrain'

pdfConfig['lambdaCPParam'] = 'lambPhi'  # default/nominal: 'lambPhi'

fastFit          = True
manualTagCatBins = False
constTagCatCoefs = True  # default: True / nominal: False
constAvgCEvenOdd = True  # default: False / nominal: True
constWTagAsyms   = 'P1'  # default/nominal: 'P1'
constCSP         = True  # default/nominal: True
constAmplitudes  = False
constLambdaCP    = ''  # default/nominal: ''

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

if not readData or manualTagCatBins :
    pdfConfig['tagCatsOS'] = [  ( 'Untagged',  0, 0.500001 )
                              , ( 'TagCat1',   1, 0.499999 )
                             ]
    pdfConfig['tagCatsSS'] = [  ( 'Untagged',  0, 0.500001 )
                              , ( 'TagCat1',   1, 0.499999 )
                             ]
    #pdfConfig['tagCatsOS'] = [  ( 'Untagged',  0, 0.500001 )
    #                          , ( 'TagCat1',   1, 0.499999 )
    #                          , ( 'TagCat2',   2, 0.40     )
    #                          , ( 'TagCat3',   3, 0.25     )
    #                         ]
    #pdfConfig['tagCatsSS'] = [  ( 'Untagged',  0, 0.500001 )
    #                          , ( 'TagCat1',   1, 0.499999 )
    #                          , ( 'TagCat2',   2, 0.30     )
    #                         ]

pdfConfig['timeEffHistFile']      = '/project/bfys/jleerdam/data/Bs2Jpsiphi/timeAcceptanceStartValues.root'\
                                    if pdfConfig['timeEffType'] == 'fit' else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504_unitAcceptance.root'
#pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['angEffMomentsFile']    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'\
                                    if not pdfConfig['nominalPdf'] and pdfConfig['transversityAngles'] else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_newTrigger'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_PHSP_Basis'

if not pdfConfig['nominalPdf'] and pdfConfig['transversityAngles'] :
    pdfConfig['angleNames'] = (  ( 'trcospsi',   'cos(#psi_{tr})'   )
                               , ( 'trcostheta', 'cos(#theta_{tr})' )
                               , ( 'trphi',      '#phi_{tr}'        )
                              )
else :
    pdfConfig['angleNames'] = (  ( 'helcosthetaK', 'cos(#theta_{K})' )
                               , ( 'helcosthetaL', 'cos(#theta_{l})' )
                               , ( 'helphi',       '#phi'            )
                              )
angleNames = pdfConfig['angleNames']

numBins = ( 50, 20, 20, 20 )
pdfConfig['numTimeBins'] = 30
numAngleBins = ( 20, 20, 20 )
pdfConfig['numAngleBins'] = ( 5, 7, 9 )


###########################################################################################################################################
## build PDF ##
###############

from P2VVLoad import RooFitOutput

# workspace
from RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

# get variables
obsSetP2VV = [ pdfBuild['observables'][obs] for obs in [ 'time', 'cpsi', 'ctheta', 'phi', 'iTagOS' ] ]
time       = obsSetP2VV[0]
angles     = obsSetP2VV[ 1 : 4 ]
iTagOS     = obsSetP2VV[4]
iTagSS     = pdfBuild['observables']['iTagSS']
BMass      = pdfBuild['observables']['BMass']
mumuMass   = pdfBuild['observables']['mumuMass']
KKMass     = pdfBuild['observables']['KKMass']
estWTagOS  = pdfBuild['observables']['estWTagOS']
timeRes    = pdfBuild['observables']['timeRes']

if not pdfConfig['SFit'] : obsSetP2VV.append(BMass)

if not pdfBuild['iTagZeroTrick'] :
    tagCatP2VVOS = pdfBuild['observables']['tagCatP2VVOS']
    tagCatP2VVSS = pdfBuild['observables']['tagCatP2VVSS']
    obsSetP2VV.append(tagCatP2VVOS)

    # tagging parameters
    numTagCats    = pdfBuild['tagCatsOS']['numTagCats']
    tagCat5Min    = pdfBuild['tagCatsOS'].traditionalCatRange(5)[0]
    taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,          numTagCats ) ] )
    tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( tagCat5Min, numTagCats ) ] )

    # tagging category ranges
    tagCatP2VVOS.setRange( 'UntaggedRange', 'Untagged'    )
    tagCatP2VVOS.setRange( 'TaggedRange',   taggedCatsStr )
    tagCatP2VVOS.setRange( 'TagCat5Range',  tagCat5Str    )

if not 'Optimize' in fitOpts or fitOpts['Optimize'] < 2 :
    # unset cache-and-track
    for par in pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )


###########################################################################################################################################
## generate data ##
###################

if generateData :
    if parameterFile :
        # read parameters from file
        pdfConfig.readParametersFromFile( filePath = parameterFile )
        pdfConfig.setParametersInPdf(pdf)

    # print parameter values
    print 120 * '='
    print 'JvLFit: observables and parameters in generation process:'
    for var in pdf.getVariables() : var.Print()
    print 120 * '='

    # generate data
    nEvents = int( pdfConfig['numEvents'] * ( pdfConfig['signalFraction'] if pdfConfig['SFit'] else 1. ) )
    print 'JvLFit: generating %d events' % nEvents
    fitData = pdf.generate( obsSetP2VV, nEvents )

    # additional observables
    if pdfConfig['nominalPdf'] or not pdfConfig['transversityAngles'] :
        from P2VVGeneralUtils import addTransversityAngles
        addTransversityAngles( fitData, 'trcospsi',          'trcostheta',        'trphi'
                                      , angles[0].GetName(), angles[1].GetName(), angles[2].GetName() )

    # write data to file
    from P2VVGeneralUtils import writeData
    writeData( dataSetFile, dataSetName, fitData )

elif pdfConfig['SFit'] :
    defData = pdfBuild['sigSWeightData']
    sigData = pdfBuild['sigSWeightData']
    bkgData = pdfBuild['bkgSWeightData']
    if corrSFitErr == 'sumWeight'\
            or ( type(corrSFitErr) != str and hasattr( corrSFitErr, '__iter__' ) and hasattr( corrSFitErr, '__getitem__' ) ) :
        from P2VVGeneralUtils import correctSWeights
        fitData = correctSWeights( pdfBuild['sigSWeightData'], 'N_bkgMass_sw'
                                  , 'KKMassCat' if pdfConfig['parameterizeKKMass'] == 'simultaneous' else ''
                                  , CorrectionFactors = None if corrSFitErr == 'sumWeight' else corrSFitErr )

    else :
        fitData = pdfBuild['sigSWeightData']

else :
    defData = pdfBuild['data']
    fitData = pdfBuild['data']
    sigData = pdfBuild['sigSWeightData']
    bkgData = pdfBuild['bkgSWeightData']

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)


###########################################################################################################################################
## fit data ##
##############

# float/fix values of some parameters
from math import sqrt
if 'lamb' in constLambdaCP.lower() :
    pdfBuild['lambdaCP'].setConstant('lambdaCPSq') if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
        else pdfBuild['lambdaCP'].setConstant('lambdaCP')
    pdfBuild['lambdaCP'].parameter('lambdaCPSq').setVal(lambCPSqVal) if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
        else pdfBuild['lambdaCP'].parameter('lambdaCP').setVal( sqrt(lambCPSqVal) )
if 'phi' in constLambdaCP.lower() :
    pdfBuild['lambdaCP'].setConstant('phiCP')
    pdfBuild['lambdaCP'].parameter('phiCP').setVal(phiCPVal)
for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
    if not pdfConfig['sameSideTagging'] :
        CEvenOdds.setConstant('avgCEven.*')
        if pdfConfig['nominalPdf'] or constAvgCEvenOdd : CEvenOdds.setConstant( 'avgCOdd.*', True )
    else :
        for CEvenOdd in CEvenOdds :
            CEvenOdd.setConstant('avgCEven.*')
            if pdfConfig['nominalPdf'] or constAvgCEvenOdd : CEvenOdd.setConstant( 'avgCOdd.*', True )

if pdfConfig['nominalPdf'] or not constTagCatCoefs : pdfBuild['taggingParams'].setConstant( 'tagCatCoef.*', False )

if not pdfConfig['nominalPdf'] and constWTagAsyms and constWTagAsyms != 'P1' :
    pdfBuild['tagCatsOS'].parameter('wTagDelP0OS').setVal(0.)
    pdfBuild['tagCatsSS'].parameter('wTagDelP0SS').setVal(0.)
    pdfBuild['tagCatsOS'].setConstant('wTagDelP0')
    pdfBuild['tagCatsSS'].setConstant('wTagDelP0')
if pdfConfig['nominalPdf'] or ( constWTagAsyms and constWTagAsyms != 'P0' ) :
    pdfBuild['tagCatsOS'].parameter('wTagDelP1OS').setVal(0.)
    pdfBuild['tagCatsSS'].parameter('wTagDelP1SS').setVal(0.)
    pdfBuild['tagCatsOS'].setConstant('wTagDelP1')
    pdfBuild['tagCatsSS'].setConstant('wTagDelP1')

if pdfConfig['parameterizeKKMass'] == 'functions' :
    for par in pdfBuild['signalKKMass'].pdf().getParameters(fitData) : par.setConstant(True)
    if not pdfConfig['SFit'] :
        for par in pdfBuild['backgroundKKMass'].pdf().getParameters(fitData) : par.setConstant(True)

if pdfConfig['nominalPdf'] or constCSP : pdfBuild['amplitudes'].setConstant('C_SP')

#if fastFit or constAmplitudes :
#    pdfBuild['amplitudes'].setConstant('A0Mag2')
#    pdfBuild['amplitudes'].setConstant('AperpMag2')
#    pdfBuild['amplitudes'].setConstant('f_S')
#    pdfBuild['amplitudes'].setConstant('AparPhase')
#    pdfBuild['amplitudes'].setConstant('AperpPhase')
#    pdfBuild['amplitudes'].setConstant('ASOddPhase')
#    pdfBuild['amplitudes'].parameter('A0Mag2').setVal(A0Mag2Val)
#    pdfBuild['amplitudes'].parameter('AperpMag2').setVal(APerpMag2Val)
#    pdfBuild['amplitudes'].parameter('f_S').setVal(f_SVal)
#    pdfBuild['amplitudes'].parameter('AparPhase').setVal(AparPhaseVal)
#    pdfBuild['amplitudes'].parameter('AperpPhase').setVal(AperpPhaseVal)
#    pdfBuild['amplitudes'].parameter('ASOddPhase').setVal(ASOddPhaseVal)

if fastFit :
    pdfBuild['lambdaCP'].setConstant('lambdaCPSq') if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
        else pdfBuild['lambdaCP'].setConstant('lambdaCP')
    pdfBuild['lambdaCP'].parameter('lambdaCPSq').setVal(lambCPSqVal) if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
        else pdfBuild['lambdaCP'].parameter('lambdaCP').setVal( sqrt(lambCPSqVal) )
    pdfBuild['lambdaCP'].setConstant('phiCP')
    pdfBuild['lambdaCP'].parameter('phiCP').setVal(phiCPVal)
    for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
        if not pdfConfig['sameSideTagging'] :
            CEvenOdds.setConstant('avgCEven.*|avgCOdd.*')
        else :
            for CEvenOdd in CEvenOdds : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')
    pdfBuild['tagCatsOS'].setConstant('.*')
    pdfBuild['tagCatsSS'].setConstant('.*')
    pdfBuild['lifetimeParams'].setConstant('dM')
    #pdfBuild['lifetimeParams'].setConstant('dGamma')
    pdfBuild['lifetimeParams'].parameter('dM').setVal(dMVal)
    #pdfBuild['lifetimeParams'].parameter('dGamma').setVal(dGammaVal)
    pdfBuild['timeResModel'].setConstant('.*')
    pdfBuild['signalBMass'].setConstant('.*')
    if not pdfConfig['SFit'] :
        pdfBuild['backgroundBMass'].setConstant('.*')
        pdfBuild['backgroundTime'].setConstant('.*')
        if hasattr( pdfBuild, '_bkgTaggingPdf' ) : pdfBuild['bkgTaggingPdf'].setConstant('.*')
    pdfBuild['amplitudes'].setConstant('C_SP')

if randomParVals :
    import random
    # give parameters random offsets
    print 'JvLFit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
          % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
    random.seed( randomParVals[1] if randomParVals[1] else None )
    for par in pdfPars :
        if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

# print parameters
print 120 * '='
print 'JvLFit: fit data:'
fitData.Print()
print 'JvLFit: observables in PDF:'
pdfObs.Print('v')
print 'JvLFit: parameters in PDF:'
pdfPars.Print('v')

if ( readData or generateData ) and doFit :
    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    RooMinPars = [ ]
    if MinosPars :
        print 'JvLFit: running Minos for parameters',
        for parName in MinosPars :
            RooMinPars.append( pdfPars.find(parName) )
            print '"%s"' % RooMinPars[-1],
        print

    if pdfConfig['SFit'] :
        fitResult = pdf.fitTo(fitData, SumW2Error = True if corrSFitErr == 'matrix' else False, Minos = RooMinPars, Save = True, **fitOpts)
    else :
        fitResult = pdf.fitTo(fitData,                                                          Minos = RooMinPars, Save = True, **fitOpts)

    # reparameterize amplitudes
    if not pdfConfig['nominalPdf'] and pdfConfig['amplitudeParam'] == 'bank' and pdfConfig['ASParam'] != 'ReIm' \
            and pdfConfig['AparParam'] == 'Mag2ReIm' :
        from ROOT import RooArgSet
        parList = fitResult.floatParsFinal()
        parSet  = RooArgSet(parList)
        parCovs = fitResult.covarianceMatrix()

        from math import pi
        from ROOT import RooRealVar, RooArgList
        from P2VVParameterizations.DecayAmplitudes import A02, Aperp2, Apar2, A0Ph, AperpPh, AparPh, f_S, AS2, ASPh
        deltaPar  = AparPh  - A0Ph
        deltaPerp = AperpPh - A0Ph
        deltaS    = ASPh    - A0Ph
        if pdfConfig['ambiguityParameters'] :
            deltaPar  = -deltaPar
            deltaPerp = pi - deltaPerp
            deltaS    = -deltaS

        # physics parameters to determine
        ampPhys = [  RooRealVar( 'A0Mag2_phys',     'A0Mag2_phys',     A02,        0.,      1.      )  # 0
                   , RooRealVar( 'AparPhase_phys',  'AparPhase_phys',  deltaPar,  -2. * pi, 2. * pi )  # 1
                   , RooRealVar( 'AperpMag2_phys',  'AperpMag2_phys',  Aperp2,     0.,      1.      )  # 2
                   , RooRealVar( 'AperpPhase_phys', 'AperpPhase_phys', deltaPerp, -2. * pi, 2. * pi )  # 3
                  ]

        if pdfConfig['parameterizeKKMass'] :
            numKKMassBins = pdfBuild['KKMassBinning'].numBins() if pdfConfig['parameterizeKKMass'] == 'functions'\
                            else pdfBuild['KKMassCat'].numTypes()
            for bin in range( numKKMassBins ) :
                ampPhys += [  RooRealVar( 'f_S_phys_bin%d' % bin,     'f_S_phys_bin%d' % bin,     f_S,     0.,      1.      )  # 4 + 2*bin
                            , RooRealVar( 'ASPhase_phys_bin%d' % bin, 'ASPhase_phys_bin%d' % bin, deltaS, -2. * pi, 2. * pi )  # 5 + 2*bin
                           ]
        else :
            ampPhys += [  RooRealVar( 'f_S_phys',     'f_S_phys',     f_S,     0.,      1.      )  # 4
                        , RooRealVar( 'ASPhase_phys', 'ASPhase_phys', deltaS, -2. * pi, 2. * pi )  # 5
                       ]

        ampPhysList = RooArgList()
        for amp in ampPhys : ampPhysList.add(amp)

        # names of parameters in likelihood fit
        ampMeasNames = [  'AparMag2'    # 0
                        , 'ReApar'      # 1
                        , 'ImApar'      # 2
                        , 'AperpMag2'   # 3
                        , 'AperpPhase'  # 4
                       ]

        if pdfConfig['parameterizeKKMass'] :
            for bin in range( numKKMassBins ) :
                ampMeasNames += [  'ASOddMag2_bin%d' % bin   # 5 + 2 * bin
                                 , 'ASOddPhase_bin%d' % bin  # 6 + 2 * bin
                                ]
        else :
            ampMeasNames += [  'ASOddMag2'   # 5
                             , 'ASOddPhase'  # 6
                            ]

        # fitted parameters in terms of physics parameters
        ampMeasFuncs = {  ampMeasNames[0] : '(1.-@0-@2)/@0'
                        , ampMeasNames[1] : 'sqrt((1.-@0-@2)/@0)*cos(@1)'
                        , ampMeasNames[2] : 'sqrt((1.-@0-@2)/@0)*sin(@1)'
                        , ampMeasNames[3] : '@2/@0'
                        , ampMeasNames[4] : '@3'
                       }

        if pdfConfig['parameterizeKKMass'] :
            for bin in range( numKKMassBins ) :
                ampMeasFuncs[ ampMeasNames[ 5 + 2 * bin ] ] = '@{0:d}/(1.-@{0:d})/@2'.format( 4 + 2 * bin )
                ampMeasFuncs[ ampMeasNames[ 6 + 2 * bin ] ] = '@{0:d}-@3{1:s}'.format( 5 + 2 * bin, '+TMath::TwoPi()'\
                                                                                       if pdfConfig['ambiguityParameters'] else '' )
        else :
            ampMeasFuncs[ ampMeasNames[5] ] = '@4/(1.-@4)/@2'
            ampMeasFuncs[ ampMeasNames[6] ] = '@5-@3'

        # create fitted parameters and functions
        from ROOT import RooFormulaVar
        ampMeasList = RooArgList()
        ampFuncsList = RooArgList()
        ampMeasInds = { }
        ampMeas = [ ]
        ampFuncs = [ ]
        for ampName in ampMeasNames :
            ampMeas.append( RooRealVar( ampName + '_meas', ampName + '_meas', parSet.getRealValue(ampName), -1e+30, 1e+30 ) )
            ampFuncs.append( RooFormulaVar( ampName + '_func', ampMeasFuncs[ampName], ampPhysList ) )
            ampMeasList.add( ampMeas[-1] )
            ampFuncsList.add( ampFuncs[-1] )
            ampMeasInds[ampName] = parList.index(ampName)
        ampMeasSet = RooArgSet(ampMeasList)

        # get covariances from likelihood fit
        from ROOT import TMatrixDSym
        ampMeasCovs = TMatrixDSym(len(ampMeasNames))
        for ampIter, ampName in enumerate(ampMeasNames) :
            for ampIter1, ampName1 in enumerate(ampMeasNames) :
                ampMeasCovs[ampIter][ampIter1] = parCovs[ampMeasInds[ampName]][ampMeasInds[ampName1]]

        # determine values of physics parameters
        from ROOT import RooDataSet, RooMultiVarGaussian
        ampsData = RooDataSet( 'ampsData', 'Measured transversity amplitudes', ampMeasSet )
        ampsGauss = RooMultiVarGaussian( 'ampsGauss', 'Gaussian for measured transversity amplitudes'
                                        , ampMeasList, ampFuncsList, ampMeasCovs )
        ampsData.add(ampMeasSet)
        ampsFitResult = ampsGauss.fitTo( ampsData, Save = True, **fitOpts )

        print '\nJvLFit: reparameterization of transversity amplitudes:'
        ampsData.Print()
        ampMeasList.Print('v')
        ampMeasCovs.Print()
        ampsFitResult.Print()
        ampsFitResult.covarianceMatrix().Print()

    print 'JvLFit: parameters:'
    fitResult.Print()
    for par in RooMinPars : par.Print()
    fitResult.covarianceMatrix().Print()
    fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

else :
    fitResult = None


###########################################################################################################################################
## make some plots ##
#####################

if ( readData or generateData ) and ( makeObservablePlots or pdfConfig['makePlots'] or makeKKMassPlots or dllPars ) :
    # import plotting tools
    #from P2VVLoad import ROOTStyle
    from P2VVLoad import LHCbStyle
    from P2VVGeneralUtils import plot, _P2VVPlotStash
    from ROOT import gStyle, TCanvas, kBlack, kBlue, kRed, kGreen, kDashed, kFullCircle, kFullSquare

    # create projection data set for conditional observables
    if pdfConfig['SFit'] :
        comps = None
    else :
        comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                }

    from RooFitWrappers import SimultaneousPdf
    projWDataSet     = [ obs for obs in pdf.ConditionalObservables() ]
    projWDataBSet    = [ obs for obs in pdf.ConditionalObservables() if obs.GetName() != iTagOS.GetName() ]
    projWDataBbarSet = [ obs for obs in pdf.ConditionalObservables() if obs.GetName() != iTagOS.GetName() ]
    if isinstance( pdf, SimultaneousPdf ) : projWDataSet.append( pdf.indexCat() )
    if projWDataSet :
        bulkData = defData.reduce( CutRange = 'Bulk' )
        projWData         = dict( ProjWData = ( defData.reduce(  ArgSet = projWDataSet     ), False ) )
        projWDataBulk     = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataSet     ), False ) )
        projWDataBBulk    = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataBSet    ), False ) )
        projWDataBbarBulk = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataBbarSet ), False ) )
    else :
        projWData         = dict()
        projWDataBulk     = dict()
        projWDataBBulk    = dict()
        projWDataBbarBulk = dict()

if pdfConfig['makePlots'] :
    # plot background time
    print 'JvLFit: plotting background lifetime distribution'
    bkgTimeCanv = TCanvas( 'bkgTimeCanv', 'Background Lifetime' )
    for ( pad, data, plotTitle, logY )\
          in zip(  bkgTimeCanv.pads( 2, 2 )
                 , 2 * [ pdfBuild['bkgRangeData'], pdfBuild['bkgSWeightData'] ]
                 , [  time.GetTitle() + ' - mass side bands - linear'
                    , time.GetTitle() + ' - mass S-weights - linear'
                    , time.GetTitle() + ' - mass side bands - logarithmic'
                    , time.GetTitle() + ' - mass S-weights - logarithmic'
                   ]
                 , 2 * [ False ] + 2 * [ True ]
                ) :
        plot(  pad, time, data, pdfBuild['bkg_t'], logy = logY
             , frameOpts  = dict( Title = plotTitle, Range = 'Bulk', Bins = pdfConfig['numTimeBins'] )
             , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2  )
            )

if makeKKMassPlots and pdfConfig['parameterizeKKMass']\
        and ( ( pdfConfig['amplitudeParam'] == 'bank' and pdfConfig['ASParam'] != 'ReIm' )\
              or ( pdfConfig['amplitudeParam'] == 'phasesSWaveFrac' and pdfConfig['ASParam'] == 'deltaPerp' ) ) :
    # create S-wave phase plots
    gStyle.SetEndErrorSize(4)
    nKKBins = pdfBuild['KKMassBinning'].numBins()

    from array import array
    KKMass1Vals    = array( 'd', [ pdfBuild['KKMassBinning'].binCenter(binIter)      + 0.35 for binIter in range(nKKBins) ] )
    KKMass1LowErr  = array( 'd', [ 0.5 * pdfBuild['KKMassBinning'].binWidth(binIter) + 0.35 for binIter in range(nKKBins) ] )
    KKMass1HighErr = array( 'd', [ 0.5 * pdfBuild['KKMassBinning'].binWidth(binIter) - 0.35 for binIter in range(nKKBins) ] )
    KKMass2Vals    = array( 'd', [ pdfBuild['KKMassBinning'].binCenter(binIter)      - 0.35 for binIter in range(nKKBins) ] )
    KKMass2LowErr  = array( 'd', [ 0.5 * pdfBuild['KKMassBinning'].binWidth(binIter) - 0.35 for binIter in range(nKKBins) ] )
    KKMass2HighErr = array( 'd', [ 0.5 * pdfBuild['KKMassBinning'].binWidth(binIter) + 0.35 for binIter in range(nKKBins) ] )

    from ROOT import TGraphAsymmErrors
    deltaS1Vals    = array( 'd', [  pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getVal()     for binIter in range(nKKBins) ] )
    deltaS1LowErr  = array( 'd', [ -pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getErrorLo() for binIter in range(nKKBins) ] )
    deltaS1HighErr = array( 'd', [  pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getErrorHi() for binIter in range(nKKBins) ] )
    deltaSGraphs   = [ TGraphAsymmErrors( len(KKMass1Vals), KKMass1Vals,                   deltaS1Vals,
                                                            KKMass1LowErr, KKMass1HighErr, deltaS1LowErr, deltaS1HighErr) ]

    from math import pi
    deltaS2Vals    = array( 'd', [ pi - pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getVal() for binIter in range(nKKBins) ] )
    deltaS2LowErr  = deltaS1HighErr
    deltaS2HighErr = deltaS1LowErr
    deltaSGraphs  += [ TGraphAsymmErrors( len(KKMass2Vals), KKMass2Vals,                   deltaS2Vals,
                                                            KKMass2LowErr, KKMass2HighErr, deltaS2LowErr, deltaS2HighErr ) ]

    if pdfConfig['ambiguityParameters'] : deltaSGraphs = [ deltaSGraphs[1], deltaSGraphs[0] ]

    deltaSGraphs[0].SetLineColor(kBlue)
    deltaSGraphs[1].SetLineColor(kBlack)

    deltaSGraphs[0].SetMarkerColor(kBlue)
    deltaSGraphs[1].SetMarkerColor(kBlack)

    deltaSGraphs[0].SetLineWidth(4)
    deltaSGraphs[1].SetLineWidth(4)

    from ROOT import kFullCircle, kFullSquare
    deltaSGraphs[0].SetMarkerStyle(kFullCircle)
    deltaSGraphs[1].SetMarkerStyle(kFullSquare)
    deltaSGraphs[0].SetMarkerSize(1.3)
    deltaSGraphs[1].SetMarkerSize(1.3)

    deltaSGraphs[0].SetMinimum( min( deltaSGraphs[0].GetYaxis().GetXmin(), deltaSGraphs[1].GetYaxis().GetXmin() )       )
    deltaSGraphs[0].SetMaximum( max( deltaSGraphs[0].GetYaxis().GetXmax(), deltaSGraphs[1].GetYaxis().GetXmax() ) + 0.3 )

    deltaSGraphs[0].GetXaxis().SetTitle('m_{KK} (MeV)')
    deltaSGraphs[0].GetYaxis().SetTitle('#delta_{S} - #delta_{#perp}    (rad)')

    deltaSGraphs[0].SetTitle('')
    deltaSGraphs[0].GetXaxis().SetTitleOffset(1.0)
    deltaSGraphs[0].GetYaxis().SetTitleOffset(0.7)

    from ROOT import TLegend
    leg = TLegend( 0.59, 0.45, 0.93, 0.63 )
    leg.AddEntry( deltaSGraphs[0], 'solution I  (#Delta#Gamma_{s} > 0)', 'LPE' )
    leg.AddEntry( deltaSGraphs[1], 'solution II (#Delta#Gamma_{s} < 0)', 'LPE' )
    leg.SetBorderSize(1)
    leg.SetFillStyle(0)

    from ROOT import TPaveText
    LHCbText = TPaveText( 0.13, 0.78, 0.51, 0.95, 'NDC' )
    LHCbText.AddText(' LHCb')
    LHCbText.AddText( '#sqrt{s} = 7 TeV, L = 0.37 fb^{-1}' if pdfConfig['dataSample'] == 'Summer2011'\
                      else '#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}' )
    LHCbText.SetShadowColor(0)
    LHCbText.SetFillStyle(0)
    LHCbText.SetBorderSize(0)
    LHCbText.SetTextAlign(12)

    from ROOT import TCanvas
    deltaSCanv = TCanvas( 'deltaSCanv', 'S-Wave Phases' )
    deltaSCanv.SetLeftMargin(0.12)
    deltaSCanv.SetRightMargin(0.04)
    deltaSCanv.SetTopMargin(0.04)
    deltaSCanv.SetBottomMargin(0.15)
    deltaSGraphs[0].Draw('AP')
    deltaSGraphs[1].Draw('P SAMES')
    leg.Draw()
    LHCbText.Draw()

else :
    deltaSCanv = None

if makeObservablePlots and not pdfBuild['iTagZeroTrick'] :
    from ROOT import RooRealVar, RooConstVar, RooCategory, RooArgSet, RooCustomizer, TPaveText

    LHCbText = TPaveText( 0.45, 0.70, 0.90, 0.90, 'NDC' )
    LHCbText.AddText('LHCb unofficial')
    LHCbText.AddText( '#sqrt{s} = 7 TeV, L = 0.37 fb^{-1}' if pdfConfig['dataSample'] == 'Summer2011'\
                      else '#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}' )
    LHCbText.SetShadowColor(0)
    LHCbText.SetFillStyle(0)
    LHCbText.SetBorderSize(0)
    LHCbText.SetTextAlign(12)
    _P2VVPlotStash.append(LHCbText)

    dummyCanv = TCanvas()
    LHCbText.Draw()
    dummyCanv.Update()

    # plot lifetime and angles
    print 'JvLFit: plotting lifetime and angular distributions'
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yScale, logY, coords )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSetP2VV[ : 4 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSetP2VV[ : 4 ] ]
                   , ( '', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                   , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                   , ( True, ) + 3 * ( False, )
                   , [ ( 0.48, 0.73, 0.93, 0.93 ), ( 0.32, 0.73, 0.77, 0.93 ), ( 0.32, 0.50, 0.77, 0.70 ), ( 0.32, 0.50, 0.77, 0.70 ) ]
                  ) :
        plot(  pad, obs, defData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                                     )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize                      )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth )
             , components = comps
            )
        _P2VVPlotStash.append( LHCbText.Clone() )
        _P2VVPlotStash[-1].SetX1NDC( coords[0] )
        _P2VVPlotStash[-1].SetY1NDC( coords[1] )
        _P2VVPlotStash[-1].SetX2NDC( coords[2] )
        _P2VVPlotStash[-1].SetY2NDC( coords[3] )
        _P2VVPlotStash[-1].Draw()

    # print canvas to file
    timeAnglesCanv.Print( plotsFile + ( '(' if deltaSCanv or pdfConfig['makePlots'] else '' ) )
    #timeCanv.Print(plotsFile)
    #timeCanv1.Print(plotsFile)
    if pdfConfig['makePlots'] :
        #anglesCanv.Print(plotsFile)
        if pdfConfig['SWeightsType'].startswith('simultaneous') and pdfConfig['parameterizeKKMass'] == 'simultaneous' :
            pdfBuild['massCanvSig'].Print(plotsFile)
            pdfBuild['massCanvLeft'].Print(plotsFile)
            pdfBuild['massCanvRight'].Print(plotsFile)
        else :
            pdfBuild['massCanv'].Print(plotsFile)
        pdfBuild['mumuMassCanv'].Print(plotsFile)
        pdfBuild['KKMassCanv'].Print(plotsFile)
        bkgTimeCanv.Print(plotsFile)
        pdfBuild['bkgAnglesSWeightCanv'].Print(plotsFile)
        pdfBuild['bkgAnglesSideBandCanv'].Print(plotsFile)
        pdfBuild['estWTagCanvOS'].Print(plotsFile)
        pdfBuild['estWTagCanvSS'].Print( plotsFile + ( '' if deltaSCanv else ')' ) )

    #else :
        #anglesCanv.Print( plotsFile + ( '' if deltaSCanv else ')' ) )

elif pdfConfig['makePlots'] :
    if pdfConfig['SWeightsType'].startswith('simultaneous') and pdfConfig['parameterizeKKMass'] == 'simultaneous' :
        pdfBuild['massCanvSig'].Print(plotsFile + '(')
        pdfBuild['massCanvLeft'].Print(plotsFile)
        pdfBuild['massCanvRight'].Print(plotsFile)
    else :
        pdfBuild['massCanv'].Print(plotsFile + '(')
    pdfBuild['mumuMassCanv'].Print(plotsFile)
    pdfBuild['KKMassCanv'].Print(plotsFile)
    bkgTimeCanv.Print(plotsFile)
    pdfBuild['bkgAnglesSWeightCanv'].Print(plotsFile)
    pdfBuild['bkgAnglesSideBandCanv'].Print(plotsFile)
    pdfBuild['estWTagCanvOS'].Print(plotsFile)
    pdfBuild['estWTagCanvSS'].Print(plotsFile + ( '' if deltaSCanv else ')' ) )

if deltaSCanv :
    deltaSCanv.Print( plotsFile + ( ')' if makeObservablePlots or pdfConfig['makePlots'] else '' ) )
