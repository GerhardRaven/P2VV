###########################################################################################################################################
## set script parameters ##
###########################

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
pdfConfig['selection']  = 'HLT1Unbiased' # 'paper2012' # 'HLT1Unbiased'
generateData            = False
doFit                   = False 
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
pdfConfig['blind']      = False
pdfConfig['nominalPdf'] = False  # nominal PDF option does not work at the moment
corrSFitErr             = 'sumWeight'     # '' / 'matrix' / 'sumWeight'
randomParVals           = ( ) # ( 2., 12345 )

plotsFile = 'plots/VSSFit.ps' if pdfConfig['SFit']\
       else 'plots/VSCFit.ps'
parameterFile = None # 'VSSFit.par' if pdfConfig['SFit'] else 'VSCFit.par'

if readData :
    pdfConfig['nTupleName'] = 'DecayTree'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20120620_MagDownMagUp.root'
    pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20120821_MagDownMagUp.root'
    pdfConfig['nominalDataSet'] = False
else :
    pdfConfig['nTupleName'] = None
    pdfConfig['nTupleFile'] = None

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'VSSFit.root' if pdfConfig['SFit'] else 'VSCFit.root'

# fit options
fitOpts = dict(  NumCPU    = 2
               , Optimize  = 2
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
              )
pdfConfig['fitOptions'] = fitOpts

# plot options
lineWidth = 2
markStyle = 8
markSize  = 0.4

# PDF options
pdfConfig['transversityAngles'] = False  # default: False | nominal: True

pdfConfig['bkgAnglePdf']          = 'hybrid'  # default/nominal: ''
pdfConfig['sigTaggingPdf']        = 'tagUntag'  # default: 'tagUntag' | nominal: 'tagCats'
pdfConfig['bkgTaggingPdf']        = 'tagUntagRelative'  # default: 'tagUntagRelative' | 'tagCatsRelative'
pdfConfig['multiplyByTagPdf']     = False
pdfConfig['multiplyByTimeEff']    = ''
pdfConfig['timeEffType']          = 'HLT1Unbiased' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['multiplyByAngEff']     = 'basis012'  # default: 'basis012'
pdfConfig['parameterizeKKMass']   = ''  # default/nominal: 'simultaneous'
pdfConfig['ambiguityParameters']  = False
pdfConfig['lifetimeRange']        = ( 0.3, 14. )
pdfConfig['SWeightsType']         = ''  # default/nominal: simultaneousFreeBkg
pdfConfig['KKMassBinBounds']      = [ 1008., 1032. ] # [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ]
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.27, 0.09), (0.05, 0.02), (0.04, 0.02), (0.17, 0.05) ]
#                                     , [ (1.4,  0.5 ), (0.5,   0.3  ), (-0.5,  0.3  ), (-0.7, 0.2 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.33, 0.09), (0.073, 0.030), (0.009, 0.012), (0.012, 0.010), (0.061, 0.027), (0.18, 0.04) ]
#                                     , [ (1.1,  0.5 ), (0.7,   0.2  ), (0.4,   0.4  ), (-0.6,  0.3  ), (-0.4, 0.2   ), (-0.7, 0.2 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ ( 0.18, 0.07 ), ( 0.07, 0.03 ), ( 0.01, 0.02 ), ( 0.02, 0.01 ), ( 0.05, 0.03 ), ( 0.15, 0.04 ) ]
#                                     , [ ( 1.4,  0.5  ), ( 0.8,  0.3  ), ( 0.3,  0.4  ), ( -0.5, 0.2  ), ( -0.5, 0.2  ), ( -0.7, 0.2  ) ] )
pdfConfig['CSPValues']            = [ 0.498 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

pdfConfig['sameSideTagging']    = False  # nominal: False
pdfConfig['conditionalTagging'] = False  # nominal: True
pdfConfig['continuousEstWTag']  = False  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 50
pdfConfig['constrainTagging']   = 'constrain'  # nominal: 'constrain'
pdfConfig['timeResType']           = '' # 'event' # 'eventNoMean'

pdfConfig['numTimeResBins']        = 100
pdfConfig['constrainTimeResScale'] = 'constrain'  # nominal: 'constrain'

pdfConfig['numEvents'] = 10000
pdfConfig['signalFraction'] = 0.45
pdfConfig['massRangeBackground'] = True

pdfConfig['amplitudeParam'] = 'phasesSWaveFrac' # default: 'bank' | nominal: 'phasesSWaveFrac'
pdfConfig['ASParam']        = 'deltaPerp'  # default/nominal: 'deltaPerp'
pdfConfig['AparParam']      = 'phase' # default: 'Mag2ReIm' | nominal: 'phase'

pdfConfig['constrainDeltaM'] = 'constrain'  # nominal: 'constrain'

pdfConfig['lambdaCPParam'] = 'lambPhi'  # default/nominal: 'lambPhi'

fastFit          = False
manualTagCatBins = False
constTagCatCoefs = True  # default: True / nominal: False
constAvgCEvenOdd = True  # default: False / nominal: True
constWTagAsyms   = True  # default/nominal: True
constCSP         = True  # default/nominal: True
constAmplitudes  = False
constLambdaCP    = ''  # default/nominal: ''

A0Mag2Val     =  0.521
APerpMag2Val  =  0.251
f_SVal        =  0.027
AparPhaseVal  =  3.34
AperpPhaseVal =  3.00
ASOddPhaseVal = -0.01

lambCPSqVal = 0.8874
phiCPVal    = 0.023

if not readData or manualTagCatBins :
    pdfConfig['tagCatsOS'] = [  ( 'Untagged',  0, 0.500001 )
                              , ( 'TagCat1',   1, 0.499999 )
                              , ( 'TagCat2',   2, 0.40     )
                              , ( 'TagCat3',   3, 0.25     )
                             ]
    pdfConfig['tagCatsSS'] = [  ( 'Untagged',  0, 0.500001 )
                              , ( 'TagCat1',   1, 0.499999 )
                              , ( 'TagCat2',   2, 0.30     )
                             ]

pdfConfig['timeEffHistFile']      = '/project/bfys/jleerdam/data/Bs2Jpsiphi/timeAcceptanceStartValues.root'\
                                    if pdfConfig['timeEffType'] == 'fit' else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504_unitAcceptance.root'
pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_40bins'
pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
#pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1ExclB_40bins'
pdfConfig['angEffMomentsFile']    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'\
                                    if not pdfConfig['nominalPdf'] and pdfConfig['transversityAngles'] else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'

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

numBins = ( 60, 30, 30, 30 )
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
    print 'VSFit: observables and parameters in generation process:'
    for var in pdf.getVariables() : var.Print()
    print 120 * '='

    # generate data
    nEvents = int( pdfConfig['numEvents'] * ( pdfConfig['signalFraction'] if pdfConfig['SFit'] else 1. ) )
    print 'VSFit: generating %d events' % nEvents
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
    if corrSFitErr == 'sumWeight' :
        from P2VVGeneralUtils import correctSWeights
        fitData = correctSWeights( pdfBuild['sigSWeightData'], 'N_bkgMass_sw'
                                  , 'KKMassCat' if pdfConfig['parameterizeKKMass'] == 'simultaneous' else '' )
    else :
        fitData = pdfBuild['sigSWeightData']

else :
    defData = pdfBuild['data']
    fitData = pdfBuild['data']
    sigData = pdfBuild['sigSWeightData']
    bkgData = pdfBuild['bkgSWeightData']

# get parameters in PDF
pdfPars = pdf.getParameters(fitData)


###########################################################################################################################################
## fit data ##
##############

# float/fix values of some parameters
if 'lamb' in constLambdaCP.lower() :
    from math import sqrt
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

if pdfConfig['nominalPdf'] or constWTagAsyms :
    pdfBuild['tagCatsOS'].parameter('wTagDelP0OS').setVal(0.)
    pdfBuild['tagCatsOS'].parameter('wTagDelP1OS').setVal(0.)
    pdfBuild['tagCatsSS'].parameter('wTagDelP0SS').setVal(0.)
    pdfBuild['tagCatsSS'].parameter('wTagDelP1SS').setVal(0.)
    pdfBuild['tagCatsOS'].setConstant('wTagDelP0')
    pdfBuild['tagCatsOS'].setConstant('wTagDelP1')
    pdfBuild['tagCatsSS'].setConstant('wTagDelP0')
    pdfBuild['tagCatsSS'].setConstant('wTagDelP1')

if pdfConfig['parameterizeKKMass'] == 'functions' :
    for par in pdfBuild['signalKKMass'].pdf().getParameters(fitData) : par.setConstant(True)
    if not pdfConfig['SFit'] :
        for par in pdfBuild['backgroundKKMass'].pdf().getParameters(fitData) : par.setConstant(True)

if pdfConfig['nominalPdf'] or constCSP : pdfBuild['amplitudes'].setConstant('C_SP')

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
    pdfBuild['lambdaCP'].setConstant('phiCP')
    for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
        if not pdfConfig['sameSideTagging'] :
            CEvenOdds.setConstant('avgCEven.*|avgCOdd.*')
        else :
            for CEvenOdd in CEvenOdds : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')
    pdfBuild['tagCatsOS'].setConstant('.*')
    pdfBuild['tagCatsSS'].setConstant('.*')
    pdfBuild['lifetimeParams'].setConstant('dM')
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
    print 'VSFit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
          % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
    random.seed( randomParVals[1] if randomParVals[1] else None )
    for par in pdfPars :
        if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

# print parameters
print 120 * '='
print 'VSFit: parameters in PDF:'
pdfPars.Print('v')

if ( readData or generateData ) and doFit :
    # fit data
    print 120 * '='
    print 'VSFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    if pdfConfig['SFit'] : fitResult = pdf.fitTo(fitData, SumW2Error = True if corrSFitErr == 'matrix' else False, Save = True, **fitOpts)
    else                 : fitResult = pdf.fitTo(fitData,                                                          Save = True, **fitOpts)

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

        print '\nVSFit: reparameterization of transversity amplitudes:'
        ampsData.Print()
        ampMeasList.Print('v')
        ampMeasCovs.Print()
        ampsFitResult.Print()
        ampsFitResult.covarianceMatrix().Print()

    print 'VSFit: parameters:'
    fitData.Print()
    fitResult.Print()
    fitResult.covarianceMatrix().Print()
    fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

else :
    fitResult = None


###########################################################################################################################################
## make some plots ##
#####################

if ( readData or generateData ) and ( makeObservablePlots or pdfConfig['makePlots'] or makeKKMassPlots ) :
    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlack, kBlue, kRed, kGreen, kDashed, kFullCircle, kFullSquare

    # create projection data set for conditional observables
    if pdfConfig['SFit'] :
        comps = None
    else :
        comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                }

    projWDataSet = []
    if   pdfConfig['continuousEstWTag']   : projWDataSet += [ tagCatP2VVOS, estWTagOS, iTagOS ]
    elif pdfConfig['conditionalTagging']  : projWDataSet += [ tagCatP2VVOS, iTagOS ]
#    if   pdfConfig['eventTimeResolution'] : projWDataSet += [ timeRes ]

    if projWDataSet :
        bulkData = data.reduce( CutRange = 'Bulk' )
        projWData     = dict( ProjWData = ( data.reduce(  ArgSet = projWDataSet ), True ) )
        projWDataBulk = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataSet ), True ) )
    else :
        projWData     = dict()
        projWDataBulk = dict()

if pdfConfig['makePlots'] :
    # plot background time
    print 'VSFit: plotting background lifetime distribution'
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

if makeKKMassPlots and pdfConfig['parameterizeKKMass'] and fitResult\
        and ( ( pdfConfig['amplitudeParam'] == 'bank' and pdfConfig['ASParam'] != 'ReIm' )\
              or ( pdfConfig['amplitudeParam'] == 'phasesSWaveFrac' and pdfConfig['ASParam'] == 'deltaPerp' ) ) :
    # create S-wave phase plots
    nKKBins = pdfBuild['KKMassBinning'].numBins()

    from array import array
    KKMassVals    = array( 'd', [ pdfBuild['KKMassBinning'].binCenter(binIter)      for binIter in range(nKKBins) ] )
    KKMassLowErr  = array( 'd', [ 0.5 * pdfBuild['KKMassBinning'].binWidth(binIter) for binIter in range(nKKBins) ] )
    KKMassHighErr = KKMassLowErr

    from ROOT import TGraphAsymmErrors
    parList = fitResult.floatParsFinal()
    deltaS1Vals    = array( 'd', [ parList.find( 'ASOddPhase_bin%d' % binIter ).getVal()   for binIter in range(nKKBins) ] )
    deltaS1LowErr  = array( 'd', [ parList.find( 'ASOddPhase_bin%d' % binIter ).getError() for binIter in range(nKKBins) ] )
    deltaS1HighErr = deltaS1LowErr
    deltaSGraphs   = [ TGraphAsymmErrors( len(KKMassVals), KKMassVals,                  deltaS1Vals,
                                                           KKMassLowErr, KKMassHighErr, deltaS1LowErr, deltaS1HighErr) ]

    from math import pi
    deltaS2Vals    = array( 'd', [ pi - parList.find( 'ASOddPhase_bin%d' % binIter ).getVal() for binIter in range(nKKBins) ] )
    deltaS2LowErr  = deltaS1LowErr
    deltaS2HighErr = deltaS2LowErr
    deltaSGraphs  += [ TGraphAsymmErrors( len(KKMassVals), KKMassVals,                  deltaS2Vals,
                                                           KKMassLowErr, KKMassHighErr, deltaS2LowErr, deltaS2HighErr ) ]

    if pdfConfig['ambiguityParameters'] : deltaSGraphs = [ deltaSGraphs[1], deltaSGraphs[0] ]

    from ROOT import kBlack, kBlue
    deltaSGraphs[0].SetLineColor(kBlue)
    deltaSGraphs[1].SetLineColor(kBlack)

    deltaSGraphs[0].SetMarkerColor(kBlue)
    deltaSGraphs[1].SetMarkerColor(kBlack)

    deltaSGraphs[0].SetLineWidth(2)
    deltaSGraphs[1].SetLineWidth(2)

    from ROOT import kFullCircle, kFullSquare
    deltaSGraphs[0].SetMarkerStyle(kFullCircle)
    deltaSGraphs[1].SetMarkerStyle(kFullSquare)

    deltaSGraphs[0].SetMinimum( min( deltaSGraphs[0].GetYaxis().GetXmin(), deltaSGraphs[1].GetYaxis().GetXmin() ) )
    deltaSGraphs[0].SetMaximum( max( deltaSGraphs[0].GetYaxis().GetXmax(), deltaSGraphs[1].GetYaxis().GetXmax() ) )

    deltaSGraphs[0].GetXaxis().SetTitle('m_{KK} (MeV)')
    deltaSGraphs[0].GetYaxis().SetTitle('#delta_{S} - #delta_{#perp}    (rad)')

    deltaSGraphs[0].SetTitle('S-Wave Phases')
    deltaSGraphs[0].GetYaxis().SetTitleOffset(1.)

    from ROOT import TLegend
    leg = TLegend( 0.52, 0.41, 0.90, 0.61 )
    leg.AddEntry( deltaSGraphs[0], 'solution I  (#Delta#Gamma_{s} > 0)', 'LPE' )
    leg.AddEntry( deltaSGraphs[1], 'solution II (#Delta#Gamma_{s} < 0)', 'LPE' )
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    from ROOT import TCanvas
    deltaSCanv = TCanvas( 'deltaSCanv', 'S-Wave Phases' )
    deltaSCanv.SetLeftMargin(0.1)
    deltaSGraphs[0].Draw('AP')
    deltaSGraphs[1].Draw('P SAMES')
    leg.Draw()

else :
    deltaSCanv = None

if makeObservablePlots and not pdfBuild['iTagZeroTrick'] :
    # plot lifetime and angles
    print 'VSFit: plotting lifetime and angular distributions'
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yScale, logY )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSetP2VV[ : 5 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSetP2VV[ : 5 ] ]
                   , ( '', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                   , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                   , ( True, ) + 3 * ( False, )
                  ) :
        plot(  pad, obs, defData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                                     )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize                      )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth )
             , components = comps
            )

    # plot lifetime
    timePlotTitles = tuple( [ time.GetTitle() + title for title in (  ' - linear'
                                                                    , ' - logarithmic'
                                                                    , ' - B/#bar{B} asymmetry'
                                                                   )
                            ] )
    timeCanv = TCanvas( 'timeCanv', 'Lifetime' )
    print 'VSFit: plotting lifetime distribution'
    for ( pad, nBins, plotTitle, yTitle, yScale, dataCuts, pdfCuts, logY )\
            in zip(  timeCanv.pads( 2, 2 )
                   , 3 * [ pdfConfig['numTimeBins'] ]
                   , timePlotTitles
                   , 2 * ( None, ) + ( 'B/#bar{B} asymmetry', )
                   , ( ( None, None ), ( 50., None ), ( None, None ) )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTagOS ), )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTagOS ), )
                   , ( False, True, False )
                  ) :
        plot(  pad, time, defData, pdf, yTitle = yTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                    )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                         )
             , pdfOpts    = dict( list( projWDataBulk.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , components = comps
            )

    # plot lifetime (tagged/untagged)
    print 'VSFit: plotting lifetime distributions tagged/untagged'
    timePlotTitles1 = tuple( [ time.GetTitle() + title for title in (  ' - untagged'
                                                                     , ' - tagging category 2'
                                                                     , ' - tagging category %d' % tagCat5Min
                                                                     , ' - B/#bar{B} asymmetry untagged'
                                                                     , ' - B/#bar{B} asymmetry tagging category 2'
                                                                     , ' - B/#bar{B} asymmetry tagging category %d' % tagCat5Min
                                                                    )
                            ] )
    timeCanv1 = TCanvas( 'timeCanv1', 'Lifetime' )
    for ( pad, nBins, plotTitle, yTitle, dataCuts, pdfCuts, logY )\
        in zip(  timeCanv1.pads( 3, 2 )
             , 6 * [ pdfConfig['numTimeBins'] ]
             , timePlotTitles1
             , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
             ,   ( dict( Cut = '%s == 0'  % ( tagCatP2VVOS.GetName()             )                     ), )
               + ( dict( Cut = '%s == 2'  % ( tagCatP2VVOS.GetName()             )                     ), )
               + ( dict( Cut = '%s == %d' % ( tagCatP2VVOS.GetName(), tagCat5Min )                     ), )
               + ( dict( Cut = '%s == 0'  % ( tagCatP2VVOS.GetName()             ), Asymmetry = iTagOS ), )
               + ( dict( Cut = '%s == 2'  % ( tagCatP2VVOS.GetName()             ), Asymmetry = iTagOS ), )
               + ( dict( Cut = '%s == %d' % ( tagCatP2VVOS.GetName(), tagCat5Min ), Asymmetry = iTagOS ), )
             ,   ( dict( Slice = ( tagCatP2VVOS, 'Untagged'              )                     ), )
               + ( dict( Slice = ( tagCatP2VVOS, 'TagCat2'               )                     ), )
               + ( dict( Slice = ( tagCatP2VVOS, 'TagCat%d' % tagCat5Min )                     ), )
               + ( dict( Slice = ( tagCatP2VVOS, 'Untagged'              ), Asymmetry = iTagOS ), )
               + ( dict( Slice = ( tagCatP2VVOS, 'TagCat2'               ), Asymmetry = iTagOS ), )
               + ( dict( Slice = ( tagCatP2VVOS, 'TagCat%d' % tagCat5Min ), Asymmetry = iTagOS ), )
             , 3 * ( False, ) + 3 * ( False, )
            ) :
        plot(  pad, time, defData, pdf, yTitle = yTitle, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                     )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , components = comps
            )

    # plot angles
    print 'VSFit: plotting angular distributions'
    if plotAnglesNoEff and pdfConfig['SFit'] and pdfConfig['multiplyByTimeEff'] not in [ 'all', 'signal' ]\
            and ( pdfConfig['nominalPdf'] or not pdfConfig['conditionalTagging'] ) :
        addPDFs = [ ws['sig_t_angles_tagCat_iTag'] ]
    else :
        addPDFs = [ ]

    anglePlotTitles =   tuple(  [ angle.GetTitle()                            for angle in angles ]\
                              + [ angle.GetTitle() + ' - B/#bar{B} asymmetry' for angle in angles ] )
    anglesCanv = TCanvas( 'anglesCanv', 'Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yTitle, dataCuts, pdfCuts )\
            in zip(  anglesCanv.pads( 3, 2 )
                   , 2 * angles
                   , 2 * numAngleBins
                   , anglePlotTitles
                   , 2 * ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                   , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
                   , 3 * ( dict( ), ) + 3 * ( dict( Asymmetry = iTagOS ), )
                   , 3 * ( dict( ), ) + 3 * ( dict( Asymmetry = iTagOS ), )
                  ) :
        plot(  pad, obs, defData, pdf, addPDFs = addPDFs, xTitle = xTitle, yTitle = yTitle
             , frameOpts   = dict( Bins = nBins, Title = plotTitle                                                )
             , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts                    )
             , pdfOpts     = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , addPDFsOpts = [ dict( list( projWData.items() ), LineColor = kRed, LineWidth = lineWidth, **pdfCuts ) ]
             , components  = comps
            )

    if not pdfConfig['SFit'] and pdfConfig['SWeightsType'].startswith('simultaneous')\
            and pdfConfig['parameterizeKKMass'] == 'simultaneous' :
        # plot signal mass
        print 'VSFit: plotting mumuKK mass distribution'
        pad = pdfBuild['massCanv'].cd(2)
        plot(  pad, BMass, defData, pdf
             , frameOpts  = dict( Range = 'Signal', Bins = pdfConfig['numBMassBins'][0], Title = BMass.GetTitle() + ' full fit - signal' )
             , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                                                      )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = 2                                            )
             , components = comps
            )

    # print canvas to file
    timeAnglesCanv.Print(plotsFile + '(')
    timeCanv.Print(plotsFile)
    timeCanv1.Print(plotsFile)
    if pdfConfig['makePlots'] :
        anglesCanv.Print(plotsFile)
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

    else :
        anglesCanv.Print( plotsFile + ( '' if deltaSCanv else ')' ) )

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


#################################__Trial_area, (Try to make a plot of the CP odd CP even components)########################



#Get the normalization of the total pdf
from P2VVGeneralUtils import getCPprojectionOFpdf, getNormOverObservables
normTotal = getNormOverObservables(pdf)

#Construct the CP_Even component of the pdf
pdfEven = getCPprojectionOFpdf(pdf, "EVEN")
normEven = getNormOverObservables(pdfEven)
f_Even = normEven /  normTotal

#Construct the CP_Odd component of the pdf
pdfOdd = getCPprojectionOFpdf(pdf, "ODD")
normOdd = getNormOverObservables(pdfOdd)
f_Odd = normOdd / normTotal

#Co## nstruct the S-wave component of the pdf
pdfSwave = getCPprojectionOFpdf(pdf, "SWAVE")
normSwave = getNormOverObservables(pdfSwave)
f_Swave = normSwave / normTotal



#Plot
from ROOT import TCanvas    
from P2VVLoad import ROOTStyle
from P2VVGeneralUtils import plot

c = TCanvas()
c.Divide(2,2)

i = 0
for i_thObs in set([angles[0], angles[1], angles[2], time]):
    plot( c.cd(i+1),i_thObs, data = defData, pdf = pdf,
          addPDFs = [pdfEven,pdfOdd,pdfSwave],
          addPDFsOpts = [dict(LineStyle = 9, Normalization = f_Even),
                         dict(LineStyle = 3, Normalization = f_Odd),
                         dict(LineStyle = 5, Normalization = f_Swave)] )    
    i += 1



#Print fractions
print "\n\nCP_Even fraction = ",f_Even
print "CP_Odd fraction = ",f_Odd
print "CP_Odd fraction = ",f_Swave
