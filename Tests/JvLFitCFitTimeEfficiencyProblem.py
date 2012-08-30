###########################################################################################################################################
## set script parameters ##
###########################

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
generateData            = False
doFit                   = True
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = False
pdfConfig['blind']      = False
pdfConfig['nominalPdf'] = False  # nominal PDF option does not work at the moment
corrSFitErr             = ''     # '' / 'matrix' / 'sumWeight'

plotsFile = 'plots/JvLSFit.ps' if pdfConfig['SFit']\
       else 'plots/JvLCFit.ps'
parameterFile = 'JvLSFit.par' if pdfConfig['SFit'] else 'JvLCFit.par'

if readData :
    pdfConfig['nTupleName'] = 'DecayTree'
    pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20120620_MagDownMagUp.root'
else :
    pdfConfig['nTupleName'] = None
    pdfConfig['nTupleFile'] = None

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLSFit.root' if pdfConfig['SFit'] else 'JvLCFit.root'

dllPars = [ ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

# fit options
fitOpts = dict(  NumCPU    = 8
               , Optimize  = 1
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
#               , Minimizer = 'Minuit2'
              )
pdfConfig['fitOptions'] = fitOpts

# plot options
lineWidth = 2
markStyle = 8
markSize  = 0.4

# PDF options
pdfConfig['transversityAngles'] = False  # default: False | nominal: True

pdfConfig['bkgAnglePdf']          = ''  # default/nominal: ''
pdfConfig['sigTaggingPdf']        = 'tagUntag'  # default: 'tagUntag' | nominal: 'tagCats'
pdfConfig['bkgTaggingPdf']        = 'tagUntagRelative'  # default: 'tagUntagRelative' | 'tagCatsRelative'
pdfConfig['multiplyByTagPdf']     = False
pdfConfig['multiplyByTimeEff']    = 'background'
pdfConfig['timeEffType']          = 'Moriond'
pdfConfig['multiplyByAngEff']     = ''  # default: 'basis012'
pdfConfig['parameterizeKKMass']   = ''  # default/nominal: ''
pdfConfig['ambiguityParameters']  = False
pdfConfig['SWeightsType']         = ''
#pdfConfig['KKMassBinBounds']      = [ 1020. - 30., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1020. + 30. ] #[ 1020. - 30., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1020. + 30. ]
#pdfConfig['SWaveAmplitudeValues'] = (  [ 0.16, 0.08, 0.02, 0.02, 0.05,  0.14 ], [ 1.5, 0.8, 0.3, -0.4, -0.5, -0.6 ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ 0.8, 0.4, 0.1, 0.1, 0.2,  0.6 ], [ 1.8, 0.6, 0.2, -0.4, -0.6, -0.6 ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ -0.12, -0.25, -0.16, -0.07, -0.18, -0.37 ], [ -0.31, -0.15, -0.10, 0.01, 0.16, 0.10 ] )
pdfConfig['CSPValues']            = [ 0.498 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]
pdfConfig['lifetimeRange']        = ( 0.3, 14. )

pdfConfig['sameSideTagging']    = False  # nominal: False
pdfConfig['conditionalTagging'] = False  # nominal: True
pdfConfig['continuousEstWTag']  = False  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 100
pdfConfig['constrainTagging']   = False  # nominal: True

pdfConfig['eventTimeResolution']   = False  # nominal: True
pdfConfig['numTimeResBins']        = 100
pdfConfig['constrainTimeResScale'] = False  # nominal: True

pdfConfig['numEvents'] = 32000
pdfConfig['signalFraction'] = 0.67
pdfConfig['massRangeBackground'] = False

pdfConfig['amplitudeParam'] = 'phasesSWaveFrac' # default: 'bank' | nominal: 'phasesSWaveFrac'
pdfConfig['ASParam']        = 'deltaPerp'  # default/nominal: 'deltaPerp'
pdfConfig['AparParam']      = 'phase' # default: 'Mag2ReIm' | nominal: 'phase'

pdfConfig['constrainDeltaM'] = False  # nominal: True

pdfConfig['lambdaCPParam'] = 'lambPhi'  # default/nominal: 'lambSqPhi'

fastFit          = True
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
                                    if pdfConfig['timeEffType'] == 'Fit' else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504_unitAcceptance.root'
pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_40bins'
pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1ExclB_40bins'
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

if ( readData or generateData ) and doFit :
    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

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

        print '\nJvLFit: reparameterization of transversity amplitudes:'
        ampsData.Print()
        ampMeasList.Print('v')
        ampMeasCovs.Print()
        ampsFitResult.Print()
        ampsFitResult.covarianceMatrix().Print()

    print 'JvLFit: parameters:'
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

if ( readData or generateData ) and ( makeObservablePlots or pdfConfig['makePlots'] or makeKKMassPlots or dllPars ) :
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
    if   pdfConfig['eventTimeResolution'] : projWDataSet += [ timeRes ]

    if projWDataSet :
        bulkData = data.reduce( CutRange = 'Bulk' )
        projWData     = dict( ProjWData = ( data.reduce(  ArgSet = projWDataSet ), True ) )
        projWDataBulk = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataSet ), True ) )
    else :
        projWData     = dict()
        projWDataBulk = dict()

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
    print 'JvLFit: plotting lifetime and angular distributions'
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
    print 'JvLFit: plotting lifetime distribution'
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
    print 'JvLFit: plotting lifetime distributions tagged/untagged'
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
    print 'JvLFit: plotting angular distributions'
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

    if not pdfConfig['SFit'] :
        # plot signal mass
        print 'JvLFit: plotting mumuKK mass distribution'
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
    pdfBuild['massCanv'].Print(plotsFile + '(')
    pdfBuild['mumuMassCanv'].Print(plotsFile)
    pdfBuild['KKMassCanv'].Print(plotsFile)
    bkgTimeCanv.Print(plotsFile)
    pdfBuild['bkgAnglesSWeightCanv'].Print(plotsFile)
    pdfBuild['bkgAnglesSideBandCanv'].Print(plotsFile)
    pdfBuild['estWTagCanvOS'].Print(plotsFile)
    pdfBuild['estWTagCanvSS'].Print(plotsFile + '' if deltaSCanv else ')')

if deltaSCanv :
    deltaSCanv.Print( plotsFile + ( ')' if makeObservablePlots or pdfConfig['makePlots'] else '' ) )


###########################################################################################################################################
## make DLL plots  ##
#####################

if dllPars :
    # make delta log-likelihood plots
    if pdfConfig['amplitudeParam'] == 'phasesSWaveFrac' :
        MparMin  =  0.21
        MparMax  =  0.25
        RparMin  = -0.480
        RparMax  = -0.455
        IparMin  = -0.15
        IparMax  = +0.15
        MperpMin =  0.23
        MperpMax =  0.27
    else :
        MparMin =  0.40
        MparMax =  0.48
        RparMin = -0.69
        RparMax = -0.61
        if pdfConfig['AparParam'] == 'Mag2ReIm' :
            IparMin = -1.
            IparMax = +0.65
        else :
            IparMin = -0.25
            IparMax = +0.25
        MperpMin =  0.44
        MperpMax =  0.52

    wsPars =\
    {  'phiCP'           : ( '#DeltaNLL #phi_{s}',                     '#phi_{s}',                     -0.22,     0.22,    1, 0.001, 0.05 )
     , 'lambdaCP'        : ( '#DeltaNLL |#lambda|',                    '|#lambda|',                     0.85,     1.0,     1, 0.001, 0.01 )
     , 'lambdaCPSq'      : ( '#DeltaNLL |#lambda|^{2}',                '|#lambda|^{2}',                 0.7,      1.0,     1, 0.001, 0.01 )
     , 'avgCOddSum'      : ( '#DeltaNLL C_{Os}^{avg}',                 'C_{Os}^{avg}',                 -0.035,    0.100,   1, 0.001, 0.01 )
     , 'avgCOddTagged'   : ( '#DeltaNLL C_{Ot}^{avg}',                 'C_{Ot}^{avg}',                 -0.100,    0.155,   1, 0.001, 0.01 )
     , 'A0Mag2'          : ( '#DeltaNLL |A_{0}|^{2}',                  '|A_{0}|^{2}',                   0.51,     0.54,    1, 0.001, 0.01 )
     , 'AparMag2'        : ( '#DeltaNLL |A_{#parallel}|^{2}',          '|A_{#parallel}|^{2}',          MparMin,  MparMax,  1, 0.001, 0.01 )
     , 'ReApar'          : ( '#DeltaNLL Re(A_{#parallel})',            'Re(A_{#parallel})',            RparMin,  RparMax,  1, 0.001, 0.01 )
     , 'ImApar'          : ( '#DeltaNLL Im(A_{#parallel})',            'Im(A_{#parallel})',            IparMin,  IparMax,  1, 0.001, 0.01 )
     , 'cosAparPhase'    : ( '#DeltaNLL cos(#delta_{#parallel})',      'cos(#delta_{#parallel})',      -1.,      -0.92,    1, 0.001, 0.01 )
     , 'AparPhase'       : ( '#DeltaNLL #delta_{#parallel}',           '#delta_{#parallel}',            2.8,      3.7,     1, 0.001, 0.01 )
     , 'AperpMag2'       : ( '#DeltaNLL |A_{#perp}|^{2}',              '|A_{#perp}|^{2}',              MperpMin, MperpMax, 1, 0.001, 0.01 )
     , 'AperpPhase'      : ( '#DeltaNLL #delta_{#perp}',               '#delta_{#perp}',                2.3,      3.7,     1, 0.001, 0.01 )
     , 'sqrtfS_Re'       : ( '#DeltaNLL #sqrt{f_{S}}^{R}',             '#sqrt{f_{S}}^{R}',             -0.22,    -0.10,    1, 0.001, 0.01 )
     , 'sqrtfS_Im'       : ( '#DeltaNLL #sqrt{f_{S}}^{I}',             '#sqrt{f_{S}}^{I}',             -0.060,    0.085,   1, 0.001, 0.01 )
     , 'ReASOdd'         : ( '#DeltaNLL Re(A_{S} / A_{#perp})',        'Re(A_{S} / A_{#perp})',         0.20,     0.44,    1, 0.001, 0.01 )
     , 'ImASOdd'         : ( '#DeltaNLL Im(A_{S} / A_{#perp})',        'Im(A_{S} / A_{#perp})',        -0.06,     0.04,    1, 0.001, 0.01 )
     , 'f_S'             : ( '#DeltaNLL f_{S}',                        'f_{S}',                         0.00,     0.045,   1, 0.001, 0.01 )
     , 'ASPhase'         : ( '#DeltaNLL #delta_{S}',                   '#delta_{S}',                    2.3,      3.7,     1, 0.001, 0.01 )
     , 'ASOddMag2'       : ( '#DeltaNLL A_{S}^{2} / A_{#perp}^{2}',    'A_{S}^{2} / A_{#perp}^{2}',     0.0,      0.18,    1, 0.001, 0.01 )
     , 'ASOddPhase'      : ( '#DeltaNLL #delta_{S} - #delta_{#perp}',  '#delta_{S} - #delta_{#perp}',  -0.22,     0.18,    1, 0.001, 0.01 )
     , 'ASOddPhase_bin0' : ( '#DeltaNLL #delta_{S0} - #delta_{#perp}', '#delta_{S0} - #delta_{#perp}', -2.,       5.,      1, 0.001, 0.01 )
     , 'ASOddPhase_bin1' : ( '#DeltaNLL #delta_{S1} - #delta_{#perp}', '#delta_{S1} - #delta_{#perp}', -2.,       5.,      1, 0.001, 0.01 )
     , 'ASOddPhase_bin2' : ( '#DeltaNLL #delta_{S2} - #delta_{#perp}', '#delta_{S2} - #delta_{#perp}', -2.,       5.,      1, 0.001, 0.01 )
     , 'ASOddPhase_bin3' : ( '#DeltaNLL #delta_{S3} - #delta_{#perp}', '#delta_{S3} - #delta_{#perp}', -2.,       5.,      1, 0.001, 0.01 )
     , 'ASOddPhase_bin4' : ( '#DeltaNLL #delta_{S4} - #delta_{#perp}', '#delta_{S4} - #delta_{#perp}', -2.,       5.,      1, 0.001, 0.01 )
     , 'ASOddPhase_bin5' : ( '#DeltaNLL #delta_{S5} - #delta_{#perp}', '#delta_{S5} - #delta_{#perp}', -2.,       5.,      1, 0.001, 0.01 )
    }

    # check DNLL parameters
    phiCPPar = False
    for par in dllPars :
        assert par[0] in wsPars, 'JvLFit - ERROR: unknown DLL parameter: "%s"' % par[0]
        assert par[0] in ws,     'JvLFit - ERROR: DLL parameter "%s" does not exist in work space' % par[0]
        if par[0] == 'phiCP' : phiCPPar = True

    # float/fix values of some parameters
    if constLambdaCP :
        pdfBuild['lambdaCP'].setConstant('lambdaCPSq') if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
            else pdfBuild['lambdaCP'].setConstant('lambdaCP')
    if constAvgCEvenOdd :
        for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')

    if 'sig_ATagBBbar' in ws : ws['sig_ATagBBbar'].setConstant()
    for bin in range(5) :
        if 'sig_ATagBBbar_bin%d' % bin in ws : ws[ 'sig_ATagBBbar_bin%d' % bin ].setConstant()
    pdfBuild['tagCatsOS'].setConstant('.*')
    pdfBuild['tagCatsSS'].setConstant('.*')
    #pdfBuild['lifetimeParams'].setConstant('dM|Gamma')
    pdfBuild['timeResModel'].setConstant('.*')
    pdfBuild['signalBMass'].setConstant('.*')
    if not pdfConfig['SFit'] :
        pdfBuild['backgroundBMass'].setConstant('.*')
        pdfBuild['backgroundTime'].setConstant('.*')

    # build NLL
    from ROOT import RooFit, RooArgSet, RooArgList, RooFormulaVar, TCanvas
    nll = pdf.createNLL( fitData, **fitOpts )

    print 120 * '='
    print 'JvLFit: parameters in NLL:'
    for par in nll.getVariables() : par.Print()
    print 120 * '='

    # create DNLL/PLL plots
    dllCanvs = [ ]
    canvFileName = plotsFile[ : -3 ] + 'DLLs.ps'
    for parIter, ( par, doDLL, doPLL, doPara ) in enumerate(dllPars) :
        rooPar = ws[par]
        parFrame = rooPar.frame(  RooFit.Range( wsPars[par][2], wsPars[par][3] )
                                , RooFit.Bins( wsPars[par][4] )
                                , RooFit.Title( wsPars[par][0] )
                               )

        if doPara and rooPar.getError() > 0. :
            parabola = RooFormulaVar(  'parabola', 'parabola'
                                     , '0.5*(@0-{0:.6f})*(@0-{0:.6f})/{1:.6f}/{1:.6f}'.format( rooPar.getVal(), rooPar.getError() )
                                     , RooArgList(rooPar)
                                    )
            parabola.plotOn( parFrame, RooFit.LineColor(RooFit.kBlack), RooFit.Precision(0.001) )

        if doDLL :
            print 'JvLFit: plotting Delta -log(L) for %s' % par
            nll.plotOn( parFrame, RooFit.ShiftToZero(), RooFit.LineColor(kBlue), RooFit.Precision( wsPars[par][5] ) )

        if doPLL :
            print 'JvLFit: plotting profile Delta -log(L) for %s' % par
            pll = nll.createProfile( RooArgSet( rooPar ) )
            pll.plotOn( parFrame, RooFit.LineColor(kRed), RooFit.Precision( wsPars[par][6] ) )

        parFrame.SetMinimum(0.)
        if parFrame.GetMaximum() > 15. : parFrame.SetMaximum(15.)
        parFrame.GetXaxis().SetTitle( wsPars[par][1] )
        parFrame.GetYaxis().SetTitle('#DeltaNLL')

        dllCanvs.append( TCanvas( 'dllCanv%d' % parIter , 'DLL canvas' ) )
        parFrame.Draw()

        for canvIter, canv in enumerate(dllCanvs) :
            if len(dllCanvs) == 1 or canvIter not in [ 0, len(dllCanvs) - 1 ] : namePF = ''
            elif canvIter == 0 : namePF = '('
            else : namePF = ')'
            canv.Print( canvFileName + namePF )
