###########################################################################################################################################
## set script parameters ##
###########################
from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
generateData            = False
doFit                   = True
fastFit                 = False
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
pdfConfig['blind']      = False
pdfConfig['nominalPdf'] = False
sumW2Error              = False

plotsFile = 'plots/JvLSFit.ps' if pdfConfig['SFit']\
       else 'plots/JvLCFit.ps'
parameterFile = 'JvLSFit.par' if pdfConfig['SFit'] else 'JvLCFit.par'

if readData :
    pdfConfig['nTupleName'] = 'DecayTree'
    pdfConfig['nTupleFile'] = '/home/raaij/data/Bs2JpsiPhi_2011_biased_unbiased.root'
else :
    pdfConfig['nTupleName'] = None
    pdfConfig['nTupleFile'] = None

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLSFit.root' if pdfConfig['SFit'] else 'JvLCFit.root'

dllPars = [ ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

# fit options
fitOpts = dict(  NumCPU              = 12
               , Optimize            = 1
               , Timer               = 1
               , Verbose             = True
#               , Minos               = True
#               , Hesse               = False
               , Minimizer           = 'Minuit2'
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
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['parameterizeKKMass']   = ''  # default/nominal: ''
pdfConfig['ambiguityParameters']  = False
pdfConfig['KKMassBinBounds']      = [ 1020. - 12., 1020. + 12. ] #[ 1020. - 30., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1020. + 30. ]
#pdfConfig['SWaveAmplitudeValues'] = (  [ 0.8, 0.4, 0.1, 0.1, 0.2,  0.6 ], [ 1.8, 0.6, 0.2, -0.4, -0.6, -0.6 ] )
pdfConfig['SWaveAmplitudeValues'] = (  [ -0.12, -0.25, -0.16, -0.07, -0.18, -0.37 ], [ -0.31, -0.15, -0.10, 0.01, 0.16, 0.10 ] )
pdfConfig['CSPValues']            = [ 0.498 ] # [ 0.4976 ] # [ 0.3263 ] # [ 0.9663, 0.9562, 0.9255, 0.9255, 0.9562, 0.9663 ]

pdfConfig['sameSideTagging']    = False  # nominal: False
pdfConfig['conditionalTagging'] = True  # nominal: True
pdfConfig['continuousEstWTag']  = False  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 100
pdfConfig['constrainTagging']   = True  # nominal: True

pdfConfig['eventTimeResolution'] = True  # nominal: True
pdfConfig['numTimeResBins']      = 100

pdfConfig['numEvents'] = 32000
pdfConfig['signalFraction'] = 0.67
pdfConfig['massRangeBackground'] = False

pdfConfig['amplitudeParam'] = 'phasesSWaveFrac' # default: 'bank' | nominal: 'phasesSWaveFrac'
pdfConfig['ASParam']        = 'deltaPerp'  # default/nominal: 'deltaPerp'
pdfConfig['AparParam']      = 'phase' # default: 'Mag2ReIm' | nominal: 'phase'

pdfConfig['constrainDeltaM'] = True  # nominal: True

pdfConfig['lambdaCPParam'] = 'lambSqPhi'  # default/nominal: 'lambSqPhi'
constLambdaCP = False  # default/nominal: False

constTagCatCoefs = True  # default: True / nominal: False
constAvgCEvenOdd = True  # default: False / nominal: True
constWTagAsyms   = True  # default/nominal: True
constCSP         = True  # default/nominal: True

if not readData :
    pdfConfig['tagCats'] = [  ( 'Untagged',  0, 0.500001, 0.500, 0.505, 0., 0.6683, 0. )
                            , ( 'TagCat1',   1, 0.499999, 0.484, 0.489, 0., 0.0107, 0. )
                            , ( 'TagCat2',   2, 0.478,    0.467, 0.470, 0., 0.0472, 0. )
                            , ( 'TagCat3',   3, 0.457,    0.447, 0.449, 0., 0.0494, 0. )
                            , ( 'TagCat4',   4, 0.437,    0.427, 0.429, 0., 0.0425, 0. )
                            , ( 'TagCat5',   5, 0.417,    0.408, 0.410, 0., 0.0359, 0. )
                            , ( 'TagCat6',   6, 0.399,    0.390, 0.391, 0., 0.0325, 0. )
                            , ( 'TagCat7',   7, 0.381,    0.372, 0.373, 0., 0.0224, 0. )
                            , ( 'TagCat8',   8, 0.363,    0.354, 0.354, 0., 0.0195, 0. )
                            , ( 'TagCat9',   9, 0.344,    0.334, 0.333, 0., 0.0143, 0. )
                            , ( 'TagCat10', 10, 0.324,    0.313, 0.311, 0., 0.0122, 0. )
                            , ( 'TagCat11', 11, 0.303,    0.291, 0.289, 0., 0.0127, 0. )
                            , ( 'TagCat12', 12, 0.280,    0.270, 0.266, 0., 0.0100, 0. )
                            , ( 'TagCat13', 13, 0.257,    0.246, 0.242, 0., 0.0086, 0. )
                            , ( 'TagCat14', 14, 0.233,    0.222, 0.217, 0., 0.0060, 0. )
                            , ( 'TagCat15', 15, 0.208,    0.195, 0.189, 0., 0.0029, 0. )
                            , ( 'TagCat16', 16, 0.181,    0.167, 0.160, 0., 0.0027, 0. )
                            , ( 'TagCat17', 17, 0.153,    0.141, 0.133, 0., 0.0019, 0. )
                            , ( 'TagCat18', 18, 0.124,    0.111, 0.102, 0., 0.0005, 0. )
                           ]

pdfConfig['timeEffHistFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root'
pdfConfig['timeEffHistName'] = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_40bins'

pdfConfig['angEffMomentsFile'] = 'effMomentsTransBasisBaseline' if not pdfConfig['nominalPdf'] and pdfConfig['transversityAngles']\
                                 else '/home/raaij/data/effMomentsHelBasisBaseline'

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
BMass      = pdfBuild['observables']['BMass']
mumuMass   = pdfBuild['observables']['mumuMass']
KKMass     = pdfBuild['observables']['KKMass']
estWTagOS  = pdfBuild['observables']['estWTagOS']
timeRes    = pdfBuild['observables']['timeRes']

if not pdfConfig['SFit'] : obsSetP2VV.append(BMass)

if not pdfBuild['iTagZeroTrick'] :
    tagCatP2VVOS = pdfBuild['observables']['tagCatP2VVOS']
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
    fitData = pdfBuild['sigSWeightData']

else :
    fitData = pdfBuild['data']


###########################################################################################################################################
## fit data ##
##############

if ( readData or generateData ) and doFit :
    # float/fix values of some parameters
    if constLambdaCP :
        pdfBuild['lambdaCP'].setConstant('lambdaCPSq') if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
            else pdfBuild['lambdaCP'].setConstant('lambdaCP')
    for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] :
        CEvenOdd.setConstant('avgCEven.*')
        if pdfConfig['nominalPdf'] or constAvgCEvenOdd : CEvenOdd.setConstant( 'avgCOdd.*', True )

    if pdfConfig['nominalPdf'] and not constTagCatCoefs : pdfBuild['taggingParams'].setConstant( 'tagCatCoef.*', False )

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

    if fastFit :
        pdfBuild['lambdaCP'].setConstant('lambdaCPSq')
        for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')
        pdfBuild['tagCatsOS'].setConstant('.*')
        pdfBuild['tagCatsSS'].setConstant('.*')
        pdfBuild['lifetimeParams'].setConstant('dM|Gamma')
        pdfBuild['timeResModel'].setConstant('.*')
        pdfBuild['signalBMass'].setConstant('.*')
        if not pdfConfig['SFit'] :
            pdfBuild['backgroundBMass'].setConstant('.*')
            pdfBuild['backgroundTime'].setConstant('.*')
        pdfBuild['amplitudes'].setConstant('C_SP')

    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    if pdfConfig['SFit'] : fitResult = pdf.fitTo( fitData, SumW2Error = sumW2Error, Save = True, **fitOpts )
    else                 : fitResult = pdf.fitTo( fitData,                          Save = True, **fitOpts )

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
        bulkData = fitData.reduce( CutRange = 'Bulk' )
        projWData     = dict( ProjWData = ( fitData.reduce(  ArgSet = projWDataSet ), True ) )
        projWDataBulk = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataSet ), True ) )
    else :
        projWData     = dict()
        projWDataBulk = dict()

if pdfConfig['makePlots'] :
    # plot background time
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

if makeKKMassPlots and pdfConfig['parameterizeKKMass'] and fitResult and pdfConfig['amplitudeParam'] == 'bank'\
        and pdfConfig['ASParam'] != 'ReIm' :
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
        plot(  pad, obs, fitData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
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
        plot(  pad, time, fitData, pdf, yTitle = yTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                    )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                         )
             , pdfOpts    = dict( list( projWDataBulk.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , components = comps
            )

    # plot lifetime (tagged/untagged)
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
        plot(  pad, time, fitData, pdf, yTitle = yTitle, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                     )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , components = comps
            )

    # plot angles
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
        plot(  pad, obs, fitData, pdf, addPDFs = addPDFs, xTitle = xTitle, yTitle = yTitle
             , frameOpts   = dict( Bins = nBins, Title = plotTitle                                                )
             , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts                    )
             , pdfOpts     = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , addPDFsOpts = [ dict( list( projWData.items() ), LineColor = kRed, LineWidth = lineWidth, **pdfCuts ) ]
             , components  = comps
            )

    if not pdfConfig['SFit'] :
        # plot signal mass
        pad = pdfBuild['massCanv'].cd(2)
        plot(  pad, BMass, fitData, pdf
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

sums = {
    'numEv'    : 0.
  , 'numOS'    : 0., 'numOSExcl'    : 0.
  , 'numSS'    : 0., 'numSSExcl'    : 0.
  , 'numComb'  : 0., 'numCombExcl'  : 0.
  , 'etaOS'    : 0., 'etaOSExcl'    : 0.
  , 'etaSS'    : 0., 'etaSSExcl'    : 0.
  , 'etaComb'  : 0., 'etaCombExcl'  : 0.
  , 'wOS'      : 0., 'wOSExcl'      : 0.
  , 'wSS'      : 0., 'wSSExcl'      : 0.
  , 'wComb'    : 0., 'wCombExcl'    : 0.
  , 'dilOS'    : 0., 'dilOSExcl'    : 0.
  , 'dilSS'    : 0., 'dilSSExcl'    : 0.
  , 'dilComb'  : 0., 'dilCombExcl'  : 0.
  , 'dil2OS'   : 0., 'dil2OSExcl'   : 0.
  , 'dil2SS'   : 0., 'dil2SSExcl'   : 0.
  , 'dil2Comb' : 0., 'dil2CombExcl' : 0.
}

for varSet in fitData :
  weight = fitData.weight()
  sums['numEv'] += weight

  if varSet.getCatIndex('tagdecision_os') != 0 :
    etaOS  = varSet.getRealValue('tagomega_os')
    wTagOS = 0.392 + 1.000 * ( etaOS - 0.392 )
    dilOS  = 1. - 2. * wTagOS

    sums['numOS']  += weight
    sums['etaOS']  += weight * etaOS
    sums['wOS']    += weight * wTagOS
    sums['dilOS']  += weight * dilOS
    sums['dil2OS'] += weight * dilOS**2

    if varSet.getCatIndex('tagdecision_ss') == 0 :
      sums['numOSExcl']  += weight
      sums['etaOSExcl']  += weight * etaOS
      sums['wOSExcl']    += weight * wTagOS
      sums['dilOSExcl']  += weight * dilOS
      sums['dil2OSExcl'] += weight * dilOS**2

      sums['numComb']  += weight
      sums['etaComb']  += weight * etaOS
      sums['wComb']    += weight * wTagOS
      sums['dilComb']  += weight * dilOS
      sums['dil2Comb'] += weight * dilOS**2

  if varSet.getCatIndex('tagdecision_ss') != 0 :
    etaSS  = varSet.getRealValue('tagomega_ss')
    wTagSS = 0.350 + 1.00 * ( etaSS - 0.350 )
    dilSS  = 1. - 2. * wTagSS

    sums['numSS']  += weight
    sums['etaSS']  += weight * etaSS
    sums['wSS']    += weight * wTagSS
    sums['dilSS']  += weight * dilSS
    sums['dil2SS'] += weight * dilSS**2

    if varSet.getCatIndex('tagdecision_os') == 0 :
      sums['numSSExcl']  += weight
      sums['etaSSExcl']  += weight * etaSS
      sums['wSSExcl']    += weight * wTagSS
      sums['dilSSExcl']  += weight * dilSS
      sums['dil2SSExcl'] += weight * dilSS**2

      sums['numComb']  += weight
      sums['etaComb']  += weight * etaSS
      sums['wComb']    += weight * wTagSS
      sums['dilComb']  += weight * dilSS
      sums['dil2Comb'] += weight * dilSS**2

  if varSet.getCatIndex('tagdecision_os') != 0 and varSet.getCatIndex('tagdecision_ss') != 0 :
    dilSign = +1. if varSet.getCatIndex('tagdecision_os') == varSet.getCatIndex('tagdecision_ss') else -1.
    dilComb = ( dilOS + dilSign * dilSS ) / ( 1. + dilSign * dilOS * dilSS )
    wTagComb = ( 1. - dilComb ) / 2.

    sums['numComb']  += weight
    sums['etaComb']  += weight * wTagComb
    sums['wComb']    += weight * wTagComb
    sums['dilComb']  += weight * dilComb
    sums['dil2Comb'] += weight * dilComb**2

    sums['numCombExcl']  += weight
    sums['etaCombExcl']  += weight * wTagComb
    sums['wCombExcl']    += weight * wTagComb
    sums['dilCombExcl']  += weight * dilComb
    sums['dil2CombExcl'] += weight * dilComb**2

print
print 'number of events:       %.4f' % sums['numEv']
print 'number of OS events:    %.4f (%.4f)'   % ( sums['numOS'],   sums['numOSExcl']   )
print 'number of SS events:    %.4f (%.4f)'   % ( sums['numSS'],   sums['numSSExcl']   )
print 'number of Comb. events: %.4f (%.4f)\n' % ( sums['numComb'], sums['numCombExcl'] )

print 'OS    eff.: %.4f (%.4f)'   % ( sums['numOS']   / sums['numEv'], sums['numOSExcl']   / sums['numEv'] )
print 'SS    eff.: %.4f (%.4f)'   % ( sums['numSS']   / sums['numEv'], sums['numSSExcl']   / sums['numEv'] )
print 'Comb. eff.: %.4f (%.4f)\n' % ( sums['numComb'] / sums['numEv'], sums['numCombExcl'] / sums['numEv'] )

print 'OS    <eta>: %.4f (%.4f)'   % ( sums['etaOS']   / sums['numOS'],   sums['etaOSExcl']   / sums['numOSExcl']   )
print 'SS    <eta>: %.4f (%.4f)'   % ( sums['etaSS']   / sums['numSS'],   sums['etaSSExcl']   / sums['numSSExcl']   )
print 'Comb. <eta>: %.4f (%.4f)\n' % ( sums['etaComb'] / sums['numComb'], sums['etaCombExcl'] / sums['numCombExcl'] )

print 'OS    <w>: %.4f (%.4f)'   % ( sums['wOS']   / sums['numOS'],   sums['wOSExcl']   / sums['numOSExcl']   )
print 'SS    <w>: %.4f (%.4f)'   % ( sums['wSS']   / sums['numSS'],   sums['wSSExcl']   / sums['numSSExcl']   )
print 'Comb. <w>: %.4f (%.4f)\n' % ( sums['wComb'] / sums['numComb'], sums['wCombExcl'] / sums['numCombExcl'] )

print 'OS    <dil>: %.4f (%.4f)'   % ( sums['dilOS']   / sums['numOS'],   sums['dilOSExcl']   / sums['numOSExcl']   )
print 'SS    <dil>: %.4f (%.4f)'   % ( sums['dilSS']   / sums['numSS'],   sums['dilSSExcl']   / sums['numSSExcl']   )
print 'Comb. <dil>: %.4f (%.4f)\n' % ( sums['dilComb'] / sums['numComb'], sums['dilCombExcl'] / sums['numCombExcl'] )

print 'OS    <dil2>: %.4f (%.4f)'   % ( sums['dil2OS']   / sums['numOS'],   sums['dil2OSExcl']   / sums['numOSExcl']   )
print 'SS    <dil2>: %.4f (%.4f)'   % ( sums['dil2SS']   / sums['numSS'],   sums['dil2SSExcl']   / sums['numSSExcl']   )
print 'Comb. <dil2>: %.4f (%.4f)\n' % ( sums['dil2Comb'] / sums['numComb'], sums['dil2CombExcl'] / sums['numCombExcl'] )

print 'OS    <eff * dil2>: %.4f%% (%.4f%%)'   % ( sums['dil2OS']   / sums['numEv'] * 100., sums['dil2OSExcl']   / sums['numEv'] * 100. )
print 'SS    <eff * dil2>: %.4f%% (%.4f%%)'   % ( sums['dil2SS']   / sums['numEv'] * 100., sums['dil2SSExcl']   / sums['numEv'] * 100. )
print 'Comb. <eff * dil2>: %.4f%% (%.4f%%)\n' % ( sums['dil2Comb'] / sums['numEv'] * 100., sums['dil2CombExcl'] / sums['numEv'] * 100. )
