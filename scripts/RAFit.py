###########################################################################################################################################
## set script parameters ##
###########################

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
pdfConfig['dataSample'] = '' # ( None, 100260, '' )  # '' / 'Summer2011' / 'runNumber % 2 == 1'
pdfConfig['selection']  = 'paper2012' # 'paper2012' # 'HLT1Unbiased' # 'timeEffFit'
generateData            = False
doFit                   = True
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
pdfConfig['blind']      = False
corrSFitErr             = 'sumWeight' # [ 1., 0.700, 0.952, 0.938, 0.764 ] # '' / 'matrix' / 'sumWeight'
randomParVals           = ( ) # ( 1., 12346 ) # ( 2., 12345 )

#plotsFile = 'plots/JvLSFit_SWavePhases.ps'
plotsFile = 'plots/JvLSFit.ps' if pdfConfig['SFit']\
       else 'plots/JvLCFit.ps'
parameterFile = None # 'JvLSFit.par' if pdfConfig['SFit'] else 'JvLCFit.par'

if readData :
    pdfConfig['nTupleName'] = 'DecayTree'
    ## pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_2012_ntupleB_20121212.root'
    pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
    #pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDown.root'
    #pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagUp.root'
    #pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Pass3-version2_Bs_050711_nocut_Phi_P2VV.root'
    #pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20121010.root'
    #pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_DGs0_MC11a_ntupleB_for_fitting_20121119.root'
    pdfConfig['nominalDataSet'] = False
else :
    pdfConfig['nTupleName'] = None
    pdfConfig['nTupleFile'] = None

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLSFit.root' if pdfConfig['SFit'] else 'JvLCFit.root'

## MinosPars = [ 'AparPhase', 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5' ]
MinosPars = [ ]
dllPars = [ ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

# fit options
fitOpts = dict(  NumCPU    = 6
               , Optimize  = 2
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
               , Offset    = True
              )
pdfConfig['fitOptions'] = fitOpts

# plot options
from ROOT import gStyle, kBlack, kBlue, kRed, kGreen, kMagenta, kSolid, kDashed, kFullCircle, kFullSquare
from P2VV.Load import RooFitOutput, LHCbStyle
lineWidth     = 3
lineColor     = kBlue
markStyle     = 8
markSize      = 0.6
markColor     = kBlack
markLineWidth = 2
gStyle.SetLineStyleString( 5, ' 40 20 10 20'  )
gStyle.SetLineStyleString( 7, ' 40 20'        )
gStyle.SetLineStyleString( 9, ' 100 20'       )

# PDF options
pdfConfig['transversityAngles'] = False  # default: False | nominal: True

pdfConfig['bkgAnglePdf']          = 'hybrid'  # default/nominal: ''
pdfConfig['sigTaggingPdf']        = 'tagUntag'  # default: 'tagUntag' | nominal: 'tagCats'
pdfConfig['bkgTaggingPdf']        = 'tagUntagRelative'  # default: 'tagUntagRelative' | 'tagCatsRelative'
pdfConfig['multiplyByTagPdf']     = False
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['timeEffType']          = 'paper2012' # 'paper2012' # 'HLT1Unbiased' # 'fit'
pdfConfig['multiplyByAngEff']     = 'weights'  # default: 'basis012'
pdfConfig['parameterizeKKMass']   = 'simultaneous'  # default/nominal: 'simultaneous'
pdfConfig['ambiguityParameters']  = False
pdfConfig['lifetimeRange']        = ( 0.3, 14. )
pdfConfig['SWeightsType']         = 'simultaneousFreeBkg'  # default/nominal: 'simultaneousFreeBkg'
pdfConfig['KKMassBinBounds']      = [ 1008., 1032. ] # [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['SWaveAmplitudeValues'] = (  [ (0.025, 0.01)]
                                      , [ (0.011,  0.16 )])
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.28, 0.11), (0.06, 0.02), (0.04, 0.02), (0.27, 0.07) ]
#                                     , [ (2.7,  0.4 ), (0.22,   0.14), (-0.11, 0.17 ), (-0.97, 0.3 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.28, 0.11), (0.05, 0.02), (0.27, 0.07) ]
#                                     , [ (2.7,  0.4 ), (0.,   0.15), (-0.97, 0.3 ) ] )
pdfConfig['CSPValues']            = [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.959, 0.498, 0.968 ]

pdfConfig['sameSideTagging']    = True  # nominal: True
pdfConfig['conditionalTagging'] = True  # nominal: True
pdfConfig['continuousEstWTag']  = True  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 50
pdfConfig['constrainTagging']   = 'constrain'  # nominal: 'constrain'

pdfConfig['timeResType']           = 'eventStQuad' # 'eventNoMean' # 'event' # 'eventNoMean' # 'eventStLinear' # 'eventStQuad'
pdfConfig['numTimeResBins']        = 400
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

pdfConfig['timeEffHistFile']      = '/bfys/raaij/p2vv/data/start_values.root'\
                                    if pdfConfig['timeEffType'] == 'fit' else\
                                    '/bfys/raaij/p2vv/data/Bs_HltPropertimeAcceptance_Data-20120816.root'
#                                    '/bfys/raaij/p2vv/data/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root'
#                                    '/bfys/raaij/p2vv/data/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504_unitAcceptance.root'
#pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['angEffMomentsFile']    = '/bfys/raaij/p2vv/data/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'\
                                    if pdfConfig['transversityAngles'] else\
                                    '/bfys/raaij/p2vv/data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
#                                    '/bfys/raaij/p2vv/data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'
#                                    '/bfys/raaij/p2vv/data/hel_UB_UT_trueTime_BkgCat050_KK30_PHSP_Basis'

if pdfConfig['transversityAngles'] :
    pdfConfig['angleNames'] = (  ( 'trcospsi',   'cos(#psi_{tr})'   )
                               , ( 'trcostheta', 'cos(#theta_{tr})' )
                               , ( 'trphi',      '#phi_{tr}'        )
                              )
else :
    pdfConfig['angleNames'] = (  ( 'helcosthetaK', 'cos(#theta_{K})'   )
                               , ( 'helcosthetaL', 'cos(#theta_{#mu})' )
                               , ( 'helphi',       '#phi_{h}'          )
                              )
angleNames = pdfConfig['angleNames']

numBins = ( 50, 21, 21, 21 )
pdfConfig['numTimeBins'] = 30
numAngleBins = ( 20, 20, 20 )
pdfConfig['numAngleBins'] = ( 5, 7, 9 )


###########################################################################################################################################
## build PDF ##
###############

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
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
    if not pdfConfig['transversityAngles'] :
        from P2VV.GeneralUtils import addTransversityAngles
        addTransversityAngles( fitData, 'trcospsi',          'trcostheta',        'trphi'
                                      , angles[0].GetName(), angles[1].GetName(), angles[2].GetName() )

    # write data to file
    from P2VV.GeneralUtils import writeData
    writeData( dataSetFile, dataSetName, fitData )

elif pdfConfig['SFit'] :
    defData = pdfBuild['sigSWeightData']
    sigData = pdfBuild['sigSWeightData']
    bkgData = pdfBuild['bkgSWeightData']
    if corrSFitErr == 'sumWeight'\
            or ( type(corrSFitErr) != str and hasattr( corrSFitErr, '__iter__' ) and hasattr( corrSFitErr, '__getitem__' ) ) :
        from P2VV.GeneralUtils import correctSWeights
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

if pdfConfig['parameterizeKKMass'] == 'functions' :
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
        if not pdfConfig['sameSideTagging'] :
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

#pdfBuild['lifetimeParams'].parameter('Gamma').setVal(0.72)

#pdfBuild['lambdaCP'].setConstant('rhoCP_m')
#pdfBuild['lambdaCP'].setConstant('rhoCP_A.*')
#pdfBuild['lambdaCP'].setConstant('phiCPRel_m')
#pdfBuild['lambdaCP'].setConstant('phiCPRel_A.*')

#ws['f_S'].setVal(0.026)
#ws['ASOddPhase'].setVal(0.03)
#ws['f_S'].setConstant()
#ws['ASOddPhase'].setConstant()

#ws['f_S_bin0'].setConstant()
#ws['f_S_bin1'].setConstant()
#ws['f_S_bin2'].setConstant()
#ws['f_S_bin3'].setConstant()
#ws['f_S_bin4'].setConstant()
#ws['f_S_bin5'].setConstant()

#ws['timeResSF'].setVal(1.4476)

#ws['ASOddPhase'].setMin(-6.)
#ws['ASOddPhase'].setMax(6.)
#ws['ASOddPhase_bin0'].setMin(-6.)
#ws['ASOddPhase_bin0'].setMax(6.)
#ws['ASOddPhase_bin1'].setMin(-6.)
#ws['ASOddPhase_bin1'].setMax(6.)
#ws['ASOddPhase_bin2'].setMin(-6.)
#ws['ASOddPhase_bin2'].setMax(6.)
#ws['ASOddPhase_bin3'].setMin(-6.)
#ws['ASOddPhase_bin3'].setMax(6.)

#ws['Gamma'].setVal(0.683)
#ws['dGamma'].setVal(0.068)
#ws['dM'].setVal(17.85)
#ws['A0Mag2'].setVal(0.6)
#ws['AperpMag2'].setVal(0.16)
#ws['f_S'].setVal(0.)
#ws['f_S'].setConstant()
#ws['AparPhase'].setVal(2.5)
#ws['AperpPhase'].setVal(-0.14)
#ws['ASOddPhase'].setVal(0.)
#ws['ASOddPhase'].setConstant()
#ws['lambdaCP'].setVal(1.)
#ws['lambdaCP'].setConstant()
#ws['phiCP'].setVal(-0.03)

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

assert(False)

from ROOT import RooArgSet
for constraint in pdfBuild._fullPdf.ExternalConstraints():
    if constraint.GetName().find('average_constraint') == -1:
        continue
    pars = constraint.getParameters(RooArgSet())
    phase = pars.find('ASOddPhase')
    if phase:
        print 'setting ASOddPhase constant for %s' % constraint.GetName()
        phase.setConstant(True)

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
    if pdfConfig['amplitudeParam'] == 'bank' and pdfConfig['ASParam'] != 'ReIm' and pdfConfig['AparParam'] == 'Mag2ReIm' :
        from ROOT import RooArgSet
        parList = fitResult.floatParsFinal()
        parSet  = RooArgSet(parList)
        parCovs = fitResult.covarianceMatrix()

        from math import pi
        from ROOT import RooRealVar, RooArgList
        from P2VV.Parameterizations.DecayAmplitudes import A02, Aperp2, Apar2, A0Ph, AperpPh, AparPh, f_S, AS2, ASPh
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
    from P2VV.GeneralUtils import plot, _P2VVPlotStash
    from ROOT import TCanvas

    # create projection data set for conditional observables
    if pdfConfig['SFit'] :
        comps = None
    else :
        comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                }

    from P2VV.RooFitWrappers import SimultaneousPdf
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
        print 'JvLFit: plot projection data set:'
        projWData['ProjWData'][0].Print()
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
    nKKBins = pdfBuild['KKMassBinning'].numBins()

    KKMassVals = [ pdfBuild['KKMassBinning'].binLow(binIter) for binIter in range(nKKBins) ]
    KKMassVals.append( pdfBuild['KKMassBinning'].highBound() )

    deltaSVals     = [  pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getVal()     for binIter in range(nKKBins) ]
    deltaSLowErrs  = [ -pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getErrorLo() for binIter in range(nKKBins) ]
    deltaSHighErrs = [  pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getErrorHi() for binIter in range(nKKBins) ]

    if pdfConfig['ambiguityParameters'] :
        from math import pi
        deltaSVals = [ pi - val for val in deltaSVals ]
        tempErrs       = deltaSLowErrs
        deltaSLowErrs  = deltaSHighErrs
        deltaSHighErrs = tempErrs

    from P2VV.GeneralUtils import plotSWavePhases
    deltaSCanv = plotSWavePhases(  MassBins         = KKMassVals
                                 , DeltaSValues     = deltaSVals
                                 , DeltaSLowErrors  = deltaSLowErrs
                                 , DeltaSHighErrors = deltaSHighErrs
                                 , LHCbTextLine1    = ' LHCb'
                                 , LHCbTextLine2    = '#sqrt{s} = 7 TeV, L = 0.37 fb^{-1}' if pdfConfig['dataSample'] == 'Summer2011'\
                                                      else '#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}'
                                )

else :
    deltaSCanv = None

if makeObservablePlots and not pdfBuild['iTagZeroTrick'] :
    from ROOT import RooRealVar, RooConstVar, RooCategory, RooArgSet, RooCustomizer, TPaveText

    # LHCb text labels
    LHCbText = TPaveText( 0.19, 0.82, 0.32, 0.91, 'NDC' )
    LHCbText.AddText('LHCb')
    LHCbText.SetShadowColor(0)
    LHCbText.SetFillStyle(0)
    LHCbText.SetBorderSize(0)
    LHCbText.SetTextAlign(12)
    _P2VVPlotStash.append(LHCbText)

    LHCbTextMore = TPaveText( 0.45, 0.70, 0.90, 0.90, 'NDC' )
    LHCbTextMore.AddText('LHCb')
    LHCbTextMore.AddText( '#sqrt{s} = 7 TeV, L = 0.37 fb^{-1}' if pdfConfig['dataSample'] == 'Summer2011'\
                      else '#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}' )
    LHCbTextMore.SetShadowColor(0)
    LHCbTextMore.SetFillStyle(0)
    LHCbTextMore.SetBorderSize(0)
    LHCbTextMore.SetTextAlign(12)
    _P2VVPlotStash.append(LHCbTextMore)

    dummyCanv = TCanvas()
    LHCbText.Draw()
    LHCbTextMore.Draw()
    dummyCanv.Update()

    # even, odd and S-wave PDFs
    AparMag2Cust = RooRealVar( 'AparMag2Cust', 'AparMag2Cust', ws['AparMag2'].getVal(), 0., 1. )
    zeroCust = RooConstVar( 'zeroCust', 'zeroCust', 1.e-6 )
    evenPdfCust = RooCustomizer( pdf._var, 'even' )
    oddPdfCust  = RooCustomizer( pdf._var, 'odd'  )
    SPdfCust    = RooCustomizer( pdf._var, 'S'    )

    evenPdfCust.replaceArg( ws['AparMag2'],  AparMag2Cust )
    evenPdfCust.replaceArg( ws['AperpMag2'], zeroCust     )
    if pdfConfig['parameterizeKKMass'] :
        for bin in range( pdfBuild['KKMassBinning'].numBins() ) : evenPdfCust.replaceArg( ws[ 'f_S_bin%d' % bin ], zeroCust )
    else :
        evenPdfCust.replaceArg( ws['f_S'], zeroCust )
    evenPdf = evenPdfCust.build()

    oddPdfCust.replaceArg( ws['AparMag2'], zeroCust )
    oddPdfCust.replaceArg( ws['A0Mag2'],   zeroCust )
    if pdfConfig['parameterizeKKMass'] :
        for bin in range( pdfBuild['KKMassBinning'].numBins() ) : oddPdfCust.replaceArg( ws[ 'f_S_bin%d' % bin ], zeroCust )
    else :
        oddPdfCust.replaceArg( ws['f_S'], zeroCust )
    oddPdf = oddPdfCust.build()

    SPdfCust.replaceArg( ws['AparMag2'],  zeroCust )
    SPdfCust.replaceArg( ws['A0Mag2'],    zeroCust )
    SPdfCust.replaceArg( ws['AperpMag2'], zeroCust )
    SPdf = SPdfCust.build()

    polSet = RooArgSet( AparMag2Cust, ws['A0Mag2'], ws['AperpMag2'], ws['AparPhase'], ws['AperpPhase'] )
    if pdfConfig['parameterizeKKMass'] :
        for bin in range( pdfBuild['KKMassBinning'].numBins() ) :
            polSet.add( ws[ 'f_S_bin%d' % bin ] )
            polSet.add( ws[ 'ASOddPhase_bin%d' % bin ] )
    else :
        polSet.add( ws['f_S'] )
        polSet.add( ws['ASOddPhase'] )

    print 'JvLFit: plot PDF polarization parameters:'
    print 'total:'
    pdf.getObservables(polSet).Print()
    print 'even:'
    evenPdf.getObservables(polSet).Print()
    print 'odd:'
    oddPdf.getObservables(polSet).Print()
    print 'S-wave:'
    SPdf.getObservables(polSet).Print()

    pdfNormSet = RooArgSet()
    for obs in pdf.Observables().difference( pdf.ConditionalObservables() ) : pdfNormSet.add(obs._var)
    pdfInt  = pdf.createIntegral(     pdfNormSet, RooArgSet() )
    evenInt = evenPdf.createIntegral( pdfNormSet, RooArgSet() )
    oddInt  = oddPdf.createIntegral(  pdfNormSet, RooArgSet() )
    SInt    = SPdf.createIntegral(    pdfNormSet, RooArgSet() )
    print 'JvLFit: plot PDF integrals:'
    pdfInt.Print()
    evenInt.Print()
    oddInt.Print()
    SInt.Print()

    if projWData :
        pdfNorm  = 0.
        evenNorm = 0.
        oddNorm  = 0.
        SNorm    = 0.
        obsSetPdf = pdf.getObservables( projWData['ProjWData'][0].get() )
        for obsSet in projWData['ProjWData'][0] :
            for obs in obsSet :
                if isinstance( obs, RooRealVar ) : obsSetPdf.setRealValue( obs.GetName(), obs.getVal()   )
                else                             : obsSetPdf.setCatIndex(  obs.GetName(), obs.getIndex() )
            pdfNorm  += projWData['ProjWData'][0].weight() * pdfInt.getVal()
            evenNorm += projWData['ProjWData'][0].weight() * evenInt.getVal()
            oddNorm  += projWData['ProjWData'][0].weight() * oddInt.getVal()
            SNorm    += projWData['ProjWData'][0].weight() * SInt.getVal()
        pdfNorm  /= projWData['ProjWData'][0].sumEntries()
        evenNorm /= projWData['ProjWData'][0].sumEntries()
        oddNorm  /= projWData['ProjWData'][0].sumEntries()
        SNorm    /= projWData['ProjWData'][0].sumEntries()

    else :
        pdfNorm  = pdfInt.getVal()
        evenNorm = evenInt.getVal()
        oddNorm  = oddInt.getVal()
        SNorm    = SInt.getVal()

    print 'JvLFit: PDF normalization factors:'
    print '    total:  %.3f (%.3f)' % ( pdfNorm, evenNorm + oddNorm + SNorm )
    print '    even:   %.3f' % evenNorm
    print '    odd:    %.3f' % oddNorm
    print '    S-wave: %.3f' % SNorm

    #pdfIntSetTagCat = RooArgSet(pdfNormSet)
    #pdfIntSetTagCat.remove(tagCatP2VVOS._var)
    #pdfIntTagCat = pdf.createIntegral( pdfIntSetTagCat, pdfNormSet )

    #pdfIntSetTimeTagCat = RooArgSet(pdfNormSet)
    #pdfIntSetTimeTagCat.remove(time._var)
    #pdfIntSetTimeTagCat.remove(tagCatP2VVOS._var)
    #pdfIntTimeTagCat = pdf.createIntegral( pdfIntSetTimeTagCat, pdfNormSet )

    #pdfIntSetTimeTagCatITag = RooArgSet(pdfIntSetTimeTagCat)
    #pdfIntSetTimeTagCatITag.remove(iTagOS._var)
    #pdfIntTimeTagCatITag = pdf.createIntegral( pdfIntSetTimeTagCatITag, pdfNormSet )

    #from math import pi
    #phiVal  = 0.5 * pi
    #lambVal = 0.9
    #pdfPars.find('phiCP').setVal(phiVal)
    #pdfPars.find('lambdaCP').setVal(lambVal)
    #pdfPars.find('avgCOddSum').setVal( -( 1. - lambVal**2 ) / ( 1. + lambVal**2 ) )
    #pdfPars.find('avgCOddOSTagged').setVal( -( 1. - lambVal**2 ) / ( 1. + lambVal**2 ) )
    #if pdfPars.find('avgCOddSSTagged') :
    #    pdfPars.find('avgCOddSSTagged').setVal( -( 1. - lambVal**2 ) / ( 1. + lambVal**2 ) )
    #    pdfPars.find('avgCOddTagged').setVal( -( 1. - lambVal**2 ) / ( 1. + lambVal**2 ) )
    #pdfPars.Print('v')

    #timeBBbarFrame = time.frame( Range = 'Bulk' )
    #timeBBbarFrame.GetYaxis().SetTitle('Decay time PDF (ps^{-1})')

    #tagCatP2VVOS.setIndex(1)
    #pdfIntTimeTagCat.plotOn( timeBBbarFrame, LineColor = kBlack, LineWidth = 4 )

    #iTagOS.setIndex(+1)
    #pdfIntTimeTagCatITag.plotOn( timeBBbarFrame, LineColor = kBlue, LineWidth = 4 )

    #iTagOS.setIndex(-1)
    #pdfIntTimeTagCatITag.plotOn( timeBBbarFrame, LineColor = kRed, LineWidth = 4 )

    #timeBBbarCanv = TCanvas( 'timeBBbarCanv', 'Decay time' )
    #timeBBbarFrame.Draw()
    #timeBBbarCanv.Print(plotsFile)

    ## plot lifetime
    #timeBBbarCanv = TCanvas( 'timeBBbarCanv', 'Lifetime' )
    #plot(  timeBBbarCanv, time, defData, pdf, addPDFs = [ pdf, pdf ], yScale = ( None, None )
    #     , frameOpts   = dict( Bins = 50, Title = '', Range = 'Bulk'                                                                )
    #     , dataOpts    = dict( MarkerStyle = kFullCircle, MarkerSize = 0.7, MarkerColor = kBlack, LineColor = kBlack, LineWidth = 3
    #                          , Cut = '%s > 0' % tagCatP2VVOS )
    #     , pdfOpts     = dict( list( projWDataBulk.items() ), LineColor = 13, LineWidth = 4, Slice = ( tagCatP2VVOS, 'TagCat1' ) )
    #     , addPDFsOpts = [  dict( list( projWDataBulk.items() ), LineColor = kBlue, LineWidth = 4
    #                             , Slices = [ ( tagCatP2VVOS, 'TagCat1' ), ( iTagOS, 'B'    ) ] )
    #                      , dict( list( projWDataBulk.items() ), LineColor = kRed,  LineWidth = 4
    #                             , Slices = [ ( tagCatP2VVOS, 'TagCat1' ), ( iTagOS, 'Bbar' ) ] )
    #                     ]
    #    )
    #LHCbText.Draw()
    #timeBBbarCanv.Print( plotsFile + '(' )

    #timeBCanv = TCanvas( 'timeBCanv', 'Lifetime' )
    #plot(  timeBCanv, time, defData, pdf, yScale = ( None, None )
    #     , frameOpts   = dict( Bins = 50, Title = '', Range = 'Bulk'                                                                )
    #     , dataOpts    = dict( MarkerStyle = kFullCircle, MarkerSize = 0.7, MarkerColor = kBlack, LineColor = kBlack, LineWidth = 3
    #                          , Cut = '%s > 0 && %s == +1' % ( tagCatP2VVOS, iTagOS ) )
    #     , pdfOpts     = dict( list( projWDataBulk.items() ), LineColor = kBlue, LineWidth = 4
    #                          , Slices = [ ( tagCatP2VVOS, 'TagCat1' ), ( iTagOS, 'B' ) ] )
    #    )
    #LHCbText.Draw()
    #timeBCanv.Print(plotsFile)

    #timeBbarCanv = TCanvas( 'timeBbarCanv', 'Lifetime' )
    #plot(  timeBbarCanv, time, defData, pdf, yScale = ( None, None )
    #     , frameOpts   = dict( Bins = 50, Title = '', Range = 'Bulk'                                                                )
    #     , dataOpts    = dict( MarkerStyle = kFullCircle, MarkerSize = 0.7, MarkerColor = kBlack, LineColor = kBlack, LineWidth = 3
    #                          , Cut = '%s > 0 && %s == -1' % ( tagCatP2VVOS, iTagOS ) )
    #     , pdfOpts     = dict( list( projWDataBulk.items() ), LineColor = kRed, LineWidth = 4
    #                          , Slices = [ ( tagCatP2VVOS, 'TagCat1' ), ( iTagOS, 'Bbar' ) ] )
    #    )
    #LHCbText.Draw()
    #timeBbarCanv.Print( plotsFile + ')' )

    # plot lifetime and angles
    print 'JvLFit: plotting time and angular distributions'
    timeCanv = TCanvas( 'timeCanv', 'Decay time'   )
    ctkCanv  = TCanvas( 'ctkCanv',  'cos(theta_K)' )
    ctlCanv  = TCanvas( 'ctlCanv',  'cos(theta_l)' )
    phiCanv  = TCanvas( 'phiCanv',  'cos(phi)'     )
    for ( pad, obs, nBins, xTitle, yScale, yScaleRel, logY, text, coords )\
            in zip(  [ timeCanv, ctkCanv, ctlCanv, phiCanv ]
                   , obsSetP2VV[ : 4 ]
                   , numBins
                   , ( '', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                   , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                   , 2 * ( ( None, None ), ) + 2 * ( ( None, 17. / 14. ), )
                   , ( True, ) + 3 * ( False, )
                   , [ LHCbTextMore ] + 3 * [ LHCbText ]
                   , [ ( 0.47, 0.73, 0.92, 0.93 ) ] + [ ( 0.26, 0.82, 0.39, 0.91 ) ] + 2 * [ ( 0.19, 0.82, 0.32, 0.91 ) ]
                  ) :
        plot(  pad, obs, defData, pdf, addPDFs = [ evenPdf, oddPdf, SPdf ]
             , xTitle = xTitle, yScale = yScale, yScaleRel = yScaleRel, logy = logY
             , frameOpts   = dict( Bins = nBins, Title = ''                                                                )
             , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize                                          )
             , pdfOpts     = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, LineStyle = kSolid )
             , addPDFsOpts = [  dict( list( projWData.items() ), LineColor = kRed, LineWidth = lineWidth, LineStyle = 9
                                     , Normalization = evenNorm / pdfNorm )
                              , dict( list( projWData.items() ), LineColor = kGreen + 3, LineWidth = lineWidth, LineStyle = 7
                                     , Normalization = oddNorm / pdfNorm )
                              , dict( list( projWData.items() ), LineColor = kMagenta + 3, LineWidth = lineWidth, LineStyle = 5
                                     , Normalization = SNorm / pdfNorm )
                             ]
             , components = comps
            )
        _P2VVPlotStash.append( text.Clone() )
        _P2VVPlotStash[-1].SetX1NDC( coords[0] )
        _P2VVPlotStash[-1].SetY1NDC( coords[1] )
        _P2VVPlotStash[-1].SetX2NDC( coords[2] )
        _P2VVPlotStash[-1].SetY2NDC( coords[3] )
        _P2VVPlotStash[-1].Draw()

    ## plot lifetime
    #timePlotTitles = tuple( [ time.GetTitle() + title for title in (  ' - linear'
    #                                                                , ' - logarithmic'
    #                                                                , ' - B/#bar{B} asymmetry'
    #                                                               )
    #                        ] )
    #timeCanv1 = TCanvas( 'timeCanv1', 'Lifetime' )
    #print 'JvLFit: plotting lifetime distribution'
    #for ( pad, nBins, plotTitle, yTitle, yScale, dataCuts, pdfCuts, logY )\
    #        in zip(  timeCanv1.pads( 2, 2 )
    #               , 3 * [ pdfConfig['numTimeBins'] ]
    #               , timePlotTitles
    #               , 2 * ( None, ) + ( 'B/#bar{B} asymmetry', )
    #               , ( ( None, None ), ( 50., None ), ( None, None ) )
    #               , 2 * ( dict(), ) + ( dict( Asymmetry = iTagOS ), )
    #               , 2 * ( dict(), ) + ( dict( Asymmetry = iTagOS ), )
    #               , ( False, True, False )
    #              ) :
    #    plot(  pad, time, defData, pdf, yTitle = yTitle, yScale = yScale, logy = logY
    #         , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                    )
    #         , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                         )
    #         , pdfOpts    = dict( list( projWDataBulk.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
    #         , components = comps
    #        )

    ## plot lifetime (tagged/untagged)
    #print 'JvLFit: plotting lifetime distributions tagged/untagged'
    #timePlotTitles1 = tuple( [ time.GetTitle() + title for title in (  ' - untagged'
    #                                                                 , ' - tagging category 2'
    #                                                                 , ' - tagging category %d' % tagCat5Min
    #                                                                 , ' - B/#bar{B} asymmetry untagged'
    #                                                                 , ' - B/#bar{B} asymmetry tagging category 2'
    #                                                                 , ' - B/#bar{B} asymmetry tagging category %d' % tagCat5Min
    #                                                                )
    #                        ] )
    #timeCanv2 = TCanvas( 'timeCanv2', 'Lifetime' )
    #for ( pad, nBins, plotTitle, yTitle, dataCuts, pdfCuts, logY )\
    #    in zip(  timeCanv2.pads( 3, 2 )
    #         , 6 * [ pdfConfig['numTimeBins'] ]
    #         , timePlotTitles1
    #         , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
    #         ,   ( dict( Cut = '%s == 0'  % ( tagCatP2VVOS.GetName()             )                     ), )
    #           + ( dict( Cut = '%s == 2'  % ( tagCatP2VVOS.GetName()             )                     ), )
    #           + ( dict( Cut = '%s == %d' % ( tagCatP2VVOS.GetName(), tagCat5Min )                     ), )
    #           + ( dict( Cut = '%s == 0'  % ( tagCatP2VVOS.GetName()             ), Asymmetry = iTagOS ), )
    #           + ( dict( Cut = '%s == 2'  % ( tagCatP2VVOS.GetName()             ), Asymmetry = iTagOS ), )
    #           + ( dict( Cut = '%s == %d' % ( tagCatP2VVOS.GetName(), tagCat5Min ), Asymmetry = iTagOS ), )
    #         ,   ( dict( Slice = ( tagCatP2VVOS, 'Untagged'              )                     ), )
    #           + ( dict( Slice = ( tagCatP2VVOS, 'TagCat2'               )                     ), )
    #           + ( dict( Slice = ( tagCatP2VVOS, 'TagCat%d' % tagCat5Min )                     ), )
    #           + ( dict( Slice = ( tagCatP2VVOS, 'Untagged'              ), Asymmetry = iTagOS ), )
    #           + ( dict( Slice = ( tagCatP2VVOS, 'TagCat2'               ), Asymmetry = iTagOS ), )
    #           + ( dict( Slice = ( tagCatP2VVOS, 'TagCat%d' % tagCat5Min ), Asymmetry = iTagOS ), )
    #         , 3 * ( False, ) + 3 * ( False, )
    #        ) :
    #    plot(  pad, time, defData, pdf, yTitle = yTitle, logy = logY
    #         , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                )
    #         , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                     )
    #         , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
    #         , components = comps
    #        )

    ## plot angles
    #print 'JvLFit: plotting angular distributions'
    #if plotAnglesNoEff and pdfConfig['SFit'] and pdfConfig['multiplyByTimeEff'] not in [ 'all', 'signal' ]\
    #        and not pdfConfig['conditionalTagging'] :
    #    addPDFs = [ ws['sig_t_angles_tagCat_iTag'] ]
    #else :
    #    addPDFs = [ ]

    #anglePlotTitles =   tuple(  [ angle.GetTitle()                            for angle in angles ]\
    #                          + [ angle.GetTitle() + ' - B/#bar{B} asymmetry' for angle in angles ] )
    #anglesCanv = TCanvas( 'anglesCanv', 'Decay Angles' )
    #for ( pad, obs, nBins, plotTitle, xTitle, yTitle, dataCuts, pdfCuts )\
    #        in zip(  anglesCanv.pads( 3, 2 )
    #               , 2 * angles
    #               , 2 * numAngleBins
    #               , anglePlotTitles
    #               , 2 * ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
    #               , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
    #               , 3 * ( dict( ), ) + 3 * ( dict( Asymmetry = iTagOS ), )
    #               , 3 * ( dict( ), ) + 3 * ( dict( Asymmetry = iTagOS ), )
    #              ) :
    #    plot(  pad, obs, defData, pdf, addPDFs = addPDFs, xTitle = xTitle, yTitle = yTitle
    #         , frameOpts   = dict( Bins = nBins, Title = plotTitle                                                )
    #         , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts                    )
    #         , pdfOpts     = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
    #         , addPDFsOpts = [ dict( list( projWData.items() ), LineColor = kRed, LineWidth = lineWidth, **pdfCuts ) ]
    #         , components  = comps
    #        )

    #if not pdfConfig['SFit'] and pdfConfig['SWeightsType'].startswith('simultaneous')\
    #        and pdfConfig['parameterizeKKMass'] == 'simultaneous' :
    #    # plot signal mass
    #    print 'JvLFit: plotting mumuKK mass distribution'
    #    pad = pdfBuild['massCanv'].cd(2)
    #    plot(  pad, BMass, defData, pdf
    #         , frameOpts  = dict( Range = 'Signal', Bins = pdfConfig['numBMassBins'][0], Title = BMass.GetTitle() + ' full fit - signal' )
    #         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                                                      )
    #         , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = 2                                            )
    #         , components = comps
    #        )

    # print canvas to file
    timeCanv.Print( plotsFile + '(' )
    ctkCanv.Print(plotsFile)
    ctlCanv.Print(plotsFile)
    phiCanv.Print( plotsFile + ( ')' if not deltaSCanv and not pdfConfig['makePlots'] else '' ) )
    #timeCanv1.Print(plotsFile)
    #timeCanv2.Print(plotsFile)
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
    gStyle.SetEndErrorSize(4)
    deltaSCanv.Update()
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


tagData = pdfBuild['sigSWeightData']

sums = {
    'numEv'       : 0.
  , 'numOS'       : 0., 'numOSExcl'       : 0.
  , 'numSS'       : 0., 'numSSExcl'       : 0.
  , 'numComb'     : 0., 'numCombExcl'     : 0.
  , 'etaOS'       : 0., 'etaOSExcl'       : 0.
  , 'etaSS'       : 0., 'etaSSExcl'       : 0.
  , 'etaComb'     : 0., 'etaCombExcl'     : 0.
  , 'wOS'         : 0., 'wOSExcl'         : 0.
  , 'wSS'         : 0., 'wSSExcl'         : 0.
  , 'wComb'       : 0., 'wCombExcl'       : 0.
  , 'dilOS'       : 0., 'dilOSExcl'       : 0.
  , 'dilSS'       : 0., 'dilSSExcl'       : 0.
  , 'dilComb'     : 0., 'dilCombExcl'     : 0.
  , 'dil2OS'      : 0., 'dil2OSExcl'      : 0.
  , 'dil2SS'      : 0., 'dil2SSExcl'      : 0.
  , 'dil2Comb'    : 0., 'dil2CombExcl'    : 0.
  , 'dilRes2Comb' : 0., 'dilRes2CombExcl' : 0.
  , 'dilTot2Comb' : 0., 'dilTot2CombExcl' : 0.
  , 'dilResAll'   : 0.
  , 'dilRes2All'  : 0.
}

from math import exp
for varSet in tagData :
  weight = tagData.weight()
  sums['numEv'] += weight

  sigmat = 1.45 * varSet.getRealValue('sigmat')
  dilResSq = exp( - 17.6**2 * sigmat**2 )
  dilRes   = sqrt(dilResSq)
  sums['dilResAll']  += weight * dilRes
  sums['dilRes2All'] += weight * dilResSq

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

      sums['numComb']     += weight
      sums['etaComb']     += weight * etaOS
      sums['wComb']       += weight * wTagOS
      sums['dilComb']     += weight * dilOS
      sums['dil2Comb']    += weight * dilOS**2
      sums['dilRes2Comb'] += weight * dilResSq
      sums['dilTot2Comb'] += weight * dilResSq * dilOS**2

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

      sums['numComb']     += weight
      sums['etaComb']     += weight * etaSS
      sums['wComb']       += weight * wTagSS
      sums['dilComb']     += weight * dilSS
      sums['dil2Comb']    += weight * dilSS**2
      sums['dilRes2Comb'] += weight * dilResSq
      sums['dilTot2Comb'] += weight * dilResSq * dilSS**2

  if varSet.getCatIndex('tagdecision_os') != 0 and varSet.getCatIndex('tagdecision_ss') != 0 :
    dilSign = +1. if varSet.getCatIndex('tagdecision_os') == varSet.getCatIndex('tagdecision_ss') else -1.
    dilComb = ( dilOS + dilSign * dilSS ) / ( 1. + dilSign * dilOS * dilSS )
    wTagComb = ( 1. - dilComb ) / 2.

    sums['numComb']     += weight
    sums['etaComb']     += weight * wTagComb
    sums['wComb']       += weight * wTagComb
    sums['dilComb']     += weight * dilComb
    sums['dil2Comb']    += weight * dilComb**2
    sums['dilRes2Comb'] += weight * dilResSq
    sums['dilTot2Comb'] += weight * dilResSq * dilComb**2

    sums['numCombExcl']     += weight
    sums['etaCombExcl']     += weight * wTagComb
    sums['wCombExcl']       += weight * wTagComb
    sums['dilCombExcl']     += weight * dilComb
    sums['dil2CombExcl']    += weight * dilComb**2
    sums['dilRes2CombExcl'] += weight * dilResSq
    sums['dilTot2CombExcl'] += weight * dilResSq * dilComb**2

print
print 'number of events:       %.4f' % sums['numEv']
print 'number of OS events:    %.4f (%.4f)'   % ( sums['numOS'],   sums['numOSExcl']   )
print 'number of SS events:    %.4f (%.4f)'   % ( sums['numSS'],   sums['numSSExcl']   )
print 'number of Comb. events: %.4f (%.4f)\n' % ( sums['numComb'], sums['numCombExcl'] )

print 'OS    eff.: %.2f%% (%.2f%%)'   % ( sums['numOS']   / sums['numEv'] * 100., sums['numOSExcl']   / sums['numEv'] * 100. )
print 'SS    eff.: %.2f%% (%.2f%%)'   % ( sums['numSS']   / sums['numEv'] * 100., sums['numSSExcl']   / sums['numEv'] * 100. )
print 'Comb. eff.: %.2f%% (%.2f%%)\n' % ( sums['numComb'] / sums['numEv'] * 100., sums['numCombExcl'] / sums['numEv'] * 100. )

print 'OS    <eta>: %.4f (%.4f)'   % ( sums['etaOS']   / sums['numOS'],   sums['etaOSExcl']   / sums['numOSExcl']   )
print 'SS    <eta>: %.4f (%.4f)'   % ( sums['etaSS']   / sums['numSS'],   sums['etaSSExcl']   / sums['numSSExcl']   )
print 'Comb. <eta>: %.4f (%.4f)\n' % ( sums['etaComb'] / sums['numComb'], sums['etaCombExcl'] / sums['numCombExcl'] )

print 'OS    <w>: %.4f (%.4f)'   % ( sums['wOS']   / sums['numOS'],   sums['wOSExcl']   / sums['numOSExcl']   )
print 'SS    <w>: %.4f (%.4f)'   % ( sums['wSS']   / sums['numSS'],   sums['wSSExcl']   / sums['numSSExcl']   )
print 'Comb. <w>: %.4f (%.4f)\n' % ( sums['wComb'] / sums['numComb'], sums['wCombExcl'] / sums['numCombExcl'] )

print 'OS    <dil>:    %.4f (%.4f)'   % ( sums['dilOS']   / sums['numOS'],   sums['dilOSExcl']   / sums['numOSExcl']   )
print 'SS    <dil>:    %.4f (%.4f)'   % ( sums['dilSS']   / sums['numSS'],   sums['dilSSExcl']   / sums['numSSExcl']   )
print 'Comb. <dil>:    %.4f (%.4f)'   % ( sums['dilComb'] / sums['numComb'], sums['dilCombExcl'] / sums['numCombExcl'] )
print 'All   <dilRes>: %.4f\n'        % ( sums['dilResAll'] / sums['numEv'] )

print 'OS    <dil2>:    %.4f (%.4f)'   % ( sums['dil2OS']   / sums['numOS'],   sums['dil2OSExcl']   / sums['numOSExcl']   )
print 'SS    <dil2>:    %.4f (%.4f)'   % ( sums['dil2SS']   / sums['numSS'],   sums['dil2SSExcl']   / sums['numSSExcl']   )
print 'Comb. <dil2>:    %.4f (%.4f)'   % ( sums['dil2Comb'] / sums['numComb'], sums['dil2CombExcl'] / sums['numCombExcl'] )
print 'Comb. <dilRes2>: %.4f (%.4f)'   % ( sums['dilRes2Comb'] / sums['numComb'], sums['dilRes2CombExcl'] / sums['numCombExcl'] )
print 'Comb. <dilTot2>: %.4f (%.4f)'   % ( sums['dilTot2Comb'] / sums['numComb'], sums['dilTot2CombExcl'] / sums['numCombExcl'] )
print 'All   <dilRes2>: %.4f\n'        % ( sums['dilRes2All'] / sums['numEv'] )

print 'OS    <eff * dil2>:    %.4f%% (%.4f%%)'   % ( sums['dil2OS']   / sums['numEv'] * 100., sums['dil2OSExcl']   / sums['numEv'] * 100. )
print 'SS    <eff * dil2>:    %.4f%% (%.4f%%)'   % ( sums['dil2SS']   / sums['numEv'] * 100., sums['dil2SSExcl']   / sums['numEv'] * 100. )
print 'Comb. <eff * dil2>:    %.4f%% (%.4f%%)'   % ( sums['dil2Comb'] / sums['numEv'] * 100., sums['dil2CombExcl'] / sums['numEv'] * 100. )
print 'Comb. <eff * dilRes2>: %.4f%% (%.4f%%)'   % (  sums['dilRes2Comb'] / sums['numEv'] * 100.
                                                    , sums['dilRes2CombExcl'] / sums['numEv'] * 100.
                                                   )
print 'Comb. <eff * dilTot2>: %.4f%% (%.4f%%)'   % (  sums['dilTot2Comb'] / sums['numEv'] * 100.
                                                    , sums['dilTot2CombExcl'] / sums['numEv'] * 100.
                                                   )
print 'All   <eff * dilRes2>: %.4f\n'            % ( sums['dilRes2All'] / sums['numEv'] )

nEv = 0.
nBB = 0.
nBbarBbar = 0.
nBbarB = 0.
nBBbar = 0.
avD_OS = 0.
avD_SS = 0.
avDD = 0.
for argSet in tagData :
    if argSet.getCatIndex('tagdecision_os') == 0 or argSet.getCatIndex('tagdecision_ss') == 0 : continue

    nEv += tagData.weight()
    if argSet.getCatIndex('tagdecision_os') == +1 and argSet.getCatIndex('tagdecision_ss') == +1 : nBB += tagData.weight()
    if argSet.getCatIndex('tagdecision_os') == -1 and argSet.getCatIndex('tagdecision_ss') == -1 : nBbarBbar += tagData.weight()
    if argSet.getCatIndex('tagdecision_os') == -1 and argSet.getCatIndex('tagdecision_ss') == +1 : nBbarB += tagData.weight()
    if argSet.getCatIndex('tagdecision_os') == +1 and argSet.getCatIndex('tagdecision_ss') == -1 : nBBbar += tagData.weight()

    D_OS = 1. - 2. * argSet.getRealValue('tagomega_os')
    D_SS = 1. - 2. * argSet.getRealValue('tagomega_ss')

    avD_OS += D_OS * tagData.weight()
    avD_SS += D_SS * tagData.weight()
    avDD   += D_OS * D_SS * tagData.weight()

avD_OS /= nEv
avD_SS /= nEv
avDD   /= nEv

Atags = ( nBB + nBbarBbar - nBbarB - nBBbar ) / nEv
print avDD, avD_OS * avD_SS, Atags