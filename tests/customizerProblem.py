###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
pdfConfig['dataSample'] = '' # ( None, 100260, '' )  # '' / 'Summer2011' / 'runNumber % 2 == 1'
pdfConfig['selection']  = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
generateData            = False
doFit                   = True
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['sFit']       = True
pdfConfig['blind']      = False
corrSFitErr             = 'sumWeight' # [ 1., 0.700, 0.952, 0.938, 0.764 ] # '' / 'matrix' / 'sumWeight'
randomParVals           = ( ) # ( 1., 12346 ) # ( 2., 12345 )

plotsFile = 'plots/paper2012_SFit.ps'
#plotsFile = 'plots/JvLSFit.ps' if pdfConfig['sFit']\
#       else 'plots/JvLCFit.ps'
parameterFile = None # 'JvLSFit.par' if pdfConfig['sFit'] else 'JvLCFit.par'

if readData :
    pdfConfig['nTupleName'] = 'DecayTree'
    pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDown.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagUp.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_rand0.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_rand1.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_BTags.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_BbarTags.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Pass3-version2_Bs_050711_nocut_Phi_P2VV.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20121010.root'
    #pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_DGs0_MC11a_ntupleB_for_fitting_20121119.root'
    pdfConfig['nominalDataSet'] = False
else :
    pdfConfig['nTupleName'] = None
    pdfConfig['nTupleFile'] = None

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLSFit.root' if pdfConfig['sFit'] else 'JvLCFit.root'

MinosPars = [  #'AparPhase'
             #, 'ASOddPhase_bin0', 'ASOddPhase_bin1',                    'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
             #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
             #, 'f_S_bin0',        'f_S_bin1',                           'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
             #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
            ]
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
pdfConfig['angularRanges']      = dict( cpsi = [ ( -1., +1. ) ], ctheta = [ ( -1., +1. ) ], phi = [ ( -pi, +pi ) ] )

pdfConfig['sigMassModel']         = '' # 'boxFixed'
pdfConfig['bkgMassModel']         = '' # 'linearConstant'
pdfConfig['bkgAnglePdf']          = 'hybrid'  # default/nominal: ''
pdfConfig['sigTaggingPdf']        = 'tagUntag'  # default: 'tagUntag' | nominal: 'tagCats'
pdfConfig['bkgTaggingPdf']        = 'tagUntagRelative'  # default: 'tagUntagRelative' | 'tagCatsRelative'
pdfConfig['multiplyByTagPdf']     = False
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['timeEffType']          = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['multiplyByAngEff']     = 'basis012'  # default: 'basis012'
pdfConfig['parameterizeKKMass']   = 'simultaneous'  # default/nominal: 'simultaneous'
pdfConfig['ambiguityParameters']  = False
pdfConfig['lifetimeRange']        = ( 0.3, 14. )
pdfConfig['SWeightsType']         = 'simultaneousFreeBkg'  # default/nominal: 'simultaneousFreeBkg'
pdfConfig['KKMassBinBounds']      = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['SWaveAmplitudeValues'] = (  [ (0.33, 0.09), (0.073, 0.030), (0.009, 0.012), (0.012, 0.010), (0.061, 0.027), (0.18, 0.04) ]
                                     , [ (1.1,  0.5 ), (0.7,   0.2  ), (0.4,   0.4  ), (-0.6,  0.3  ), (-0.4, 0.2   ), (-0.7, 0.2 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.23, 0.08), (0.069, 0.029), (0.009, 0.013), (0.016, 0.011), (0.056, 0.027), (0.17, 0.04) ]
#                                     , [ (1.3,  0.6 ), (0.8,   0.3  ), (0.5,   0.5  ), (-0.5,  0.3  ), (-0.5, 0.2   ), (-0.7, 0.2 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.23, 0.08), (0.068, 0.029), (0.013, 0.008), (0.013, 0.008), (0.055, 0.026), (0.17,  0.04) ]
#                                     , [ (1.3,  0.7 ), (0.76,   0.27), (0.39,  0.25 ), (-0.57, 0.28 ), (-0.46, 0.21 ), (-0.64, 0.20) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.22, 0.07),  (0.057, 0.028), (0.003, 0.007), (0.012, 0.009), (0.046, 0.026), (0.16,  0.04) ]
#                                     , [ (1.43, 0.85 ), (0.83,  0.34 ), (0.8,   1.4  ), (-0.62, 0.33 ), (-0.52, 0.25 ), (-0.67, 0.21) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.22, 0.08), (0.058, 0.028), (0.009, 0.006), (0.009, 0.006), (0.047, 0.026), (0.16,  0.04) ]
#                                     , [ (1.4,  0.8 ), (0.81,  0.32 ), (0.42,  0.30 ), (-0.72, 0.38 ), (-0.51, 0.24 ), (-0.66, 0.20) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.28, 0.11), (0.06, 0.02), (0.04, 0.02), (0.27, 0.07) ]
#                                     , [ (2.7,  0.4 ), (0.22,   0.14), (-0.11, 0.17 ), (-0.97, 0.3 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.22, 0.07), (0.017, 0.013), (0.020, 0.010), (0.16,  0.04) ]
#                                     , [ (1.4,  0.9 ), (0.7,   0.4),   (-0.64, 0.25),  (-0.66, 0.20) ] )
pdfConfig['CSPValues']            = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.964, 0.770, 0.824, 0.968 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.959, 0.498, 0.968 ]

pdfConfig['sameSideTagging']    = True  # nominal: True
pdfConfig['conditionalTagging'] = True  # nominal: True
pdfConfig['continuousEstWTag']  = True  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 50
pdfConfig['constrainTagging']   = 'constrain'  # nominal: 'constrain'

pdfConfig['timeResType']           = 'eventNoMean' # 'event' # 'eventNoMean'
pdfConfig['numTimeResBins']        = 25
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

pdfConfig['timeEffHistFiles'] = dict(  file      = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
                                     , hlt1UB    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                     , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                    )
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504_unitAcceptance.root'
#pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
#pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['angEffMomsFiles']    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'\
                                    if pdfConfig['transversityAngles'] else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_newTrigger'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_newTrigger_simpleForPHSPWeights'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_newTrigger_simpleForCorrPHSPWeights'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_PHSP_Basis'

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

numBins = ( 50, 20, 20, 20 )
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

if not pdfConfig['sFit'] : obsSetP2VV.append(BMass)

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
    nEvents = int( pdfConfig['numEvents'] * ( pdfConfig['signalFraction'] if pdfConfig['sFit'] else 1. ) )
    print 'JvLFit: generating %d events' % nEvents
    fitData = pdf.generate( obsSetP2VV, nEvents )

    # additional observables
    if not pdfConfig['transversityAngles'] :
        from P2VV.Utilities.DataHandling import addTransversityAngles
        addTransversityAngles( fitData, 'trcospsi',          'trcostheta',        'trphi'
                                      , angles[0].GetName(), angles[1].GetName(), angles[2].GetName() )

    # write data to file
    from P2VV.Utilities.DataHandling import writeData
    writeData( dataSetFile, dataSetName, fitData )

elif pdfConfig['sFit'] :
    defData = pdfBuild['sigSWeightData']
    sigData = pdfBuild['sigSWeightData']
    bkgData = pdfBuild['bkgSWeightData']
    if corrSFitErr == 'sumWeight'\
            or ( type(corrSFitErr) != str and hasattr( corrSFitErr, '__iter__' ) and hasattr( corrSFitErr, '__getitem__' ) ) :
        from P2VV.Utilities.DataHandling import correctWeights
        fitData = correctWeights( pdfBuild['sigSWeightData'], 'KKMassCat' if pdfConfig['parameterizeKKMass'] == 'simultaneous' else ''
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

# print parameters
print 120 * '='
print 'JvLFit: fit data:'
fitData.Print()
print 'JvLFit: observables in PDF:'
pdfObs.Print('v')
print 'JvLFit: parameters in PDF:'
pdfPars.Print('v')

###########################################################################################################################################
## customize PDF ##
###################

from P2VV.RooFitWrappers import RealVar, Customizer
newTime  = RealVar( 'newTime',  Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. ) )
newTime1 = RealVar( 'newTime1', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. ) )
custPdf  = Customizer( Pdf = pdf, OriginalArgs = [ time ], SubstituteArgs = [ newTime ],  NameSuffix = 'newTime'  )

print '\noriginal PDF:'
pdf.Print()

print '\ncustomized PDF:'
custPdf.Print()

from ROOT import RooArgSet
intSet       = RooArgSet()
normSetOrig  = RooArgSet( ws['time'],     ws['helcosthetaK'], ws['helcosthetaL'], ws['helphi'] )
normSetCust  = RooArgSet( ws['newTime'],  ws['helcosthetaK'], ws['helcosthetaL'], ws['helphi'] )

pdfInt   = pdf.createIntegral(      intSet, normSetOrig  )
custInt  = custPdf.createIntegral(  intSet, normSetCust  )

print '\noriginal PDF integral:'
pdfInt.Print()

print '\ncustomized PDF integral:'
custInt.Print()

custPdf1 = Customizer( Pdf = pdf, OriginalArgs = [ time ], SubstituteArgs = [ newTime1 ], NameSuffix = 'newTime1' )
print '\nanother customized PDF:'
custPdf1.Print()

normSetCust1 = RooArgSet( ws['newTime1'], ws['helcosthetaK'], ws['helcosthetaL'], ws['helphi'] )
custInt1 = custPdf1.createIntegral( intSet, normSetCust1 )
print '\nanother customized PDF integral:'
custInt1.Print()

print '\ndo the time observables have the original ranges?'
print 'time',     ws['time'].hasRange('R0_time_Paper2012_FitAcceptance_hlt1_excl_biased_dec_notExclB_shape_default_I_time')
print 'newTime',  ws['newTime'].hasRange('R0_time_Paper2012_FitAcceptance_hlt1_excl_biased_dec_notExclB_shape_default_I_time')
print 'newTime1', ws['newTime1'].hasRange('R0_time_Paper2012_FitAcceptance_hlt1_excl_biased_dec_notExclB_shape_default_I_time')

print '\ndo the time observables have their own ranges?'
print 'time',     ws['time'].hasRange('R0_time_Paper2012_FitAcceptance_hlt1_excl_biased_dec_notExclB_shape_default_I_time')
print 'newTime',  ws['newTime'].hasRange('R0_newTime_Paper2012_FitAcceptance_hlt1_excl_biased_dec_notExclB_shape_newTime_default_I_newTime')
print 'newTime1', ws['newTime1'].hasRange('R0_newTime1_Paper2012_FitAcceptance_hlt1_excl_biased_dec_notExclB_shape_newTime1_default_I_newTime1')
