



###########################################################################################################################################
# Naive fit of the KK mass peak, (Intended to make a nice looking plot of the KK mass.)
###########################################################################################################################################







###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
pdfConfig['dataSample'] = '' # ( None, 100260, '' )  # '' / 'Summer2011' / 'runNumber % 2 == 1'
pdfConfig['selection']  = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
generateData            = False
doFit                   = False
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = False
pdfConfig['blind']      = False
pdfConfig['nominalPdf'] = False  # nominal PDF option does not work at the moment
corrSFitErr             = 'sumWeight' # [ 1., 0.700, 0.952, 0.938, 0.764 ] # '' / 'matrix' / 'sumWeight'
randomParVals           = ( ) # ( 1., 12346 ) # ( 2., 12345 )

plotsFile = 'plots/paper2012_SFit.ps'
#plotsFile = 'plots/JvLSFit.ps' if pdfConfig['SFit']\
#       else 'plots/JvLCFit.ps'
parameterFile = None # 'JvLSFit.par' if pdfConfig['SFit'] else 'JvLCFit.par'

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
    dataSetFile = 'JvLSFit.root' if pdfConfig['SFit'] else 'JvLCFit.root'

MinosPars = [  #'AparPhase'
             #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
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
from P2VVLoad import RooFitOutput, LHCbStyle
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
pdfConfig['angularRanges']      = dict( cpsi = ( -1., +1. ), ctheta = ( -1., +1. ), phi = ( -pi, +pi ) )

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
pdfConfig['KKMassBinBounds']      = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['SWaveAmplitudeValues'] = (  [ (0.33, 0.09), (0.073, 0.030), (0.009, 0.012), (0.012, 0.010), (0.061, 0.027), (0.18, 0.04) ]
                                     , [ (1.1,  0.5 ), (0.7,   0.2  ), (0.4,   0.4  ), (-0.6,  0.3  ), (-0.4, 0.2   ), (-0.7, 0.2 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.34, 0.11),  (0.057, 0.031), (0.003, 0.006), (0.032, 0.016), (0.069, 0.032), (0.16,  0.05) ]
#                                     , [ (0.60, 0.33 ), (0.81,  0.34 ), (0.8,   1.3  ), (-0.48, 0.19 ), (-0.47, 0.20 ), (-0.61, 0.21) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.28, 0.11), (0.06, 0.02), (0.04, 0.02), (0.27, 0.07) ]
#                                     , [ (2.7,  0.4 ), (0.22,   0.14), (-0.11, 0.17 ), (-0.97, 0.3 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.28, 0.11), (0.05, 0.02), (0.27, 0.07) ]
#                                     , [ (2.7,  0.4 ), (0.,   0.15), (-0.97, 0.3 ) ] )
pdfConfig['CSPValues']            = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.959, 0.498, 0.968 ]

pdfConfig['sameSideTagging']    = True  # nominal: True
pdfConfig['conditionalTagging'] = True  # nominal: True
pdfConfig['continuousEstWTag']  = True  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 50
pdfConfig['constrainTagging']   = 'constrain'  # nominal: 'constrain'

pdfConfig['timeResType']           = 'eventNoMean' # 'event' # 'eventNoMean'
pdfConfig['numTimeResBins']        = 50
pdfConfig['constrainTimeResScale'] = 'fixed'  # nominal: 'constrain'

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
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_newTrigger_simpleForWeights'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_newTrigger_simpleForCorrWeights'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'
#                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_PHSP_Basis'

if not pdfConfig['nominalPdf'] and pdfConfig['transversityAngles'] :
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
## get data ##
###################
#This is the unweighted dataset, after full selection and trigger
defData = pdfBuild['data']




assert(False)

##################################################
## Make the B mass plots
##################################################
from P2VVGeneralUtils import plot
from ROOT import TCanvas
from P2VVLoad import LHCbStyle

# JpsiPhi Pdf
JpsiPhiPdf = pdfBuild._massPdf

#Pdf Components
comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
         , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
          }

##!!!!!!!!!!!!!! TODO: MAKE PLOT STYLE IDENTICAL TO THAT OF THE PAPER AND APLY THE INSTRUCTION FOR THE PAPER LAYOUT
plot(  TCanvas('JpsiPhiMass'), BMass, defData, JpsiPhiPdf, components = comps)


assert(False)




##################################################
## Make the mumu mass plots
##################################################
from P2VVParameterizations.MassPDFs import CB_Signal_Mass, Linear_Background_Mass

mumuSig = CB_Signal_Mass         ( dict( Name = 'sig_mumu', mass = mumuMass ))
mumuBkg = Linear_Background_Mass ( dict( Name = 'bkg_mumu', mass = mumuMass ))

from RooFitWrappers import buildPdf
mumuMassPdf = buildPdf( [ mumuSig, mumuBkg ], Observables = [ mumuMass ], Name = 'mumuMass' )



##!!!!!!!!!!! defData DO NOT CONTAIN mdau1!!!!!!!!!!!!!

assert(False)






##################################################
## Make the KK mass plots
##################################################
#Build KK mas pdf
from ROOT import RooRelBreitWigner, RooConstVar, RooFFTConvPdf, RooGaussModel

#Observable
#KKmass = defData.get().find('mdau2')

#Parameters
mean   = RooRealVar ( 'mean'  , 'mean'  , 1020, 1015, 1020 ) 
width  = RooRealVar ( 'width' , 'width' , 2   ,    1,    7 ) 
spin   = RooConstVar( 'spin'  , 'spin'  , 1                )  
radius = RooRealVar ( 'radius', 'radius', 1   ,    0,    5 )
K1mass = RooConstVar( 'K1mass', 'K1mass', 493.7            )
K2mass = RooConstVar( 'K2mass', 'K2mass', 493.7            )


#Phi Signal Mass pdf 
PhiMassPdf  = RooRelBreitWigner('PhiMassPdf', 'PhiMassPdf', KKmass, mean, width, spin, radius, K1mass, K2mass )
 #Resolution Model
resMean    = RooRealVar   ( 'resMean'   , 'resMean'   , 0, -2, 2  ) 
resSigma   = RooRealVar   ( 'resSigma'  , 'resSigma'  , 2,  1, 10 ) 
GaussModel = RooGaussModel('GaussModel', 'GaussModel' , KKmass, resMean,resSigma )
 #Phi Mass Signal Total Pdf 
PhiMassConvPdf = RooFFTConvPdf('PhiMassConvPdf', 'PhiMassConvPdf', KKmass, PhiMassPdf, GaussModel )


##KK Bkg mass pdf



assert(False)



from ROOT import TCanvas
c2 = TCanvas()

fr = KKmass.frame()

PhiMassConvPdf.plotOn(fr)
fr.Draw()

data.plotOn(fr)


data = defData.reduce(ArgSet = [KKmass] )
