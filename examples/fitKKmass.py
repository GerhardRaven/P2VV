
###########################################################################################################################################
## This script makes some plots for the paper             ##
## THe plots are: Bmass, mumuMass, Flavour Tagging plots  ##
###########################



###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                  = True
pdfConfig['dataSample']   = '' # ( None, 100260, '' )  # '' / 'Summer2011' / 'runNumber % 2 == 1'
pdfConfig['selection']    = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
generateData              = False
doFit                     = False
makeObservablePlots       = False
makeKKMassPlots           = False
plotAnglesNoEff           = False
pdfConfig['makePlots']    = True
pdfConfig['AllTagPlots']  = True
pdfConfig['SFit']         = True
pdfConfig['blind']        = False
pdfConfig['nominalPdf']   = False  # nominal PDF option does not work at the moment
corrSFitErr               = 'sumWeight' # [ 1., 0.700, 0.952, 0.938, 0.764 ] # '' / 'matrix' / 'sumWeight'
randomParVals             = ( ) # ( 1., 12346 ) # ( 2., 12345 )

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
             #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',                           'f_S_bin4',        'f_S_bin5'
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
pdfConfig['KKMassBinBounds']      = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.33, 0.09), (0.073, 0.030), (0.009, 0.012), (0.012, 0.010), (0.061, 0.027), (0.18, 0.04) ]
#                                     , [ (1.1,  0.5 ), (0.7,   0.2  ), (0.4,   0.4  ), (-0.6,  0.3  ), (-0.4, 0.2   ), (-0.7, 0.2 ) ] )
pdfConfig['SWaveAmplitudeValues'] = (  [ (0.23, 0.08), (0.068, 0.029), (0.013, 0.008), (0.013, 0.008), (0.055, 0.026), (0.17,  0.04) ]
                                     , [ (1.3,  0.7 ), (0.76,   0.27), (0.39,  0.25 ), (-0.57, 0.28 ), (-0.46, 0.21 ), (-0.64, 0.20) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.22, 0.07),  (0.057, 0.028), (0.003, 0.007), (0.012, 0.009), (0.046, 0.026), (0.16,  0.04) ]
#                                     , [ (1.43, 0.85 ), (0.83,  0.34 ), (0.8,   1.4  ), (-0.62, 0.33 ), (-0.52, 0.25 ), (-0.67, 0.21) ] )
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
                               , ( 'trphi',      '#varphi_{tr}'        )
                              )
else :
    pdfConfig['angleNames'] = (  ( 'helcosthetaK', 'cos(#theta_{K})'   )
                               , ( 'helcosthetaL', 'cos(#theta_{#mu})' )
                               , ( 'helphi',       '#varphi_{h}'          )
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

#Nominal-Unweighted data
nomData = pdfBuild['data']
#sWeighted data
sigSWeigtData = pdfBuild['sigSWeightData']

#Import ploting stuff
from P2VVGeneralUtils import plot
from ROOT import TCanvas, TPaveText, gStyle
from P2VVLoad import LHCbStyle

#########################
#Markers and lines style
#########################
lineWidth = LHCbStyle.lhcbWidth
markStyle = LHCbStyle.lhcbStyle.GetMarkerStyle()
markSize  = LHCbStyle.lhcbStyle.GetMarkerSize()

###################
#LHCbLabel
###################
lhcbName = TPaveText( 0.24, 0.733, 0.391, 0.803, 'BRNDC')
lhcbName.AddText("LHCb")
lhcbName.SetFillColor(0)
lhcbName.SetTextAlign(12)
lhcbName.SetBorderSize(0)

############################################################
## Make the estimated OS and SS tagging probabilities plots
############################################################
from P2VVGeneralUtils import _P2VVPlotStash
tagPlots = []
canvList = []
for plot in _P2VVPlotStash:
    if plot.GetName() and plot.GetName().startswith('tagomega'):
        oldTitle = plot.GetYaxis().GetTitle()
        try: measure  = round(float(oldTitle[11:20]), 6)
        except ValueError:
            'Cannot round the y axzis measure, set measure number starting point.'
            'Setting default measure value (it might not be the correct one)'
            measure = 0.01

        plot.SetYTitle( 'Candidates / ' +  str(measure) ) 
        plot.SetMinimum(0)

        tagPlots.append(plot)
        canvList.append( TCanvas(plot.GetName(),plot.GetTitle()) )

for ( canv, plot) in zip (canvList , tagPlots):
    canv.cd()
    plot.Draw()
    lhcbName.Draw()


##################################################
## Make the B mass plots
##################################################

for plot in _P2VVPlotStash:
   if plot.GetName():
        name = plot.GetName()
        if plot.GetName() == 'mass fit - signal': BMassPlot = plot 

BMassPlot.SetXTitle( 'M(J/#psiK^{+}K^{-}) [MeV/c^{2}]' )
OldTitle = BMassPlot.GetYaxis().GetTitle()

BMassPlot.SetYTitle('Candidates /  (' + OldTitle[11:15] + ' MeV/c^{2})' )

BMassCanv = TCanvas('BMassCanv')
BMassPlot.Draw()
lhcbName.Draw()
assert(False)
##################################################
## Make the mumu mass plots
##################################################
from P2VVParameterizations.MassPDFs import Signal_PsiMass, Linear_Background_Mass

mumuSig = Signal_PsiMass(         Name = 'sig_mumu', mass = mumuMass )
mumuBkg = Linear_Background_Mass( Name = 'bkg_mumu', mass = mumuMass )

#Initialize Components
from RooFitWrappers import buildPdf, Component
mumuMassSigComp = Component('mumuSig', ( mumuSig.pdf(), ), Yield = (21000,10000,50000) )
mumuMassBkgComp = Component('mumuBkg', ( mumuBkg.pdf(), ), Yield = (10000,5000,15000 ) )
#Build mumu Mass Pdf
mumuMassPdf = buildPdf(Components = ( mumuMassSigComp, mumuMassBkgComp ), Observables = (mumuMass, ), Name='mumuMassPdf'  )

#Fit
ws.var('m_bkg_arg').setRange(-1e-4,.5)
mumuMassPdf.fitTo(nomData)

#Plot
mumuMassPlot = mumuMass.frame()
mumuMassPlot.SetXTitle('m(#mu^{+}#mu^{-})' + '[MeV/c^{2}]')

sigArgs = {'Components':'sig_mumu', 'LineColor':kRed,     'LineStyle':10, 'LineWidth':3}
bkgArgs = {'Components':'bkg_mumu', 'LineColor':kGreen+3, 'LineStyle': 2, 'LineWidth':3}

nomData.plotOn(mumuMassPlot)
mumuMassPdf.plotOn(mumuMassPlot)
mumuMassPdf.plotOn(mumuMassPlot, **sigArgs)
mumuMassPdf.plotOn(mumuMassPlot, **bkgArgs)

OldTitle = mumuMassPlot.GetYaxis().GetTitle()
mumuMassPlot.SetYTitle('Candidates /  ' + OldTitle[11:15] + ' MeV/c^{2}'  )

mumuMassCanv = TCanvas('mumuMassCanv')
mumuMassPlot.Draw()
lhcbName.Draw()


##################################################
## Make the KK mass plots
##################################################
#Build KK mas pdf
from ROOT import RooRealVar, RooRelBreitWigner, RooConstVar, RooFFTConvPdf, RooGaussModel

#Mass Resolution Model
KKMassVar = KKMass._target_()
resMean    = RooRealVar   ( 'resMean'   , 'resMean'   , -0.3, -1.5, 1.2) 
resSigma   = RooRealVar   ( 'resSigma'  , 'resSigma'  ,  0.1, 0.05, 0.2) 
GaussModel = RooGaussModel('GaussModel', 'GaussModel' , KKMassVar, resMean,resSigma )

#Build Phi Mass Pdf
  #Phi Mass parameters
mean   = RooRealVar ( 'mean'  , 'mean'  , 1021, 1014, 1026 ) 
width  = RooRealVar ( 'width' , 'width' ,  5.3,    3,  6.5 )
spin   = RooConstVar( 'spin'  , 'spin'  ,    1             )  
radius = RooRealVar ( 'radius', 'radius', -3.8, -4.5,    0 )
K1mass = RooConstVar( 'K1mass', 'K1mass', 493.7            )
K2mass = RooConstVar( 'K2mass', 'K2mass', 493.7            )
 #Phi Mass pdf 
_PhiMassPdf  = RooRelBreitWigner('_PhiMassPdf', '_PhiMassPdf', KKMassVar, mean, width, spin, radius, K1mass, K2mass )
PhiMassPdf   = RooFFTConvPdf    ('PhiMassPdf' , 'PhiMassPdf' , KKMassVar, _PhiMassPdf, GaussModel )

#s-Wave ( Background ) Pdf
  #s-Wave Parameters
from ROOT import RooDstD0BG as PhaseSpaceFunc
dm0   = RooRealVar('dm0'  , 'dm0'  ,  990,          )
par_A = RooRealVar('par_A', 'par_A',  4,   1,  6)
par_B = RooRealVar('par_B', 'par_B',  8,    7,  11)
par_C = RooRealVar('par_C', 'par_C',  14,     13, 16)
  #s-Wave Pdf
_KKbkgPdf = PhaseSpaceFunc('_KKbkgPdf', '_KKbkgPdf', KKMassVar, dm0, par_A, par_B, par_C )
KKbkgPdf  = RooFFTConvPdf ('KKbkgPdf' , 'KKbkgPdf' , KKMassVar, _KKbkgPdf , GaussModel )

#Total Pdf
from ROOT import RooAddPdf
relCoef         = RooRealVar('relCoef', 'relCoef', 0.5, 0.25, 0.75) 
KKMassPdf  = RooAddPdf ('totPdf', 'totPdf', PhiMassPdf, KKbkgPdf, relCoef )


assert(False)


#Fit
KKMassPdf.fitTo(nomData)

#Plot
from ROOT import TCanvas
KKmassCanv = TCanvas()
comps = { 'PhiMassPdf':dict( LineColor=kRed,      LineStyle=10, LineWidth=lineWidth ), 
          'KKbkgPdf'  :dict( LineColor=kGreen +3, LineStyle= 2, LineWidth=lineWidth )
            }
plot(KKmassCanv, KKMass, nomData, KKMassPdf
     , components = comps
     , frameOpts = {'Bins':50}

     )

KKmassCanv.cd()
lhcbName.Draw()


assert(False)



##################################################
## Save plots
#################################################

#path = '~/../../project/bfys/vsyropou/PhD/plots/paperPlots/'
#BMassCanv.Print(path + 'BMass.pdf')
#mumuMassCanv.Print(path + 'mumuMass.pdf')
#for ( canv, plot) in zip (canvList , tagPlots):
#    canv.Print( path + plot.GetName() + '.pdf')



