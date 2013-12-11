import sys
import os
from P2VV.ToyMCUtils import Toy

toy = Toy()
parser = toy.parser()
(options, args) = toy.configure()

###########################################################################################################################################
## set script parameters ##
###########################

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
pdfConfig['selection']  = 'HLT1Unbiased' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['makePlots']  = False
pdfConfig['sFit']       = True
pdfConfig['blind']      = False
corrSFitErr             = 'sumWeight'     # '' / 'matrix' / 'sumWeight'

parameterFile = 'CFit_unbiased.pars'
pdfConfig['nTupleName'] = 'DecayTree'
pdfConfig['nTupleFile'] = '/bfys/raaij/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120821_MagDownMagUp.root'

dataSetName = 'JpsiphiData'
dataSetFile = 'SFit_generate.root'

dllPars = [ ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

# fit options
fitOpts = dict(  NumCPU    = 4
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
pdfConfig['SWeightsType']         = ''
pdfConfig['KKMassBinBounds']      = [ 1008., 1032. ] # [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ]
#pdfConfig['SWaveAmplitudeValues'] = (  [ (0.27, 0.09), (0.079, 0.032), (0.016, 0.015), (0.011, 0.010), (0.057, 0.029), (0.17, 0.05) ]
#                                     , [ (1.4,  0.5 ), (0.8,   0.3  ), (0.3,   0.3  ), (-0.5,  0.3  ), (-0.5, 0.2   ), (-0.7, 0.2 ) ] )
#pdfConfig['SWaveAmplitudeValues'] = (  [ ( 0.18, 0.07 ), ( 0.07, 0.03 ), ( 0.01, 0.02 ), ( 0.02, 0.01 ), ( 0.05, 0.03 ), ( 0.15, 0.04 ) ]
#                                     , [ ( 1.4,  0.5  ), ( 0.8,  0.3  ), ( 0.3,  0.4  ), ( -0.5, 0.2  ), ( -0.5, 0.2  ), ( -0.7, 0.2  ) ] )
pdfConfig['CSPValues']            = [ 0.498 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

pdfConfig['sameSideTagging']    = True  # nominal: False
pdfConfig['conditionalTagging'] = True  # nominal: True
pdfConfig['continuousEstWTag']  = True  # default: False | nominal: True
pdfConfig['numEstWTagBins']     = 100
pdfConfig['constrainTagging']   = 'constraint'  # nominal: 'constraint'

pdfConfig['numTimeResBins']        = 100
pdfConfig['constrainTimeResScale'] = 'constraint'  # nominal: 'constraint'

pdfConfig['massRangeBackground'] = True

pdfConfig['amplitudeParam'] = 'phasesSWaveFrac' # default: 'bank' | nominal: 'phasesSWaveFrac'
pdfConfig['ASParam']        = 'deltaPerp'  # default/nominal: 'deltaPerp'
pdfConfig['AparParam']      = 'phase' # default: 'Mag2ReIm' | nominal: 'phase'

pdfConfig['constrainDeltaM'] = 'constraint'  # nominal: 'constraint'

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

pdfConfig['timeEffHistFiles'] = dict(  file      = '/bfys/raaij/p2vv/data/Bs_HltPropertimeAcceptance_Data-20120816.root'
                                     , hlt1UB    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                     , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                    )
pdfConfig['angEffMomsFiles']    = '/bfys/raaij/p2vv/data/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'\
                                    if pdfConfig['transversityAngles'] else\
                                    '/bfys/raaij/p2vv/data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'

if pdfConfig['transversityAngles'] :
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

from P2VV.Load import RooFitOutput

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

# get variables
obsSetP2VV = [pdfBuild['observables'][obs] for obs in ['time', 'cpsi', 'ctheta', 'phi']]
time       = obsSetP2VV[0]
angles     = obsSetP2VV[ 1 : 4 ]
BMass      = pdfBuild['observables']['BMass']
mumuMass   = pdfBuild['observables']['mumuMass']
KKMass     = pdfBuild['observables']['KKMass']
estWTagOS  = pdfBuild['observables']['estWTagOS']
timeRes    = pdfBuild['observables']['timeRes']

gen_observables = [time] + angles

if not pdfConfig['sFit'] : obsSetP2VV.append(BMass)

if not pdfBuild['iTagZeroTrick'] :
    tagCatP2VVOS = pdfBuild['observables']['tagCatP2VVOS']
    tagCatP2VVSS = pdfBuild['observables']['tagCatP2VVSS']

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

## Read data to use protodata
from P2VV.Utilities.DataHandling import readData
proto_data = pdfBuild['data']
from ROOT import RooArgSet
conditionals = RooArgSet()
for c in pdf.ConditionalObservables():
    r = conditionals.add(c._target_())
proto_data = proto_data.reduce(conditionals)

# read parameters from file
pdfConfig.readParametersFromFile( filePath = parameterFile )
pdfConfig.setParametersInPdf(pdf)

# Hack to set maximum value
pdf.setMaxVal(0.15)

# Wrong PV shape
from P2VV.Parameterizations.WrongPV import ShapeBuilder
from P2VV.RooFitWrappers import Component, Projection
wpv = ShapeBuilder(time, {'B' : BMass}, UseKeysPdf = True, Weights = 'B',
                   InputFile = "/bfys/raaij/p2vv/data/Bs2JpsiPhiPrescaled_2011.root")
wpv_signal = wpv.shape('B')

# Projection over time of signal PDF
projection = Projection(Name = 'pdf_angles', Original = pdf, ProjectVars = [time])
signal_wpv = Component('signal_wpv', (projection, wpv_signal), Yield = (100, 50, 10000))

signal = pdfBuild._signalComps
signal.setYield(10000, 5000, 50000)

from P2VV.RooFitWrappers import buildPdf
genPdf = buildPdf(Components = [signal_wpv, signal], Observables = angles + [time], Name = 'genPdf')

# run the toy
genPdf.fixAddCoefNormalization(RooArgSet(*gen_observables))

toy.set_fit_opts(**dict(Verbose = False))
toy.run(Observables = gen_observables, Pdf = pdf, GenPdf = genPdf, ProtoData = proto_data)

toy.write_output()
