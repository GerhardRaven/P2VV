###########################################################################################################################################
## set script parameters ##
###########################

from math import pi

runPeriods = '3fb' # '2011' # '2012' # '3fb'
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig( RunPeriods = runPeriods )

# job parameters
from P2VV.Parameterizations.FullPDFs import SimulCatSettings
generateData        = False
doFit               = True #'NLL'
makeObservablePlots = False
makeKKMassPlots     = False
plotAnglesNoEff     = False
corrSFitErr         = ( 'sumWeight', [ 'runPeriod', 'KKMassCat' ] )
randomParVals       = ( ) #( 0.2, 12345 )
dataPath            = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'

plotsFile     = 'temp.ps'   #'/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/20112012Reco14_angEffSimple_timeLin.ps'
plotsROOTFile = 'temp.root' #'/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/20112012Reco14_angEffSimple_timeLin.root'
parFileIn     = '20112012Reco14DataFitValues_6KKMassBins.par' # '20112012Reco14DataFitValues_6KKMassBins.par' # '20112012Reco14DataFitValues_6KKMassBins_CPVDecay.par'
parFileOut    = '' # ( 'polDep_3fb.par', dict( Format = 'common' ) )

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'paper2012_SFit.root' if pdfConfig['sFit'] else 'paper2012_CFit.root'
elif pdfConfig['sFit'] :
    dataSetName = 'JpsiKK_sigSWeight'
    if pdfConfig['runPeriods'] == [ 2011 ] :
        dataSetFile = dataPath + 'P2VVDataSets2011Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
    elif pdfConfig['runPeriods'] == [ 2012 ] :
        dataSetFile = dataPath + 'P2VVDataSets2012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
    else :
        dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'
else :
    dataSetName = 'JpsiKK'
    dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats.root'

MinosPars = [#  '__phiCP__', '__dGamma__', 'lambdaCP', 'Gamma', 'dM'
             #, 'A0Mag2', 'AperpMag2', 'AparPhase', 'AperpPhase'
             #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
             #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
             #, 'betaTimeEff_p2011', 'betaTimeEff_p2012'
            ]
dllPars = [ ] # [ ( 'ASOddPhase_bin2', True, True, True ) ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

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
manualTagCatBins  = False
constTagCatCoefs  = True
constAvgCEvenOdd  = True
constWTagAsyms    = ''
constCSP          = True
constAmplitudes   = False
constLambdaCP     = ''  # 'lamb'
equalAbsLambdaCPs = False

# PDF options
pdfConfig['timeResType']           = 'event3fb' # 'eventNoMean'
pdfConfig['timeEffType']           = 'paper2012' # 'fit_uniformUB' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['timeEffConstraintType'] = 'poisson'
pdfConfig['timeEffParameters']     = { } # dict(RandomBinOrder = False) # dict( Parameterization = 'Spline', Fit = False ) # dict( Fit = False, RandomBinOrder = False )
pdfConfig['constrainBeta']         = ''  # '' / 'constrain' / 'fixed' / 'noBeta'
timeResSystType                    = ''  # 'mean_param' / 'quadratic_no_offset' / 'linear'

timeEffFile2011 = dataPath + 'timeAcceptanceFit_2011.root' # 'Bs_HltPropertimeAcceptance_Data_2011_40bins_TOS.root'
timeEffFile2012 = dataPath + 'timeAcceptanceFit_2012.root' # 'Bs_HltPropertimeAcceptance_Data_2012_40bins_TOS.root'
timeEffHist2011UB = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached' # 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'

from collections import defaultdict
splitConstr = defaultdict(dict)
if timeResSystType == 'mean_param':
    pdfConfig['timeResType'] += '_mean_param'
    ## splitConstr['timeResMu']['2011']        = ( -0.00259 - 0.00075, 0. )
    ## splitConstr['timeResMu']['2012']        = ( -0.00333 - 0.0012 , 0. )
    splitConstr['mu_placeholder']['2011']  = (  0.0350, 0. )
    ##                                           data       sig MC    prompt MC
    splitConstr['timeResMu_offset']['2011'] = ( -0.00255 - (0.00174 - 0.00056), 0. )
    splitConstr['timeResMu_slope']['2011']  = (  0       - (0.135   - 0      ), 0. )
    splitConstr['timeResMu_quad']['2011']   = (  0       - (3.50    - 0      ), 0. )
    
    splitConstr['mu_placeholder']['2012']  =  (  0.0349, 0. )
    ##                                           data       sig MC    prompt MC
    splitConstr['timeResMu_offset']['2012'] = ( -0.00322 - (0.00169 - 0.00072), 0. )
    splitConstr['timeResMu_slope']['2012']  = (  0       - (0.127   - 0      ), 0. )
    splitConstr['timeResMu_quad']['2012']   = (  0       - (2.29    - 0      ), 0. )
    split_runPeriod = set(pdfConfig['splitParams']['runPeriod'])
    for pn in ('timeResMu',):
        pdfConfig['externalConstr'].pop(pn)
        pdfConfig['splitParams']['runPeriod'].remove(pn)
    for par, vals in splitConstr.iteritems() :
        if par not in split_runPeriod:
            pdfConfig['splitParams']['runPeriod'].append(par)

if pdfConfig['timeResType'].endswith('linear'):
    from collections import defaultdict
    splitConstr = defaultdict(dict)
    splitConstr['sf_placeholder']['2011']  = (  0.0350, 0. )
    splitConstr['sf_mean_offset']['2011']  = ( 0.04956, 0. )
    splitConstr['sf_mean_slope']['2011']   = (   1.313, 0. )
    splitConstr['timeResFrac2']['2011']    = (  0.2449, 0. )
    splitConstr['sf_sigma_offset']['2011'] = ( 0.01302, 0. )
    splitConstr['sf_sigma_slope']['2011']  = (  0.2668, 0. )

    splitConstr['sf_placeholder']['2012']  =  (  0.0349, 0. )
    splitConstr['sf_mean_offset']['2012']  =  ( 0.05139, 0. )
    splitConstr['sf_mean_slope']['2012']   =  (   1.304, 0. )
    splitConstr['timeResFrac2']['2012']    =  (   0.244, 0. )
    splitConstr['sf_sigma_offset']['2012'] =  ( 0.01385, 0. )
    splitConstr['sf_sigma_slope']['2012']  =  (   0.242, 0. )
    for pn in ('sf_mean_quad', 'sf_sigma_quad'):
        pdfConfig['externalConstr'].pop(pn)
        pdfConfig['splitParams']['runPeriod'].remove(pn)
    split_runPeriod = set(pdfConfig['splitParams']['runPeriod'])
    for par, vals in splitConstr.iteritems() :
        if par not in split_runPeriod:
            pdfConfig['splitParams']['runPeriod'].append(par)

if timeResSystType == 'quadratic_no_offset':
    from collections import defaultdict
    splitConstr = defaultdict(dict)
    splitConstr['sf_placeholder']['2011']  = (  0.0350,  0. )
    splitConstr['sf_mean_slope']['2011']   = (  1.427, 0. )
    splitConstr['sf_mean_quad']['2011']    = ( -7.802, 0. )
    splitConstr['timeResFrac2']['2011']    = ( 0.1572, 0. )
    splitConstr['sf_sigma_slope']['2011']  = ( 0.3301, 0. )
    splitConstr['sf_sigma_quad']['2011']   = (  2.755, 0. )

    splitConstr['sf_placeholder']['2012']  =  (  0.0349,  0. )
    splitConstr['sf_mean_slope']['2012']   =  (  1.435, 0. )
    splitConstr['sf_mean_quad']['2012']    =  (  5.233, 0. )
    splitConstr['timeResFrac2']['2012']    =  (  0.312, 0. )
    splitConstr['sf_sigma_slope']['2012']  =  ( 0.2775, 0. )
    splitConstr['sf_sigma_quad']['2012']   =  (  1.846, 0. )
elif timeResSystType == 'linear':

    splitConstr['sf_placeholder']['2011']  = (  0.0349, 0. )
    splitConstr['sf_mean_offset']['2011']  = (  0.04499, 0. )
    splitConstr['sf_mean_slope']['2011']   = (    1.285, 0. )
    splitConstr['timeResFrac2']['2011']    = (   0.3254, 0. )
    splitConstr['sf_sigma_offset']['2011'] = ( 0.009135, 0. )
    splitConstr['sf_sigma_slope']['2011']  = (    0.261, 0. )
    
    splitConstr['sf_placeholder']['2012']  =  ( 0.0350, 0. )
    splitConstr['sf_mean_offset']['2012']  =  (  0.04475, 0. )
    splitConstr['sf_mean_slope']['2012']   =  (    1.282, 0. )
    splitConstr['timeResFrac2']['2012']    =  (   0.4581, 0. )
    splitConstr['sf_sigma_offset']['2012'] =  ( 0.006587, 0. )
    splitConstr['sf_sigma_slope']['2012']  =  (   0.1887, 0. )

for par, vals in splitConstr.iteritems() :
    constr = SimulCatSettings( '%sConstr' % par )
    constr.addSettings( [ 'runPeriod' ], [ [ 'p2011' ] ], vals['2011'] )
    constr.addSettings( [ 'runPeriod' ], [ [ 'p2012' ] ], vals['2012'] )
    pdfConfig['externalConstr'][par] = constr

if pdfConfig['timeEffType'].startswith('fit') :
    pdfConfig['timeEffData']['file'] = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140309.root'
    pdfConfig['externalConstr']['acceptance'] = SimulCatSettings('acceptanceConstr')
    pdfConfig['externalConstr']['acceptance'].addSettings(['runPeriod'], [['p2011']],
                                                          {('hlt1_excl_biased_dec', 'exclB') : (0.65, 0.01),
                                                           ('hlt2_biased', 'B') : (0.65, 0.01)})
    pdfConfig['externalConstr']['acceptance'].addSettings(['runPeriod'], [['p2012']],
                                                          {('hlt1_excl_biased_dec', 'exclB') : (0.65, 0.01),
                                                           ('hlt2_biased', 'B') : (0.65, 0.01)})

if pdfConfig['runPeriods'] == [ 2011 ] :
    pdfConfig['timeEffHistFiles']['file'] = timeEffFile2011
    pdfConfig['timeEffHistFiles']['hlt1UB'] = timeEffHist2011UB
elif pdfConfig['runPeriods'] == [ 2012 ] :
    pdfConfig['timeEffHistFiles']['file'] = timeEffFile2012
elif pdfConfig['runPeriods'] == [ 2011, 2012 ] :
    pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeEffFile2011
    pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeEffFile2012
    pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['hlt1UB'] = timeEffHist2011UB
    #pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['hlt1UB']\
    #        = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'

pdfConfig['anglesEffType'] = 'weights' # 'basis012' # 'basisSig6'
pdfConfig['constAngEffCoefs'] = True
pdfConfig['angEffMomsFiles'] = dataPath + 'angEffNominalRew_moms.par' # 'angEffNominalRew_5thOrder.par' # 'angEffNominalRew_moms.par' # 'angEffNominalNoRew_moms.par' # 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
#angEffMomsFiles = SimulCatSettings('angEffMomsFiles')
#angEffMomsFiles.addSettings( [ 'KKMassCat' ], [ [ 'bin0', 'bin1', 'bin2' ] ]
#                            , dataPath + 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK09901020_Phys_moms_norm'
#                           )
#angEffMomsFiles.addSettings( [ 'KKMassCat' ], [ [ 'bin3', 'bin4', 'bin5' ] ]
#                            , dataPath + 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK10201050_Phys_moms_norm'
#                           )
#angEffMomsFiles.addSettings( [ 'tagCatP2VV' ], [ [ 'Untagged' ] ]
#                            , dataPath + 'Sim08_20112012_hel_UB_U_trueTime_BkgCat050_KK30_Phys_moms_norm'
#                           )
#angEffMomsFiles.addSettings( [ 'tagCatP2VV' ], [ [ 'Tagged' ] ]
#                            , dataPath + 'Sim08_20112012_hel_UB_T_trueTime_BkgCat050_KK30_Phys_moms_norm'
#                           )
#pdfConfig['angEffMomsFiles'] = angEffMomsFiles

pdfConfig['SSTagging']        = True
pdfConfig['condTagging']      = True
pdfConfig['contEstWTag']      = True
pdfConfig['constrainTagging'] = ''  # 'constrain'
if generateData or manualTagCatBins :
    pdfConfig['tagCatsOS'] = [  ( 'Untagged', 0, 0.5000001, 0.5,   0.5,   0.0, 0.669, 0.0 )
                              , ( 'Tagged',   1, 0.4999999, 0.392, 0.392, 0.0, 0.331, 0.0 )
                             ]
    pdfConfig['tagCatsSS'] = [  ( 'Untagged', 0, 0.5000001, 0.5,   0.5,   0.0, 0.896, 0.0 )
                              , ('Tagged',    1, 0.4999999, 0.359, 0.359, 0.0, 0.104, 0.0 )
                             ]

pdfConfig['paramKKMass']     = 'simultaneous'  # 'simultaneous'
pdfConfig['KKMassBinBounds'] = [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['CSPValues']       = [ 0.9178, 0.9022, 0.8619, 0.8875, 0.9360, 0.9641 ] # [ 0.9152, 0.8797, 0.8357, 0.8599, 0.9207, 0.9624 ] # [ 0.9586, 0.9110, 0.8618, 0.8828, 0.9227, 0.9110 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.326 ]
KKMassPars = pdfConfig['obsDict']['KKMass']
pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                  , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )

pdfConfig['lambdaCPParam'] = 'lambPhi' # 'lambPhi' # 'observables_CPVDecay' # 'lambPhi_CPVDecay' # 'lambPhiRel_CPVDecay' # 'lambPhi_CPVDecay_PSWaves'
if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    pdfConfig['splitParams']['KKMassCat'] = [ 'av' + par if par == 'f_S' else par for par in pdfConfig['splitParams']['KKMassCat'] ]

if 'Parameterization' in pdfConfig['timeEffParameters'] and pdfConfig['timeEffParameters']['Parameterization'] == 'Spline' :
    pdfConfig['splitParams']['hlt1_excl_biased_dec'] = [ 'tagCatCoef0_1' ]

#pdfConfig['externalConstr'].pop('dM')
#pdfConfig['externalConstr']['wTagP0OS'] = ( 0.3853, None )
#pdfConfig['externalConstr']['wTagP1OS'] = ( 1.000,  None )
#pdfConfig['externalConstr']['wTagP0SS'] = ( 0.450,  None )
#pdfConfig['externalConstr']['wTagP1SS'] = ( 1.00,   None )
#pdfConfig['externalConstr']['wTagDelP0OS'] = (  0.0140, 0. )
#pdfConfig['externalConstr']['wTagDelP1OS'] = (  0.066,  0. )
#pdfConfig['externalConstr']['wTagDelP0SS'] = ( -0.016,  0. )
#pdfConfig['externalConstr']['wTagDelP1SS'] = (  0.007,  0. )
#pdfConfig['externalConstr']['dM']          = (  17.768, 0.024  )
#pdfConfig['externalConstr'].pop('betaTimeEff')
#pdfConfig['externalConstr']['betaTimeEff'] = ( 0., 0. )
#pdfConfig['externalConstr']['betaTimeEff'][1] = ( dict( runPeriod = [ 'p2012' ] ), ( -0.0135, 0.004 ) )
#pdfConfig['externalConstr']['betaTimeEff'] = SimulCatSettings('betaConstr')
#pdfConfig['externalConstr']['betaTimeEff'].addSettings( [ 'runPeriod' ], [ [ 'p2011' ] ], ( -0.0083, 0. ) )
#pdfConfig['externalConstr']['betaTimeEff'].addSettings( [ 'runPeriod' ], [ [ 'p2012' ] ], ( -0.0083, 0. ) )
#pdfConfig['splitParams']['runPeriod'].remove('betaTimeEff')
#pdfConfig['splitParams']['runPeriod'].append('betaTimeEff')
#pdfConfig['splitParams']['runPeriod'].append('Gamma')
#pdfConfig['splitParams']['runPeriod'].remove('Gamma')

for par in [ 'wTagP0OS', 'wTagP0SS', 'wTagP1OS', 'wTagP1SS', 'wTagDelP0OS', 'wTagDelP0SS', 'wTagDelP1OS', 'wTagDelP1SS' ] :
    pdfConfig['externalConstr'].pop(par)

#for par in [ 'tres_placeholder', 'timeResMu', 'timeResFrac2', 'sf_mean_offset', 'sf_mean_slope', 'sf_sigma_offset', 'sf_sigma_slope' ] :
#    pdfConfig['splitParams']['runPeriod'].remove(par)
#    pdfConfig['externalConstr'].pop(par)

#pdfConfig['splitParams']['hlt1_excl_biased_dec'] = [ 'sf_mean_offset' ]

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

# plot options
from ROOT import gStyle, kBlack, kBlue, kRed, kGreen, kMagenta, kSolid, kDashed, kFullCircle, kFullSquare, kFullDotLarge
from P2VV.Load import RooFitOutput, LHCbStyle
#gStyle.SetEndErrorSize(3)
gStyle.SetPalette(1)
gStyle.SetNumberContours(50)
lineWidth     = 3
lineColor     = kBlue
markStyle     = kFullDotLarge
markSize      = 0.7
markColor     = kBlack
markLineWidth = 3
gStyle.SetLineStyleString( 5, ' 40 20 10 20'  )
gStyle.SetLineStyleString( 7, ' 40 20'        )
gStyle.SetLineStyleString( 9, ' 100 20'       )

if pdfConfig['transAngles'] :
    angleNames = (  ( 'trcospsi',   'cos#kern[0.1]{#psi_{tr}}'   )
                  , ( 'trcostheta', 'cos#kern[0.1]{#theta_{tr}}' )
                  , ( 'trphi',      '#varphi_{tr} [rad]'         )
                 )
else :
    angleNames = (  ( 'helcosthetaK', 'cos#kern[0.1]{#theta_{K}}'   )
                  , ( 'helcosthetaL', 'cos#kern[0.1]{#theta_{#mu}}' )
                  , ( 'helphi',       '#varphi_{h} [rad]'           )
                 )

obsNames = {  'timeLin' : ( 'time', 'Decay time [ps]' )
            , 'timeLog' : ( 'time', 'Decay time [ps]' )
            , 'ctk'     : angleNames[0]
            , 'ctl'     : angleNames[1]
            , 'phi'     : angleNames[2]
           }
obsYTitles = {  'timeLin' : 'Candidates / (%.4f ps)'
              , 'timeLog' : 'Candidates / (%.3f ps)'
              , 'ctk'     : 'Candidates / %.3f'
              , 'ctl'     : 'Candidates / %.3f'
              , 'phi'     : 'Candidates / (%.3f#kern[0.3]{#pi} rad)'
             }
yTitleOffsets = {  'timeLin' : 1.15
                 , 'timeLog' : 1.00
                 , 'ctk'     : 1.15
                 , 'ctl'     : 1.15
                 , 'phi'     : 1.15
                }
LHCbCoords = {  'timeLin' : ( 0.75, 0.81, 0.88, 0.89 )
              , 'timeLog' : ( 0.75, 0.81, 0.88, 0.89 )
              , 'ctk'     : ( 0.67, 0.81, 0.80, 0.89 )
              , 'ctl'     : ( 0.75, 0.81, 0.88, 0.89 )
              , 'phi'     : ( 0.75, 0.81, 0.88, 0.89 )
             }

doTimeAnglePlots = { 'timeLin' : True, 'timeLog' : False, 'ctk' : False, 'ctl' : False, 'phi' : False }
plotComponents   = { 'phi' : False, 'even' : False, 'odd' : False, 'S' : False }
numBins = { 'timeLin' : 40, 'timeLog' : 50, 'ctk' : 30, 'ctl' : 30, 'phi' : 30 }


###########################################################################################################################################
## build PDF ##
###############

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

if not generateData :
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
KKMass     = pdfBuild['observables']['KKMass']
estWTagOS  = pdfBuild['observables']['wTagOS']
timeRes    = pdfBuild['observables']['timeRes']

if not pdfConfig['sFit'] : obsSetP2VV.append(BMass)

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


###########################################################################################################################################
## generate data ##
###################

if generateData :
    # print parameter values
    print 120 * '='
    print 'JvLFit: observables and parameters in generation process:'
    for var in pdf.getVariables() : var.Print()
    print 120 * '='

    # generate data
    nEvents = int( pdfConfig['numEvents'] * ( pdfConfig['sigFrac'] if pdfConfig['sFit'] else 1. ) )
    print 'JvLFit: generating %d events' % nEvents
    import P2VV.Load.MultiCatGen
    from ROOT import RooMsgService,RooFit
    RooMsgService.instance().addStream(RooFit.INFO, Topic =  RooFit.Generation)
    #from ROOT import RooAbsPdf
    #cfg = RooAbsPdf.defaultGeneratorConfig()
    #cfg.getConfigSection("RooFoamGenerator").setRealValue("chatLevel",1)
    #cfg.Print()
    dataSet = pdf.generate( obsSetP2VV, nEvents )

    # additional observables
    if not pdfConfig['transAngles'] :
        from P2VV.Utilities.DataHandling import addTransversityAngles
        addTransversityAngles( dataSet, 'trcospsi',          'trcostheta',        'trphi'
                                      , angles[0].GetName(), angles[1].GetName(), angles[2].GetName() )

    # write data to file
    from P2VV.Utilities.DataHandling import writeData
    writeData( dataSetFile, dataSetName, dataSet )
    fitData = dataSet

elif pdfConfig['sFit'] :
    if corrSFitErr[0] == 'sumWeight'\
            or ( type(corrSFitErr[0]) != str and hasattr( corrSFitErr[0], '__iter__' ) and hasattr( corrSFitErr[0], '__getitem__' ) ) :
        from P2VV.Utilities.DataHandling import correctWeights
        fitData = correctWeights( dataSet, corrSFitErr[1], CorrectionFactors = None if corrSFitErr[0] == 'sumWeight' else corrSFitErr[0] )

    else :
        fitData = dataSet

else :
    fitData = dataSet

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)


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
    for par in pdfBuild['signalKKMass'].pdf().getParameters(fitData) : par.setConstant(True)
    if not pdfConfig['sFit'] :
        for par in pdfBuild['backgroundKKMass'].pdf().getParameters(fitData) : par.setConstant(True)

for par in pdfPars :
    if par.GetName().startswith('C_SP') : par.setConstant( True if constCSP else False )

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
    if not pdfConfig['sFit'] :
        pdfBuild['backgroundBMass'].setConstant('.*')
        pdfBuild['backgroundTime'].setConstant('.*')
        if hasattr( pdfBuild, '_bkgTaggingPdf' ) : pdfBuild['bkgTaggingPdf'].setConstant('.*')
    pdfBuild['amplitudes'].setConstant('C_SP')

if pdfConfig['lambdaCPParam'].startswith('lambPhi_CPVDecay') :
    if equalAbsLambdaCPs : pdfBuild['lambdaCP'].setConstant('rhoCP_A.*')
    else :                 pdfBuild['lambdaCP'].setConstant('rhoCP_m')

## if pdfConfig['timeEffParameters']['Fit'] :
##     for period in [ 'p2011', 'p2012' ] :
##         for cat, coefs in pdfBuild['timeResModels'][ '{bin0;%s}' % period ].shapes()[0].coefficients().iteritems() :
##             for coef in coefs :
##                 if coef.GetName().endswith('001') :
##                     coef.setConstant(True)
##                 if cat[2] == 'notExclB' :
##                     coef.setVal( 0.5 * coef.getVal() )

#ws['timeResMu_p2011'].setVal(0.)
#ws['timeResMu_p2012'].setVal(0.)

#pdfBuild['lifetimeParams'].parameter('Gamma').setVal(0.72)

#ws['f_S'].setVal(0.)
#ws['ASOddPhase'].setVal(0.)
#ws['f_S'].setConstant()
#ws['ASOddPhase'].setConstant()

#for it in range(6) :
#    ws[ 'avf_S_bin%d' % it ].setVal(0.)
#    ws[ 'avf_S_bin%d' % it ].setConstant(True)
#    ws[ 'ASOddPhase_bin%d' % it ].setConstant(True)

#ws['f_S_bin2'].setVal(0.)
#ws['f_S_bin2'].setConstant()
#ws['ASOddPhase_bin2'].setVal(0.)
#ws['ASOddPhase_bin2'].setConstant()

#ws['timeResSigmaSF'].setVal(1.4)
#ws['timeResSigmaSF'].setConstant(False)

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
#ws['dGamma'].setVal(0.)
#ws['dGamma'].setConstant(True)
#ws['dM'].setVal(17.500)   #  17.600  17.667  17.700  17.768  17.800
#ws['A0Mag2'].setVal(0.534)
#ws['AperpMag2'].setVal(0.238)
#ws['f_S'].setVal(0.)
#ws['f_S'].setConstant()
#ws['AparPhase'].setVal(pi)
#ws['AperpPhase'].setVal(pi)
#ws['AparPhase'].setConstant(True)
#ws['AperpPhase'].setConstant(True)
#ws['ASOddPhase'].setVal(0.)
#ws['ASOddPhase'].setConstant()
#ws['lambdaCP'].setVal(1.)
#ws['lambdaCP'].setConstant()
#ws['phiCP'].setVal(-0.03)

#ws['rhoCP_P'].setConstant()
#ws['rhoCP_S'].setConstant()

if randomParVals :
    import random
    # give parameters random offsets
    print 'JvLFit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
          % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
    random.seed( randomParVals[1] if randomParVals[1] else None )
    for par in pdfPars :
        if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

## set S-wave fractions for the two middle bins equal
#print 'JvLFit: setting S-wave fractions for bin 2 and bin 3 equal'
#from P2VV.RooFitWrappers import Customizer
#pdf = Customizer( Pdf = pdf, OriginalArgs = [ ws['f_S_bin3'] ], SubstituteArgs = [ ws['f_S_bin2'] ] )
#pdfObs  = pdf.getObservables(fitData)
#pdfPars = pdf.getParameters(fitData)

# print parameters
print 120 * '='
print 'JvLFit: fit data:'
fitData.Print()
print 'JvLFit: observables in PDF:'
pdfObs.Print('v')
print 'JvLFit: parameters in PDF:'
pdfPars.Print('v')
print 'JvLFit: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()

if doFit :
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

    if doFit == 'NLL' :
        opts = fitOpts.copy()
        minimizer = opts.pop( 'Minimizer', 'Minuit2' )
        verbose   = opts.pop( 'Verbose',   False     )
        fitTime   = opts.pop( 'Timer',     False     )
        hesse     = opts.pop( 'Hesse',     True      )
        minos     = opts.pop( 'Minos',     False     )

        from ROOT import RooMinimizer
        nll = pdf.createNLL( fitData, **opts )
        minNLL = RooMinimizer(nll)
        minNLL.setMinimizerType(minimizer)
        minNLL.setVerbose(verbose)
        minNLL.setProfile(fitTime)

        minNLL.minimize(minimizer)
        if hesse : minNLL.hesse()

        if RooMinPars :
            from ROOT import RooArgSet
            minosSet = RooArgSet( par for par in RooMinPars )
            minNLL.minos(minosSet)

        elif minos :
            minNLL.minos()

        fitResult = minNLL.save(  'fitresult_%s_%s' % ( pdf.GetName(), fitData.GetName() )
                                , 'Result of fit of p.d.f. %s to dataset %s' % ( pdf.GetName(), fitData.GetName() )
                               )

    else :
        if pdfConfig['sFit'] :
            fitResult = pdf.fitTo( fitData, SumW2Error = True if corrSFitErr == 'matrix' else False
                                  , Minos = RooMinPars, Save = True, Range = fitRange
                                  , **fitOpts
                                 )
        else :
            fitResult = pdf.fitTo( fitData
                                  , Minos = RooMinPars, Save = True, Range = fitRange
                                  , **fitOpts
                                 )

    # reparameterize amplitudes
    if pdfConfig['amplitudeParam'] == 'bank' and pdfConfig['ASParam'] != 'ReIm' and pdfConfig['AparParam'] == 'Mag2ReIm' :
        from ROOT import RooArgSet
        parList = fitResult.floatParsFinal()
        parSet  = RooArgSet(parList)
        parCovs = fitResult.covarianceMatrix()

        from ROOT import RooRealVar, RooArgList
        from P2VV.Parameterizations.DecayAmplitudes import A02, Aperp2, Apar2, A0Ph, AperpPh, AparPh, f_S, AS2, ASPh
        deltaPar  = AparPh  - A0Ph
        deltaPerp = AperpPh - A0Ph
        deltaS    = ASPh    - A0Ph
        if pdfConfig['ambiguityPars'] :
            deltaPar  = -deltaPar
            deltaPerp = pi - deltaPerp
            deltaS    = -deltaS

        # physics parameters to determine
        ampPhys = [  RooRealVar( 'A0Mag2_phys',     'A0Mag2_phys',     A02,        0.,      1.      )  # 0
                   , RooRealVar( 'AparPhase_phys',  'AparPhase_phys',  deltaPar,  -2. * pi, 2. * pi )  # 1
                   , RooRealVar( 'AperpMag2_phys',  'AperpMag2_phys',  Aperp2,     0.,      1.      )  # 2
                   , RooRealVar( 'AperpPhase_phys', 'AperpPhase_phys', deltaPerp, -2. * pi, 2. * pi )  # 3
                  ]

        if pdfConfig['paramKKMass'] :
            numKKMassBins = pdfBuild['KKMassBinning'].numBins() if pdfConfig['paramKKMass'] == 'functions'\
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

        if pdfConfig['paramKKMass'] :
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

        if pdfConfig['paramKKMass'] :
            for bin in range( numKKMassBins ) :
                ampMeasFuncs[ ampMeasNames[ 5 + 2 * bin ] ] = '@{0:d}/(1.-@{0:d})/@2'.format( 4 + 2 * bin )
                ampMeasFuncs[ ampMeasNames[ 6 + 2 * bin ] ] = '@{0:d}-@3{1:s}'.format( 5 + 2 * bin, '+TMath::TwoPi()'\
                                                                                       if pdfConfig['ambiguityPars'] else '' )
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

    if 'CPVDecay' in pdfConfig['lambdaCPParam'] :
        from P2VV.Imports import parNames, parValuesCPVDecay as parValues
    elif 'lamb' in constLambdaCP.lower() :
        from P2VV.Imports import parNames, parValuesFixLamb as parValues
    else :
        from P2VV.Imports import parNames, parValues
    print 'JvLFit: parameters:'
    fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    print 'JvLFit: correlation matrix:'
    fitResult.correlationMatrix().Print()
    print 'JvLFit: covariance matrix:'
    fitResult.covarianceMatrix().Print()

    print 120 * '=' + '\n'

else :
    fitResult = None

if parFileOut :
    # write parameters to file
    filePath = parFileOut[0] if type(parFileOut) != str else parFileOut
    fileOpts = parFileOut[1] if type(parFileOut) != str else { }
    pdfConfig.getParametersFromPdf( pdf, fitData )
    pdfConfig.writeParametersToFile(  filePath = filePath
                                    , FitStatus = ( fitResult.status(), fitResult.minNll(), fitResult.edm() )\
                                                  if fitResult else ( -1, 0., 0. )
                                    , **fileOpts
                                   )


###########################################################################################################################################
## make some plots ##
#####################

if makeObservablePlots or makeKKMassPlots or dllPars :
    # import plotting tools
    from P2VV.Utilities.Plotting import plot, _P2VVPlotStash
    from ROOT import TCanvas

    if type(makeObservablePlots) == str :
        # reduce data set
        dataSetOld = dataSet
        dataSet = dataSet.reduce(makeObservablePlots)

    print 'JvLFit: plot data set:'
    dataSet.Print()

    # create projection data set for conditional observables
    if pdfConfig['sFit'] :
        comps = None
    else :
        comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                }

    from P2VV.RooFitWrappers import SimultaneousPdf
    projWDataSet     = [ obs for obs in pdf.ConditionalObservables() ]
    if isinstance( pdf, SimultaneousPdf ) :
        if pdf.indexCat().isFundamental() :
            projWDataSet += [ pdf.indexCat() ]
        else :
            projWDataSet += [ cat for cat in pdf.indexCat().inputCatList() ]

    if projWDataSet :
        #bulkData = dataSet.reduce( CutRange = 'Bulk' )
        projWData         = dict( ProjWData = ( dataSet.reduce(  ArgSet = projWDataSet     ), False ) )
        #projWDataBulk     = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataSet     ), False ) )
        print 'JvLFit: plot projection data sets:'
        projWData['ProjWData'][0].Print()
        #projWDataBulk['ProjWData'][0].Print()
    else :
        projWData         = dict()
        #projWDataBulk     = dict()

if makeKKMassPlots and pdfConfig['paramKKMass']\
        and ( ( pdfConfig['amplitudeParam'] == 'bank' and pdfConfig['ASParam'] != 'ReIm' )\
              or ( pdfConfig['amplitudeParam'] == 'phasesSWaveFrac' and pdfConfig['ASParam'] == 'deltaPerp' ) ) :
    # create S-wave phase plots
    nKKBins = pdfBuild['KKMassBinning'].numBins()

    KKMassVals = [ pdfBuild['KKMassBinning'].binLow(binIter) for binIter in range(nKKBins) ]
    KKMassVals.append( pdfBuild['KKMassBinning'].highBound() )

    deltaSVals     = [  pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getVal()     for binIter in range(nKKBins) ]
    deltaSLowErrs  = [ -pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getErrorLo() for binIter in range(nKKBins) ]
    deltaSHighErrs = [  pdfPars.find( 'ASOddPhase_bin%d' % binIter ).getErrorHi() for binIter in range(nKKBins) ]

    if pdfConfig['ambiguityPars'] :
        from math import pi
        deltaSVals = [ pi - val for val in deltaSVals ]
        tempErrs       = deltaSLowErrs
        deltaSLowErrs  = deltaSHighErrs
        deltaSHighErrs = tempErrs

    from P2VV.Utilities.Plotting import plotSWavePhases
    deltaSCanv = plotSWavePhases(  MassBins         = KKMassVals
                                 , DeltaSValues     = deltaSVals
                                 , DeltaSLowErrors  = deltaSLowErrs
                                 , DeltaSHighErrors = deltaSHighErrs
                                 , LHCbTextLine1    = ' LHCb'
                                 , LHCbTextLine2    = ''
                                )

else :
    deltaSCanv = None

plotCanvs = [ ]
plotsRFile = None
if makeObservablePlots :
    from ROOT import RooRealVar, RooConstVar, RooCategory, RooArgSet, RooCustomizer, RooExplicitNormPdf, TPaveText

    # LHCb text labels
    LHCbText = TPaveText( 0.24, 0.81, 0.37, 0.89, 'BRNDC')
    LHCbText.SetName('LHCbLabel')
    LHCbText.AddText('LHCb')
    LHCbText.SetFillStyle(0)
    LHCbText.SetFillColor(0)
    LHCbText.SetShadowColor(0)
    LHCbText.SetBorderSize(0)
    LHCbText.SetTextAlign(12)
    LHCbText.SetTextSize(0.072)
    _P2VVPlotStash.append(LHCbText)

    dummyCanv = TCanvas()
    LHCbText.Draw()
    dummyCanv.Update()

    # phi, even, odd and S-wave PDFs
    AparMag2Cust = RooRealVar( 'AparMag2Cust', 'AparMag2Cust', ws['AparMag2'].getVal(), 0., 1. )
    zeroCust = RooConstVar( 'zeroCust', 'zeroCust', 1.e-6 )
    phiPdfCust  = RooCustomizer( pdf._var, 'phi'  )
    evenPdfCust = RooCustomizer( pdf._var, 'even' )
    oddPdfCust  = RooCustomizer( pdf._var, 'odd'  )
    SPdfCust    = RooCustomizer( pdf._var, 'S'    )

    if pdfConfig['paramKKMass'] :
        for bin in range( pdfBuild['KKMassBinning'].numBins() ) : phiPdfCust.replaceArg( ws[ 'f_S_bin%d' % bin ], zeroCust )
    else :
        phiPdfCust.replaceArg( ws['f_S'], zeroCust )
    phiPdf = phiPdfCust.build()

    evenPdfCust.replaceArg( ws['AparMag2'],  AparMag2Cust )
    evenPdfCust.replaceArg( ws['AperpMag2'], zeroCust     )
    if pdfConfig['paramKKMass'] :
        for bin in range( pdfBuild['KKMassBinning'].numBins() ) : evenPdfCust.replaceArg( ws[ 'f_S_bin%d' % bin ], zeroCust )
    else :
        evenPdfCust.replaceArg( ws['f_S'], zeroCust )
    evenPdf = evenPdfCust.build()

    oddPdfCust.replaceArg( ws['AparMag2'], zeroCust )
    oddPdfCust.replaceArg( ws['A0Mag2'],   zeroCust )
    if pdfConfig['paramKKMass'] :
        for bin in range( pdfBuild['KKMassBinning'].numBins() ) : oddPdfCust.replaceArg( ws[ 'f_S_bin%d' % bin ], zeroCust )
    else :
        oddPdfCust.replaceArg( ws['f_S'], zeroCust )
    oddPdf = oddPdfCust.build()

    SPdfCust.replaceArg( ws['AparMag2'],  zeroCust )
    SPdfCust.replaceArg( ws['A0Mag2'],    zeroCust )
    SPdfCust.replaceArg( ws['AperpMag2'], zeroCust )
    SPdf = SPdfCust.build()

    polSet = RooArgSet( zeroCust, AparMag2Cust, ws['A0Mag2'], ws['AperpMag2'] )
    if pdfConfig['paramKKMass'] :
        for bin in range( pdfBuild['KKMassBinning'].numBins() ) : polSet.add( ws[ 'f_S_bin%d' % bin ] )
    else :
        polSet.add( ws['f_S'] )

    print 'JvLFit: plot PDF polarization parameters:'
    print 'total:'
    pdf.getObservables(polSet).Print()
    print 'phi:'
    phiPdf.getObservables(polSet).Print()
    print 'even:'
    evenPdf.getObservables(polSet).Print()
    print 'odd:'
    oddPdf.getObservables(polSet).Print()
    print 'S-wave:'
    SPdf.getObservables(polSet).Print()

    pdfNormSet = RooArgSet()
    for obs in pdf.Observables().difference( pdf.ConditionalObservables() ) : pdfNormSet.add(obs._var)
    timePdfSet = RooArgSet(pdfNormSet)
    ctkPdfSet  = RooArgSet(pdfNormSet)
    ctlPdfSet  = RooArgSet(pdfNormSet)
    phiPdfSet  = RooArgSet(pdfNormSet)
    timePdfSet.remove( pdfNormSet.find( time.GetName()      ) )
    ctkPdfSet.remove(  pdfNormSet.find( angles[0].GetName() ) )
    ctlPdfSet.remove(  pdfNormSet.find( angles[1].GetName() ) )
    phiPdfSet.remove(  pdfNormSet.find( angles[2].GetName() ) )
    timeSet = RooArgSet( pdfNormSet.find( time.GetName()      ) )
    ctkSet  = RooArgSet( pdfNormSet.find( angles[0].GetName() ) )
    ctlSet  = RooArgSet( pdfNormSet.find( angles[1].GetName() ) )
    phiSet  = RooArgSet( pdfNormSet.find( angles[2].GetName() ) )

    obsToPlot    = [ ]
    totPdfs      = [ ]
    phiPdfs      = [ ]
    evenPdfs     = [ ]
    oddPdfs      = [ ]
    SPdfs        = [ ]

    timeNormBulk = ( time.getMax('Bulk') - time.getMin('Bulk') ) * dataSet.sumEntries() / float( numBins['timeLin'] )
    timeNorm = ( time.getMax() - time.getMin() ) * dataSet.sumEntries() / float( numBins['timeLog'] )
    if doTimeAnglePlots['timeLin'] :
        # linear time PDFs
        obsToPlot.append('timeLin')
        plotCanvs.append( TCanvas( 'timeLinCanv', 'Decay time' ) )
        if projWData :
            totPdfs.append( RooExplicitNormPdf( 'timePdfLin', 'timePdf', timeSet, timePdfSet, pdf._var, pdf._var, timeNormBulk
                                                , projWData['ProjWData'][0] ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf( 'timePhiLin', 'timePhi',  timeSet, timePdfSet, phiPdf,  pdf._var, timeNormBulk
                                                    , projWData['ProjWData'][0] ) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf( 'timeEvenLin','timeEven', timeSet, timePdfSet, evenPdf, pdf._var, timeNormBulk
                                                    , projWData['ProjWData'][0] ) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf( 'timeOddLin', 'timeOdd',  timeSet, timePdfSet, oddPdf,  pdf._var, timeNormBulk
                                                    , projWData['ProjWData'][0] ) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf( 'timeSLin',   'timeS',    timeSet, timePdfSet, SPdf,    pdf._var, timeNormBulk
                                                    , projWData['ProjWData'][0] ) )
        else :
            totPdfs.append( RooExplicitNormPdf( 'timePdfLin',  'timePdf',  timeSet, timePdfSet, pdf._var, pdf._var, timeNormBulk ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf( 'timePhiLin',  'timePhi',  timeSet, timePdfSet, phiPdf,  pdf._var, timeNormBulk ) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf( 'timeEvenLin', 'timeEven', timeSet, timePdfSet, evenPdf, pdf._var, timeNormBulk ) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf( 'timeOddLin',  'timeOdd',  timeSet, timePdfSet, oddPdf,  pdf._var, timeNormBulk ) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf( 'timeSLin',    'timeS',    timeSet, timePdfSet, SPdf,    pdf._var, timeNormBulk ) )

    if doTimeAnglePlots['timeLog'] :
        # logarithmic time PDFs
        obsToPlot.append('timeLog')
        plotCanvs.append( TCanvas( 'timeLogCanv', 'Decay time' ) )
        if projWData :
            totPdfs.append( RooExplicitNormPdf( 'timePdfLog', 'timePdf', timeSet, timePdfSet, pdf._var, pdf._var, timeNorm
                                               , projWData['ProjWData'][0] ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf( 'timePhiLog',  'timePhi',  timeSet, timePdfSet, phiPdf,  pdf._var, timeNorm
                                                    , projWData['ProjWData'][0] ) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf( 'timeEvenLog', 'timeEven', timeSet, timePdfSet, evenPdf, pdf._var, timeNorm
                                                    , projWData['ProjWData'][0] ) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf( 'timeOddLog',  'timeOdd',  timeSet, timePdfSet, oddPdf,  pdf._var, timeNorm
                                                    , projWData['ProjWData'][0] ) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf( 'timeSLog',    'timeS',    timeSet, timePdfSet, SPdf,    pdf._var, timeNorm
                                                    , projWData['ProjWData'][0] ) )
        else :
            totPdfs.append( RooExplicitNormPdf( 'timePdfLog',  'timePdf',  timeSet, timePdfSet, pdf._var, pdf._var, timeNorm ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf( 'timePhiLog',  'timePhi',  timeSet, timePdfSet, phiPdf,  pdf._var, timeNorm ) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf( 'timeEvenLog', 'timeEven', timeSet, timePdfSet, evenPdf, pdf._var, timeNorm ) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf( 'timeOddLog',  'timeOdd',  timeSet, timePdfSet, oddPdf,  pdf._var, timeNorm ) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf( 'timeSLog',    'timeS',    timeSet, timePdfSet, SPdf,    pdf._var, timeNorm ) )

    if doTimeAnglePlots['ctk'] :
        # cos(theta_K) PDFs
        obsToPlot.append('ctk')
        plotCanvs.append( TCanvas( 'ctkCanv', 'Decay time' ) )
        ctkNorm = ( angles[0].getMax() - angles[0].getMin() ) * dataSet.sumEntries() / float( numBins['ctk'] )
        if projWData :
            totPdfs.append( RooExplicitNormPdf( 'ctkPdf',  'ctkPdf', ctkSet, ctkPdfSet, pdf._var, pdf._var, ctkNorm, projWData['ProjWData'][0] ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf('ctkPhi',  'ctkPhi',  ctkSet, ctkPdfSet, phiPdf,  pdf._var, ctkNorm, projWData['ProjWData'][0]) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf('ctkEven', 'ctkEven', ctkSet, ctkPdfSet, evenPdf, pdf._var, ctkNorm, projWData['ProjWData'][0]) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf('ctkOdd',  'ctkOdd',  ctkSet, ctkPdfSet, oddPdf,  pdf._var, ctkNorm, projWData['ProjWData'][0]) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf('ctkS',    'ctkS',    ctkSet, ctkPdfSet, SPdf,    pdf._var, ctkNorm, projWData['ProjWData'][0]) )
        else :
            totPdfs.append(  RooExplicitNormPdf( 'ctkPdf',  'ctkPdf',  ctkSet, ctkPdfSet, pdf._var, pdf._var, ctkNorm ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf( 'ctkPhi',  'ctkPhi',  ctkSet, ctkPdfSet, phiPdf,  pdf._var, ctkNorm ) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf( 'ctkEven', 'ctkEven', ctkSet, ctkPdfSet, evenPdf, pdf._var, ctkNorm ) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf( 'ctkOdd',  'ctkOdd',  ctkSet, ctkPdfSet, oddPdf,  pdf._var, ctkNorm ) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf( 'ctkS',    'ctkS',    ctkSet, ctkPdfSet, SPdf,    pdf._var, ctkNorm ) )

    if doTimeAnglePlots['ctl'] :
        # cos(theta_l) PDFs
        obsToPlot.append('ctl')
        plotCanvs.append( TCanvas( 'ctlCanv', 'Decay time' ) )
        ctlNorm = ( angles[1].getMax() - angles[1].getMin() ) * dataSet.sumEntries() / float( numBins['ctl'] )
        if projWData :
            totPdfs.append( RooExplicitNormPdf( 'ctlPdf', 'ctlPdf', ctlSet, ctlPdfSet, pdf._var, pdf._var, ctlNorm, projWData['ProjWData'][0] ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf('ctlPhi',  'ctlPhi',  ctlSet, ctlPdfSet, phiPdf,  pdf._var, ctlNorm, projWData['ProjWData'][0]) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf('ctlEven', 'ctlEven', ctlSet, ctlPdfSet, evenPdf, pdf._var, ctlNorm, projWData['ProjWData'][0]) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf('ctlOdd',  'ctlOdd',  ctlSet, ctlPdfSet, oddPdf,  pdf._var, ctlNorm, projWData['ProjWData'][0]) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf('ctlS',    'ctlS',    ctlSet, ctlPdfSet, SPdf,    pdf._var, ctlNorm, projWData['ProjWData'][0]) )
        else :
            totPdfs.append( RooExplicitNormPdf( 'ctlPdf', 'ctlPdf', ctlSet, ctlPdfSet, pdf._var, pdf._var, ctlNorm ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf( 'ctlPhi',  'ctlPhi',  ctlSet, ctlPdfSet, phiPdf,  pdf._var, ctlNorm ) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf( 'ctlEven', 'ctlEven', ctlSet, ctlPdfSet, evenPdf, pdf._var, ctlNorm ) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf( 'ctlOdd',  'ctlOdd',  ctlSet, ctlPdfSet, oddPdf,  pdf._var, ctlNorm ) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf( 'ctlS',    'ctlS',    ctlSet, ctlPdfSet, SPdf,    pdf._var, ctlNorm ) )

    if doTimeAnglePlots['phi'] :
        # phi PDFs
        obsToPlot.append('phi')
        plotCanvs.append( TCanvas( 'phiCanv', 'Decay time' ) )
        phiNorm = ( angles[2].getMax() - angles[2].getMin() ) * dataSet.sumEntries() / float( numBins['phi'] )
        if projWData :
            totPdfs.append( RooExplicitNormPdf( 'phiPdf', 'phiPdf', phiSet, phiPdfSet, pdf._var, pdf._var, phiNorm, projWData['ProjWData'][0] ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf('phiPhi',  'phiPhi',  phiSet, phiPdfSet, phiPdf,  pdf._var, phiNorm, projWData['ProjWData'][0]) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf('phiEven', 'phiEven', phiSet, phiPdfSet, evenPdf, pdf._var, phiNorm, projWData['ProjWData'][0]) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf('phiOdd',  'phiOdd',  phiSet, phiPdfSet, oddPdf,  pdf._var, phiNorm, projWData['ProjWData'][0]) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf('phiS',    'phiS',    phiSet, phiPdfSet, SPdf,    pdf._var, phiNorm, projWData['ProjWData'][0]) )
        else :
            totPdfs.append( RooExplicitNormPdf( 'phiPdf', 'phiPdf',  phiSet, phiPdfSet, pdf._var, pdf._var, phiNorm ) )
            if plotComponents['phi'] :
                phiPdfs.append(  RooExplicitNormPdf( 'phiPhi',  'phiPhi',  phiSet, phiPdfSet, phiPdf,  pdf._var, phiNorm ) )
            if plotComponents['even'] :
                evenPdfs.append( RooExplicitNormPdf( 'phiEven', 'phiEven', phiSet, phiPdfSet, evenPdf, pdf._var, phiNorm ) )
            if plotComponents['odd'] :
                oddPdfs.append(  RooExplicitNormPdf( 'phiOdd',  'phiOdd',  phiSet, phiPdfSet, oddPdf,  pdf._var, phiNorm ) )
            if plotComponents['S'] :
                SPdfs.append(    RooExplicitNormPdf( 'phiS',    'phiS',    phiSet, phiPdfSet, SPdf,    pdf._var, phiNorm ) )

    # plot lifetime and angles
    print 'JvLFit: plotting time and angular distributions'
    from ROOT import TFile
    plotsRFile  = TFile.Open( plotsROOTFile, 'RECREATE' )
    addPDFsOpts = [ ]
    if phiPdfs  : addPDFsOpts.append( dict( LineColor = kRed,         LineWidth = lineWidth, LineStyle = 7 ) )
    if evenPdfs : addPDFsOpts.append( dict( LineColor = kRed,         LineWidth = lineWidth, LineStyle = 9 ) )
    if oddPdfs  : addPDFsOpts.append( dict( LineColor = kGreen + 3,   LineWidth = lineWidth, LineStyle = 7 ) )
    if SPdfs    : addPDFsOpts.append( dict( LineColor = kMagenta + 3, LineWidth = lineWidth, LineStyle = 5 ) )

    for iter, ( obsName, pad, totPdf ) in enumerate( zip(  obsToPlot, plotCanvs, totPdfs ) ) :
        pad.SetLeftMargin(0.18)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.18)
        pad.SetTopMargin(0.05)

        obs = obsSetP2VV[ 1 if obsName == 'ctk' else 2 if obsName == 'ctl' else 3 if obsName == 'phi' else 0 ]
        plotRange = 'Bulk' if obsName == 'timeLin' else ''
        binWidth = ( obs.getMax(plotRange) - obs.getMin(plotRange) ) / float(numBins[obsName]) / ( pi if obsName == 'phi' else 1. )
        plots = plot(  pad, obs, dataSet, totPdf, addPDFs = [ pdfList[iter] for pdfList in [ phiPdfs, evenPdfs, oddPdfs, SPdfs ] if pdfList ]
                     , xTitle       = obsNames[obsName][1]
                     , yTitle       = obsYTitles[obsName] % binWidth, logy = True if obsName == 'timeLog' else False
                     , xTitleOffset = 1.10
                     , yTitleOffset = yTitleOffsets[obsName]
                     , frameOpts    = dict( Range = plotRange, Bins = numBins[obsName], Title = ''       )
                     , dataOpts     = dict( MarkerStyle = markStyle, MarkerSize = 0.7, LineWidth = 3     )
                     , pdfOpts      = dict( LineColor = kBlue, LineWidth = lineWidth, LineStyle = kSolid )
                     , addPDFsOpts  = addPDFsOpts
                     , components   = comps
                    )

        _P2VVPlotStash.append( LHCbText.Clone() )
        _P2VVPlotStash[-1].SetName( 'LHCbLabel%d_%s' % ( iter, obs.GetName() ) )
        _P2VVPlotStash[-1].SetX1NDC( LHCbCoords[obsName][0] )
        _P2VVPlotStash[-1].SetY1NDC( LHCbCoords[obsName][1] )
        _P2VVPlotStash[-1].SetX2NDC( LHCbCoords[obsName][2] )
        _P2VVPlotStash[-1].SetY2NDC( LHCbCoords[obsName][3] )
        _P2VVPlotStash[-1].Draw()

        plotsRFile.Add(pad)
        plotsRFile.Add( _P2VVPlotStash[-1] )

    plotsRFile.Write()

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
    #plot(  timeBBbarCanv, time, dataSet, pdf, addPDFs = [ pdf, pdf ], yScale = ( None, None )
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
    #plot(  timeBCanv, time, dataSet, pdf, yScale = ( None, None )
    #     , frameOpts   = dict( Bins = 50, Title = '', Range = 'Bulk'                                                                )
    #     , dataOpts    = dict( MarkerStyle = kFullCircle, MarkerSize = 0.7, MarkerColor = kBlack, LineColor = kBlack, LineWidth = 3
    #                          , Cut = '%s > 0 && %s == +1' % ( tagCatP2VVOS, iTagOS ) )
    #     , pdfOpts     = dict( list( projWDataBulk.items() ), LineColor = kBlue, LineWidth = 4
    #                          , Slices = [ ( tagCatP2VVOS, 'TagCat1' ), ( iTagOS, 'B' ) ] )
    #    )
    #LHCbText.Draw()
    #timeBCanv.Print(plotsFile)

    #timeBbarCanv = TCanvas( 'timeBbarCanv', 'Lifetime' )
    #plot(  timeBbarCanv, time, dataSet, pdf, yScale = ( None, None )
    #     , frameOpts   = dict( Bins = 50, Title = '', Range = 'Bulk'                                                                )
    #     , dataOpts    = dict( MarkerStyle = kFullCircle, MarkerSize = 0.7, MarkerColor = kBlack, LineColor = kBlack, LineWidth = 3
    #                          , Cut = '%s > 0 && %s == -1' % ( tagCatP2VVOS, iTagOS ) )
    #     , pdfOpts     = dict( list( projWDataBulk.items() ), LineColor = kRed, LineWidth = 4
    #                          , Slices = [ ( tagCatP2VVOS, 'TagCat1' ), ( iTagOS, 'Bbar' ) ] )
    #    )
    #LHCbText.Draw()
    #timeBbarCanv.Print( plotsFile + ')' )

    ## plot angles
    #print 'JvLFit: plotting angular distributions'
    #updateAngPlots = False

    #from ROOT import TFile
    #if not updateAngPlots :
    #    angPlotsFileRead = None
    #    angPlotsFile     = TFile.Open( 'anglePlotsFile.root',       'RECREATE' )
    #else :
    #    angPlotsFileRead = TFile.Open( 'anglePlotsFile.root',       'READ'     )
    #    angPlotsFile     = TFile.Open( 'anglePlotsFileUpdate.root', 'RECREATE' )

    #timeCanv = TCanvas( 'timeCanv', 'Decay time'   )
    #ctkCanv  = TCanvas( 'ctkCanv',  'cos(theta_K)' )
    #ctlCanv  = TCanvas( 'ctlCanv',  'cos(theta_l)' )
    #phiCanv  = TCanvas( 'phiCanv',  'cos(phi)'     )
    #for ( canv, obs, nBins, xTitle, binWidthStr, textXCoords )\
    #        in zip(  [ ctkCanv, ctlCanv, phiCanv ]
    #               , obsSetP2VV[ 1 : 4 ]
    #               , numBins[ 1 : 4 ]
    #               , ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
    #               , ( 'Candidates / %.3f', 'Candidates / %.3f', 'Candidates / (%.2f rad)' )
    #               , ( ( 0.54, 0.69 ), ( 0.32, 0.47 ), ( 0.32, 0.47 ) )
    #              ) :
    #    canv.cd()
    #    canv.SetLeftMargin(0.28)

    #    if not updateAngPlots :
    #        obsFrame = obs.frame( Bins = nBins, Title = obs.GetTitle() )
    #        xAxis = obsFrame.GetXaxis()
    #        yAxis = obsFrame.GetYaxis()
    #        _P2VVPlotStash.append(obsFrame)

    #        obsFrame.SetMinimum(0.)
    #        obsFrame.SetMaximum(2100.)
    #        xAxis.SetTitle(xTitle)
    #        yAxis.SetTitle( binWidthStr % ( ( obs.getMax() - obs.getMin() ) / float(nBins) ) )
    #        yAxis.SetTitleOffset(1.0)

    #        dataSet.plotOn( obsFrame, Name = 'data', MarkerStyle = markStyle, MarkerSize = markSize, MarkerColor = markColor
    #                    , LineColor = markColor, LineWidth = markLineWidth )
    #        pdf.plotOn( obsFrame, Name = 'pdf', LineColor = kBlue, LineWidth = lineWidth, **projWData )
    #        obsFrame.drawAfter( 'pdf', 'data' )

    #    else :
    #        obsFrame = None
    #        keysList = angPlotsFileRead.GetListOfKeys()
    #        for key in keysList :
    #            if key.GetName().startswith( 'frame_' + obs.GetName() ) : obsFrame = angPlotsFileRead.Get( key.GetName() )

    #        if obsFrame :
    #            angPlotsFile.Add(obsFrame)
    #            pdf.plotOn( obsFrame, Name = 'pdf1', LineColor = kRed, LineStyle = 7, LineWidth = lineWidth, **projWData )
    #            obsFrame.drawAfter( 'pdf1', 'data' )
    #            obsFrame.drawAfter( 'pdf',  'pdf1' )

    #    if obsFrame : obsFrame.Draw()

    #    _P2VVPlotStash.append( LHCbText.Clone() )
    #    _P2VVPlotStash[-1].SetX1NDC( textXCoords[0] )
    #    _P2VVPlotStash[-1].SetY1NDC(0.82)
    #    _P2VVPlotStash[-1].SetX2NDC( textXCoords[1] )
    #    _P2VVPlotStash[-1].SetY2NDC(0.91)
    #    #_P2VVPlotStash[-1].SetBorderSize(1)
    #    _P2VVPlotStash[-1].Draw()

    #    angPlotsFile.Add(canv)

    #angPlotsFile.Write()

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
    #    plot(  pad, time, dataSet, pdf, yTitle = yTitle, yScale = yScale, logy = logY
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
    #    plot(  pad, time, dataSet, pdf, yTitle = yTitle, logy = logY
    #         , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                )
    #         , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                     )
    #         , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
    #         , components = comps
    #        )

    ## plot angles
    #print 'JvLFit: plotting angular distributions'
    #if plotAnglesNoEff and pdfConfig['sFit'] and pdfConfig['multiplyByTimeEff'] not in [ 'all', 'signal' ]\
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
    #    plot(  pad, obs, dataSet, pdf, addPDFs = addPDFs, xTitle = xTitle, yTitle = yTitle
    #         , frameOpts   = dict( Bins = nBins, Title = plotTitle                                                )
    #         , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts                    )
    #         , pdfOpts     = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
    #         , addPDFsOpts = [ dict( list( projWData.items() ), LineColor = kRed, LineWidth = lineWidth, **pdfCuts ) ]
    #         , components  = comps
    #        )

    ## plot background time
    #print 'JvLFit: plotting background lifetime distribution'
    #bkgTimeCanv = TCanvas( 'bkgTimeCanv', 'Background Lifetime' )
    #for ( pad, data, plotTitle, logY )\
    #      in zip(  bkgTimeCanv.pads( 2, 2 )
    #             , 2 * [ pdfBuild['bkgRangeData'], pdfBuild['bkgSWeightData'] ]
    #             , [  time.GetTitle() + ' - mass side bands - linear'
    #                , time.GetTitle() + ' - mass S-weights - linear'
    #                , time.GetTitle() + ' - mass side bands - logarithmic'
    #                , time.GetTitle() + ' - mass S-weights - logarithmic'
    #               ]
    #             , 2 * [ False ] + 2 * [ True ]
    #            ) :
    #    plot(  pad, time, data, pdfBuild['bkg_t'], logy = logY
    #         , frameOpts  = dict( Title = plotTitle, Range = 'Bulk', Bins = pdfConfig['numTimeBins'] )
    #         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
    #         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3  )
    #        )

if deltaSCanv :
    gStyle.SetEndErrorSize(4)
    deltaSCanv.Update()
    plotsCanvs.append(deltaSCanv)

if plotCanvs :
    # print canvases to file
    for iter, canv in enumerate(plotCanvs) :
        canv.Print( plotsFile + ( '(' if iter == 0 and len(plotCanvs) > 1 else ')' if iter == len(plotCanvs) - 1 and len(plotCanvs) > 1 else '' ) )
    if plotsRFile : plotsRFile.Close()


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
     , 'timeResSigmaSF'  : ( '#DeltaNLL #sigma_{t} s.f.',              '#sigma_{t} s.f.',               1.27,     1.63,    1, 0.001, 0.01 )
    }

    # check DNLL parameters
    phiCPPar = False
    for par in dllPars :
        assert par[0] in wsPars, 'JvLFit - ERROR: unknown DLL parameter: "%s"' % par[0]
        assert par[0] in ws,     'JvLFit - ERROR: DLL parameter "%s" does not exist in work space' % par[0]
        if par[0] == 'phiCP' : phiCPPar = True

    # float/fix values of some parameters
    #if constLambdaCP :
    #    pdfBuild['lambdaCP'].setConstant('lambdaCPSq') if pdfConfig['lambdaCPParam'] == 'lambSqPhi'\
    #        else pdfBuild['lambdaCP'].setConstant('lambdaCP')
    #if constAvgCEvenOdd :
    #    for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')

    #if 'sig_ATagBBbar' in ws : ws['sig_ATagBBbar'].setConstant()
    #for bin in range(5) :
    #    if 'sig_ATagBBbar_bin%d' % bin in ws : ws[ 'sig_ATagBBbar_bin%d' % bin ].setConstant()
    #pdfBuild['tagCatsOS'].setConstant('.*')
    #pdfBuild['tagCatsSS'].setConstant('.*')
    ##pdfBuild['lifetimeParams'].setConstant('dM|Gamma')
    #pdfBuild['timeResModel'].setConstant('.*')
    #pdfBuild['signalBMass'].setConstant('.*')
    #if not pdfConfig['sFit'] :
    #    pdfBuild['backgroundBMass'].setConstant('.*')
    #    pdfBuild['backgroundTime'].setConstant('.*')

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
            parabola.plotOn( parFrame, RooFit.LineColor(RooFit.kBlack), RooFit.Precision(0.0001) )

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


tagData = dataSet

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
DelmVal  = 17.768
etaOSVal = 0.380 #0.382 #0.392
p0OSVal  = 0.380 #0.382 #0.392
p1OSVal  = 1.000 #0.981 #1.000
etaSSVal = 0.437 #0.438 #0.35
p0SSVal  = 0.437 #0.439 #0.35
p1SSVal  = 1.000  #1.04  #1.0
print
print 'Delta m_s = %.3f' % DelmVal
print 'w_OS = %.3f + %.3f * (eta - %.3f)' % ( p0OSVal, p1OSVal, etaOSVal )
print 'w_SS = %.3f + %.3f * (eta - %.3f)' % ( p0SSVal, p1SSVal, etaSSVal )
for varSet in tagData :
  weight = tagData.weight()
  sums['numEv'] += weight

  sigmat = 1.45 * varSet.getRealValue('sigmat')
  dilResSq = exp( -DelmVal**2 * sigmat**2 )
  dilRes   = sqrt(dilResSq)
  sums['dilResAll']  += weight * dilRes
  sums['dilRes2All'] += weight * dilResSq

  if varSet.getCatIndex('tagCatP2VVOS') != 0 :
    etaOS  = varSet.getRealValue('tagomega_os_cb')
    wTagOS = p0OSVal + p1OSVal * ( etaOS - etaOSVal )
    dilOS  = 1. - 2. * wTagOS

    sums['numOS']  += weight
    sums['etaOS']  += weight * etaOS
    sums['wOS']    += weight * wTagOS
    sums['dilOS']  += weight * dilOS
    sums['dil2OS'] += weight * dilOS**2

    if varSet.getCatIndex('tagCatP2VVSS') == 0 :
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

  if varSet.getCatIndex('tagCatP2VVSS') != 0 :
    etaSS  = varSet.getRealValue('tagomega_ss_nn')
    wTagSS = p0SSVal + p1SSVal * ( etaSS - etaSSVal )
    dilSS  = 1. - 2. * wTagSS

    sums['numSS']  += weight
    sums['etaSS']  += weight * etaSS
    sums['wSS']    += weight * wTagSS
    sums['dilSS']  += weight * dilSS
    sums['dil2SS'] += weight * dilSS**2

    if varSet.getCatIndex('tagCatP2VVOS') == 0 :
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

  if varSet.getCatIndex('tagCatP2VVOS') != 0 and varSet.getCatIndex('tagCatP2VVSS') != 0 :
    dilSign = +1. if varSet.getCatIndex('iTagOS') == varSet.getCatIndex('iTagSS') else -1.
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

if sums['numOS'] :       print 'OS    <eta>: %.4f (%.4f)'   % ( sums['etaOS']   / sums['numOS'],   sums['etaOSExcl']   / sums['numOSExcl']   )
if sums['numSS'] :       print 'SS    <eta>: %.4f (%.4f)'   % ( sums['etaSS']   / sums['numSS'],   sums['etaSSExcl']   / sums['numSSExcl']   )
if sums['numCombExcl'] : print 'Comb. <eta>: %.4f (%.4f)\n' % ( sums['etaComb'] / sums['numComb'], sums['etaCombExcl'] / sums['numCombExcl'] )

if sums['numOS'] :       print 'OS    <w>: %.4f (%.4f)'   % ( sums['wOS']   / sums['numOS'],   sums['wOSExcl']   / sums['numOSExcl']   )
if sums['numSS'] :       print 'SS    <w>: %.4f (%.4f)'   % ( sums['wSS']   / sums['numSS'],   sums['wSSExcl']   / sums['numSSExcl']   )
if sums['numCombExcl'] : print 'Comb. <w>: %.4f (%.4f)\n' % ( sums['wComb'] / sums['numComb'], sums['wCombExcl'] / sums['numCombExcl'] )

if sums['numOS'] :       print 'OS    <dil>:    %.4f (%.4f)'   % ( sums['dilOS']   / sums['numOS'],   sums['dilOSExcl']   / sums['numOSExcl']   )
if sums['numSS'] :       print 'SS    <dil>:    %.4f (%.4f)'   % ( sums['dilSS']   / sums['numSS'],   sums['dilSSExcl']   / sums['numSSExcl']   )
if sums['numCombExcl'] : print 'Comb. <dil>:    %.4f (%.4f)'   % ( sums['dilComb'] / sums['numComb'], sums['dilCombExcl'] / sums['numCombExcl'] )
print 'All   <dilRes>: %.4f\n'        % ( sums['dilResAll'] / sums['numEv'] )

if sums['numOS'] :       print 'OS    <dil2>:    %.4f (%.4f)'   % ( sums['dil2OS']   / sums['numOS'],   sums['dil2OSExcl']   / sums['numOSExcl']   )
if sums['numSS'] :       print 'SS    <dil2>:    %.4f (%.4f)'   % ( sums['dil2SS']   / sums['numSS'],   sums['dil2SSExcl']   / sums['numSSExcl']   )
if sums['numCombExcl'] : print 'Comb. <dil2>:    %.4f (%.4f)'   % ( sums['dil2Comb'] / sums['numComb'], sums['dil2CombExcl'] / sums['numCombExcl'] )
if sums['numCombExcl'] : print 'Comb. <dilRes2>: %.4f (%.4f)'   % ( sums['dilRes2Comb'] / sums['numComb'], sums['dilRes2CombExcl'] / sums['numCombExcl'] )
if sums['numCombExcl'] : print 'Comb. <dilTot2>: %.4f (%.4f)'   % ( sums['dilTot2Comb'] / sums['numComb'], sums['dilTot2CombExcl'] / sums['numCombExcl'] )
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

from math import log
print 'effective decay time resolution: %.4f\n' % (  sqrt( -log( sums['dilRes2All'] / sums['numEv'] ) / DelmVal**2 ) )

nEv = 0.
nBB = 0.
nBbarBbar = 0.
nBbarB = 0.
nBBbar = 0.
avD_OS = 0.
avD_SS = 0.
avDD = 0.
for argSet in tagData :
    if argSet.getCatIndex('tagCatP2VVOS') == 0 or argSet.getCatIndex('tagCatP2VVSS') == 0 : continue

    nEv += tagData.weight()
    if argSet.getCatIndex('iTagOS') == +1 and argSet.getCatIndex('iTagSS') == +1 : nBB += tagData.weight()
    if argSet.getCatIndex('iTagOS') == -1 and argSet.getCatIndex('iTagSS') == -1 : nBbarBbar += tagData.weight()
    if argSet.getCatIndex('iTagOS') == -1 and argSet.getCatIndex('iTagSS') == +1 : nBbarB += tagData.weight()
    if argSet.getCatIndex('iTagOS') == +1 and argSet.getCatIndex('iTagSS') == -1 : nBBbar += tagData.weight()

    D_OS = 1. - 2. * argSet.getRealValue('tagomega_os_cb')
    D_SS = 1. - 2. * argSet.getRealValue('tagomega_ss_nn')

    avD_OS += D_OS * tagData.weight()
    avD_SS += D_SS * tagData.weight()
    avDD   += D_OS * D_SS * tagData.weight()

avD_OS /= nEv
avD_SS /= nEv
avDD   /= nEv

Atags = ( nBB + nBbarBbar - nBbarB - nBBbar ) / nEv
print avDD, avD_OS * avD_SS, Atags

########################################################################
#plotsFile = 'timeAngleData_PID_NextVert_GL0020_KKPT750_DTF6.ps'
#
#plotVars = [ time ] + angles
#plots = [ time.frame(30), angles[0].frame(20), angles[1].frame(40), angles[2].frame(20) ]
#plots[0].SetXTitle('Decay time [ps]')
#for iter, plot in enumerate( plots[ 1 : ] ) : plot.SetXTitle( angleNames[iter][1] )
#
#from ROOT import kFullCircle, kBlue, kRed
#for plot in plots :
#    bkgData.plotOn( plot, MarkerStyle = kFullCircle, MarkerSize = 0.6, MarkerColor = kRed,  LineWidth = 3, LineColor = kRed  )
#    sigData.plotOn( plot, MarkerStyle = kFullCircle, MarkerSize = 0.6, MarkerColor = kBlue, LineWidth = 3, LineColor = kBlue )
#
#plots[0].SetMinimum(0.1)
##plots[0].SetMaximum(4.e4)
##plots[1].SetMinimum(0.)
##plots[1].SetMaximum(1000.)
##plots[2].SetMinimum(0.)
##plots[2].SetMaximum(700.)
##plots[3].SetMinimum(0.)
##plots[3].SetMaximum(1000.)
#
#from ROOT import TCanvas
#motherCanv = TCanvas()
#motherCanv.Divide( 2, 2 )
#canvs = [ motherCanv.cd(iter) for iter in range( 1, 5 ) ]
##canvs = [ TCanvas( var.GetName() ) for var in plotVars ]
#canvs[0].SetLogy()
#
#for iter, ( plot, canv ) in enumerate( zip( plots, canvs ) ) :
#    canv.cd()
#    plot.Draw()
#    #canv.Print( plotsFile + ( '(' if iter == 0 else ')' if iter == len(plots) - 1 else '' ) )
#motherCanv.Print(plotsFile)

#execfile('plotNLL.py')
