###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
pdfConfig['dataSample'] = '' #'sel_cleantail==1' # ( None, 100260, '' )  # '' / 'Summer2011' / 'runNumber % 2 == 1'
pdfConfig['selection']  = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
generateData            = False
doFit                   = True #'NLL'
makeObservablePlots     = False
makeKKMassPlots         = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
corrSFitErr             = 'sumWeight' # '' / 'sumWeight' / ( 0.887, [ 0.566, 0.863, 0.956, 0.948, 0.855, 0.662 ] ) / 'matrix'
randomParVals           = ( ) #( 0.2, 12345 )
pdfConfig['blind']      = {#  'phiCP'  : ( 'UnblindUniform', 'BsPhis2013EPS',  0.2  )
                           #, 'dGamma' : ( 'UnblindUniform', 'BsDGs2013EPS',   0.02 )
                          }

plotsFile     = 'plots/2011Data_SFit.ps'
plotsROOTFile = '2011Data_SFit_plots.root'
parFileIn     = '2011DataFitValues.par'
parFileOut    = ''

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'paper2012_SFit.root' if pdfConfig['SFit'] else 'paper2012_CFit.root'
elif pdfConfig['SFit'] :
    dataSetName = 'JpsiKK_sigSWeight'
    dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_narrowKKMass_noKKMassBins_noTriggerBins_2TagCats.root'
else :
    dataSetName = 'JpsiKK'
    dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_4KKMassBins_2TagCats.root'

MinosPars = [#  'phiCP', 'lambdaCP'
             #, 'AparPhase', 'AperpPhase'
             #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3'#, 'ASOddPhase_bin4', 'ASOddPhase_bin5'
             #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3'#,        'f_S_bin4',        'f_S_bin5'
             #, 'timeResSigmaSF'
            ]
dllPars = [ ] # [ ( 'ImApar', True, True, True ) ] / [ ( 'phiCP', True, True, True ) ]

# fit options
fitRange = ''
fitOpts = dict(  NumCPU    = 6
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

# plot options
from ROOT import gStyle, kBlack, kBlue, kRed, kGreen, kMagenta, kSolid, kDashed, kFullCircle, kFullSquare, kFullDotLarge
from P2VV.Load import RooFitOutput, LHCbStyle
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

# PDF options
pdfConfig['transversityAngles'] = False
pdfConfig['angularRanges']      = dict( ctheta = [  ( 'ctlBin0', -1.,   -0.70 )
                                                  , ( 'ctlBin1', -0.70, -0.25 )
                                                  , ( 'ctlBin2', -0.25, +0.25 )
                                                  , ( 'ctlBin3', +0.25, +0.70 )
                                                  , ( 'ctlBin4', +0.70, +1.   )
                                                 ]
                                      )
pdfConfig['sigMassModel']         = '' # 'boxFixed'
pdfConfig['bkgMassModel']         = '' # 'linearConstant'
pdfConfig['bkgAnglePdf']          = 'hybrid'
pdfConfig['sigTaggingPdf']        = 'tagUntag'
pdfConfig['bkgTaggingPdf']        = 'tagUntagRelative'  # 'tagUntagRelative' # 'tagCatsRelative'
pdfConfig['multiplyByTagPdf']     = False
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['timeEffType']          = 'paper2012_multi' # 'paper2012' # 'paper2012_multi' # 'HLT1Unbiased'
pdfConfig['multiplyByAngEff']     = 'weights'  # 'weights' # 'basis012' # 'basisSig4'
pdfConfig['parameterizeKKMass']   = ''  # 'simultaneous'
pdfConfig['ambiguityParameters']  = False
pdfConfig['KKMassBinBounds']      = [ 1008., 1032. ] # [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 988., 1020. - 12., 1020., 1020. + 12., 1050. ]
pdfConfig['SWaveAmplitudeValues'] = ( )
pdfConfig['CSPValues']            = [ 0.498 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.959, 0.770, 0.824, 0.968 ] # [ 0.959, 0.498, 0.968 ]

KKMassPars = pdfConfig['obsDict']['KKMass']
pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                  , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )

pdfConfig['sameSideTagging']    = False
pdfConfig['conditionalTagging'] = True
pdfConfig['continuousEstWTag']  = True
pdfConfig['numEstWTagBins']     = 50
pdfConfig['constrainTagging']   = 'constrain'  # 'constrain'

pdfConfig['timeResType']           = 'eventNoMean' # 'eventDoubleGaussConstantFixedMean' # 'event' # 'eventNoMean'
pdfConfig['numTimeResBins']        = 40
pdfConfig['constrainTimeResScale'] = 'constrain' # 'fixed'  # 'constrain'

pdfConfig['numEvents'] = 54755
pdfConfig['signalFraction'] = 0.504

pdfConfig['amplitudeParam'] = 'phasesSWaveFrac' # 'bank' # 'phasesSWaveFrac'
pdfConfig['ASParam']        = 'deltaPerp'  # 'deltaPerp'
pdfConfig['AparParam']      = 'phase' # 'Mag2ReIm' # 'phase'

pdfConfig['constrainDeltaM'] = 'constrain'  # 'constrain' # fixed

pdfConfig['lambdaCPParam'] = 'lambPhi' # 'lambPhi_CPVDecay_PSWaves'  # 'lambPhi'

from P2VV.Imports import extConstraintValues
extConstraintValues.setVal( 'DM',      (  17.63, 0.11 ) )

fastFit           = False
manualTagCatBins  = False
constTagCatCoefs  = True
constAvgCEvenOdd  = True
constWTagAsyms    = 'P1'
constCSP          = True
constAmplitudes   = False
constLambdaCP     = ''  # 'lamb'
equalAbsLambdaCPs = False

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

if generateData or manualTagCatBins :
    pdfConfig['tagCatsOS'] = [  ( 'Untagged',  0, 0.500001 )
                              , ( 'TagCat1',   1, 0.499999 )
                             ]
    pdfConfig['tagCatsSS'] = [  ( 'Untagged',  0, 0.500001 )
                              , ( 'TagCat1',   1, 0.499999 )
                             ]

pdfConfig['timeEffHistFiles'] = dict(  file      = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
                                     , hlt1UB    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                     , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                    )
pdfConfig['angEffMomsFiles']    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'\
                                    if pdfConfig['transversityAngles'] else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

if pdfConfig['transversityAngles'] :
    pdfConfig['angleNames'] = (  ( 'trcospsi',   'cos#kern[0.1]{#psi_{tr}}'   )
                               , ( 'trcostheta', 'cos#kern[0.1]{#theta_{tr}}' )
                               , ( 'trphi',      '#varphi_{tr} [rad]'         )
                              )
else :
    pdfConfig['angleNames'] = (  ( 'helcosthetaK', 'cos#kern[0.1]{#theta_{K}}'   )
                               , ( 'helcosthetaL', 'cos#kern[0.1]{#theta_{#mu}}' )
                               , ( 'helphi',       '#varphi_{h} [rad]'           )
                              )
obsNames = {  'timeLin' : ( 'time', 'Decay time [ps]' )
            , 'timeLog' : ( 'time', 'Decay time [ps]' )
            , 'ctk'     : pdfConfig['angleNames'][0]
            , 'ctl'     : pdfConfig['angleNames'][1]
            , 'phi'     : pdfConfig['angleNames'][2]
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
plotComponents   = { 'phi' : False, 'even' : True, 'odd' : True, 'S' : True }
numBins = { 'timeLin' : 40, 'timeLog' : 50, 'ctk' : 30, 'ctl' : 30, 'phi' : 30 }
pdfConfig['numTimeBins'] = 30
numAngleBins = ( 20, 20, 20 )
pdfConfig['numAngleBins'] = ( 5, 7, 9 )


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
    pdfConfig['dataSet'] = dataSet

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
mumuMass   = pdfBuild['observables']['mumuMass']
KKMass     = pdfBuild['observables']['KKMass']
estWTagOS  = pdfBuild['observables']['wTagOS']
timeRes    = pdfBuild['observables']['timeRes']

if not pdfConfig['SFit'] : obsSetP2VV.append(BMass)

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
    nEvents = int( pdfConfig['numEvents'] * ( pdfConfig['signalFraction'] if pdfConfig['SFit'] else 1. ) )
    print 'JvLFit: generating %d events' % nEvents
    dataSet = pdf.generate( obsSetP2VV, nEvents )

    # additional observables
    if not pdfConfig['transversityAngles'] :
        from P2VV.Utilities.DataHandling import addTransversityAngles
        addTransversityAngles( dataSet, 'trcospsi',          'trcostheta',        'trphi'
                                      , angles[0].GetName(), angles[1].GetName(), angles[2].GetName() )

    # write data to file
    from P2VV.Utilities.DataHandling import writeData
    writeData( dataSetFile, dataSetName, dataSet )

elif pdfConfig['SFit'] :
    if corrSFitErr == 'sumWeight'\
            or ( type(corrSFitErr) != str and hasattr( corrSFitErr, '__iter__' ) and hasattr( corrSFitErr, '__getitem__' ) ) :
        from P2VV.Utilities.DataHandling import correctSWeights
        fitData = correctSWeights( dataSet, 'N_cbkgMass_sw'
                                  , 'KKMassCat' if pdfConfig['parameterizeKKMass'] == 'simultaneous' else ''
                                  , CorrectionFactors = None if corrSFitErr == 'sumWeight' else corrSFitErr )

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

if pdfConfig['parameterizeKKMass'] == 'amplitudes' :
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

if pdfConfig['lambdaCPParam'].startswith('lambPhi_CPVDecay') :
    if equalAbsLambdaCPs : pdfBuild['lambdaCP'].setConstant('rhoCP_A.*')
    else :                 pdfBuild['lambdaCP'].setConstant('rhoCP_m')

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
        if pdfConfig['SFit'] :
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

    from P2VV.Imports import parNames, parValues4KKBins as parValues
    print 'JvLFit: parameters:'
    fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    fitResult.covarianceMatrix().Print()
    fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

else :
    fitResult = None

if parFileOut :
    # write parameters to file
    pdfConfig.getParametersFromPdf( pdf, fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut )
