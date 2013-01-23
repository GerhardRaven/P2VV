###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
readMoments = False
multPdfEff  = False
makePlots   = False
transAngles = False
normPDF     = True
tResModel   = ''
trigger     = ''
timeInt     = False

momentsFile = '%s_UB_UT_trueTime_BkgCat050_KK30' % ( 'trans' if transAngles else 'hel' )
plotsFile   = '%s_UB_UT_trueTime_BkgCat050_KK30' % ( 'trans' if transAngles else 'hel' ) + '.ps'

dataSetName = 'DecayTree'
dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20121010.root'
#dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_DGs0_MC11a_ntupleB_for_fitting_20121119.root'
#dataSetFile = '/data/bfys/jleerdam/Bs2Jpsiphi/Bs2JpsiPhi_PHSP_MC11a_ntupleB_for_fitting_20121031.root'
#dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/angRes/angRes.root'

# transversity amplitudes
A0Mag2Val    = 0.60
AperpMag2Val = 0.16
AparMag2Val  = 1. - A0Mag2Val - AperpMag2Val

A0PhVal      =  0.
AperpPhVal   = -0.17
AparPhVal    =  2.50

# CP violation parameters
phiCPVal      = -0.04

# B lifetime parameters
GammaVal  = 0.679
dGammaVal = 0.060
dMVal     = 17.8
tResSigma = 0.045

# plot options
if transAngles : angleNames = ( 'cos(#psi_{tr})',  'cos(#theta_{tr})', '#phi_{tr}' )
else           : angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{l})',  '#phi'      )
numBins         = ( 60, 60, 60, 60 )
lineWidth       = 2
markStyle       = 8
markSize        = 0.4


###########################################################################################################################################
## create variables and read data ##
####################################

# import RooFit wrappers
from P2VVLoad import RooFitOutput

# workspace
from RooFitWrappers import RooObject
ws = RooObject(workspace = 'ws').ws()

# angular functions
if transAngles :
    from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )
    #angleFuncs = AngleFuncs( cpsi = 'ctkRecoLHCb', ctheta = 'ctlRecoLHCb', phi = 'phiRecoLHCb' )
    #angleFuncs = AngleFuncs( cpsi = 'ctkTrueCalc', ctheta = 'ctlTrueCalc', phi = 'phiTrueCalc' )

# variables in PDF
from RooFitWrappers import RealVar, Category
time     = RealVar(  'time',     Title = 'Decay time',      Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. ) )
trueTime = RealVar(  'truetime', Title = 'True decay time', Unit = 'ps', Observable = True, Value = 0.,  MinMax = ( 0.,  20. ) )
iTag     = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'Untagged' : 0 } )
angles   = [ angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] ]

obsSet = [ time if tResModel in [ 'Gauss', '3Gauss' ] else trueTime ] + angles

# read ntuple
bkgcatCut      = '(bkgcat == 0 || bkgcat == 50)'
trackChiSqCuts = 'muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4. && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
massCuts       = 'mass > 5200. && mass < 5550. && mdau1 > 3030. && mdau1 < 3150. && mdau2 > 990. && mdau2 < 1050.'
timeCuts       = 'time > 0.3 && time < 14. && sigmat < 0.12'
tagCuts        = '(tagdecision == 0 || tagdecision == -1 || tagdecision == +1)'

from P2VVGeneralUtils import readData
cuts = bkgcatCut + ' && ' + trackChiSqCuts + ' && ' + massCuts + ' && ' + timeCuts + ' && ' + tagCuts
if trigger == 'ExclBiased' :
    cuts  = 'sel == 1 && hlt1_excl_biased_dec == 1 && hlt2_biased == 1 && ' + cuts
    data = readData(  dataSetFile, dataSetName = dataSetName, NTuple = True, observables = obsSet, Rename = 'DecayTreeB'
                    , ntupleCuts = cuts )

elif trigger == 'Unbiased' :
    cuts = 'sel == 1 && hlt1_unbiased_dec == 1 && hlt2_biased == 1 && ' + cuts
    data = readData(  dataSetFile, dataSetName = dataSetName, NTuple = True, observables = obsSet, Rename = 'DecayTreeUB'
                    , ntupleCuts = cuts )

else :
    cuts = 'sel == 1 && (hlt1_unbiased_dec == 1 || hlt1_biased == 1) && hlt2_biased == 1 && ' + cuts
    data = readData(  dataSetFile, dataSetName = dataSetName, NTuple = True, observables = obsSet, Rename = 'DecayTreeUBB'
                    , ntupleCuts = cuts )


###########################################################################################################################################
## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
#####################################################################

if normPDF :
    # transversity amplitudes
    from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
    amplitudes = Amplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                            , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                            , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                            , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                            , ReAS    = 0.
                            , ImAS    = 0.
                           )

    # B lifetime
    from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
    lifetimeParams = LifetimeParams( Gamma = GammaVal, dGamma = dGammaVal, dM = dMVal )

    tResArgs = { }
    if tResModel == 'Gauss' :
        from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
        tResArgs['time']         = time
        tResArgs['timeResSigma'] = tResSigma
    elif tResModel == '3Gauss' :
        from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
        tResArgs['time'] = time
    else :
        from P2VVParameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
        tResArgs['time'] = trueTime
    timeResModel = TimeResolution( **tResArgs )

    # CP violation parameters
    from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
    lambdaCP = CPParam( lambdaCPSq = 1., phiCP = phiCPVal )

    # tagging parameters
    from P2VVParameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
    taggingParams = TaggingParams()

    # coefficients for time functions
    from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
    timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp' ] )

    # build signal PDF
    args = dict(  time                   = time if tResModel in [ 'Gauss', '3Gauss' ] else trueTime
                , iTag                   = iTag
                , tau                    = lifetimeParams['MeanLifetime']
                , dGamma                 = lifetimeParams['dGamma']
                , dm                     = lifetimeParams['dM']
                , dilution               = taggingParams['dilution']
                , ADilWTag               = taggingParams['ADilWTag']
                , avgCEven               = taggingParams['avgCEven']
                , avgCOdd                = taggingParams['avgCOdd']
                , coshCoef               = timeBasisCoefs['cosh']
                , sinhCoef               = timeBasisCoefs['sinh']
                , cosCoef                = timeBasisCoefs['cos']
                , sinCoef                = timeBasisCoefs['sin']
                , resolutionModel        = timeResModel['model']
               )

    from RooFitWrappers import BTagDecay
    pdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

    intSet  = [ ]
    normSet = angleFuncs.angles.values()
    if timeInt :
        # integrate over time, normalize over time + angles
        intSet.append( time if tResModel in [ 'Gauss', '3Gauss' ] else trueTime )
        normSet.append( intSet[-1] )


###########################################################################################################################################
## compute angular efficiency moments and multiply PDF with efficiency function ##
##################################################################################

# print PDF, data, variables and parameters
print '\nData set:'
data.Print()

if normPDF :
    print '\nPDF:'
    pdf.Print()
    print '\nIntegration set for PDF:'
    for var in intSet : print var.GetName(),
    print '\n\nNormalization set for PDF:'
    for var in normSet : print var.GetName(),
    print '\n\nLifetime resolution model:'
    timeResModel['model'].Print()
    print '\nVariables in PDF:'
    for var in pdf.getObservables(data) : var.Print()
    print '\nParameters in PDF:'
    for par in pdf.getParameters(data) : par.Print()
    print

else :
    print '\nUsing uniform angular PDF'
    print

# moments builder with angular functions from physics PDF
from P2VVGeneralUtils import RealMomentsBuilder
if normPDF :
    from RooFitWrappers import RealEffMoment
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( func, 1, pdf, intSet, normSet )\
                                                  for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func
                                                )
                                    )

else :
    from RooFitWrappers import RealMoment
    physMoments = RealMomentsBuilder( Moments = ( RealMoment( func, 1 )\
                                                  for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func
                                                )
                                    )


# moments builder with angular basis functions
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(4) for YIndex0 in range(4) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(5) for YIndex0 in range(5) for YIndex1 in range( -YIndex0, YIndex0 + 1 )\
#             if PIndex == 4 or YIndex0 == 4 ]
#indices = [ ( PIndex, 2, YIndex1 ) for PIndex in range(40) for YIndex1 in [ +1, -1 ] ]

basisMoments = RealMomentsBuilder()
if normPDF :
    basisMoments.appendPYList( angleFuncs.angles, indices, PDF = pdf, IntSet = intSet, NormSet = normSet )
else :
    basisMoments.appendPYList( angleFuncs.angles, indices )

PDFInt = 1. if normPDF else 8. * pi
if readMoments :
    # read moments from file
    physMoments.read(  momentsFile + '_Phys'  )
    basisMoments.read( momentsFile + '_Basis' )
else :
    # compute moments from data set
    physMoments.compute(data)
    basisMoments.compute(data)

    physMoments.write(  momentsFile + '_Phys',  Scale = PDFInt / 16. / sqrt(pi) )
    basisMoments.write( momentsFile + '_Basis', Scale = PDFInt                  )

# print moments to screen
physMoments.Print(  Scale = PDFInt / 16. / sqrt(pi)                       )
basisMoments.Print( Scale = PDFInt /  2. / sqrt(pi)                       )
basisMoments.Print( Scale = PDFInt /  2. / sqrt(pi), MinSignificance = 5. )


###########################################################################################################################################
## build efficiency function ##
###############################

effTerms = basisMoments.buildPDFTerms()
effFunc = effTerms.buildAddition('effFunc')

from ROOT import RooArgSet
ctkSet = RooArgSet( angles[1], angles[2] )
ctlSet = RooArgSet( angles[0], angles[2] )
phiSet = RooArgSet( angles[0], angles[1] )
effFuncCtk = effFunc.createIntegral( ctkSet, RooArgSet() )
effFuncCtl = effFunc.createIntegral( ctlSet, RooArgSet() )
effFuncPhi = effFunc.createIntegral( phiSet, RooArgSet() )


###########################################################################################################################################
## multiply PDF with angular efficiency ##
##########################################

if multPdfEff and normPDF :
  effPdf = basisMoments * pdf

  basisMomentsSignif = RealMomentsBuilder()
  basisMomentsSignif.appendPYList( angleFuncs.angles, [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if not transAngles \
                                                 else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
                                 )
  basisMomentsSignif.read(momentsFile + '_Basis')
  basisMomentsSignif.Print( Scale = 1. / 2. / sqrt(pi) )

  effSignifPdf = basisMomentsSignif.multiplyPDFWithEff( pdf, Name = 'sig_t_angles_tagCat_iTag_x_EffSignif', EffName = 'effSignif' )


###########################################################################################################################################
## make some plots ##
#####################

if multPdfEff and makePlots and normPDF :
    # import plotting tools
    from P2VVLoad import LHCbStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen

    # plot efficiency
    effCanv = TCanvas( 'effCanv', 'Efficiency' )
    for ( pad, obs, func, norm )\
            in zip(  effCanv.pads( 2, 2 )
                   , obsSet[ 1 : 5 ]
                   , [ effFuncCtk, effFuncCtl, effFuncPhi ]
                   , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
                  ) :
        plot(  pad, obs, None, func, addPDFs = [ effFunc ]
             , yScale      = ( 0.85, 1.15 )
             , pdfOpts     = dict( LineColor = kBlue, LineWidth = 2, Normalization = norm )
             , addPDFsOpts = [ dict( LineColor = kRed,  LineWidth = 2 ) ]
            )

    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, logY )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSet[ : 5 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSet[ : 5 ] ]
                   , ( '', ) + angleNames
                   , ( True, False, False, False )
                  ) :
        plot(  pad, obs, data, effSignifPdf, addPDFs = [ pdf, effPdf ], xTitle = xTitle, logy = logY
             , frameOpts   = dict( Bins = nBins, Title = plotTitle                )
             , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize )
             , pdfOpts     = dict( LineColor = kBlue, LineWidth = lineWidth       )
             , addPDFsOpts = [  dict( LineColor = kRed,       LineWidth = lineWidth )
                              , dict( LineColor = kGreen + 2, LineWidth = lineWidth )
                             ]
            )

    # print canvases to file
    effCanv.Print( plotsFile + '(' )
    timeAnglesCanv.Print( plotsFile + ')' )
