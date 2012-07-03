###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
readMoments = False
multPdfEff  = True
makePlots   = True
transAngles = False
tResModel   = ''
trigger     = ''

momentsFile = 'effMoments' + ( 'Trans' if transAngles else 'Hel' )
plotsFile   = 'effMoments' + ( 'Trans' if transAngles else 'Hel' ) + '.ps'

dataSetName = 'DecayTree'
dataSetFile = '/data/bfys/jleerdam/Bs2Jpsiphi/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20120606.root'

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
from RooFitWrappers import *
from P2VVLoad import RooFitOutput

# workspace
ws = RooObject(workspace = 'ws').ws()

# angular functions
if transAngles :
    from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

# variables in PDF
time     = RealVar(  'time',     Title = 'Decay time',      Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. ) )
trueTime = RealVar(  'truetime', Title = 'True decay time', Unit = 'ns', Observable = True, Value = 0.,  MinMax = ( 0.,  20. ) )
iTag     = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'Untagged' : 0 } )
angles   = [ angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] ]

# ntuple variables
BMass    = RealVar( 'mass',   Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = ( 5200., 5550. ), Value = 5368. )
mumuMass = RealVar( 'mdau1',  Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. )    )
KKMass   = RealVar( 'mdau2',  Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = ( 1020. - 30., 1020. + 30. )    )
timeRes  = RealVar( 'sigmat', Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = ( 0.0, 0.12 )                   )

tagDecision = Category( 'tagdecision', Title = 'Tag decision', Observable = True, States = { 'Untagged' : 0 } )#{ 'Untagged' : 0, 'B' : +1, 'Bbar' : -1 } )

sel    = Category( 'sel',                       Title = 'Selection',            Observable = True, States = { 'selected' : +1 } )
trigUB = Category( 'triggerDecisionUnbiased',   Title = 'Trigger Unbiased',     Observable = True, States = { 'selected' : +1 } )
trigB  = Category( 'triggerDecisionBiasedExcl', Title = 'Trigger Excl. Biased', Observable = True, States = { 'selected' : +1 } )
bkgcat = Category( 'bkgcat',                    Title = 'Background category',  Observable = True, States = { 'cat0': 0, 'cat10': 10 } )

muPlusTrackChi2 = RealVar( 'muplus_track_chi2ndof',  Title = 'mu+ track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
muMinTrackChi2  = RealVar( 'muminus_track_chi2ndof', Title = 'mu- track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
KPlusTrackChi2  = RealVar( 'Kplus_track_chi2ndof',   Title = 'K+ track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )
KMinTrackChi2   = RealVar( 'Kminus_track_chi2ndof',  Title = 'K- track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )

obsSetUB = [ time, trueTime ] + angles +  [ BMass, mumuMass, KKMass, timeRes ] + [tagDecision] \
           + [ sel, trigUB, bkgcat, muPlusTrackChi2, muMinTrackChi2, KPlusTrackChi2, KMinTrackChi2 ]
obsSetB  = [ time, trueTime ] + angles +  [ BMass, mumuMass, KKMass, timeRes ] + [tagDecision] \
           + [ sel, trigB, bkgcat, muPlusTrackChi2, muMinTrackChi2, KPlusTrackChi2, KMinTrackChi2 ]

# read ntuple
from P2VVGeneralUtils import readData
dataUB = readData( dataSetFile, dataSetName = dataSetName, NTuple = True, observables = obsSetUB, Rename = 'DecayTreeUB' )
dataB  = readData( dataSetFile, dataSetName = dataSetName, NTuple = True, observables = obsSetB,  Rename = 'DecayTreeB'  )

if trigger == 'ExclBiased' :
    data = dataB
elif trigger == 'Unbiased' :
    data = dataUB
else :
    dataSetVars = dataUB.get()
    dataSetVars.remove(dataSetVars.find('triggerDecisionUnbiased'))

    from ROOT import RooDataSet, RooFit
    data = RooDataSet( 'DecayTreeUBB', 'Decay Tree UB+B', dataSetVars, RooFit.Import(dataUB) )
    data.append(dataB)

obsSetP2VV = [ time if tResModel in [ 'Gauss', '3Gauss' ] else trueTime ] + angles


###########################################################################################################################################
## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
#####################################################################

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

pdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )


###########################################################################################################################################
## compute angular efficiency moments and multiply PDF with efficiency function ##
##################################################################################

# print PDF, data, variables and parameters
print '\nData set:'
data.Print()
print '\nPDF:'
pdf.Print()
print '\nLifetime resolution model:'
timeResModel['model'].Print()
print '\nVariables in PDF:'
for var in pdf.getObservables(data) : var.Print()
print '\nParameters in PDF:'
for par in pdf.getParameters(data) : par.Print()
print

normSet = angleFuncs.angles.values()

# moments builder with angular functions from physics PDF
from P2VVGeneralUtils import RealMomentsBuilder
physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( func, 1, pdf, normSet )\
                                              for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func
                                            )
                                )

# moments builder with angular basis functions
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#indices = [ ( PIndex, 2, YIndex1 ) for PIndex in range(40) for YIndex1 in [ +1, -1 ] ]
#indices = [ ( PIndex, 2, YIndex1 ) for PIndex in range(40) for YIndex1 in [ -2, 1 ] ]

basisMoments = RealMomentsBuilder()
basisMoments.appendPYList( angleFuncs.angles, indices, PDF = pdf, NormSet = normSet )

if readMoments :
    # read moments from file
    physMoments.read(  momentsFile + 'Phys'  )
    basisMoments.read( momentsFile + 'Basis' )
else :
    # compute moments from data set
    physMoments.compute(data)
    basisMoments.compute(data)

    physMoments.write(  momentsFile + 'Phys'  )
    basisMoments.write( momentsFile + 'Basis' )

# print moments to screen
physMoments.Print(  Scale = 1. / 16. / sqrt(pi)                       )
basisMoments.Print( Scale = 1. /  2. / sqrt(pi), MinSignificance = 5. )


###########################################################################################################################################
## multiply PDF with angular efficiency ##
##########################################

if multPdfEff :
  effPdf = basisMoments * pdf

  basisMomentsSignif = RealMomentsBuilder()
  basisMomentsSignif.appendPYList( angleFuncs.angles, [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if not transAngles \
                                                 else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
                                 )
  basisMomentsSignif.read(momentsFile + 'Basis')
  basisMomentsSignif.Print( Scale = 1. / 2. / sqrt(pi) )

  effSignifPdf = basisMomentsSignif.multiplyPDFWithEff( pdf, Name = 'sig_t_angles_tagCat_iTag_x_EffSignif', EffName = 'effSignif' )


###########################################################################################################################################
## make some plots ##
#####################

if multPdfEff and makePlots :
    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen

    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, logY )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSetP2VV[ : 5 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSetP2VV[ : 5 ] ]
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

    # print canvas to file
    timeAnglesCanv.Print(plotsFile)
