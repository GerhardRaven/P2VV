###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
readMoments = False
makePlots   = True
transAngles = True

momentsFile = 'effMoments' + ( 'Trans' if transAngles else 'Hel' )
plotsFile   = 'effMoments' + ( 'Trans' if transAngles else 'Hel' ) + '.ps'

dataSetName = 'DecayTree'
#dataSetFile = '/data/bfys/dveijk/MC/2012/Bs2JpsiPhi_MC11a_ntupleB_for_fitting_20120109.root'
dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_MC11a_ntupleB_for_fitting_20120209.root'

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
ws = RooObject(workspace = 'ws')

# angular functions
if transAngles :
    from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

# variables in PDF
time     = RealVar(  'time',     Title = 'Decay time',      Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14.    ) )
trueTime = RealVar(  'truetime', Title = 'True decay time', Unit = 'ns', Observable = True, Value = 0.,  MinMax = ( 0.,   0.020 ) )
iTag     = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'Untagged' : 0 } )
angles   = [ angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] ]

# ntuple variables
BMass   = RealVar( 'mass',   Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = ( 5200., 5550. ), Value = 5368. )
mpsi    = RealVar( 'mdau1',  Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. )    )
mphi    = RealVar( 'mdau2',  Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = ( 1020. - 12., 1020. + 12. )    )
timeRes = RealVar( 'sigmat', Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = ( 0.0, 0.12 )                   )

sel    = Category( 'sel',             Title = 'Selection',           Observable = True, States = { 'selected' : +1 } )
trig   = Category( 'triggerDecision', Title = 'Trigger Decision',    Observable = True, States = { 'selected' : +1 } )
bkgcat = Category( 'bkgcat',          Title = 'Background category', Observable = True, States = { 'cat0': 0, 'cat10': 10 } )

obsSetNTuple = [ time, trueTime ] + angles +  [ BMass, mpsi, mphi, timeRes ] + [ sel, trig, bkgcat ]

# read ntuple
from P2VVGeneralUtils import readData
data = readData( dataSetFile, dataSetName = dataSetName, NTuple = True, observables = obsSetNTuple )

# add true lifetime to data set
trueTimeP2VV = RealVar(  'trueTimeP2VV', Title = 'True decay time', Unit = 'ps', Observable = True
                       , Formula = ( '1000. * @0', [ trueTime ], data )
                      )

obsSetP2VV = [ trueTimeP2VV ] + angles


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

from P2VVParameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
timeResModel = TimeResolution( time = trueTimeP2VV )

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
args = dict(  time                   = trueTimeP2VV
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

# print variables and parameters
print 'Variables in PDF:'
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
#indices += [ ( PIndex, 2, YIndex1 ) for PIndex in range( 3, 10 ) for YIndex1 in [ -2, 1 ] ]

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
basisMoments.Print( Scale = 1. /  2. / sqrt(pi), MinSignificance = 3. )


###########################################################################################################################################
## multiply PDF with angular efficiency ##
##########################################

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

if makePlots :
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

