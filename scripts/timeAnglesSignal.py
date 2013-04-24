###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
generateData = True
fitData      = True
makePlots    = True

nEvents = 10000
dataSetName = 'JpsiphiData'
dataSetFile = 'timeAnglesSignal.root'
NTuple = False

# transversity amplitudes
A0Mag2Val    =  0.4
A0PhVal      =  0.
AparMag2Val  =  0.3
AparPhVal    = -2.4
AperpMag2Val =  0.3
AperpPhVal   = -0.79
ASMag2Val    =  0.3
ASPhVal      =  2.4

# CP violation parameters
carthLambdaCP = False
phiCPVal      = -pi / 4.
lambdaCPSqVal = 0.6

# B lifetime parameters
GammaVal        = 0.68
dGammaVal       = 0.05
dmVal           = 17.8
timeResSigmaVal = 0.05

# asymmetries
AProdVal =  0.4
ANormVal = -( 1. - lambdaCPSqVal ) / ( 1. + lambdaCPSqVal )

# tagging parameters
WTagVal  = ( 0.1 + 0.2 ) / 2.
AWTagVal = ( 0.1 - 0.2 ) / ( 0.1 + 0.2 )

# plot options
angleNames   = ( 'cos(#theta_{K})', 'cos(#theta_{l})', '#phi' )
numBins      = ( 50, 50, 50, 50 )
numTimeBins  = ( 30, 30, 30 )
numAngleBins = ( 30, 30, 30 )
lineWidth    = 2
markStyle    = 8
markSize     = 0.4


###########################################################################################################################################
## build the PDF ##
###################

# import RooFit wrappers
from P2VV.RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# angular functions
from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles
angleFuncs = JpsiphiHelicityAngles( cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi' )

# variables in PDF
time = RealVar(  't',          Title = 'Decay time', Unit = 'ps',   Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
iTag = Category( 'tagInitial', Title = 'Initial state flavour tag', Observable = True, States = {'B': +1, 'Bbar': -1} )#, 'Untagged' : 0} )

angles      = ( angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] )
observables = [ time ] + list(angles) + [ iTag ]

# transversity amplitudes
from P2VV.Parameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet
transAmps = JpsiVCarthesian_AmplitudeSet(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                         , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                         , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                         , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                         , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                                         , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                                        )

# B lifetime
from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = GammaVal, dGamma = dGammaVal, dM = dmVal )

from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution
timeResModel = Gaussian_TimeResolution( time = time, timeResSigma = timeResSigmaVal )

# CP violation parameters
if carthLambdaCP :
  # carthesian lambda
  from P2VV.Parameterizations.CPVParams import LambdaCarth_CPParam
  lambdaCP = LambdaCarth_CPParam( ReLambdaCP = sqrt(lambdaCPSqVal) * cos(-phiCPVal), ImLambdaCP = sqrt(lambdaCPSqVal) * sin(-phiCPVal) )

else :
  # polar lambda
  from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam
  lambdaCP = LambdaSqArg_CPParam( lambdaCPSq = lambdaCPSqVal, phiCP = phiCPVal )

# tagging parameters
from P2VV.Parameterizations.FlavourTagging import WTagsCoefAsyms_TaggingParams
taggingParams = WTagsCoefAsyms_TaggingParams( WTag = WTagVal, AWTag = AWTagVal, AProd = AProdVal, ANorm = ANormVal )

# coefficients for time functions
from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
timeBasisCoefs = JpsiphiBTagDecayBasisCoefficients( angleFuncs.functions, transAmps, lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] ) 

# build the B_s -> J/psi phi signal PDF
args = {
    'time'            : time
  , 'iTag'            : iTag
  , 'tau'             : lifetimeParams['MeanLifetime']
  , 'dGamma'          : lifetimeParams['dGamma']
  , 'dm'              : lifetimeParams['dM']
  , 'dilution'        : taggingParams['dilution']
  , 'ADilWTag'        : taggingParams['ADilWTag']
  , 'avgCEven'        : taggingParams['avgCEven']
  , 'avgCOdd'         : taggingParams['avgCOdd']
  , 'coshCoef'        : timeBasisCoefs['cosh']
  , 'sinhCoef'        : timeBasisCoefs['sinh']
  , 'cosCoef'         : timeBasisCoefs['cos']
  , 'sinCoef'         : timeBasisCoefs['sin']
  , 'resolutionModel' : timeResModel['model']
}

pdf = BTagDecay('JpsiphiPDF', **args)


###########################################################################################################################################
## generate/read data and fit ##
################################

# generate data
from P2VV.Load import RooFitOutput
if generateData :
  print 'timeAnglesSignal: generating %d events' % nEvents
  data = pdf.generate( observables, nEvents )

  from P2VV.GeneralUtils import writeData
  writeData( dataSetFile, dataSetName, data, NTuple )

else :
  from P2VV.GeneralUtils import readData
  data = readData( dataSetFile, dataSetName = dataSetName, NTuple = NTuple, observables = observables )

  # TODO: a trick to change the observables in a data set ( iTag( +1, -1 ) <--> iTag( +1, -1, 0 ) )
  data = RooDataSet( dataSetName + '1', '', data, ( obs._var for obs in observables ) )

if fitData :
  # fix values of some parameters
  #ANuissance.setConstant('avgCOdd')
  #lambdaCP.setConstant('phiCP')
  #taggingParams.setConstant('AWTag')

  # fit data
  print 'timeAnglesSignal: fitting %d events' % data.numEntries()
  pdf.fitTo( data, NumCPU = 2, Timer = 1 )#, ConditionalObservables = [iTag] )


###########################################################################################################################################
## make some plots ##
#####################

if makePlots :
    # import plotting tools
    from P2VV.Load import LHCbStyle
    from P2VV.GeneralUtils import plot
    from ROOT import TCanvas

    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , observables[ : -1 ]
                   , numBins
                   , [ var.GetTitle() for var in observables[ : -1 ] ]
                   , ( '', ) + angleNames
                  ) :
        plot(  pad, obs, data, pdf, xTitle = xTitle
             , frameOpts = { 'Bins' : nBins, 'Title' : plotTitle }
             , dataOpts  = { 'MarkerStyle' : markStyle, 'MarkerSize' : markSize, 'XErrorSize' : 0 }
             , pdfOpts   = { 'LineWidth' : lineWidth }
            )

    # plot lifetime
    timePlotTitles = tuple( [ time.GetTitle() + str for str in (  ' - B (linear)'
                                                                , ' - B (logarithmic)'
                                                                , ' - #bar{B} (linear)'
                                                                , ' - #bar{B} (logarithmic)'
                                                               )
                            ] )
    timeCanv = TCanvas( 'timeCanv', 'Lifetime' )
    for ( pad, nBins, plotTitle, dataCuts, pdfCuts, logY )\
            in zip(  timeCanv.pads( 2, 2 )
                   , 2 * numTimeBins
                   , timePlotTitles
                   , 2 * ( { 'Cut' : iTag.GetName() + ' == +1' }, ) + 2 * ( { 'Cut' : iTag.GetName() + ' == -1' }, )
                   , 2 * ( { 'Slice' : ( iTag, 'B' ) }, ) + 2 * ( { 'Slice' : ( iTag, 'Bbar' ) }, )
                   , 2 * ( False, True )
                  ) :
        plot(  pad, time, data, pdf, logy = logY
             , frameOpts = dict( Bins = nBins, Title = plotTitle )
             , dataOpts  = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts   = dict( LineWidth = lineWidth, **pdfCuts )
            )

    # set Y-axis maximum for lifetime plots
    timeYMax = max( frame.GetMaximum() for frame in timeCanv.frameHists() )
    map( lambda obj : obj.SetMaximum(timeYMax), ( frame for frame in timeCanv.frameHists() ) )
    for pad in timeCanv.pads() : pad.Draw()

    # plot angles
    anglePlotTitles =   tuple(  [ angle.GetTitle() + ' - B'       for angle in angles ]\
                              + [ angle.GetTitle() + ' - #bar{B}' for angle in angles ] )
    anglesCanv = TCanvas( 'anglesCanv', 'Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, dataCuts, pdfCuts )\
            in zip(  anglesCanv.pads( 3, 2 )
                   , 2 * angles
                   , 2 * numAngleBins
                   , anglePlotTitles
                   , 2 * angleNames
                   , 3 * ( { 'Cut' : iTag.GetName() + ' == +1' }, ) + 3 * ( { 'Cut' : iTag.GetName() + ' == -1' }, )
                   , 3 * ( { 'Slice' : ( iTag, 'B' ) }, ) + 3 * ( { 'Slice' : ( iTag, 'Bbar' ) }, )
                  ) :
        plot(  pad, obs, data, pdf, xTitle = xTitle
             , frameOpts = dict( Bins = nBins, Title = plotTitle )
             , dataOpts  = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts )
             , pdfOpts   = dict( LineWidth = lineWidth, **pdfCuts )
            )

    # set Y-axis maximum for angles plots
    from collections import defaultdict
    maxVal = defaultdict(float)
    for frame in anglesCanv.frameHists() :
        maxVal[ frame.GetXaxis().GetTitle() ] = max( frame.GetMaximum(), maxVal[ frame.GetXaxis().GetTitle() ] )
    for frame in anglesCanv.frameHists() :
        frame.SetMaximum( maxVal[ frame.GetXaxis().GetTitle() ] )
    for pad  in anglesCanv.pads() : pad.Draw()
