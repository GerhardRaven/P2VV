###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
generateData = True
fitData      = False
makePlots    = True

nEvents = 200000
dataSetName = 'JpsiphiData'
dataSetFile = 'JvLFitTagCats.root'
#dataSetFile = '/data/bfys/jleerdam/Bs2Jpsiphi/testSample.root'
NTuple = False
plotsFile = 'JvLFitTagCats.ps'

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
AProdVal = 0.4

# tagging parameters
numTagCats    = 11
cat5Min       = 9
taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,       numTagCats ) ] )
tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( cat5Min, numTagCats ) ] )
tagCatCoefs   = [ 0.090, 0.060, 0.045, 0.025, 0.016, 0.013, 0.010, 0.008, 0.006, 0.004 ]
ATagEffs      = ( numTagCats - 1 ) * [ 0. ]
wTags         = [ 0.42, 0.38, 0.35, 0.32, 0.27, 0.25, 0.24, 0.20, 0.15, 0.10 ]
wTagBars      = wTags

# plot options
angleNames   = ( 'cos(#theta_{K})', 'cos(#theta_{l})', '#phi' )
numBins      = ( 60, 30, 30, 30 )
numTimeBins  = ( 60, 60 )
numAngleBins = ( 30, 30, 30 )
lineWidth    = 2
markStyle    = 8
markSize     = 0.4


###########################################################################################################################################
## build the PDF ##
###################

# import RooFit wrappers
from RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# constants
zero = ConstVar( Name = 'zero', Value = 0 )

# angular functions
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles
angleFuncs = JpsiphiHelicityAngles( cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi' )

# variables in PDF
time   = RealVar(  't',          Title = 'Decay time', Unit = 'ps',   Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
iTag   = Category( 'tagInitial', Title = 'Initial state flavour tag', Observable = True, States = {'B' : +1, 'Bbar' : -1} )
tagCat = Category( 'tagCat'    , Title = 'Tagging Category',          Observable = True
                  , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, numTagCats ) ]
                 )

angles      = ( angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] )
observables = [ time ] + list(angles) + [ iTag, tagCat ]

# variable ranges
tagCat.setRange( 'UntaggedRange', 'Untagged'    )
tagCat.setRange( 'TaggedRange',   taggedCatsStr )
tagCat.setRange( 'TagCat5Range',  tagCat5Str    )

# transversity amplitudes
from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesianAmplitudes
transAmps = JpsiVCarthesianAmplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                      , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                      , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                      , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                      , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                                      , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                                     )

# B lifetime
from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = GammaVal, deltaGamma = dGammaVal, deltaM = dmVal )

from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution
timeResModel = Gaussian_TimeResolution( time = time, timeResSigma = timeResSigmaVal )

# CP violation parameters
from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
lambdaCP = LambdaSqArg_CPParam( lambdaCPSq = lambdaCPSqVal, phiCP = phiCPVal )

# tagging parameters
from P2VVParameterizations.FlavourTagging import WTagsCoefAsyms_TaggingParams
tagParamVals  = dict(  [ ( 'tagCatCoef%d' % ( cat + 1 ), coef ) for cat, coef in enumerate(tagCatCoefs) ]
                     + [ ( 'ATagEff%d'    % ( cat + 1 ), asym ) for cat, asym in enumerate(ATagEffs)    ]
                     + [ ( 'wTag%d'       % ( cat + 1 ), wTag ) for cat, wTag in enumerate(wTags)       ]
                     + [ ( 'wTagBar%d'    % ( cat + 1 ), wTag ) for cat, wTag in enumerate(wTagBars)    ]
                    )
taggingParams = WTagsCoefAsyms_TaggingParams(  NumTagCats = tagCat.numTypes(), AProd = AProdVal, ANorm = -lambdaCP['C'].getVal()
                                             , **tagParamVals
                                            )

# coefficients for time functions
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
timeBasisCoefs = JpsiphiBTagDecayBasisCoefficients( angleFuncs.functions, transAmps, lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] ) 

# build the B_s -> J/psi phi signal PDF
args = {
    'time'            : time
  , 'iTag'            : iTag
  , 'tagCat'          : tagCat
  , 'tau'             : lifetimeParams['MeanLifetime']
  , 'dGamma'          : lifetimeParams['deltaGamma']
  , 'dm'              : lifetimeParams['deltaM']
  , 'dilutions'       : taggingParams['dilutions']
  , 'ADilWTags'       : taggingParams['ADilWTags']
  , 'avgCEvens'       : taggingParams['avgCEvens']
  , 'avgCOdds'        : taggingParams['avgCOdds']
  , 'tagCatCoefs'     : taggingParams['tagCatCoefs']
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
from P2VVLoad import RooFitOutput
if generateData :
  print 'JvLFit: generating %d events' % nEvents
  data = pdf.generate( observables, nEvents )

  from P2VVGeneralUtils import writeData
  writeData( dataSetFile, dataSetName, data, NTuple )

else :
  from P2VVGeneralUtils import readData
  data = readData( dataSetFile, dataSetName, NTuple, observables )

if fitData :
  # fix values of some parameters
  for CEvenOdd in taggingParams['CEvenOdds'][ 1 : ] :
      CEvenOdd.setConstant('avgCEven.*')
      CEvenOdd.setConstant('avgCOdd.*')
  taggingParams.setConstant('wTag.*')

  # fit data
  print 'JvLFit: fitting %d events' % data.numEntries()
  pdf.fitTo( data, NumCPU = 12, Timer = 1 )


###########################################################################################################################################
## make some plots ##
#####################

if makePlots :
    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
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
             , dataOpts  = { 'MarkerStyle' : markStyle, 'MarkerSize' : markSize }
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

    # plot lifetime (tagged/untagged)
    timePlotTitles1 = tuple( [ time.GetTitle() + str for str in (  ' - B (untagged)'
                                                                 , ' - B (tagged)'
                                                                 , ' - B (tagging category 5)'
                                                                 , ' - #bar{B} (untagged)'
                                                                 , ' - #bar{B} (tagged)'
                                                                 , ' - #bar{B} (tagging category 5)'
                                                                )
                            ] )
    timeCanv1 = TCanvas( 'timeCanv1', 'Lifetime' )
    for ( pad, nBins, plotTitle, dataCuts, pdfCuts )\
        in zip(  timeCanv1.pads( 3, 2 )
               , 3 * numTimeBins
               , timePlotTitles1
               ,   ( { 'Cut' : '{tag} == +1 && {cat} == 0'.format(    tag = iTag.GetName(), cat = tagCat.GetName()                 ) }, )
                 + ( { 'Cut' : '{tag} == +1 && {cat} >  0'.format(    tag = iTag.GetName(), cat = tagCat.GetName()                 ) }, )
                 + ( { 'Cut' : '{tag} == +1 && {cat} >= {c5:d}'.format( tag = iTag.GetName(), cat = tagCat.GetName(), c5 = cat5Min ) }, )
                 + ( { 'Cut' : '{tag} == -1 && {cat} == 0'.format(    tag = iTag.GetName(), cat = tagCat.GetName()                 ) }, )
                 + ( { 'Cut' : '{tag} == -1 && {cat} >  0'.format(    tag = iTag.GetName(), cat = tagCat.GetName()                 ) }, )
                 + ( { 'Cut' : '{tag} == -1 && {cat} >= {c5:d}'.format( tag = iTag.GetName(), cat = tagCat.GetName(), c5 = cat5Min ) }, )
               ,   ( { 'Slice' : ( iTag, 'B'    ), 'ProjectionRange' : 'UntaggedRange' }, )
                 + ( { 'Slice' : ( iTag, 'B'    ), 'ProjectionRange' : 'TaggedRange'   }, )
                 + ( { 'Slice' : ( iTag, 'B'    ), 'ProjectionRange' : 'TagCat5Range'  }, )
                 + ( { 'Slice' : ( iTag, 'Bbar' ), 'ProjectionRange' : 'UntaggedRange' }, )
                 + ( { 'Slice' : ( iTag, 'Bbar' ), 'ProjectionRange' : 'TaggedRange'   }, )
                 + ( { 'Slice' : ( iTag, 'Bbar' ), 'ProjectionRange' : 'TagCat5Range'  }, )
              ) :
        plot(  pad, time, data, pdf
             , frameOpts = dict( Bins = nBins, Title = plotTitle )
             , dataOpts  = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts   = dict( LineWidth = lineWidth, **pdfCuts )
            )

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

    # print canvas to file
    timeAnglesCanv.Print(plotsFile + '(')
    timeCanv.Print(plotsFile)
    timeCanv1.Print(plotsFile)
    anglesCanv.Print(plotsFile + ')')

