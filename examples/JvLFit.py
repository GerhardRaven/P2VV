###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
generateData   = False
addTaggingVars = True
fitData        = True
makePlots      = True

plotsFile = 'JvLFitTagCats.ps'

#dataSetName = 'JpsiphiData'
#dataSetFile = 'JvLFitTagCats.root'
#realData = False

dataSetName = 'DecayTree'
dataSetFile = '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
realData = True

nEvents    = 200000
sigFrac    = 0.8

# PDF options
nominalFit = False
components = '' # 'signal' # 'background'
bkgAngles  = ''#'histPdf'

# transversity amplitudes
amplitudeParam = 'phasesSWaveFrac'

A0Mag2Val    = 0.522
AperpMag2Val = 0.249
Apar2MagVal  = 1. - A0Mag2Val - AperpMag2Val

A0PhVal      = 0.
AperpPhVal   = 2.77
AparPhVal    = 3.32

fSVal        = 0.016
ASMag2Val    = fSVal / ( 1. - fSVal )
ASPhVal      = 2.74

# CP violation parameters
carthLambdaCP = False
phiCPVal      = 0.17
lambdaCPSqVal = 1.

# B lifetime parameters
GammaVal        = 0.667
dGammaVal       = 0.12

# asymmetries
AProdVal = 0.

# fit options
fitOpts = dict(  NumCPU              = 10
               , Timer               = 1
               , Minos               = False
               , Hesse               = False
               , Save                = True
              )

# plot options
if nominalFit : angleNames = ( 'cos(#psi_{tr})',  'cos(#theta_{tr})', '#phi_{tr}' )
else          : angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{l})',  '#phi'      )
numBins         = ( 60, 30, 30, 30 )
numTimeBins     = ( 60, 60 )
numAngleBins    = ( 20, 40, 20 )
numBkgAngleBins = ( 5, 40, 5 )
lineWidth       = 2
markStyle       = 8
markSize        = 0.4


###########################################################################################################################################
## create variables and read real data ##
#########################################

# import RooFit wrappers
from RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# angular functions
if nominalFit :
    from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

# variables in PDF (except for tagging category)
time   = RealVar(  'time',         Title = 'Decay time', Unit = 'ps',   Observable = True, Value = 0.5,   MinMax = ( 0.3, 14. ) )
iTag   = Category( 'iTag',         Title = 'Initial state flavour tag', Observable = True, States = { 'B' : +1, 'Bbar' : -1 } )
angles = [ angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] ]

BMass = RealVar( 'mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, Value = 5368., MinMax = ( 5200., 5550. ), nBins = 48
                       ,  Ranges =  {  'LeftSideBand'  : ( None,  5330. )
                                     , 'Signal'        : ( 5330., 5410. )
                                     , 'RightSideBand' : ( 5410., None  )
                                    }
               )

obsSetP2VV = [ time ] + angles + [ iTag, BMass ]

# tagging categories
if not generateData and realData :
    # ntuple variables
    mpsi = RealVar( 'mdau1', Title = 'M(#mu#mu)', Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. ), nBins =  32 )
    mphi = RealVar( 'mdau2', Title = 'M(KK)',     Unit = 'MeV', Observable = True, MinMax = ( 1020. - 12., 1020. + 12. ), nBins =  16 )

    timeRes = RealVar( 'sigmat', Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = ( 0.0, 0.12 ), nBins =  50 )

    tagDecision = Category( 'tagdecision_os', Title = 'Tag decision', Observable = True
                           , States = { 'B' : +1, 'Bbar' : -1 , 'Untagged' : 0 }
                          )
    tagOmega = RealVar( 'tagomega_os', Title = 'Estimated wrong tag', Observable = True, Value = 0.25, MinMax = ( 0., 0.50001 ) )
    tagCat = Category( 'tagcat_os',   Title = 'Tagging Category', Observable = True
                      , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, 6 ) ]
                     )

    sel   = Category( 'sel',             Title = 'Selection',        Observable = True, States = { 'selected' : +1 } )
    trig  = Category( 'triggerDecision', Title = 'Trigger Decision', Observable = True, States = { 'selected' : +1 } )

    obsSetNTuple = [ time ] + angles +  [ BMass, mpsi, mphi ] + [ tagDecision, tagOmega, tagCat ] + [ sel, trig ]

    # read ntuple and add P2VV tagging variables
    from P2VVGeneralUtils import readData
    data = readData( dataSetFile, dataSetName = dataSetName, NTuple = True, observables = obsSetNTuple )

else :
    data = None

from P2VVParameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
tagCats = TaggingCategories( tagCat = 'tagCatP2VV', DataSet = data, wTagP0Constraint = True, wTagP1Constraint = True )
tagCatP2VV = tagCats['tagCat']
obsSetP2VV.append( tagCatP2VV )

# tagging parameters
numTagCats    = tagCats['numTagCats']
cat5Min       = 1 + 4 * int( float( numTagCats - 1 ) / 5. )
taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,       numTagCats ) ] )
tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( cat5Min, numTagCats ) ] )

# tagging category ranges
tagCatP2VV.setRange( 'UntaggedRange', 'Untagged'    )
tagCatP2VV.setRange( 'TaggedRange',   taggedCatsStr )
tagCatP2VV.setRange( 'TagCat5Range',  tagCat5Str    )


###########################################################################################################################################
## initialize PDF component objects ##
######################################

nSignal     = data.numEntries() * sigFrac          if data else nEvents * sigFrac
nBackground = data.numEntries() * ( 1. - sigFrac ) if data else nEvents * ( 1. - sigFrac )

signalComps     = Component('signal', [ ], Yield = ( nSignal,     0., 1.1 * nSignal     ) )
backgroundComps = Component('bkg'   , [ ], Yield = ( nBackground, 0., 2.0 * nBackground ) )


###########################################################################################################################################
## build mass PDFs ##
#####################

# build the signal and background mass PDFs
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass, LP2011_Background_Mass as BackgroundBMass
signalBMass     = SignalBMass(     Name = 'sig_m', mass = BMass )
backgroundBMass = BackgroundBMass( Name = 'bkg_m', mass = BMass )

signalComps     += signalBMass.pdf()
backgroundComps += backgroundBMass.pdf()


###########################################################################################################################################
## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
#####################################################################

# transversity amplitudes
if amplitudeParam == 'phasesSWaveFrac' :
    from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
    amplitudes = Amplitudes(  A0Mag2    = A0Mag2Val
                            , A0Phase   = A0PhVal
                            , AperpMag2 = AperpMag2Val
                            , AparPhase = AparPhVal
                            , f_S_Re    = fSVal * cos(ASPhVal)
                            , f_S_Im    = fSVal * sin(ASPhVal)
                           )

else :
    from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
    amplitudes = Amplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                            , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                            , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                            , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                            , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                            , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                           )

# B lifetime
from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
lifetimeParams = LifetimeParams(  Gamma = dict(Value = GammaVal)
                                , deltaGamma = dict(  Name = 'dGamma', Value = dGammaVal
#                                                    , Blind = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
                                                   )
                                , deltaMConstraint = True
                               )

from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
timeResModel = TimeResolution( time = time, timeResSFConstraint = True )

# CP violation parameters
from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
lambdaCP = CPParam(  lambdaCPSq = lambdaCPSqVal
                   , phiCP = dict(  Value = phiCPVal
#                                  , Blind = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
                                 )
                  )

# tagging parameters
from P2VVParameterizations.FlavourTagging import WTagCatsCoefAsyms_TaggingParams as TaggingParams
taggingParams = TaggingParams( AProd = AProdVal, ANorm = -lambdaCP['C'].getVal(), **tagCats.tagCatsDict() )

# coefficients for time functions
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] ) 

# build signal PDF
args = dict(  time                   = time
            , iTag                   = iTag
            , tagCat                 = tagCatP2VV
            , tau                    = lifetimeParams['MeanLifetime']
            , dGamma                 = lifetimeParams['deltaGamma']
            , dm                     = lifetimeParams['deltaM']
            , dilutions              = taggingParams['dilutions']
            , ADilWTags              = taggingParams['ADilWTags']
            , avgCEvens              = taggingParams['avgCEvens']
            , avgCOdds               = taggingParams['avgCOdds']
            , tagCatCoefs            = taggingParams['tagCatCoefs']
            , coshCoef               = timeBasisCoefs['cosh']
            , sinhCoef               = timeBasisCoefs['sinh']
            , cosCoef                = timeBasisCoefs['cos']
            , sinCoef                = timeBasisCoefs['sin']
            , resolutionModel        = timeResModel['model']
            , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                       + timeResModel.externalConstraints()\
                                       + tagCats.externalConstraints()
#            , ConditionalObservables = timeResModel.ConditionalObservables()
           )

sig_t_angles_tagCat_iTag = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )
signalComps += sig_t_angles_tagCat_iTag


###########################################################################################################################################
## build background time PDF ##
###############################

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
backgroundTime = BackgroundTime( Name = 'bkg_t', time = time, resolutionModel = timeResModel['model'] )
backgroundComps += backgroundTime.pdf()


###########################################################################################################################################
## build background angular PDF ##
##################################

sideBandData = data.reduce( CutRange = 'LeftSideBand' )
sideBandData.append( data.reduce( CutRange = 'RightSideBand' ) )

if nominalFit :
    bkg_angles = HistPdf(  Name = 'bkg_angles'
                         , Observables = angles
                         , Binning =  {  angleFuncs.angles['cpsi']   : 5
                                       , angleFuncs.angles['ctheta'] : 7
                                       , angleFuncs.angles['phi' ]   : 9
                                      }
                         , Data = sideBandData
                        )
elif bkgAngles == 'histPdf' :
    bkg_angles = HistPdf(  Name = 'bkg_angles'
                         , Observables = angles
                         , Binning =  {  angleFuncs.angles['cpsi']   : 5
                                       , angleFuncs.angles['ctheta'] : 40
                                       , angleFuncs.angles['phi' ]   : 5
                                      }
                         , Data = sideBandData
                        )

else :
    from array import array
    from ROOT import RooBinning
    ctKBinBounds = array( 'd', [ -1.,                      -0.6,      -0.2,      0.2,      0.6,                   1. ] )
    ctlBinBounds = array( 'd', [ -1., -0.95, -0.90, -0.85, -0.6,      -0.2,      0.2,      0.6, 0.85, 0.90, 0.95, 1. ] )
    phiBinBounds = array( 'd', [ -pi,                      -0.6 * pi, -0.2 * pi, 0.2 * pi, 0.6 * pi,              pi ] )
    ctKBins = RooBinning( len(ctKBinBounds) - 1, ctKBinBounds, 'bkg_ctKBins' )
    ctlBins = RooBinning( len(ctlBinBounds) - 1, ctlBinBounds, 'bkg_ctlBins' )
    phiBins = RooBinning( len(phiBinBounds) - 1, phiBinBounds, 'bkg_phiBins' )
    angles[0].setBinning( ctKBins, 'bkg_ctKBins' )
    angles[1].setBinning( ctlBins, 'bkg_ctlBins' )
    angles[2].setBinning( phiBins, 'bkg_phiBins' )
    bkg_angCoefs = [ RealVar(  'bkg_angBin_%d_%d_%d' % ( bin0, bin1, bin2 )
                             , Title = 'Background angles bin %d-%d-%d' % ( bin0, bin1, bin2 )
                             , Value = 1. / ( len(ctKBinBounds) - 1 ) / ( len(ctlBinBounds) - 1 ) / ( len(phiBinBounds) - 1 )
                             , MinMax = ( 0., 1. )
                            )\
                     if bin0 != 0 or bin1 != 0 or bin2 != 0 else None\
                     for bin2 in range( len(phiBinBounds) - 1 )\
                     for bin1 in range( len(ctlBinBounds) - 1 )\
                     for bin0 in range( len(ctKBinBounds) - 1 )
                   ]
    del bkg_angCoefs[0]

    bkg_angles = BinnedPdf(  Name = 'bkg_angles'
                           , Observables = angles
                           , Binnings = [ ctKBins, ctlBins, phiBins ]
                           , Coefficients = bkg_angCoefs
                           , BinIntegralCoefs = True
                          )

    sideBandDataTree = sideBandData.tree()
    for bin, coef in enumerate(bkg_angCoefs) :
        ctKBin, ctlBin, phiBin = coef.GetName()[ 11 : ].split('_')
        selection = dict(  ctK = angles[0].GetName(), ctl = angles[1].GetName(), phi = angles[2].GetName()
                         , ctKMin = ctKBinBounds[ int(ctKBin) ], ctKMax = ctKBinBounds[ int(ctKBin) + 1 ]
                         , ctlMin = ctlBinBounds[ int(ctlBin) ], ctlMax = ctlBinBounds[ int(ctlBin) + 1 ]
                         , phiMin = phiBinBounds[ int(phiBin) ], phiMax = phiBinBounds[ int(phiBin) + 1 ]
                        )

        value = float( sideBandDataTree.GetEntries( (     '%(ctK)s >= %(ctKMin)f && %(ctK)s < %(ctKMax)f'\
                                                     + '&& %(ctl)s >= %(ctlMin)f && %(ctl)s < %(ctlMax)f'\
                                                     + '&& %(phi)s >= %(phiMin)f && %(phi)s < %(phiMax)f'
                                                    ) % selection
                                                  )
                     ) / sideBandDataTree.GetEntries()
        coef.setVal(value)
        coef.setError( 0.1 * value )
        coef.setConstant()

backgroundComps += bkg_angles

if makePlots :
    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed

    # plot background angles
    bkgAnglesCanv = TCanvas( 'bkgAnglesCanv', 'Background Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle )\
            in zip( bkgAnglesCanv.pads( 3, 2 ), angles, numBkgAngleBins, [ angle.GetTitle() for angle in angles ], angleNames ) :
        plot(  pad, obs, sideBandData, bkg_angles, xTitle = xTitle
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth       )
            )

###########################################################################################################################################
## build background tagging PDF ##
##################################

backgroundComps += BinnedPdf(  'bkg_tagCat_iTag'
                             , Categories   = ( tagCatP2VV, iTag )
                             , Coefficients = [  taggingParams['tagCatCoefs']
                                               , [  RealVar( 'bkg_BbarFrac'
                                                            , Title    = 'Anti-B fraction in background'
                                                            , Value    = 0.5
                                                            , MinMax = ( 0., 1. )
                                                            , Constant = True
                                                           )
                                                 ]
                                              ]
                            )


###########################################################################################################################################
## build full PDF ##
####################

if components == 'signal' :
    pdf = buildPdf( [ signalComps                  ], Observables = obsSetP2VV, Name = 'JpsiphiPDF' )
elif components == 'background' :
    pdf = buildPdf( [ backgroundComps              ], Observables = obsSetP2VV, Name = 'JpsiphiPDF' )
else :
    pdf = buildPdf( [ signalComps, backgroundComps ], Observables = obsSetP2VV, Name = 'JpsiphiPDF' )


###########################################################################################################################################
## generate/read data and fit ##
################################

# generate data
from P2VVLoad import RooFitOutput
if generateData :
    print 'JvLFit: generating %d events' % nEvents
    data = pdf.generate(obsSetP2VV)

    from P2VVGeneralUtils import writeData
    writeData( dataSetFile, dataSetName, data )

else :
    from P2VVGeneralUtils import readData
    if realData :
        from P2VVGeneralUtils import addTaggingObservables
        addTaggingObservables(data, iTag.GetName(), tagCatP2VV.GetName(), tagDecision.GetName(), tagOmega.GetName(), tagCats['tagCats'])

    else :
        data = readData( dataSetFile, dataSetName = dataSetName, observables = obsSetP2VV )

if fitData :
    # fix values of some parameters
    lifetimeParams.setConstant('deltaM')
    lifetimeParams.setConstant('dGamma')
    lambdaCP.setConstant('lambdaCPSq')
    lambdaCP.setConstant('phiCP')

    for CEvenOdd in taggingParams['CEvenOdds'] :
        CEvenOdd.setConstant('avgCEven.*')
        CEvenOdd.setConstant('avgCOdd.*')
    taggingParams.setConstant('tagCatCoef.*')
    tagCats.setConstant('wTag.*')

    #amplitudes.setConstant('f_S|ASPhase')
    #amplitudes.setConstant('.*')

    #timeResModel.setConstant('timeResSF')
    #backgroundTime.setConstant('.*')
    #signalBMass.setConstant('.*')
    #backgroundBMass.setConstant('.*')

    # fit data
    print 'JvLFit: fitting %d events' % data.numEntries()
    fitResult = pdf.fitTo( data, **fitOpts )


###########################################################################################################################################
## make some plots ##
#####################

if makePlots :
    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed

    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSetP2VV[ : 5 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSetP2VV[ : 5 ] ]
                   , ( '', ) + angleNames
                  ) :
        plot(  pad, obs, data, pdf, xTitle = xTitle
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth       )
             , components = { 'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                            , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                            }
            )

    # plot lifetime
    timePlotTitles = tuple( [ time.GetTitle() + title for title in (  ' - B (linear)'
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
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                            )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts        )
             , components = { 'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                            , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                            }
            )

    # set Y-axis maximum for lifetime plots
    timeYMax = max( frame.GetMaximum() for frame in timeCanv.frameHists() )
    map( lambda obj : obj.SetMaximum(timeYMax), ( frame for frame in timeCanv.frameHists() ) )
    for pad in timeCanv.pads() : pad.Draw()

    # plot lifetime (tagged/untagged)
    timePlotTitles1 = tuple( [ time.GetTitle() + title for title in (  ' - B (untagged)'
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
             ,   ( { 'Cut' : '{tag} == +1 && {cat} == 0'.format(    tag = iTag.GetName(), cat = tagCatP2VV.GetName()                 ) }, )
               + ( { 'Cut' : '{tag} == +1 && {cat} >  0'.format(    tag = iTag.GetName(), cat = tagCatP2VV.GetName()                 ) }, )
               + ( { 'Cut' : '{tag} == +1 && {cat} >= {c5:d}'.format( tag = iTag.GetName(), cat = tagCatP2VV.GetName(), c5 = cat5Min ) }, )
               + ( { 'Cut' : '{tag} == -1 && {cat} == 0'.format(    tag = iTag.GetName(), cat = tagCatP2VV.GetName()                 ) }, )
               + ( { 'Cut' : '{tag} == -1 && {cat} >  0'.format(    tag = iTag.GetName(), cat = tagCatP2VV.GetName()                 ) }, )
               + ( { 'Cut' : '{tag} == -1 && {cat} >= {c5:d}'.format( tag = iTag.GetName(), cat = tagCatP2VV.GetName(), c5 = cat5Min ) }, )
             ,   ( { 'Slice' : ( iTag, 'B'    ), 'ProjectionRange' : 'UntaggedRange' }, )
               + ( { 'Slice' : ( iTag, 'B'    ), 'ProjectionRange' : 'TaggedRange'   }, )
               + ( { 'Slice' : ( iTag, 'B'    ), 'ProjectionRange' : 'TagCat5Range'  }, )
               + ( { 'Slice' : ( iTag, 'Bbar' ), 'ProjectionRange' : 'UntaggedRange' }, )
               + ( { 'Slice' : ( iTag, 'Bbar' ), 'ProjectionRange' : 'TaggedRange'   }, )
               + ( { 'Slice' : ( iTag, 'Bbar' ), 'ProjectionRange' : 'TagCat5Range'  }, )
            ) :
        plot(  pad, time, data, pdf
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                            )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts      )
             , components = { 'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                            , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                            }
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
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                             )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts       )
             , components = { 'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                            , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                            }
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
    anglesCanv.Print(plotsFile)
    bkgAnglesCanv.Print(plotsFile + ')')

