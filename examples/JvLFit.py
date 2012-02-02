###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
generateData   = False
doFit          = True
fitType        = 'SFit' # 'SFit' / 'CFit'
makePlots      = True
blind          = True
nominalPdf     = True
addTaggingVars = True

plotsFile = 'JvLSFit.ps' if fitType == 'SFit' else 'JvLCFit.ps'
angEffMomentsFile = 'effMomentsTransBasis' if nominalPdf else 'effMomentsHelBasis'

nTupleName = 'DecayTree'
nTupleFile = '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20120120.root'
if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLFit.root'

nEvents    = 30000
sigFrac    = 0.67

# PDF options
components      = 'all' # 'all' / 'signal' / 'background'
bkgAngles       = '' # '' / 'histPdf'
perEventTimeRes = False
multiplyTimeEff = '' # 'all' / 'signal'

# transversity amplitudes
amplitudeParam = 'phasesSWaveFrac'

A0Mag2Val    = 0.50
AperpMag2Val = 0.25
AparMag2Val  = 1. - A0Mag2Val - AperpMag2Val

A0PhVal      = 0.
AperpPhVal   = 3.
AparPhVal    = 3.

fSVal        = 0.02
ASMag2Val    = fSVal / ( 1. - fSVal )
ASPhVal      = 3.

# CP violation parameters
carthLambdaCP = False
phiCPVal      = 0.2 if blind else 0.
lambdaCPSqVal = 1.

# B lifetime parameters
GammaVal        = 0.67
dGammaVal       = 0.1

# asymmetries
AProdVal = 0.

# fit options
fitOpts = dict(  NumCPU              = 1
               , Timer               = 1
#               , Minos               = False
#               , Hesse               = False
#               , Save                = True
              )

# plot options
if nominalPdf : angleNames = ( 'cos(#psi_{tr})',  'cos(#theta_{tr})', '#phi_{tr}' )
else          : angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{l})',  '#phi'      )
numBins         = ( 60, 30, 30, 30 )
numTimeBins     = ( 30, 30, 30 )
numAngleBins    = ( 20, 20, 20 )
numBkgAngleBins = ( 5, 7, 9 )
lineWidth       = 2
markStyle       = 8
markSize        = 0.4


###########################################################################################################################################
## create variables (except for tagging category) and read real data ##
#######################################################################

# import RooFit wrappers
from RooFitWrappers import *
from P2VVLoad import RooFitOutput

# workspace
ws = RooObject(workspace = 'ws')

# angular functions
if nominalPdf :
    from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

# variables in PDF (except for tagging category)
time   = RealVar(  'time', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. )
                 , Ranges = dict( Bulk = ( None, 5. ) )
                )
iTag   = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'B' : +1, 'Bbar' : -1 } )
angles = [ angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] ]

BMass = RealVar( 'mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, Value = 5368., MinMax = ( 5200., 5550. ), nBins = 48
                ,  Ranges = dict(  LeftSideBand  = ( None,  5330. )
                                 , Signal        = ( 5330., 5410. )
                                 , RightSideBand = ( 5410., None  )
                                )
               )

obsSetP2VV = [ time ] + angles + [ iTag ]
if fitType != 'SFit' : obsSetP2VV += [ BMass ]

# ntuple variables
mpsi = RealVar( 'mdau1', Title = 'M(#mu#mu)', Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. ), nBins =  32 )
mphi = RealVar( 'mdau2', Title = 'M(KK)',     Unit = 'MeV', Observable = True, MinMax = ( 1020. - 12., 1020. + 12. ), nBins =  16 )

timeRes = RealVar( 'sigmat', Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = ( 0., 0.12 ), nBins =  50 )

tagDecision = Category( 'tagdecision_os', Title = 'Tag decision', Observable = True
                       , States = { 'B' : +1, 'Bbar' : -1 , 'Untagged' : 0 }
                      )
tagOmega = RealVar( 'tagomega_os', Title = 'Estimated wrong tag', Observable = True, Value = 0.25, MinMax = ( 0., 0.50001 ) )
tagCat = Category( 'tagcat_os',   Title = 'Tagging Category', Observable = True
                  , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, 6 ) ]
                 )

sel   = Category( 'sel',             Title = 'Selection',        Observable = True, States = { 'selected' : +1 } )
trig  = Category( 'triggerDecision', Title = 'Trigger Decision', Observable = True, States = { 'selected' : +1 } )

obsSetNTuple = [ time ] + angles +  [ BMass, mpsi, mphi, timeRes ] + [ tagDecision, tagOmega, tagCat ] + [ sel, trig ]

# read ntuple
from P2VVGeneralUtils import readData
data = readData( filePath = nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSetNTuple, Rename = 'JpsiphiData' )


###########################################################################################################################################
## build tagging categories ##
##############################

# build tagging categories
from P2VVParameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
tagCats = TaggingCategories(  tagCat = 'tagCatP2VV', DataSet = data, estWTagName = tagOmega.GetName(), TagCats = [ ], NumSigmaTagBins = 1.
                            , wTagP0Constraint = True, wTagP1Constraint = True
                           )
tagCatP2VV = tagCats['tagCat']
obsSetP2VV.append( tagCatP2VV )

# add tagging category to data set
from P2VVGeneralUtils import addTaggingObservables
addTaggingObservables(data, iTag.GetName(), tagCatP2VV.GetName(), tagDecision.GetName(), tagOmega.GetName(), tagCats['tagCats'])

# tagging parameters
numTagCats    = tagCats['numTagCats']
tagCat5Min    = tagCats.traditionalCatRange(5)[0]
taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,          numTagCats ) ] )
tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( tagCat5Min, numTagCats ) ] )

# tagging category ranges
tagCatP2VV.setRange( 'UntaggedRange', 'Untagged'    )
tagCatP2VV.setRange( 'TaggedRange',   taggedCatsStr )
tagCatP2VV.setRange( 'TagCat5Range',  tagCat5Str    )


###########################################################################################################################################
## initialize PDF component objects ##
######################################

nSignal     = data.numEntries() * sigFrac          if data else nEvents * sigFrac
nBackground = data.numEntries() * ( 1. - sigFrac ) if data else nEvents * ( 1. - sigFrac )

if fitType == 'SFit' :
    signalComps  = Component( 'signal',  [ ]                                                 )
    sigMassComps = Component( 'sigMass', [ ], Yield = ( nSignal,     0., 1.1 * nSignal     ) )
    bkgMassComps = Component( 'bkgMass', [ ], Yield = ( nBackground, 0., 2.0 * nBackground ) )
else :
    signalComps     = Component( 'signal', [ ], Yield = ( nSignal,     0., 1.1 * nSignal     ) )
    backgroundComps = Component( 'bkg'   , [ ], Yield = ( nBackground, 0., 2.0 * nBackground ) )


###########################################################################################################################################
## build mass PDFs ##
#####################

# build the signal and background mass PDFs
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass, LP2011_Background_Mass as BackgroundBMass
signalBMass     = SignalBMass(     Name = 'sig_m', mass = BMass )
backgroundBMass = BackgroundBMass( Name = 'bkg_m', mass = BMass )

if fitType == 'SFit' :
    sigMassComps += signalBMass.pdf()
    bkgMassComps += backgroundBMass.pdf()
    massPdf = buildPdf( [ sigMassComps, bkgMassComps ], Observables = [ BMass ],  Name = 'JpsiphiMass' )

else :
    signalComps     += signalBMass.pdf()
    backgroundComps += backgroundBMass.pdf()


###########################################################################################################################################
## compute S-weights ##
#######################

if fitType == 'SFit' :
    print 120 * '='
    print 'JvLFit: computing S-weights'

    from P2VVGeneralUtils import SData, splot
    massPdf.fitTo( data, **fitOpts )
    for par in massPdf.Parameters() : par.setConstant( not par.getAttribute('Yield') )
    SWeightData = SData( Pdf = massPdf, Data = data, Name = 'massSData' )
    sigData = SWeightData.data('sigMass')

    print 120 * '=' + '\n'


###########################################################################################################################################
## time acceptance function ##
##############################

if multiplyTimeEff in [ 'all', 'signal', 'background' ] :
    from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance as TimeAcceptance
    timeAcceptance = TimeAcceptance(  time = time
                                    , Input = '/data/bfys/dveijk/DataJpsiPhi/2012/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root'
                                    , Histogram = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_20bins'
                                   )


###########################################################################################################################################
## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
#####################################################################

# transversity amplitudes
if nominalPdf :
    from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
    amplitudes = Amplitudes(  A0Mag2    = A0Mag2Val
                            , A0Phase   = A0PhVal
                            , AperpMag2 = AperpMag2Val
                            , AparPhase = AparPhVal
                            , f_S       = fSVal
                            , ASPhase   = ASPhVal
                           )

elif amplitudeParam == 'phasesSWaveFrac' :
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
dGamma = dict( Name = 'dGamma', Value = dGammaVal )
if blind : dGamma['Blind'] = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
lifetimeParams = LifetimeParams( Gamma = dict(Value = GammaVal), deltaGamma = dGamma, deltaMConstraint = True )

from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
timeResModel = TimeResolution( time = time, timeResSFConstraint = True )

# CP violation parameters
from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
phiCP = dict( Value = phiCPVal )
if blind: phiCP['Blind'] = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
lambdaCP = CPParam( lambdaCPSq = lambdaCPSqVal, phiCP = phiCP )

# tagging parameters
from P2VVParameterizations.FlavourTagging import WTagCatsCoefAsyms_TaggingParams as TaggingParams
taggingParams = TaggingParams( AProd = AProdVal, ANorm = -lambdaCP['C'].getVal(), **tagCats.tagCatsDict() )

# coefficients for time functions
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] ) 

# build signal PDF
args = dict(  time                = time
            , iTag                = iTag
            , tagCat              = tagCatP2VV
            , tau                 = lifetimeParams['MeanLifetime']
            , dGamma              = lifetimeParams['deltaGamma']
            , dm                  = lifetimeParams['deltaM']
            , dilutions           = taggingParams['dilutions']
            , ADilWTags           = taggingParams['ADilWTags']
            , avgCEvens           = taggingParams['avgCEvens']
            , avgCOdds            = taggingParams['avgCOdds']
            , tagCatCoefs         = taggingParams['tagCatCoefs']
            , coshCoef            = timeBasisCoefs['cosh']
            , sinhCoef            = timeBasisCoefs['sinh']
            , cosCoef             = timeBasisCoefs['cos']
            , sinCoef             = timeBasisCoefs['sin']
            , resolutionModel     = timeResModel['model']
            , ExternalConstraints = lifetimeParams.externalConstraints()\
                                    + timeResModel.externalConstraints()\
                                    + tagCats.externalConstraints()
           )

sig_t_angles_tagCat_iTag = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

if angEffMomentsFile :
    # multiply signal PDF with angular efficiency
    from P2VVGeneralUtils import RealMomentsBuilder
    moments = RealMomentsBuilder()
    moments.appendPYList( angleFuncs.angles, [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if not nominalPdf \
                                        else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
                        )
    moments.read(angEffMomentsFile)
    moments.Print()
    sig_t_angles_tagCat_iTag = moments * sig_t_angles_tagCat_iTag

if multiplyTimeEff in [ 'all', 'signal' ] :
    # multiply signal PDF with time acceptance
    sig_t_angles_tagCat_iTag = timeAcceptance * sig_t_angles_tagCat_iTag

signalComps += sig_t_angles_tagCat_iTag


###########################################################################################################################################
## build background time PDF ##
###############################

if not fitType == 'SFit' :
    from P2VVParameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
    backgroundTime = BackgroundTime( Name = 'bkg_t', time = time, resolutionModel = timeResModel['model'] )
    bkg_t = backgroundTime.pdf()

    if multiplyTimeEff in [ 'all', 'background' ] :
        # multiply background time PDF with time acceptance
        bkg_t = timeAcceptance * bkg_t

    backgroundComps += bkg_t


###########################################################################################################################################
## build background angular PDF ##
##################################

if not fitType == 'SFit' :
    sideBandData = data.reduce( CutRange = 'LeftSideBand' )
    sideBandData.append( data.reduce( CutRange = 'RightSideBand' ) )

    if bkgAngles == 'histPdf' :
        bkg_angles = HistPdf(  Name = 'bkg_angles'
                             , Observables = angles
                             , Binning =  {  angleFuncs.angles['cpsi']   : 5
                                           , angleFuncs.angles['ctheta'] : 7 if nominalPdf else 40
                                           , angleFuncs.angles['phi' ]   : 9 if nominalPdf else 5
                                          }
                             , Data = sideBandData
                            )

    else :
        from array import array
        if nominalPdf :
            cpsBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 5. * float(i)   for i in range( 1, 5 ) ] + [ 1. ] )
            cthBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 7. * float(i)   for i in range( 1, 7 ) ] + [ 1. ] )
            phiBinBounds = array( 'd', [ -pi ] + [ pi * ( -1. + 2. / 9. * float(i) ) for i in range( 1, 9 ) ] + [ pi ] )
        else :
            cpsBinBounds = array( 'd', [ -1.,                      -0.6,      -0.2,      0.2,      0.6,                   1. ] )
            cthBinBounds = array( 'd', [ -1., -0.95, -0.90, -0.85, -0.6,      -0.2,      0.2,      0.6, 0.85, 0.90, 0.95, 1. ] )
            phiBinBounds = array( 'd', [ -pi,                      -0.6 * pi, -0.2 * pi, 0.2 * pi, 0.6 * pi,              pi ] )

        from ROOT import RooBinning
        cpsBins = RooBinning( len(cpsBinBounds) - 1, cpsBinBounds, 'bkg_cpsBins' )
        cthBins = RooBinning( len(cthBinBounds) - 1, cthBinBounds, 'bkg_cthBins' )
        phiBins = RooBinning( len(phiBinBounds) - 1, phiBinBounds, 'bkg_phiBins' )
        angles[0].setBinning( cpsBins, 'bkg_cpsBins' )
        angles[1].setBinning( cthBins, 'bkg_cthBins' )
        angles[2].setBinning( phiBins, 'bkg_phiBins' )
        bkg_angCoefs = [ RealVar(  'bkg_angBin_%d_%d_%d' % ( bin0, bin1, bin2 )
                                 , Title = 'Background angles bin %d-%d-%d' % ( bin0, bin1, bin2 )
                                 , Value = 1. / ( len(cpsBinBounds) - 1 ) / ( len(cthBinBounds) - 1 ) / ( len(phiBinBounds) - 1 )
                                 , MinMax = ( 0., 1. )
                                )\
                         if bin0 != 0 or bin1 != 0 or bin2 != 0 else None\
                         for bin2 in range( len(phiBinBounds) - 1 )\
                         for bin1 in range( len(cthBinBounds) - 1 )\
                         for bin0 in range( len(cpsBinBounds) - 1 )
                       ]
        del bkg_angCoefs[0]

        bkg_angles = BinnedPdf(  Name = 'bkg_angles'
                               , Observables = angles
                               , Binnings = [ cpsBins, cthBins, phiBins ]
                               , Coefficients = bkg_angCoefs
                               , BinIntegralCoefs = True
                              )

        sideBandDataTree = sideBandData.buildTree( ','.join( angle.GetName() for angle in angles ) )
        for bin, coef in enumerate(bkg_angCoefs) :
            cpsBin, cthBin, phiBin = coef.GetName()[ 11 : ].split('_')
            selection = dict(  cps = angles[0].GetName(), cth = angles[1].GetName(), phi = angles[2].GetName()
                             , cpsMin = cpsBinBounds[ int(cpsBin) ], cpsMax = cpsBinBounds[ int(cpsBin) + 1 ]
                             , cthMin = cthBinBounds[ int(cthBin) ], cthMax = cthBinBounds[ int(cthBin) + 1 ]
                             , phiMin = phiBinBounds[ int(phiBin) ], phiMax = phiBinBounds[ int(phiBin) + 1 ]
                            )

            value = float( sideBandDataTree.GetEntries( (     '%(cps)s >= %(cpsMin)f && %(cps)s < %(cpsMax)f'\
                                                         + '&& %(cth)s >= %(cthMin)f && %(cth)s < %(cthMax)f'\
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

if not fitType == 'SFit' :
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

if fitType == 'SFit' :
    pdf = buildPdf( [ signalComps ], Observables = obsSetP2VV, Name = 'Jpsiphi' )

else :
    if components == 'signal' :
        pdf = buildPdf( [ signalComps                  ], Observables = obsSetP2VV, Name = 'JpsiphiSig' )
    elif components == 'background' :
        pdf = buildPdf( [ backgroundComps              ], Observables = obsSetP2VV, Name = 'JpsiphiBkg' )
    else :
        pdf = buildPdf( [ signalComps, backgroundComps ], Observables = obsSetP2VV, Name = 'Jpsiphi'    )


###########################################################################################################################################
## generate/read data and fit ##
################################

if generateData :
    # generate data
    print 'JvLFit: generating %d events' % nEvents
    fitData = pdf.generate(obsSetP2VV)

    from P2VVGeneralUtils import writeData
    writeData( dataSetFile, dataSetName, fitData )
elif fitType == 'SFit' :
    fitData = sigData
else :
    fitData = data

if doFit :
    # fix values of some parameters
    #lifetimeParams.setConstant('deltaM')
    #lifetimeParams.setConstant('dGamma')
    lambdaCP.setConstant('lambdaCPSq')
    #lambdaCP.setConstant('phiCP')

    for CEvenOdd in taggingParams['CEvenOdds'] :
        CEvenOdd.setConstant('avgCEven.*')
        CEvenOdd.setConstant('avgCOdd.*')
    taggingParams.setConstant('tagCatCoef.*')
    #tagCats.setConstant('wTag.*')

    #amplitudes.setConstant('f_S|ASPhase')
    #amplitudes.setConstant('.*')

    #timeResModel.setConstant('timeResSF')
    #backgroundTime.setConstant('.*')
    #signalBMass.setConstant('.*')
    #backgroundBMass.setConstant('.*')

    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    if fitType == 'SFit' : fitResult = pdf.fitTo( fitData, SumW2Error = True, **fitOpts )
    else                 : fitResult = pdf.fitTo( fitData,                    **fitOpts )

    print 120 * '=' + '\n'


###########################################################################################################################################
## make some plots ##
#####################

if makePlots :
    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed

    if fitType == 'SFit' :
        comps = None
    else :
        comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                }

    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yScale, logY )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSetP2VV[ : 5 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSetP2VV[ : 5 ] ]
                   , ( '', ) + angleNames
                   , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                   , ( True, ) + 3 * ( False, )
                  ) :
        plot(  pad, obs, fitData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth       )
             , components = comps
            )

    # plot lifetime
    timePlotTitles = tuple( [ time.GetTitle() + title for title in (  ' - linear'
                                                                    , ' - logarithmic'
                                                                    , ' - B/#bar{B} asymmetry'
                                                                   )
                            ] )
    timeCanv = TCanvas( 'timeCanv', 'Lifetime' )
    for ( pad, nBins, plotTitle, yTitle, yScale, dataCuts, pdfCuts, logY )\
            in zip(  timeCanv.pads( 2, 2 )
                   , numTimeBins
                   , timePlotTitles
                   , 2 * ( None, ) + ( 'B/#bar{B} asymmetry', )
                   , ( ( None, None ), ( 50., None ), ( None, None ) )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTag ), )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTag ), )
                   , ( False, True, False )
                  ) :
        plot(  pad, time, fitData, pdf, yTitle = yTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'            )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts        )
             , components = comps
            )

    # plot lifetime (tagged/untagged)
    timePlotTitles1 = tuple( [ time.GetTitle() + title for title in (  ' - untagged'
                                                                     , ' - tagging category 2'
                                                                     , ' - tagging category %d' % tagCat5Min
                                                                     , ' - B/#bar{B} asymmetry untagged'
                                                                     , ' - B/#bar{B} asymmetry tagging category 2'
                                                                     , ' - B/#bar{B} asymmetry tagging category %d' % tagCat5Min
                                                                    )
                            ] )
    timeCanv1 = TCanvas( 'timeCanv1', 'Lifetime' )
    for ( pad, nBins, plotTitle, yTitle, dataCuts, pdfCuts, logY )\
        in zip(  timeCanv1.pads( 3, 2 )
             , 3 * numTimeBins
             , timePlotTitles1
             , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
             ,   ( dict( Cut = '%s == 0'  % ( tagCatP2VV.GetName()             )                   ), )
               + ( dict( Cut = '%s == 2'  % ( tagCatP2VV.GetName()             )                   ), )
               + ( dict( Cut = '%s == %d' % ( tagCatP2VV.GetName(), tagCat5Min )                   ), )
               + ( dict( Cut = '%s == 0'  % ( tagCatP2VV.GetName()             ), Asymmetry = iTag ), )
               + ( dict( Cut = '%s == 2'  % ( tagCatP2VV.GetName()             ), Asymmetry = iTag ), )
               + ( dict( Cut = '%s == %d' % ( tagCatP2VV.GetName(), tagCat5Min ), Asymmetry = iTag ), )
             ,   ( dict( Slice = ( tagCatP2VV, 'Untagged'              )                   ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat2'               )                   ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat%d' % tagCat5Min )                   ), )
               + ( dict( Slice = ( tagCatP2VV, 'Untagged'              ), Asymmetry = iTag ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat2'               ), Asymmetry = iTag ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat%d' % tagCat5Min ), Asymmetry = iTag ), )
             , 3 * ( False, ) + 3 * ( False, )
            ) :
        plot(  pad, time, fitData, pdf, yTitle = yTitle, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'            )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts      )
             , components = comps
            )

    # plot angles
    anglePlotTitles =   tuple(  [ angle.GetTitle()                            for angle in angles ]\
                              + [ angle.GetTitle() + ' - B/#bar{B} asymmetry' for angle in angles ] )
    anglesCanv = TCanvas( 'anglesCanv', 'Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yTitle, dataCuts, pdfCuts )\
            in zip(  anglesCanv.pads( 3, 2 )
                   , 2 * angles
                   , 2 * numAngleBins
                   , anglePlotTitles
                   , 2 * angleNames
                   , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
                   , 3 * ( dict(), ) + 3 * ( dict( Asymmetry = iTag ), )
                   , 3 * ( dict(), ) + 3 * ( dict( Asymmetry = iTag ), )
                  ) :
        plot(  pad, obs, fitData, pdf, xTitle = xTitle, yTitle = yTitle
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                             )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts       )
             , components = comps
            )

    # print canvas to file
    timeAnglesCanv.Print(plotsFile + '(')
    timeCanv.Print(plotsFile)
    timeCanv1.Print(plotsFile)
    if fitType == 'SFit' :
      anglesCanv.Print(plotsFile + ')')
    else :
      anglesCanv.Print(plotsFile)
      bkgAnglesCanv.Print(plotsFile + ')')

