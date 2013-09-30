###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
makePlots   = True
physPdf     = True
tResModel   = ''
trigger     = ''
timeInt     = False

momentsFile = 'hel_UB_UT_trueTime_BkgCat050_KK30'
plotsFile   = 'hel_UB_UT_trueTime_BkgCat050_KK30.ps'
dataSetFile = 'hel_UB_UT_trueTime_BkgCat050_KK30.root'

nTupleName = 'DecayTree'
nTupleFile = 'data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20121010.root'

# transversity amplitudes
A0Mag2Val    = 0.60
AperpMag2Val = 0.16
AparMag2Val  = 1. - A0Mag2Val - AperpMag2Val

A0PhVal    =  0.
AperpPhVal = -0.17
AparPhVal  =  2.50

# CP violation parameters
phiCPVal      = -0.04

# B lifetime parameters
GammaVal  = 0.679
dGammaVal = 0.060
dMVal     = 17.8
tResSigma = 0.045

# plot options
numEffBins   = ( 20, 20, 20 )
numDistrBins = ( 40, 40, 40 )
lineWidth    = 3
markStyle    = 8
markSize     = 0.4

angleNames = ( 'cos#kern[0.1]{#theta_{K}}', 'cos#kern[0.1]{#theta_{l}}', '#varphi [rad]' )
effLabels  = (  '#int d_{}cos#theta_{#mu} d#varphi #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
              , '#int d_{}cos#theta_{K} d#varphi #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
              , '#int d_{}cos#theta_{K} dcos#theta_{#mu} #varepsilon_{#Omega}(#Omega) / (4 #LT#varepsilon_{#Omega}#GT)'
             )

###########################################################################################################################################
## create variables and read data ##
####################################

# import RooFit wrappers
from P2VV.Load import RooFitOutput

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'angEff' ).ws()

# angular functions
from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

# variables in PDF
from P2VV.RooFitWrappers import RealVar, Category
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

from P2VV.Utilities.DataHandling import readData
cuts = bkgcatCut + ' && ' + trackChiSqCuts + ' && ' + massCuts + ' && ' + timeCuts + ' && ' + tagCuts
if trigger == 'ExclBiased' :
    cuts  = 'sel == 1 && sel_cleantail==1 && hlt1_excl_biased_dec == 1 && hlt2_biased == 1 && ' + cuts
    data = readData(  nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSet, ntupleCuts = cuts )

elif trigger == 'Unbiased' :
    cuts = 'sel == 1 && sel_cleantail==1 && hlt1_unbiased_dec == 1 && hlt2_biased == 1 && ' + cuts
    data = readData(  nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSet, ntupleCuts = cuts )

else :
    cuts = 'sel == 1 && sel_cleantail==1 && (hlt1_unbiased_dec == 1 || hlt1_biased == 1) && hlt2_biased == 1 && ' + cuts
    data = readData(  nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSet, ntupleCuts = cuts )


###########################################################################################################################################
## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
#####################################################################

if physPdf :
    # transversity amplitudes
    from P2VV.Parameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
    amplitudes = Amplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                            , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                            , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                            , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                            , ReAS    = 0.
                            , ImAS    = 0.
                           )

    # B lifetime
    from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
    lifetimeParams = LifetimeParams( Gamma = GammaVal, dGamma = dGammaVal, dM = dMVal )

    tResArgs = { }
    if tResModel == 'Gauss' :
        from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
        tResArgs['time']         = time
        tResArgs['timeResSigma'] = tResSigma
    elif tResModel == '3Gauss' :
        from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
        tResArgs['time'] = time
    else :
        from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
        tResArgs['time'] = trueTime
    timeResModel = TimeResolution( **tResArgs )

    # CP violation parameters
    from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
    lambdaCP = CPParam( lambdaCPSq = 1., phiCP = phiCPVal )

    # tagging parameters
    from P2VV.Parameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
    taggingParams = TaggingParams()

    # coefficients for time functions
    from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
    timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp' ] )

    # build underlying physics PDF
    args = dict(  time            = time if tResModel in [ 'Gauss', '3Gauss' ] else trueTime
                , iTag            = iTag
                , tau             = lifetimeParams['MeanLifetime']
                , dGamma          = lifetimeParams['dGamma']
                , dm              = lifetimeParams['dM']
                , dilution        = taggingParams['dilution']
                , ADilWTag        = taggingParams['ADilWTag']
                , avgCEven        = taggingParams['avgCEven']
                , avgCOdd         = taggingParams['avgCOdd']
                , coshCoef        = timeBasisCoefs['cosh']
                , sinhCoef        = timeBasisCoefs['sinh']
                , cosCoef         = timeBasisCoefs['cos']
                , sinCoef         = timeBasisCoefs['sin']
                , resolutionModel = timeResModel['model']
               )

    from P2VV.RooFitWrappers import BTagDecay
    pdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )


###########################################################################################################################################
## compute angular efficiency moments ##
########################################

# print PDF, data, variables and parameters
print '\nData set:'
data.Print()

if physPdf :
    print '\nPDF:'
    pdf.Print()
    print '\nNormalization set for PDF:'
    for var in angles : print var.GetName(),
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
from P2VV.Utilities.DataMoments import RealMomentsBuilder
if physPdf :
    from P2VV.RooFitWrappers import RealEffMoment
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                Norm = 1., PDF = pdf, IntSet = [ ], NormSet = angles )\
                                                  for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func
                                                )
                                    )

else :
    from P2VV.RooFitWrappers import RealMoment
    physMoments = RealMomentsBuilder( Moments = ( RealMoment( Name = func.GetName(), BasisFunc = func, Norm = 1. )\
                                                  for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func
                                                )
                                    )

physMoments.initCovariances()

# moments builder with angular basis functions
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
indices += [ ( 0, 4, 0 ), ( 0, 4, 2 ), ( 0, 4, 4 ) ]

basisMoments = RealMomentsBuilder()
if physPdf :
    basisMoments.appendPYList( angleFuncs.angles, indices, PDF = pdf, IntSet = [ ], NormSet = angles )
else :
    basisMoments.appendPYList( angleFuncs.angles, indices )

PDFInt = 1. if physPdf else 8. * pi

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
## add efficiency weights column to data set ##
###############################################

if makePlots and physPdf :
    # add column with efficiency weights to data set
    from ROOT import RooArgSet, RooArgList, RooFormulaVar, RooConstVar
    rooIntSet  = RooArgSet()
    rooNormSet = RooArgSet( var._var for var in angles )
    pdfInt = pdf.createIntegral( rooIntSet, rooNormSet )
    nEvents = RooConstVar( 'nEvents', 'Number of events', data.sumEntries() )
    effWeightList = RooArgList( nEvents, pdfInt )
    effWeight = RooFormulaVar( 'effWeight', 'Efficiency weight', '1./@0/@1', effWeightList )
    data.addColumn(effWeight)

    # apply efficiency weights to events in data set
    from ROOT import RooDataSet
    weightedData = RooDataSet( 'effWeightData', 'effWeightData', data.get(), Import = data, WeightVar = ( 'effWeight', True ) )


###########################################################################################################################################
## write data set to file ##
############################

if dataSetFile :
    from ROOT import TFile
    dataFile = TFile.Open( dataSetFile, 'RECREATE' )
    dataFile.Add(data)
    dataFile.Write()
    dataFile.Close()


###########################################################################################################################################
## make some plots ##
#####################

if makePlots :
    # build efficiency function
    effTerms = basisMoments.buildPDFTerms()
    effFunc = effTerms.buildAddition('effFunc')

    from ROOT import RooArgSet
    ctkSet = RooArgSet( angles[1], angles[2] )
    ctlSet = RooArgSet( angles[0], angles[2] )
    phiSet = RooArgSet( angles[0], angles[1] )
    effInts = (  effFunc.createIntegral( ctkSet, RooArgSet() )
               , effFunc.createIntegral( ctlSet, RooArgSet() )
               , effFunc.createIntegral( phiSet, RooArgSet() )
              )

    # import plotting tools
    from P2VV.Load import LHCbStyle
    from P2VV.Utilities.Plotting import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kFullDotLarge
    canvs = [ ]

    # LHCb label
    from ROOT import TPaveText
    LHCbText = TPaveText( 0.31, 0.83, 0.67, 0.92, 'NDC' )
    LHCbText.AddText('LHCb simulation')
    LHCbText.SetShadowColor(0)
    LHCbText.SetFillStyle(0)
    LHCbText.SetBorderSize(0)
    LHCbText.SetTextAlign(12)

    # plot (binned) efficiency function integrals
    canvs += [  TCanvas( 'cpsiDataIntCanv',   'Angular Efficiency' )
              , TCanvas( 'cthetaDataIntCanv', 'Angular Efficiency' )
              , TCanvas( 'phiDataIntCanv',    'Angular Efficiency' )
             ]
    for ( pad, func, angle, xTitle, yTitle, yScale, norm, nBins )\
        in zip(  canvs[ : 3 ]
               , effInts
               , angles
               , angleNames
               , effLabels
               , [ ( 0.88, 1.12 ), ( 0.9328, 1.1872 ), ( 0.88, 1.12 ) ]
               , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
               , numEffBins
              ) :
        pad.SetLeftMargin(0.28)
        plot(  pad, angle, weightedData, func
             , xTitle       = xTitle
             , yTitle       = yTitle
             , yScale       = yScale
             , yTitleOffset = 1.0
             , frameOpts    = dict( Title = angle.GetTitle(), Bins = nBins )
             , dataOpts     = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.9, LineWidth = 3
                                   , Rescale = norm * float(nBins) / ( angle.getMax() - angle.getMin() ) )
             , pdfOpts      = dict( LineColor = kBlue, LineWidth = 4, Normalization = norm )
            )
        LHCbText.Draw()

    if physPdf :
        # plot lifetime and angles
        canvs.append( TCanvas( 'ctkCanv',  'cos(theta_K)' ) )
        canvs.append( TCanvas( 'ctlCanv',  'cos(theta_l)' ) )
        canvs.append( TCanvas( 'phiCanv',  'phi'          ) )
        for ( pad, obs, nBins, plotTitle, xTitle )\
                in zip(  canvs[ -3 : ]
                       , obsSet[ 1 : 5 ]
                       , numDistrBins
                       , [ var.GetTitle() for var in obsSet[ 1 : 5 ] ]
                       , angleNames
                      ) :
            plot(  pad, obs, data, pdf, xTitle = xTitle
                 , frameOpts   = dict( Bins = nBins, Title = plotTitle                )
                 , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize )
                 , pdfOpts     = dict( LineColor = kBlue, LineWidth = lineWidth       )
                )

        # multiply PDF with angular efficiency
        basisMoments * pdf

        # plot efficiency PDF
        from P2VV.Utilities.Plotting import _P2VVPlotStash
        global _P2VVPlotStash
        for it in range(3) :
            pdf.plotOn( _P2VVPlotStash[ -it - 1 ], Name = 'effPdf', LineColor = kRed, LineWidth = lineWidth )
            canvs[ -it - 1 ].cd()
            _P2VVPlotStash[ -it - 1 ].Draw()

    # print canvases to file
    for canvIt, canv in enumerate(canvs) : canv.Print( plotsFile + ( '(' if canvIt == 0 else ')' if canvIt == len(canvs) - 1 else '' ) )
