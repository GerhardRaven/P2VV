###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
MCProd      = 'Sim08_2012'
trueAngles  = True
readMoments = False
makePlots   = False
transAngles = False
incSWave    = False
normPdf     = True
tResModel   = '' # '3fb'
trigger     = ''
timeInt     = False
addInvPdf   = False
weightVar   = 'wKin' # 'sWeights_ipatia'
doSelection = False
blind       = { } # { 'phiCP' : ( 'UnblindUniform', 'BsPhisCombination', 0.2 ), 'dGamma' : ( 'UnblindUniform', 'BsDGsCombination', 0.02 ) }
parFileIn   = '' # '../it6/fitPars.par'
customSel   = True
timeRange   = (0., 25.) if customSel else (0.3, 14.)
sigmatCut   = False
blind       = { } # { 'phiCP' : ( 'UnblindUniform', 'BsPhisCombination', 0.2 ), 'dGamma' : ( 'UnblindUniform', 'BsDGsCombination', 0.02 ) }
parFileIn   = '' # '../it6/fitPars.par'
trueTimeVarName = 'truetime'


momentsFile = '%s_%s_UB_UT_%s_BkgCat050_KK30'      % ( MCProd, 'trans' if transAngles else 'hel', tResModel if tResModel else 'trueTime' )
plotsFile   = '%s_%s_UB_UT_%s_BkgCat050_KK30.ps'   % ( MCProd, 'trans' if transAngles else 'hel', tResModel if tResModel else 'trueTime' )
dataSetFile = '' # '%s_%s_UB_UT_%s_BkgCat050_KK30.root' % ( MCProd, 'trans' if transAngles else 'hel', tResModel if tResModel else 'trueTime' )

nTupleName = 'DecayTree'
if MCProd == 'real' :
    nTupleFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/nTupleC_w0_add.root'
elif MCProd == 'Sim08' :
    #nTupleFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/Bs2JpsiPhi_20112012_Sim08_ntupleB_201309_add.root'
    nTupleFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/reweighted/Bs2JpsiPhi_Sim08a_reWeighted.root'
elif MCProd == 'Sim08_2011' :
    nTupleFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_angEff.root' if not trueAngles else\
        '/project/bfys/vsyropou/PhD/workdir/jpsiphiAngAccChecks/DTT_stripingAndBkg60/MC2011_all_BsJpsiPhi_DTT_20131006_minimal_truth_roel.root'    
    # '/project/bfys/vsyropou/data/Bs2JpsiPhi/dataSets/mcWithTrueAngles/Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_trueAngles.root'
elif MCProd == 'Sim08_2012' :
    nTupleFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/Bs2JpsiPhi_MC2012_ntupleB_20130904_angEff.root' if not trueAngles else\
        '/project/bfys/vsyropou/PhD/workdir/jpsiphiAngAccChecks/DTT_stripingAndBkg60/MC2012_all_BsJpsiPhi_DTT_20130709_minimal_truth_roel_v2.root'
        #'/project/bfys/vsyropou/data/Bs2JpsiPhi/dataSets/mcWithTrueAngles/Bs2JpsiPhi_MC2012_ntupleB_20130909_trueAngles.root'
else :
    nTupleFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20121010.root'

# transversity amplitudes
A0Mag2Val    = 0.722**2 / (0.722**2 + 0.480**2 + 0.499**2) if MCProd.startswith('Sim08') else 0.5208 if MCProd == 'real'\
               else 0.775**2 / (0.775**2 + 0.490**2 + 0.400**2)
AperpMag2Val = 0.499**2 / (0.722**2 + 0.480**2 + 0.499**2) if MCProd.startswith('Sim08') else 0.2535 if MCProd == 'real'\
               else 0.400**2 / (0.775**2 + 0.490**2 + 0.400**2)
AparMag2Val  = 0.480**2 / (0.722**2 + 0.480**2 + 0.499**2) if MCProd.startswith('Sim08') else 1. - 0.5208 - 0.2535 if MCProd == 'real'\
               else 0.490**2 / (0.775**2 + 0.490**2 + 0.400**2)

A0PhVal    = 0.
AperpPhVal = 3.07 if MCProd.startswith('Sim08') else 3.19 if MCProd == 'real' else -0.17
AparPhVal  = 3.30 if MCProd.startswith('Sim08') else 3.26 if MCProd == 'real' else  2.50

KKMassStates = dict( [ ('bin%d' % i, i) for i in range(6) ] )
if incSWave :
    SWaveAmps = dict(  f_S = dict(  bin0 = 0.445  if MCProd == 'real' else 0.
                                  , bin1 = 0.064  if MCProd == 'real' else 0.
                                  , bin2 = 0.0087 if MCProd == 'real' else 0.
                                  , bin3 = 0.0089 if MCProd == 'real' else 0.
                                  , bin4 = 0.045  if MCProd == 'real' else 0.
                                  , bin5 = 0.210  if MCProd == 'real' else 0.
                                 )
                     , ASOddPhase = dict(  bin0 = +0.80 if MCProd == 'real' else 0.
                                         , bin1 = +2.32 if MCProd == 'real' else 0.
                                         , bin2 = +0.46 if MCProd == 'real' else 0.
                                         , bin3 = -0.36 if MCProd == 'real' else 0.
                                         , bin4 = -0.68 if MCProd == 'real' else 0.
                                         , bin5 = -0.82 if MCProd == 'real' else 0.
                                        )
                     , C_SP = dict(  bin0 = 0.966 if MCProd == 'real' else 1.
                                   , bin1 = 0.956 if MCProd == 'real' else 1.
                                   , bin2 = 0.926 if MCProd == 'real' else 1.
                                   , bin3 = 0.926 if MCProd == 'real' else 1.
                                   , bin4 = 0.956 if MCProd == 'real' else 1.
                                   , bin5 = 0.966 if MCProd == 'real' else 1.
                                  )
                    )
else :
    SWaveAmps = { }

# CP violation parameters
phiCPVal  = +0.07 if MCProd.startswith('Sim08') else +0.0827 if MCProd == 'real' else -0.04
lambCPVal = 0.964 if MCProd == 'real' else 1.

# B lifetime parameters
GammaVal  = 1. / 1.503 if MCProd.startswith('Sim08') else 0.6606 if MCProd == 'real' else 0.679
dGammaVal = 1. / 1.406 - 1. / 1.614 if MCProd.startswith('Sim08') else 0.087 if MCProd == 'real' else 0.060
dMVal     = 17.8
tResSigma = 0.045

# resolution model parameters
runPeriodStates = dict( p2011 = 2011, p2012 = 2012 )
if tResModel == '3fb' :
    tResPars = dict(  tres_placeholder = dict( p2011 =  0.0349,  p2012 =  0.0347  )
                    , timeResMu        = dict( p2011 = -0.00259, p2012 = -0.00333 )
                    , timeResFrac2     = dict( p2011 =  0.242,   p2012 =  0.239   )
                    , sf_mean_offset   = dict( p2011 =  1.4273,  p2012 =  1.4887  )
                    , sf_mean_slope    = dict( p2011 = -1.93,    p2012 = -3.88    )
                    , sf_sigma_offset  = dict( p2011 =  0.3857,  p2012 =  0.4143  )
                    , sf_sigma_slope   = dict( p2011 = -0.63,    p2012 = -2.80    )
                   )
else :
    tResPars = { }

# plot options
if transAngles : angleNames = ( 'cos(#psi_{tr})',  'cos(#theta_{tr})',  '#phi_{tr}' )
else           : angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{#mu})', '#phi_{h}'  )
numBins         = ( 60, 60, 60 )
lineWidth       = 2
markStyle       = 8
markSize        = 0.4


###########################################################################################################################################
## create variables and read data ##
####################################

# import RooFit wrappers
from P2VV.Load import RooFitOutput

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'ws').ws()

# angular functions
if transAngles :
    from P2VV.Parameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )
    #angleFuncs = AngleFuncs( cpsi = 'ctkRecoLHCb', ctheta = 'ctlRecoLHCb', phi = 'phiRecoLHCb' )
    #angleFuncs = AngleFuncs( cpsi = 'ctkTrueCalc', ctheta = 'ctlTrueCalc', phi = 'phiTrueCalc' )

# variables in PDF
from P2VV.RooFitWrappers import RealVar, Category
time      = RealVar(  'time',      Title = 'Decay time',      Unit = 'ps', Observable = True, Value = 0.5,  MinMax = ( 0.3,    14.  ) )
trueTime  = RealVar(  'truetime',  Title = 'True decay time', Unit = 'ps', Observable = True, Value = 0.,   MinMax = ( 0.,     20.  ) )
sigmat    = RealVar(  'sigmat',    Title = 'Time resolution', Unit = 'ps', Observable = True, Value = 0.01, MinMax = ( 0.0001, 0.12 ) )
iTag      = Category( 'iTag',      Title = 'Initial state flavour tag', Observable = True, States = { 'Untagged' : 0 } )
runPeriod = Category( 'runPeriod', Title = 'Running period',            Observable = True, States = runPeriodStates )
KKMassCat = Category( 'KKMassCat', Title = 'KK-mass category',          Observable = True, States = KKMassStates )
angles    = [ angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] ]
if weightVar :
    weight = RealVar( weightVar, Title = 'Signal weight', Observable = True, Value = 1. )

sigmat.setBins( 40, 'cache' )

obsSet = [ time if tResModel in [ 'Gauss', '3fb' ] else trueTime ] + angles + [ runPeriod, KKMassCat ]
if tResModel == '3fb' : obsSet.append(sigmat)
if tResPars : obsSet.append(runPeriod)
if SWaveAmps : obsSet.append(KKMassCat)
if weightVar : obsSet.append(weight)

# read ntuple
bkgcatCut      = '(bkgcat == 0 || bkgcat == 50)'
trackChiSqCuts = 'muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4. && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
massCuts       = 'mass > 5200. && mass < 5550. && mdau1 > 3030. && mdau1 < 3150. && mdau2 > 990. && mdau2 < 1050.'
trueTimeCut    = 'truetime > 0.'
timeCuts       = 'time > 0.3 && time < 14. && sigmat < 0.12'
tagCuts        = '(tagdecision == 0 || tagdecision == -1 || tagdecision == +1)'

cuts = ''
if doSelection :
    cuts = trackChiSqCuts + ' && ' + massCuts + ' && ' + timeCuts + ' && ' + tagCuts
    if MCProd != 'real' :
        cuts += ' && ' + bkgcatCut + ' && ' + trueTimeCut
    if trigger == 'ExclBiased' :
        #cuts = 'sel == 1 && sel_cleantail==1 && hlt1_excl_biased_dec == 1 && hlt2_biased == 1 && ' + cuts
        cuts = 'sel == 1 && sel_cleantail==1 && hlt1_excl_biased == 1 && hlt2_biased == 1 && ' + cuts
    elif trigger == 'Unbiased' :
        #cuts = 'sel == 1 && sel_cleantail==1 && hlt1_unbiased_dec == 1 && hlt2_biased == 1 && ' + cuts
        cuts = 'sel == 1 && sel_cleantail==1 && hlt1_unbiased == 1 && hlt2_biased == 1 && ' + cuts
    else :
        #cuts = 'sel == 1 && sel_cleantail==1 && (hlt1_unbiased_dec == 1 || hlt1_biased == 1) && hlt2_biased == 1 && ' + cuts
        cuts = 'sel == 1 && sel_cleantail==1 && (hlt1_unbiased == 1 || hlt1_biased == 1) && hlt2_biased == 1 && ' + cuts

# customise selection
if customSel:
    cuts = ''
    cuts += 'abs(B_s0_TRUEID)==531'
    cuts += ' && truetime>%s && truetime<%s'%(timeRange[0],timeRange[1])

readDataOpts = { }
if cuts      : readDataOpts['ntupleCuts'] = cuts
if weightVar : readDataOpts['WeightVar']  = ( weightVar, True )

from P2VV.Utilities.DataHandling import readData
data = readData( nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSet, **readDataOpts )


###########################################################################################################################################
## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
#####################################################################

if normPdf :
    # transversity amplitudes
    from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
    amplitudes = Amplitudes(  A0Mag2 = A0Mag2Val, AperpMag2 = AperpMag2Val
                            , AparPhase = AparPhVal - A0PhVal, AperpPhase = AperpPhVal - A0PhVal
                            , f_S = 0., ASOddPhase = 0., C_SP = 1.
                           )

    # B lifetime
    dGammaVar = dict( Value = dGammaVal )
    if blind and 'dGamma' in blind : dGammaVar['Blind'] = blind['dGamma']
    from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
    lifetimeParams = LifetimeParams( Gamma = GammaVal, dGamma = dGammaVar, dM = dMVal )

    tResArgs = { }
    if tResModel == 'Gauss' :
        from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
        tResArgs['time']         = time
        tResArgs['timeResSigma'] = tResSigma
    elif tResModel == '3fb' :
        from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
        tResArgs['time'] = time
        tResArgs['sigmat'] = sigmat
        tResArgs['Cache'] = True
        tResArgs['PerEventError'] = True
        tResArgs['Parameterise'] = 'RMS'
        tResArgs['TimeResSFParam'] = 'linear'
    else :
        from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
        tResArgs['time'] = trueTime
    timeResModel = TimeResolution( **tResArgs )

    # CP violation parameters
    phiCPVar = dict( Value = phiCPVal )
    if blind and 'phiCP' in blind : phiCPVar['Blind'] = blind['phiCP']
    from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
    lambdaCP = CPParam( lambdaCP = lambCPVal, phiCP = phiCPVar )

    # tagging parameters
    from P2VV.Parameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
    taggingParams = TaggingParams()

    # coefficients for time functions
    from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
    timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] )

    # build signal PDF
    assert all( timeBasisCoefs[comp][ind] and not timeBasisCoefs[comp][abs(ind - 1)]\
           for comp, ind in [ ( 'cosh', 0 ), ( 'sinh', 0 ), ( 'cos', 1 ), ( 'sin', 1 ) ] )\
           , 'in the standard RooBTagDecay it is assumed that the "cosh" and "sinh" coefficients don\'t flip sign between B and Bbar and the "cos" and "sin" coefficients do'
    args = dict(  time                   = time if tResModel in [ 'Gauss', '3fb' ] else trueTime
                , iTag                   = iTag
                , tau                    = lifetimeParams['MeanLifetime']
                , dGamma                 = lifetimeParams['dGamma']
                , dm                     = lifetimeParams['dM']
                , dilution               = taggingParams['dilution']
                , ADilWTag               = taggingParams['ADilWTag']
                , avgCEven               = taggingParams['avgCEven']
                , avgCOdd                = taggingParams['avgCOdd']
                , coshCoef               = timeBasisCoefs['cosh'][0]
                , sinhCoef               = timeBasisCoefs['sinh'][0]
                , cosCoef                = timeBasisCoefs['cos'][1]
                , sinCoef                = timeBasisCoefs['sin'][1]
                , resolutionModel        = timeResModel['model']
               )

    from P2VV.RooFitWrappers import BTagDecay
    protoPdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

    intSet  = [ ]
    normSet = angleFuncs.angles.values()
    if timeInt :
        # integrate over time, normalize over time + angles
        intSet.append( time if tResModel in [ 'Gauss', '3fb' ] else trueTime )
        normSet.append( intSet[-1] )

    splitCats = [ ]
    splitPars = [ ]
    if tResPars :
        splitCats.append( [ runPeriod ] )
        splitPars.append( [ ws[par] for par in tResPars.iterkeys() ] )
    if SWaveAmps :
        splitCats.append( [ KKMassCat ] )
        splitPars.append( [ ws[par] for par in SWaveAmps.iterkeys() ] )
    if splitCats :
        # split PDF for run period and KK-mass bins
        from P2VV.RooFitWrappers import SimultaneousPdf
        pdf = SimultaneousPdf( protoPdf.GetName() + '_simul', MasterPdf = protoPdf
                              , SplitCategories = splitCats, SplitParameters = splitPars )

        # set values of split parameters
        for par, vals in tResPars.iteritems() :
            for cat, val in vals.iteritems() : ws[ par + '_' + cat ].setVal(val)
        for par, vals in SWaveAmps.iteritems() :
            for cat, val in vals.iteritems() : ws[ par + '_' + cat ].setVal(val)

    else :
        pdf = protoPdf

if parFileIn :
    # read parameter values from file
    from P2VV.Parameterizations.FullPDFs import PdfConfiguration
    pdfConfig = PdfConfiguration()
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)


###########################################################################################################################################
## compute angular efficiency moments ##
########################################

# print PDF, data, variables and parameters
print '\nData set:'
data.Print()

if normPdf :
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
from P2VV.Utilities.DataMoments import RealMomentsBuilder
if normPdf :
    from P2VV.RooFitWrappers import RealEffMoment
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,Norm = 1., PDF = pdf
                                                                , IntSet = intSet, NormSet = normSet )\
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
#indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(4) for YIndex0 in range(4) for YIndex1 in range( -YIndex0, YIndex0 + 1 )\
#             if PIndex == 3 or YIndex0 == 3 ]
#indices  = [ ( 0, YIndex0, YIndex1 ) for YIndex0 in range(6) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#indices = [ ( PIndex, 2, YIndex1 ) for PIndex in range(40) for YIndex1 in [ +1, -1 ] ]

basisMoments = RealMomentsBuilder()
if normPdf :
    basisMoments.appendPYList( angleFuncs.angles, indices, PDF = pdf, IntSet = intSet, NormSet = normSet )
else :
    basisMoments.appendPYList( angleFuncs.angles, indices )

PDFInt = 1. if normPdf else 8. * pi
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
physMoments.convertEffWeightsToMoments( OutputFilePath = momentsFile + '_Basis_weights', Scale = PDFInt / 16. / sqrt(pi) )
basisMoments.Print( Scale = PDFInt /  2. / sqrt(pi)                       )
basisMoments.Print( Scale = PDFInt /  2. / sqrt(pi), MinSignificance = 5. )


###########################################################################################################################################
## add efficiency weights column to data set ##
###############################################

if addInvPdf and normPdf :
    print 'adding efficiency weights to dataset'
    from ROOT import RooArgSet, RooArgList, RooFormulaVar, RooConstVar
    rooIntSet  = RooArgSet( var._var for var in intSet  )
    rooNormSet = RooArgSet( var._var for var in normSet )
    pdfInt = pdf.createIntegral( rooIntSet, rooNormSet )
    nEvents = RooConstVar( 'nEvents', 'Number of events', data.sumEntries() )
    effWeightList = RooArgList( nEvents, pdfInt )
    effWeight = RooFormulaVar( 'effWeight', 'Efficiency weight', '1./@0/@1', effWeightList )
    data.addColumn(effWeight)


###########################################################################################################################################
## write data set to file ##
############################

if dataSetFile :
    print 'converting dataset to n-tuple'
    from ROOT import TFile, TObject
    dataFile = TFile.Open( dataSetFile, 'RECREATE' )
    nTuple = data.buildTree( Name = nTupleName, Title = nTupleName, WeightName = weightVar, RooFitFormat = False )
    print 'writing n-tuple to file %s' % dataSetFile
    dataFile.Write( dataSetFile, TObject.kOverwrite )
    dataFile.Close()


###########################################################################################################################################
## make some plots ##
#####################

if makePlots and normPdf :
    # build efficiency function
    effTerms = basisMoments.buildPDFTerms()
    effFunc = effTerms.buildAddition('effFunc')

    from ROOT import RooArgSet
    ctkSet = RooArgSet( angles[1], angles[2] )
    ctlSet = RooArgSet( angles[0], angles[2] )
    phiSet = RooArgSet( angles[0], angles[1] )
    effFuncCtk = effFunc.createIntegral( ctkSet, RooArgSet() )
    effFuncCtl = effFunc.createIntegral( ctlSet, RooArgSet() )
    effFuncPhi = effFunc.createIntegral( phiSet, RooArgSet() )

    # import plotting tools
    from P2VV.Load import LHCbStyle
    from P2VV.Utilities.Plotting import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen

    # plot efficiency
    effCanv = TCanvas( 'effCanv', 'Efficiency' )
    for ( pad, obs, func, norm, xTitle )\
            in zip(  effCanv.pads( 2, 2 )
                   , obsSet[ 1 : 4 ]
                   , [ effFuncCtk, effFuncCtl, effFuncPhi ]
                   , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
                   , angleNames
                  ) :
        plot(  pad, obs, None, func, addPDFs = [ effFunc ]
             , yScale      = ( 0.85, 1.15 )
             , pdfOpts     = dict( LineColor = kBlue, LineWidth = 2, Normalization = norm )
             , addPDFsOpts = [ dict( LineColor = kRed,  LineWidth = 2 ) ]
            )

    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSet[ 1 : 4 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSet[ 1 : 4 ] ]
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
        timeAnglesCanv.cd( 3 - it )
        _P2VVPlotStash[ -it - 1 ].Draw()

    # print canvases to file
    effCanv.Print( plotsFile + '(' )
    timeAnglesCanv.Print( plotsFile + ')' )
