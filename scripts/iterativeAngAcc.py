########################################################################################################################
## Configuration ##
################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-R', '--KKmomRew',  dest='KKmomRew',   default = 'vertical',            help='KK momentum reweighting approach (vertical/horizontal)')
parser.add_option('-o', '--rewOrder',  dest='rewOrder',   default = 'Bmom mkk phys KKmom', help='reweghting steps order')
parser.add_option('-s', '--MCProd',    dest='MCProd',     default = '2011_2012',           help='choose mc sample ( 2011/2012/2011+2012 )')
parser.add_option('-n', '--numIters',  dest='numIters',   default = 1, type=int,           help='number of iterations')
parser.add_option('-a', '--nomAngAcc', dest='nomAngAcc',  default = '',                    help='nominal angular acceptance')
parser.add_option('-d', '--physPars',  dest='physPars',   default = '',                    help='physics aprameters')
(options, args) = parser.parse_args()

# reweightng flow control
initialFitOnData      = False
NumbOfIterations      = options.numIters
RewApproach           = options.KKmomRew
MCProd                = options.MCProd
reweightBmomentum     = True if 'Bmom'  in options.rewOrder else False
reweightMkk           = True if 'mkk'   in options.rewOrder else False
reweightPhysics       = True if 'phys'  in options.rewOrder else False
reweightKKmom         = True if 'KKmom' in options.rewOrder else False
OneDverticalRewNbins  = 1000
TwoDverticalRewNbins  = 50
EqualStatsBins        = True
physWeightName        = 'phys'
mKKWeightsName        = 'mKK'
KmomentaWeightsName   = 'KKmom'
BmomentumWeightsName  = 'Bmom'

# fit configuration
doFit       = True
nCPU        = 8

# plotig configuration
makePlots          = False
plotAtTheseSteps   = [ NumbOfIterations ] # plot only at the last iteration
plotFinalPdfonData = False

# source distribution
dataSetsPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/angEff/'
monteCarloPath  = dataSetsPath
mcTupleName     = 'DecayTree'
monteCarloData  = {'2011' : dataSetsPath + 'Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_angEff.root',
                   '2012' : dataSetsPath + 'Bs2JpsiPhi_MC2012_ntupleB_20130904_angEff.root'
                   }

# target distribution
from ROOT import TFile 
data2011Path = dataSetsPath + 'P2VVDataSets2011Reco14_I2Mass_6KKMassBins_2TagCats_kinematics_HLT2B.root'
data2012Path = dataSetsPath + 'P2VVDataSets2012Reco14_I2Mass_6KKMassBins_2TagCats_kinematics_HLT2B.root'
fitDataPath  = dataSetsPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
sDataName    = 'JpsiKK_sigSWeight'
sWeightsName = 'sWeights_ipatia'
data2011     = TFile.Open(data2011Path).Get(sDataName)
data2012     = TFile.Open(data2012Path).Get(sDataName)

# nominal angluar acceptance 
#nomAngEffMomentsFile     = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/uncorrecteEffMoments/'
nomAngEffMomentsFile = options.nomAngAcc if options.nomAngAcc \
    else '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
outputEffMomentsFileName = 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys'

# source generating physics  parameters
## TODO:: set mc pars from a file just like setdata pars
from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as monteCarloParameters
#/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/nominalFitResults/monteCarloSim08parameterValues.par

# target physics parameters
dataParameters = options.physPars if options.physPars\
    else'/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/nominalFitResults/realDataNominalFitResult.par'

###########################################################################################################################
## Begin iterative procedure  ##
################################
if RewApproach == 'vertical':  
    from P2VV.Utilities.MCReweighting import TwoDimentionalVerticalReweighting, OneDimentionalVerticalReweighting
elif RewApproach == 'horizontal': 
    from P2VV.Utilities.MCReweighting import MatchWeightedDistributions
from P2VV.Utilities.MCReweighting import MatchPhysics, BuildBs2JpsiKKFit
from P2VV.Utilities.MCReweighting import compareDistributions, cleanP2VVPlotStash, WeightedDataSetsManager
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from P2VV import RooFitDecorators
from math import pi, sqrt
import gc

# define a workspace
worksp = RooObject( workspace = 'iterativeProcedure' ).ws()

# build data pdf and prepare the sFit ( This pdf will not be multiplied by the angular acceptance !! ).
Bs2JpsiKKFit = BuildBs2JpsiKKFit( dataSetPath=fitDataPath, dataSetName=sDataName, Ncpu=nCPU )

if initialFitOnData:
    Bs2JpsiKKFit.doFit( angAccFile=nomAngEffMomentsFile )
    assert False

# initialise physics matching class and build MC pdf
PhysicsReweight = MatchPhysics( monteCarloData, mcTupleName , MonteCarloProduction=MCProd, BlindPdf=Bs2JpsiKKFit.getBlindString() )

# keep track of the weights (avoid creating too many datasets)
PhysicsReweight.selectDataSet('2011')
mcDataMngr2011 = WeightedDataSetsManager( source = PhysicsReweight.getDataSet() )
PhysicsReweight.selectDataSet('2012')
mcDataMngr2012 = WeightedDataSetsManager( source = PhysicsReweight.getDataSet() )

# get observables
angles     = [worksp[o] for o in ['helcosthetaK','helcosthetaL','helphi']]
time       = [worksp['time']]
truetime   = [worksp['truetime']]
muMomenta  = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['muplus','muminus'] for comp in ('P','PX','PY','PZ') ]    ]
Kmomenta   = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['Kplus','Kminus']   for comp in ('P','PX','PY','PZ') ]    ]
Bmomenta   = [ worksp['B_P'], worksp['B_Pt'] ]
KKMass     = [worksp['mdau2']]
KKMassCat  = Bs2JpsiKKFit.getPdf().indexCat()

# initialise horizontal kinematic reweighting classs
if RewApproach == 'horizontal':
    KinematicReweight = MatchWeightedDistributions( outTree        = Bs2JpsiKKFit.getDataSet(), # Target: Distribution to be matched with
                                                    reweightVars   = ['Kminus_P'],              # Variables that enter the transformation
                                                    inWeightName   = physWeightName,
                                                    outWeightName  = sWeightsName,
                                                    observables    = angles + time + truetime,
                                                    nonObsVars     = muMomenta + Kmomenta + Bmomenta + KKMass,
                                                    nBins          = 1000                       # preceision of the transformation
                                                    )


# match B momentum
if reweightBmomentum:
    for source, target, dataMngr in zip( [mcDataMngr2011.getDataSet(),mcDataMngr2012.getDataSet()], 
                                         [data2011, data2012], 
                                         [mcDataMngr2011,mcDataMngr2012] 
                                         ):
        BmomentumWeights = TwoDimentionalVerticalReweighting( source, target, TwoDverticalRewNbins, ['B_P','B_Pt'],
                                                              equalStatsBins = EqualStatsBins,
                                                              )
        dataMngr['saveIntermediateDatasets'] = True if makePlots else False
        dataMngr.appendWeights( BmomentumWeightsName, BmomentumWeights, permanetnWeigts=True, permanentDataSet=False )
assert False

# start looping.
for iterNumb in range( 1, NumbOfIterations + 1 ):
    print 'P2VV - INFO: Iteratitive procedure, begining of iteration %s.'%str(iterNumb)
    dataMngr['iterationNumber'] = iterNumb
    dataMngr['saveIntermediateDatasets'] = True if iterNumb in plotAtTheseSteps and makePlots else False

    # match mKK
    if reweightMkk:
        mKKweights = OneDimentionalVerticalReweighting( dataMngr.getDataSet(),      # source distribution
                                                        Bs2JpsiKKFit.getDataSet(),  # target distribution
                                                        OneDverticalRewNbins, 'mdau2', iterationNumber = iterNumb, # nBins, variable
                                                        equalStatsBins = EqualStatsBins
                                                        )
        dataMngr.appendWeights( mKKWeightsName, mKKweights )

    # match physics
    if reweightPhysics:
        physWeights = PhysicsReweight.calculateWeights( iterNumb, dataParameters )
        dataMngr.appendWeights( physWeightName, physWeights )
    
        assert False
    
    # match KK momenta
    if reweightKKmom and RewApproach == 'vertical':
        KKMomWeights = TwoDimentionalVerticalReweighting(dataMngr.getDataSet(),      # source distribution
                                                         Bs2JpsiKKFit.getDataSet(),  # target distribution
                                                         TwoDverticalRewNbins, ['Kplus_P','Kminus_P'], iterationNumber = iterNumb ,
                                                         # number of bins per dimention, variables
                                                         equalStatsBins = EqualStatsBins
                                                         )
        dataMngr.appendWeights( KmomentaWeightsName, KKMomWeights )
    elif reweightKKmom and RewApproach == 'horizontal':
        KinematicReweight.reweight( iterNumb, PhysicsReweight.getDataSet(weighted=True) )
        reweightedData = KinematicReweight.getDataSet()
   
    # comparition plots
    if makePlots and iterNumb in plotAtTheseSteps: # plot data after each reweighting step
        compPlots = compareDistributions( mcData          = dataMngr.getDataSet('initSource'),
                                          mcDataPhysRew   = dataMngr.getDataSet(physWeightName),
                                          MomRewData      = dataMngr.getDataSet(KmomentaWeightsName),
                                          BmomRewData     = dataMngr.getDataSet(BmomentumWeightsName) if reweightBmomentum else '', 
                                          mkkRewData      = dataMngr.getDataSet(mKKWeightsName) if reweightMkk else '',
                                          sData           = Bs2JpsiKKFit.getDataSet(),
                                          obsSet          = angles + time + muMomenta + Kmomenta + Bmomenta + KKMass,
                                          itNumb          = iterNumb
                                          )
        
        # plot weights
        plotWeightsList = [physWeightName, KmomentaWeightsName]
        if reweightBmomentum: plotWeightsList += [ BmomentumWeightsName ]
        if reweightMkk:       plotWeightsList += [ mKKWeightsName ]
        for wList in plotWeightsList: dataMngr.plotWeights(wList)

    # compute angular efficiency moments from the new reweighted mc dataset.    
    PhysicsReweight.setDataFitParameters(dataParameters) # set data pars to pdf (reweighted data has the data physics now)
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = PhysicsReweight.getPdf(), IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in PhysicsReweight.getAngleFunctions().functions.itervalues() for func in complexFunc if func )
                                      )
    
    scaleFactor = 1 / 16. / sqrt(pi)
    physMoments.initCovariances()
    physMoments.compute(dataMngr.getDataSet()) 
    physMoments.write( outputEffMomentsFileName + '_Phys_%s_Iteration'%iterNumb, Scale=scaleFactor )
    physMoments.convertEffWeightsToMoments( OutputFilePath    = outputEffMomentsFileName + '_weights_%s_Iteration'%iterNumb, 
                                            Scale             = scaleFactor,
                                            WeightNamesPrefix = PhysicsReweight.getParNamePrefix()
                                            )
    
    # perform sFit on data using the new angular acceptance and update the data physics parameters
    if doFit:
        Bs2JpsiKKFit.doFit( itNum=iterNumb, angAccFile= outputEffMomentsFileName + '_weights_%s_Iteration'%iterNumb )
        Bs2JpsiKKFit.updateDataParameters( dataParameters, itNum=iterNumb ) 
    
    # save memory
    del physMoments
    cleanP2VVPlotStash()
    dataMngr.clear()
    gc.collect() 

# observables plot with the corrected angular acceptance
if makePlots and plotFinalPdfonData:
    from P2VV.Utilities.Plotting import plot
    from ROOT import RooAbsData, TCanvas, RooArgSet
    c = TCanvas( 'sFit: CorrAngAcc', 'sFit: CorrAngAcc' )
    c.Divide(2,2)

    pdf  = Bs2JpsiKKFit.getPdf()
    data = Bs2JpsiKKFit.getDataSet()
    projectionArgSet = RooArgSet( list( pdf.ConditionalObservables() ) + [ pdf.indexCat() ] )

    for can, obs, Logy in zip( [ c.cd(i) for i in xrange(1,5) ],  
                               angles + time, 
                               3*[False] + [True]
                               ):
        plot( can, obs, data, pdf, plotResidHist=True, logy=Logy, 
              pdfOpts = dict( ProjWData = ( data.reduce(projectionArgSet), False ),
                              LineWidth = 1 
                              ),
              dataOpts = dict( MarkerSize = .5, 
                               DataError  = RooAbsData.SumW2,
                               XErrorSize = 0)
              )
    c.Print( 'anglesTimeCorrAngAcc.pdf' )
