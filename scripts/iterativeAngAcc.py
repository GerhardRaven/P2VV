########################################################################################################################
## Configuration ##
################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-R', '--KKmomRew', dest='KKmomRew',  default = 'vertical',       help='KK momentum reweighting approach (vertical/horizontal)')
parser.add_option('-c', '--MCProd',   dest='MCProd',    default = 'Sim08_reduced',  help='MC simulation conditions (Sim08_2011/Sim08_2012/Sim08)')
parser.add_option('-n', '--numIters', dest='numIters',  default = 7, type=int,      help='number of iterations')
parser.add_option('-f', '--InitFit',  dest='InitFit',   default = False,            help='Initial fit on data')
(options, args) = parser.parse_args()

# reweightng flow control
NumbOfIterations      = options.numIters
kinematicRewApproach  = options.KKmomRew
MCProd                = options.MCProd
initialFitOnData      = options.InitFit
combinedFit           = True
reweightBmomentum     = True
OneDverticalRewNbins  = 1000
TwoDverticalRewNbins  = 50
physWeightName        = 'phys'
mKKWeightsName        = 'mKK'
KmomentaWeightsName   = 'KKmom'
BmomentumWeightsName  = 'Bmom'

# plotig control
makePlots          = True
plotAtTheseSteps   = [ NumbOfIterations ]  # [ i for i in xrange(1,NumbOfIterations+1) ]
plotFinalPdfonData = False

# specify datasets, sim. conditions, nominal accceptance weights and data physics parameters
# source distribution
mcTuplePath = '/project/bfys/vsyropou/data/iterativeProcedure/Bs2JpsiPhi_20112012_Sim08_ntupleB_201309_add_afterFullSel_trackMom_BMom.root'
mcTupleName = 'DecayTree'

# target distribution
sDataPath    = '/project/bfys/vsyropou/data/iterativeProcedure/'
sDataName    = 'JpsiKK_sigSWeight'
sWeightsName = 'sWeights_ipatia'

# nominal angluar acceptance path 
nomAngEffMomentsFile     = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/uncorrecteEffMoments/'
outputEffMomentsFileName = 'hel_UB_UT_trueTime_BkgCat050_KK30' 

# source generating parameters
from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as monteCarloParameters

if  '2011' in MCProd:
    from P2VV.Utilities.MCReweighting import parValues6KKmassBins2011 as dataParameters
    sDataPath += 'P2VVDataSets2011Reco14_I2DiegoMass_6KKMassBins_2TagCats_trackMom_BMom.root'
    if 'reduced' in MCProd: nomAngEffMomentsFile += 'MC11_Sim08_reduced/Sim08_2011_hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
    else:                   nomAngEffMomentsFile += 'MC11_Sim08/Sim08_2011_hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
elif '2012' in MCProd:
    from P2VV.Utilities.MCReweighting import parValues6KKmassBins2012 as dataParameters
    sDataPath += 'P2VVDataSets2012Reco14_I2DiegoMass_6KKMassBins_2TagCats_trackMom_BMom.root'
    nomAngEffMomentsFile += 'MC12_Sim08/Sim08_2012_hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
else:
    sDataPath += 'P2VVDataSets20112012Reco14_I2DiegoMass_6KKMassBins_2TagCats_trackMom_BMom.root'
    from P2VV.Utilities.MCReweighting import parValues6KKmassBins20112012 as dataParameters
    nomAngEffMomentsFile += 'MC20112012_Sim08/Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
if combinedFit: sDataPath = '/project/bfys/vsyropou/data/iterativeProcedure/P2VVDataSets20112012Reco14_I2DiegoMass_6KKMassBins_2TagCats_trackMom_BMom.root'

###########################################################################################################################
## Begin iterative procedure  ##
################################
if kinematicRewApproach == 'vertical':  
    from P2VV.Utilities.MCReweighting import TwoDimentionalVerticalReweighting, OneDimentionalVerticalReweighting
elif kinematicRewApproach == 'horizontal': 
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
Bs2JpsiKKFit = BuildBs2JpsiKKFit( dataSetPath=sDataPath, dataSetName=sDataName, weightsName=sWeightsName) #, MonteCarloProduction=MCProd  )
if initialFitOnData:
    Bs2JpsiKKFit.doFit( angAccFile=nomAngEffMomentsFile )
    assert False

# initialise physics matching class and build MC pdf
PhysicsReweight = MatchPhysics( mcTuplePath, mcTupleName , MonteCarloProduction=MCProd, BlindPdf=Bs2JpsiKKFit.getBlindString() )
   
# keep track of the weights (avoid creating too many datasets)
dataMngr = WeightedDataSetsManager( source = PhysicsReweight.getDataSet() )

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
if kinematicRewApproach == 'horizontal':
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
    BmomentumWeights = OneDimentionalVerticalReweighting( dataMngr.getDataSet(),      # source distribution
                                                          Bs2JpsiKKFit.getDataSet(),  # target distribution
                                                          OneDverticalRewNbins, 'B_P', # nBins, variable
                                                          #xCheckPlots = True
                                                          )
    dataMngr.appendWeights( BmomentumWeightsName, BmomentumWeights, permanetnWeigts=True )

# start looping.
for iterNumb in range( 1, NumbOfIterations + 1 ):
    print 'P2VV - INFO: Iteratitive procedure, begining of iteration %s.'%str(iterNumb)
    dataMngr['iterationNumber'] = iterNumb
    dataMngr['saveIntermediateDatasets'] = True if iterNumb in plotAtTheseSteps else False
            
    # match physics
    physWeights = PhysicsReweight.calculateWeights( iterNumb, dataParameters )
    dataMngr.appendWeights( physWeightName, physWeights )

    # match mKK 
    mKKweights = OneDimentionalVerticalReweighting( dataMngr.getDataSet(),      # source distribution
                                                    Bs2JpsiKKFit.getDataSet(),  # target distribution
                                                    OneDverticalRewNbins, 'mdau2', iterationNumber = iterNumb, # nBins, variable
                                                    )
    dataMngr.appendWeights( mKKWeightsName, mKKweights )
    
    # match KK momenta
    if kinematicRewApproach == 'vertical':
        KKMomWeights = TwoDimentionalVerticalReweighting(dataMngr.getDataSet(),      # source distribution
                                                         Bs2JpsiKKFit.getDataSet(),  # target distribution
                                                         TwoDverticalRewNbins, ['Kplus_P','Kminus_P'], iterationNumber = iterNumb ,
                                                         # number of bins per dimention, variables
                                                         )
        dataMngr.appendWeights( KmomentaWeightsName, KKMomWeights )
    elif kinematicRewApproach == 'horizontal':
        KinematicReweight.reweight( iterNumb, PhysicsReweight.getDataSet(weighted=True) )
        reweightedData = KinematicReweight.getDataSet()
   
    # comparition plots
    if makePlots and iterNumb in plotAtTheseSteps: # plot data after each reweighting step
        compPlots = compareDistributions( mcData          = dataMngr.getDataSet('initSource'),
                                          mcDataPhysRew   = dataMngr.getDataSet(physWeightName),
                                          mkkRewData      = dataMngr.getDataSet(mKKWeightsName),
                                          MomRewData      = dataMngr.getDataSet(KmomentaWeightsName),
                                          BmomRewData     = dataMngr.getDataSet(BmomentumWeightsName) if reweightBmomentum else '', 
                                          sData           = Bs2JpsiKKFit.getDataSet(),
                                          obsSet          = angles + time + muMomenta + Kmomenta + Bmomenta + KKMass,
                                          itNumb          = iterNumb,
                                          )
        # plot weights
        plotWeightsList = [physWeightName, mKKWeightsName, KmomentaWeightsName]
        if reweightBmomentum: plotWeightsList += [ BmomentumWeightsName ]
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
    physMoments.write( outputEffMomentsFileName + '_Phys_%s_Iteration'%iterNumb , Scale=scaleFactor )
    physMoments.convertEffWeightsToMoments( OutputFilePath    = outputEffMomentsFileName + '_weights_%s_Iteration'%iterNumb, 
                                            Scale             = scaleFactor,
                                            WeightNamesPrefix = PhysicsReweight.getParNamePrefix()
                                            )

    # perform sFit on data using the new angular acceptance and update the data physics parameters
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
