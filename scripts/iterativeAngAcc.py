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

NumbOfIterations      = options.numIters
kinematicRewApproach  = options.KKmomRew
MCProd                = options.MCProd
initialFitOnData      = options.InitFit
physWeightName        = 'weightPhys'

# plotig control
makePlots        = True
plotAtTheseSteps = [ NumbOfIterations ]  # [ i for i in xrange(1,NumbOfIterations+1) ]

# specify datasets, sim. conditions, nominal accceptance weights and data physics parameters
from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as monteCarloParameters

# source distribution
mcTuplePath = '/project/bfys/vsyropou/data/iterativeProcedure/Bs2JpsiPhi_20112012_Sim08_ntupleB_201309_add_afterFullSel_trackMom_BMom.root'
mcTupleName = 'DecayTree'

# target distribution
sDataPath    = '/project/bfys/vsyropou/data/iterativeProcedure/'
sDataName    = 'JpsiKK_sigSWeight'
sWeightsName = 'sWeights_ipatia'

# nominal angluar acceptance path 
nomAngEffMomentsFile = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/uncorrecteEffMoments/'
outputEffMomentsFileName = 'hel_UB_UT_trueTime_BkgCat050_KK30' 

if   '2011' in MCProd:
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

###########################################################################################################################
## Begin iterative procedure  ##
################################
if kinematicRewApproach == 'vertical':  
    from P2VV.Utilities.MCReweighting import TwoDimentionalVerticalReweighting, OneDimentionalVerticalReweighting
elif kinematicRewApproach == 'horizontal': 
    from P2VV.Utilities.MCReweighting import MatchWeightedDistributions
from P2VV.Utilities.MCReweighting import MatchPhysics, compareDistributions, BuildBs2JpsiKKFit, cleanP2VVPlotStash #, destroyRootObject
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.Utilities.Plotting import plot
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from P2VV import RooFitDecorators
from ROOT import RooArgSet, TFile, TCanvas
from math import pi, sqrt
import gc

# define a workspace
worksp = RooObject( workspace = 'iterativeProcedure' ).ws()

# build data pdf and prepare the sFit ( This pdf will not be multiplied by the angular acceptance ).
Bs2JpsiKKFit = BuildBs2JpsiKKFit( dataSetPath = sDataPath, dataSetName = sDataName, weightsName = sWeightsName, MonteCarloProduction = MCProd  )
if initialFitOnData:
    Bs2JpsiKKFit.doFit( angAccFile=nomAngEffMomentsFile )
    assert False

# initialise physics matching class and build MC pdf
PhysicsReweight = MatchPhysics( mcTuplePath,  mcTupleName, MonteCarloProduction = MCProd, BlindPdf = Bs2JpsiKKFit.getBlindString() )

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

# start looping.
for iterNumb in range( 1, NumbOfIterations + 1 ):
    print 'P2VV - INFO: Iteratitive procedure, begining of iteration %s.'%str(iterNumb)

    # save memory
    cleanP2VVPlotStash()
    #if iterNumb > 1: destroyRootObject(PhysicsReweight.getDataSet(weighted=True))

    # match mc physics to sData
    PhysicsReweight.calculateWeights( iterNumb, dataParameters )
    PhysicsReweight.writeWeights( weightsName=physWeightName )
    
    # reweight track momenta
    if kinematicRewApproach == 'vertical':
        reweightedData = TwoDimentionalVerticalReweighting(PhysicsReweight.getDataSet(),      # source distribution
                                                           Bs2JpsiKKFit.getDataSet(),         # target distribution
                                                           50,                                # number of bins per dimention
                                                           ['Kplus_P','Kminus_P'],            # variables to reweight
                                                           'MomRew',                          # weights name
                                                           SourceWeightName = physWeightName, # weight name of the source if any
                                                           TargetWeightName = sWeightsName,
                                                           iterationNumber  = iterNumb        
                                                           )
    elif kinematicRewApproach == 'horizontal':
        KinematicReweight.reweight( iterNumb, PhysicsReweight.getDataSet(weighted=True) )
        reweightedData = KinematicReweight.getDataSet()

    if makePlots and iterNumb in plotAtTheseSteps: # plot data after each reweighting step
        compPlots = compareDistributions( mcData          = PhysicsReweight.getDataSet(),
                                          mcDataPhysRew   = PhysicsReweight.getDataSet(weighted=True),
                                          MomRewData      = reweightedData,
                                          sData           = Bs2JpsiKKFit.getDataSet(),
                                          obsSet          = angles + time + muMomenta + Kmomenta + Bmomenta + KKMass,
                                          itNumb          = iterNumb,
                                          )
        # plot physics matching weights
        PhysicsReweight.plotWeights()
    
    # compute angular efficiency moments for the new reweighted MC sample.    
    PhysicsReweight.setDataFitParameters(dataParameters) # set data pars to pdf (reweighted data has the data physics now)
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = PhysicsReweight.getPdf(), IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in PhysicsReweight.getAngleFunctions().functions.itervalues() for func in complexFunc if func )
                                      )
    
    scaleFactor = 1 / 16. / sqrt(pi) # scale all efficiency weights 
    physMoments.initCovariances()
    physMoments.compute(reweightedData) 
    physMoments.write( outputEffMomentsFileName + '_Phys_%s_Iteration'%iterNumb , Scale=scaleFactor )
    physMoments.convertEffWeightsToMoments( OutputFilePath    = outputEffMomentsFileName + '_weights_%s_Iteration'%iterNumb, 
                                            Scale             = scaleFactor,
                                            WeightNamesPrefix = PhysicsReweight.getParNamePrefix() # eff moments in mcPdf have a prefix 
                                            )

    # perform sFit on data using the new angular acceptance and update the data physics parameters
    Bs2JpsiKKFit.doFit( itNum=iterNumb, angAccFile= outputEffMomentsFileName + '_weights_%s_Iteration'%iterNumb )
    Bs2JpsiKKFit.updateDataParameters( dataParameters, itNum=iterNumb ) 
    
    # save memory
    del PhysicsReweight._weightedData
    del reweightedData
    del physMoments
    gc.collect() 

assert False
# observables plot with the corrected angular acceptance
if makePlots:
    from ROOT import RooAbsData 
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
