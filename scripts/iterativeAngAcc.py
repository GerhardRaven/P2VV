########################################################################################################################
## Configuration ##
################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-r', '--KKmomRew',  dest='KKmomRew',   default = 'vertical',            help='KK momentum reweighting approach (vertical/horizontal)')
parser.add_option('-o', '--rewSteps',  dest='rewSteps',   default = 'BmommkkphysKKmom',    help='reweghting steps order')
parser.add_option('-b', '--Bmom2DRew', dest='Bmom2DRew',  default = 'True',                help='2 dimentional Bmom reweighting switch')
parser.add_option('-s', '--MCProd',    dest='MCProd',     default = '2011',                help='choose mc sample ( 2011,2012 )')
parser.add_option('-n', '--iterNum',   dest='iterNum',    default = 1, type=int,           help='iteration number')
parser.add_option('-a', '--nomAngAcc', dest='nomAngAcc',  default = '',                    help='nominal angular acceptance')
parser.add_option('-d', '--physPars',  dest='physPars',   default = '',                    help='physics aprameters')
parser.add_option('-f', '--fit',       dest='fit',        default = 'False',               help='switch on/off fitting')
parser.add_option('-w', '--writeData', dest='writeData',  default = 'False',               help='save mc datasets to file')
parser.add_option('-p', '--makePlots', dest='makePlots',  default = 'False',               help='switch on/off plotting')
parser.add_option('-c', '--combMoms',  dest='combMoms',   default = 'False',               help='combine 2011,2012 moments')
parser.add_option('-R', '--reduced',   dest='reduced',    default = 'False',                help='apply a mass cut for a reduced sample')
(options, args) = parser.parse_args()

# reweightng flow control
iterNumb              = options.iterNum
RewApproach           = options.KKmomRew
MCProd                = options.MCProd
reweightBmomentum     = True if 'Bmom'  in options.rewSteps else False
reweightMkk           = True if 'mkk'   in options.rewSteps else False
reweightPhysics       = True if 'phys'  in options.rewSteps else False
reweightKKmom         = True if 'KKmom' in options.rewSteps else False
twoDimensionalBmomRew = True if 'True'  in options.Bmom2DRew else False
EqualStatsBins        = True
OneDverticalRewNbins  = 1000
TwoDverticalRewNbins  = 50
physWeightName        = 'phys'
mKKWeightsName        = 'mKK'
KmomentaWeightsName   = 'KKmom'
BmomentumWeightsName  = 'Bmom'
writeWeightedData     = True if 'True' in options.writeData else False
combineEffMoments     = True if 'True' in options.combMoms else False
reduced               = True if 'True' in options.reduced else False

# fit configuration
doFit            = True if 'True' in options.fit else False
nCPU             = 8
initialFitOnData = False

# plotig configuration
makePlots = True if 'True' in options.makePlots else False

# source distribution
dataSetsPath     = '/project/bfys/jleerdam/data/Bs2Jpsiphi/angEff/'
mcData11FileName = 'Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_angEff.root'
mcData12FileName = 'Bs2JpsiPhi_MC2012_ntupleB_20130904_angEff.root'
mcTupleName      = 'DecayTree'
if MCProd == '2011': monteCarloData = dataSetsPath + mcData11FileName
if MCProd == '2012': monteCarloData = dataSetsPath + mcData12FileName

# target distribution
dataPath     = dataSetsPath + 'P2VVDataSets%sReco14_I2Mass_6KKMassBins_2TagCats_kinematics_HLT2B.root'%MCProd
fitDataPath  = dataSetsPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root' if doFit else None
sDataName    = 'JpsiKK_sigSWeight'
sWeightsName = 'sWeights_ipatia'

# nominal angluar acceptance 
nomAngEffMomentsFile = options.nomAngAcc if options.nomAngAcc \
    else '/project/bfys/vsyropou/data/uncorrecteEffMoments/MC20112012_Sim08/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_weights'
outputEffMomentsBaselineName = 'hel_UB_UT_trueTime_BkgCat050_KK30'

# source generating physics  parameters
## TODO:: set mc pars from a file just like setdatapars
from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as monteCarloParameters
# /project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/nominalFitResults/monteCarloSim08parameterValues.par

# target physics parameters
dataParameters = options.physPars if options.physPars\
    else '/project/bfys/vsyropou/data/nominalFitResults/20112012Reco14DataFitValues_6KKMassBins.par'

###########################################################################################################################
## Begin iterative procedure  ##
################################
if RewApproach == 'horizontal': from P2VV.Utilities.MCReweighting import MatchWeightedDistributions
if doFit: from P2VV.Utilities.MCReweighting import  BuildBs2JpsiKKFit
from P2VV.Utilities.MCReweighting import TwoDimentionalVerticalReweighting, OneDimentionalVerticalReweighting
from P2VV.Utilities.MCReweighting import MatchPhysics, combineMoments
from P2VV.Utilities.MCReweighting import compareDistributions, cleanP2VVPlotStash, WeightedDataSetsManager
from P2VV.Utilities.DataMoments import RealMomentsBuilder, normalizeMoments, convertEffWeightsToMoments, readMoments
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from P2VV import RooFitDecorators
from ROOT import TFile 
from math import pi, sqrt
import os, gc

# define a workspace
worksp = RooObject( workspace = 'iterativeProcedure' ).ws()

# build data pdf and prepare the sFit ( This pdf will not be multiplied by the angular acceptance !! ).
Bs2JpsiKKFit = BuildBs2JpsiKKFit( dataSetPath=fitDataPath, dataSetName=sDataName, Ncpu=nCPU ) if doFit else None

if initialFitOnData and doFit:
    Bs2JpsiKKFit.doFit( angAccFile=nomAngEffMomentsFile )
    assert False

# initialise physics matching class and build MC pdf
print 'P2VV - INFO: Iteration Number %s. Running reweighting procedure in sample %s'%(iterNumb, mcData11FileName if MCProd=='2011' else mcData12FileName)
PhysicsReweight = MatchPhysics( monteCarloData, mcTupleName , MonteCarloProduction = MCProd, Reduced = reduced )

# manage the weights (avoid creating too many datasets)
mcDataMngr = WeightedDataSetsManager( source = PhysicsReweight.getDataSet() )
mcDataMngr['saveIntermediateDatasets'] = True if makePlots else False

# get observables
angles     = [worksp[o] for o in ['helcosthetaK','helcosthetaL','helphi']]
time       = [worksp['time']]
truetime   = [worksp['truetime']]
muMomenta  = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['muplus','muminus'] for comp in ('P','PX','PY','PZ') ]    ]
Kmomenta   = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['Kplus','Kminus']   for comp in ('P','PX','PY','PZ') ]    ]
Bmomenta   = [ worksp['B_P'], worksp['B_Pt'] ]
KKMass     = [worksp['mdau2']]
KKMassCat  = PhysicsReweight.getPdf().indexCat()

# begin reweighting procedure
print 'P2VV - INFO: Start reweighting mc data_%s.'%MCProd
print 'P2VV - INFO:\nSource distribution file: %s. \nTarget distribution file: %s.'%(monteCarloData,dataPath)
source = lambda: mcDataMngr.getDataSet()
target = TFile.Open(dataPath).Get(sDataName)
mcDataMngr['iterationNumber'] = iterNumb

# match B momentum
if reweightBmomentum:
    BmomentumWeights = TwoDimentionalVerticalReweighting( source(), target, TwoDverticalRewNbins, ['B_P','B_Pt'],
                                                          equalStatsBins = EqualStatsBins
                                                          ) if twoDimensionalBmomRew else \
                       OneDimentionalVerticalReweighting( source(), target, TwoDverticalRewNbins, 'B_P',
                                                          equalStatsBins = EqualStatsBins
                                                          )
    mcDataMngr.appendWeights( BmomentumWeightsName, BmomentumWeights, permanetnWeigts=True )

# match mKK
if reweightMkk:
    mKKweights = OneDimentionalVerticalReweighting( source(), target, OneDverticalRewNbins, 'mdau2', 
                                                    equalStatsBins = EqualStatsBins
                                                    )
    mcDataMngr.appendWeights( mKKWeightsName, mKKweights, permanetnWeigts=True )

# match physics
if reweightPhysics:
    PhysicsReweight.setDataSet( source() )
    physWeights = PhysicsReweight.calculateWeights( iterNumb, dataParameters )
    mcDataMngr.appendWeights( physWeightName, physWeights )
          
# match KK momenta
if reweightKKmom and RewApproach == 'vertical':
    KKMomWeights = TwoDimentionalVerticalReweighting(source(), target, TwoDverticalRewNbins, ['Kplus_P','Kminus_P'], 
                                                     iterationNumber = iterNumb, equalStatsBins = EqualStatsBins
                                                     )
    mcDataMngr.appendWeights( KmomentaWeightsName, KKMomWeights )
elif reweightKKmom and RewApproach == 'horizontal':
    KinematicReweight = MatchWeightedDistributions( outTree        = target, # Target: Distribution to be matched with
                                                    reweightVars   = ['Kminus_P'],  # Variables that enter the transformation
                                                    inWeightName   = mcDataMngr.getWeightName(),
                                                    outWeightName  = sWeightsName,
                                                    observables    = angles + time + truetime,
                                                    nonObsVars     = muMomenta + Kmomenta + Bmomenta + KKMass,
                                                    nBins          = 1000           # preceision of the transformation
                                                    )
    
    KinematicReweight.reweight( iterNumb, source() )
    KmomentaWeightsName = 'hor' + KmomentaWeightsName
    mcDataMngr.setDataSet( KinematicReweight.getDataSet(),  KmomentaWeightsName )
      
# compute angular efficiency moments from the new reweighted mc dataset.
if reweightPhysics: # set data pars to pdf (reweighted data has the data physics now)
    PhysicsReweight.setDataFitParameters(dataParameters) 
else:
    PhysicsReweight.setMonteCarloParameters()
physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), 
                                                             BasisFunc = func,
                                                             Norm = 1., 
                                                             PDF = PhysicsReweight.getPdf(), 
                                                             IntSet = [ ], 
                                                             NormSet = angles 
                                                             )\
       for complexFunc in PhysicsReweight.getAngleFunctions().functions.itervalues() for func in complexFunc if func )
                                  )
  
scaleFactor = 1 / 16. / sqrt(pi)
physMoments.initCovariances()
physMoments.compute(mcDataMngr.getDataSet()) 
physMoments.write( 'Sim08_{0}_{1}_Phys_{2}'.format(MCProd,outputEffMomentsBaselineName,iterNumb), Scale=scaleFactor )
  
# normalize effciency moments
normalizeMoments( 'Sim08_{0}_{1}_Phys_{2}'.format(MCProd,outputEffMomentsBaselineName,iterNumb),
                  'Sim08_{0}_{1}_Phys_norm_{2}'.format(MCProd,outputEffMomentsBaselineName,iterNumb),
                  normMoment = 'mc_Re_ang_A0_A0',
                  printMoms  = False
                  )
os.remove('Sim08_{0}_{1}_Phys_{2}'.format(MCProd,outputEffMomentsBaselineName,iterNumb) )

# combine 2011,2012 acceptances
if combineEffMoments:
    angAcc2011  = 'Sim08_2011_%s_Phys_norm_%s'%(outputEffMomentsBaselineName,iterNumb) 
    angAcc2012  = 'Sim08_2012_%s_Phys_norm_%s'%(outputEffMomentsBaselineName,iterNumb) 
    combAccName = 'Sim08_20112012_{0}_Phys_norm_{1}'.format(outputEffMomentsBaselineName,iterNumb)
    combineMoments( angAcc2011, angAcc2012, combAccName, Prefix = 'mc' )
    
# convert effyciency weights to efficiency moments
    moments, correlations = {}, {}
    readMoments( combAccName, BasisFuncNames = [], Moments = moments, Correlations = correlations, ProcessAll = True )
    convertEffWeightsToMoments( moments, OutputFilePath = combAccName.replace('Phys','weights').replace('_norm',''),
                                WeightNamesPrefix = PhysicsReweight.getParNamePrefix(),
                                PrintMoments      = False
                                )

# perform sFit on data using the new angular acceptance and update the data physics parameters
if doFit:
    ParFileOut = '20112012Reco14DataFitValues_6KKMassBins_angEff_%s.par'%iterNumb 
    Bs2JpsiKKFit.doFit( angAccFile = combAccName.replace('Phys','weights'), 
                        parFileOut = ParFileOut, 
                        itNum      = iterNumb 
                        )
    dataParameters = ParFileOut
  
# plot data after each reweighting step
if makePlots: 
    from P2VV.Utilities.MCReweighting import plotingScenarios, weightNamesDataKeysMap
    plot = True
    try: plotingScenarios[options.rewSteps]
    except KeyError: plot = False         
    if plot:
        for keys in plotingScenarios[options.rewSteps]:
            dataSets = { keys[0]: mcDataMngr.getDataSet(weightNamesDataKeysMap[ keys[0] ]),
                         keys[1]: mcDataMngr.getDataSet(weightNamesDataKeysMap[ keys[1] ])
                         }
            compPlots = compareDistributions( sData    = target,
                                              obsSet   = angles + time + muMomenta + Kmomenta + Bmomenta + KKMass,
                                              itNumb   = iterNumb,
                                              prodData = MCProd,
                                              **dataSets
                                              )

        # plot weights
        plotWeightsList = []
        if (reweightKKmom,RewApproach)==(True,'vertical'):  plotWeightsList += [ KmomentaWeightsName ]    
        if reweightPhysics:                                 plotWeightsList += [ physWeightName ]
        if reweightMkk:                                     plotWeightsList += [ mKKWeightsName ]
        if reweightBmomentum:                               plotWeightsList += [ BmomentumWeightsName ]
        for wList in plotWeightsList: mcDataMngr.plotWeights(wList)
    else: print 'P2VV - WARNING: There is no point in ploting with this reweighting configuration, skipping plotting.'

# write weighted mc data to a file
if writeWeightedData: # TODO: got errors fix this
    weightedMcFileName = monteCarloData.partition(dataSetsPath)[2].partition('.root')[0] + '_' + str(iterNumb) + '.root'
    wMcFile = TFile.Open( weightedMcFileName, 'recreate')
    mcDataMngr.getDataSet().Write()
    wMcFile.Close()
    del wMcFile
    print 'P2VV -INFO: Wrote reweighted MC dataset to file %s'%weightedMcFileName

# save memory
del physMoments
cleanP2VVPlotStash()
mcDataMngr.clear()
gc.collect()
