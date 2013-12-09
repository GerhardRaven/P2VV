########################################################################################################################
## Specify paths and paramters ##
################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-R', '--KKmomRew', dest='KKmomRew', default = 'vertical', type='str', help='KK momentum reweighting approach (vertical/horizontal)')
(options, args) = parser.parse_args()

MCProd = 'Sim08'

# varius flags / names
NumbOfIterations     = 7
physWeightName       = 'weightPhys'
kinematicRewApproach = options.KKmomRew

initialFitOnData = False
makePlots        = True
plotAfterFitting = False
canvs            = {}

# nominal angular acceptance file
angEffMomentsFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'

# set dataset paths and physics and physics parameters parameters for data and mc gen
if MCProd == 'Sim08':
    # source distribution
    mcTuplePath = '/project/bfys/vsyropou/data/iterativeProcedure/Bs2JpsiPhi_20112012_Sim08_ntupleB_201309_add_afterFullSel_trackMom_BMom.root'
    mcTupleName = 'DecayTree'

    # target distribution
    sDataPath    = '/project/bfys/vsyropou/data/iterativeProcedure/P2VVDataSets20112012Reco14_I2DiegoMass_6KKMassBins_2TagCats_trackMom_BMom.root'
    sDataName    = 'JpsiKK_sigSWeight'
    sWeightsName = 'sWeights_ipatia'
    
    from P2VV.Utilities.MCReweighting import parValues6KKmassBins20112012 as dataParameters
    from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as monteCarloParameters

elif MCProd == 'Sim06':
    mcTuplePath = [ '/project/bfys/vsyropou/data/iterativeProcedure/P2VVDataSetsMC11a_noKKMassBins_2TagCats_forReweighting_part%s.root'%n for n in [1,2] ]
    mcTupleName = 'JpsiKK'

    # sData input file
    sDataPath    = '/project/bfys/vsyropou/data/iterativeProcedure/P2VVDataSets2011Reco12_wideKKMass_noKKMassBins_2TagCats_forReweighting.root'
    sDataName    = 'JpsiKK_sigSWeight'
    sWeightsName = 'N_sigMass_sw'

    from P2VV.Utilities.MCReweighting import parValuesNoKKBinsWideKKWindow as dataParameters
    from P2VV.Utilities.MCReweighting import parValuesMc2011Gen as monteCarloParameters

###########################################################################################################################
## Begin iterative procedure  ##
################################
if kinematicRewApproach == 'vertical':  
    from P2VV.Utilities.MCReweighting import TwoDimentionalVerticalReweighting, OneDimentionalVerticalReweighting
elif kinematicRewApproach == 'horizontal': 
    from P2VV.Utilities.MCReweighting import MatchWeightedDistributions
from P2VV.Utilities.MCReweighting import MatchPhysics, compareDistributions, BuildBs2JpsiKKFit, cleanP2VVPlotStash, destroyRootObject
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.Utilities.Plotting import plot
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from P2VV import RooFitDecorators
from ROOT import RooArgSet, TFile, TCanvas
from math import pi, sqrt

# define a workspace
worksp = RooObject( workspace = 'iterativeProcedure' ).ws()

# build data pdf and prepare the sFit ( This pdf will not be multiplied by the angular acceptance ).
Bs2JpsiKKFit = BuildBs2JpsiKKFit( dataSetPath = sDataPath, dataSetName = sDataName, weightsName = sWeightsName )

# initialise physics matching class and build MC pdf
PhysicsReweight = MatchPhysics( mcTuplePath,  mcTupleName, MonteCarloProduction = MCProd )

# get observables
angles     = [worksp[o] for o in ['helcosthetaK','helcosthetaK','helphi']]
time       = [worksp['time']]
truetime   = [worksp['truetime']]
muMomenta  = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['muplus','muminus'] for comp in ('P','PX','PY','PZ') ]    ]
Kmomenta   = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['Kplus','Kminus']   for comp in ('P','PX','PY','PZ') ]    ]
Bmomenta   = [ worksp['B_P'], worksp['B_Pt'] ]
KKMass     = [worksp['mdau2']]
KKMassCat  = Bs2JpsiKKFit.getPdf().indexCat()

# initialise kinematic reweighting classs
if kinematicRewApproach == 'horizontal':
    KinematicReweight = MatchWeightedDistributions( outTree        = Bs2JpsiKKFit.getDataSet(), # Target: Distribution to be matched with
                                                    reweightVars   = ['Kminus_P'],              # Variables that enter the transformation
                                                    inWeightName   = physWeightName,
                                                    outWeightName  = sWeightsName,
                                                    observables    = angles + time + truetime,
                                                    nonObsVars     = muMomenta + Kmomenta + Bmomenta + KKMass,
                                                    nBins          = 1000                       # preceision of the transformation
                                                    ) 

# perform initial fit on data with the nominal angular acceptance
if initialFitOnData:
    Bs2JpsiKKFit.doFit( angAccFile=angEffMomentsFile )
    if makePlots and plotAfterFitting:
        condObsSet = Bs2JpsiKKFit.getPdf().ConditionalObservables().union( set([KKMassCat]) )
        projDataSet = Bs2JpsiKKFit.getDataSet().reduce( RooArgSet(condObsSet) )
        canvs['nomFit'] = TCanvas( 'initFit', 'intFit' )
        canvs['nomFit'].Divide(2,2)
        for can, obs, Logy in zip( [ canvs['nomFit'].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
            plot( can, obs,Bs2JpsiKKFit.getDataSet(), Bs2JpsiKKFit.getPdf(), plotResidHist=True, logy=Logy,  
                  pdfOpts=dict( ProjWData=projDataSet ) )
            canvs['nomFit'].Print('sFit_nom.pdf')
    # Update the data parameter values with the ones obtained from the fit 
    Bs2JpsiKKFit.updateDataParameters( dataParameters )

# start looping.
for iterNumb in range( 1, NumbOfIterations + 1 ):
    print 'P2VV - INFO: Iteratitive procedure, begining of iteration %s.'%str(iterNumb)

    # save memory
    cleanP2VVPlotStash()
    if iterNumb > 1: destroyRootObject(PhysicsReweight.getDataSet(weighted=True))

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

    if makePlots: # plot data after each reweighting step
        compPlots = compareDistributions( mcData          = PhysicsReweight.getDataSet(),
                                          mcDataPhysRew   = PhysicsReweight.getDataSet(weighted=True),
                                          MomRewData      = reweightedData,
                                          sData           = Bs2JpsiKKFit.getDataSet(),
                                          obsSet          = angles + muMomenta + Kmomenta + Bmomenta + KKMass,
                                          itNumb          = iterNumb,
                                          )
        # plot physics matching weights
        PhysicsReweight.plotWeights()
    assert False
    # compute angular efficiency moments for the new reweighted MC sample.##
    nominalEffMoms = 'hel_UB_UT_trueTime_BkgCat050_KK30'              # standard eff moments output file 
    effMomentsFile = nominalEffMoms + '_Phys_%s_Iteration'%iterNumb   # o
    effWeightsFile = nominalEffMoms + '_weights_%s_Iteration'%iterNumb
    
    # build and write effciency moments.
    PhysicsReweight.setDataFitParameters(dataParameters) # set data pars to pdf (reweighted data has the data physics.)
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = PhysicsReweight.getPdf(), IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in PhysicsReweight.getAngleFunctions().functions.itervalues() for func in complexFunc if func )
                                      )
    
    scaleFactor = 1 / 16. / sqrt(pi) # scale all efficiency weights 
    physMoments.initCovariances()
    physMoments.compute(reweightedData) 
    physMoments.write( effMomentsFile , Scale=scaleFactor )
    physMoments.convertEffWeightsToMoments( OutputFilePath    = effWeightsFile, 
                                            Scale             = scaleFactor,
                                            WeightNamesPrefix = PhysicsReweight.getParNamePrefix() # eff moments in mcPdf have a prefix 
                                            )

    # perform sFit on data using the new angular acceptance
    Bs2JpsiKKFit.doFit( itNum=iterNumb, angAccFile=effWeightsFile )
    
    if makePlots and plotAfterFitting: # plot
        canvs['%siter_sFit'%iterNumb] = TCanvas( 'sFit, %s AngAccCorr'%iterNumb, 'sFit, %s AngAccCorr'%iterNumb )
        canvs['%siter_sFit'%iterNumb].Divide(2,2)
        for can, obs, Logy in zip( [ canvs['%siter_sFit'%iterNumb].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
            plot( can, obs, Bs2JpsiKKFit.getDataSet(), Bs2JpsiKKFit.getPdf(), plotResidHist=True, logy=Logy, 
                  pdfOpts=dict( ProjWData=projDataSet ) )
        canvs['%siter_sFit'%iterNumb].Print( 'sFit_%siter.pdf'%iterNumb )

    # update data physics parameters dictionary
    Bs2JpsiKKFit.updateDataParameters( dataParameters, itNum=iterNumb ) 
    
    # collect garbage
    import gc
    gc.collect()
    




# Question
# C_SP factor in MC pdf ???

# Improvement ideas:
# Unify the two classes into 1.

# plot stuff
# PhysicsReweight.setMonteCarloParameters()
# c3 = TCanvas('el','skase')
# c3.Divide(2,2)
# angles = [worksp[o] for o in ['helcosthetaL','helcosthetaK','helphi']]
# for o, canv in zip(angles + [worksp['truetime']], [c3.cd(i) for i in [1,2,3,4]] ): plot(canv, o, PhysicsReweight.getDataSet(), PhysicsReweight.getPdf())

