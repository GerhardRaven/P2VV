
########################################################################################################################
## Specify paths and paramters ##
################################

# mc input File    
mcTuplePath = [ '/project/bfys/vsyropou/data/P2VVDataSetsMC11a_noKKMassBins_2TagCats_forReweighting_part%s.root'%n for n in [1,2] ]
mcTupleName = 'JpsiKK'

# sData input file
sDataPath    = '/project/bfys/vsyropou/data/P2VVDataSets2011Reco12_wideKKMass_noKKMassBins_2TagCats_forReweighting.root'
sDataName    = 'JpsiKK_sigSWeight' # 'DecayTree' # 'JpsiKK_sigSWeight', # JpsiKK
sWeightsName = 'N_sigMass_sw'      # 'sWeight'   # 'N_sigMass_sw'       # 'weightVar'

# time acceptance
from P2VV.Parameterizations.FullPDFs import SimulCatSettings
timeEffType = 'paper2012' 
timeEffHistFiles = SimulCatSettings('timeEffHistFiles')
timeEffPath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/'
timeEffHistFiles.addSettings( '', '', dict(  file      = timeEffPath + 'Bs_HltPropertimeAcceptance_Data-20120816.root'
                                           , hlt1UB    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                           , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached')
                              )

# nominal angular acceptance file
angEffMomentsFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
anglesEffType = 'weights' 

# varius flags / names
NumbOfIterations     = 7
physWeightName       = 'weightPhys'
kinematicRewApproach = 'horizontal' # 'horizontal'

initialFitOnData = False
makePlots        = True
plotAfterFitting = False
canvs            = {}

# parameters / flags for varius xChecks
devel        = False
convergeTest = False
untaggedFit  = False  
convergeTestDataFile = '/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting.root'
convergeTestSourceDataName = 'source'
convergeTestTargetDataName = 'target'
convergeTestWeights        = 'trivialWeights'

# import data parameters obtained from an initial fit to data
if not convergeTest: 
    from P2VV.Utilities.MCReweighting import parValuesNoKKBinsWideKKWindow as dataParameters,\
                                             parValuesMc2011Gen as monteCarloParameters
else: 
    from P2VV.Utilities.MCReweighting import parValuesMc2012Fit as dataParameters, \
                                             parValuesMc2012Fit as monteCarloParameters
    mcTuplePath  = convergeTestDataFile
    mcTupleName  = convergeTestSourceDataName
    sDataPath    = convergeTestDataFile
    sDataName    = convergeTestTargetDataName
    sWeightsName = convergeTestWeights 

###########################################################################################################################
## Begin iterative procedure  ##
################################
# import stuff  initialize objects.
from P2VV.Utilities.MCReweighting import MatchPhysics, MatchWeightedDistributions, compareDistributions,\
                                         BuildBs2JpsiKKFit, TwoDimentionalVerticalReweighting, cleanP2VVPlotStash
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.Utilities.Plotting import plot
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from P2VV import RooFitDecorators
from ROOT import RooArgSet, TFile, TCanvas
from math import pi, sqrt

# define a workspace
worksp = RooObject( workspace = 'iterativeProcedure' ).ws()

# build data pdf and prepare the sFit ( This pdf is not multiplied by the angular acceptance ).
Bs2JpsiKKFit = BuildBs2JpsiKKFit( dataSetPath     = sDataPath,
                                  dataSetName     = sDataName,
                                  weightsName     = sWeightsName,
                                  timeEffHistFile = timeEffHistFiles, 
                                  KKmassBins      = None,  
                                  doUntaggedFit   = True if untaggedFit else False,
                                  doNullTest      = True if convergeTest else False
                                  )

# initialise physics matching class and build MC pdf
PhysicsReweight = MatchPhysics( mcTuplePath, 
                                nTupleName          = mcTupleName,
                                mcParameters        = monteCarloParameters, 
                                timeEffType         = timeEffType,
                                timeEffHistFiles    = timeEffHistFiles,
                                anglesEffType       = anglesEffType,
                                angEffMomsFiles     = angEffMomentsFile,
                                monteCarloParams    = monteCarloParameters
                                )

# get observables
angles     = Bs2JpsiKKFit.getObservables('angles') 
time       = Bs2JpsiKKFit.getObservables('time') # reconstructed decay time
mcTime     = PhysicsReweight.getMcTime()         # true / reco decay time
muMomenta  = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['muplus','muminus'] for comp in ('P','PX','PY','PZ') ]    ]
Kmomenta   = [worksp[o] for o in [ '%s_%s' % ( part, comp ) for part in ['Kplus','Kminus']   for comp in ('P','PX','PY','PZ') ]    ]
Bmomenta   = [ worksp['B_P'], worksp['B_Pt'] ]
KKMass     = [worksp['mdau2']]
KKMassCat  = Bs2JpsiKKFit.getPdf().indexCat()
condObsSet = Bs2JpsiKKFit.getPdf().ConditionalObservables().union( set([KKMassCat]) )

# get projection dataset for ploting.
projDataSet = Bs2JpsiKKFit.getDataSet().reduce( RooArgSet(condObsSet) )

# initialise kinematic reweighting classs
KinematicReweight = MatchWeightedDistributions( inTree         = None,                      # Source: Distribution to be altered 
                                                outTree        = Bs2JpsiKKFit.getDataSet(), # Target: Distribution to be matched with
                                                reweightVars   = ['Kminus_P'],              # Variables that enter the transformation
                                                inWeightName   = physWeightName,
                                                outWeightName  = sWeightsName if not convergeTest else convergeTestWeights,
                                                observables    = angles + time + mcTime,
                                                nonObsVars     = muMomenta + Kmomenta + Bmomenta + KKMass,  
                                                nBins          = 1000                       # preceision of the transformation
                                                ) 

# perform initial fit on data with the nominal angular acceptance
if initialFitOnData:
    Bs2JpsiKKFit.doFit( angAccFile=angEffMomentsFile )
    if makePlots and plotAfterFitting:
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

    # match mc physics to sData
    PhysicsReweight.calculateWeights( iterNumb, dataParameters )
    PhysicsReweight.writeWeights( weightsName=physWeightName)
    
    # clean plot stash to save memory
    cleanP2VVPlotStash()

    # reweight track momenta
    if kinematicRewApproach == 'vertical':
        reweightedData = TwoDimentionalVerticalReweighting(PhysicsReweight.getDataSet(),      # source distribution
                                                           Bs2JpsiKKFit.getDataSet(),         # target distribution
                                                           50,                                # number of bins per dimention
                                                           ['Kplus_P','Kminus_P'],            # variables to reweight
                                                           'MomRew',                          # weights name
                                                           SourceWeightName = physWeightName, # weight name of the source if any
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
                                          obsSet          = angles + time + muMomenta + Kmomenta + Bmomenta + KKMass,
                                          itNumb          = iterNumb,
                                         ## physWeightsName = physWeightName,
                                          nullTest        = True if convergeTest else False 
                                          )
        # plot physics matching weights
        PhysicsReweight.plotWeights()
    assert False
    # compute angular efficiency moments for the new reweighted MC sample.##
    nominalEffMoms = 'hel_UB_UT_trueTime_BkgCat050_KK30' # efficeincy moments output file 
    effMomentsFile = nominalEffMoms + '_Phys_%s_Iteration'%iterNumb 
    effWeightsFile = nominalEffMoms + '_weights_%s_Iteration'%iterNumb
    
    angleFuncs  = PhysicsReweight.getAngleFunctions()    # grab angular functions from mc pdf
    PhysicsReweight.setDataFitParameters(dataParameters) # set data pars to pdf (reweighted data has the data physics.)

    # build and write effciency moments.
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = PhysicsReweight.getPdf(), IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func )
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

