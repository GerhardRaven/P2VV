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
timeEffPath       = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
timeEffUBName     = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
timeEffExclBName  = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
angEffMomentsFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

# varius flags / names
NumbOfIterations = 10
physWeightName   = 'weightPhys'
initialFitOnData = False
makePlots, plotAfterFitting, canvs = True, False, {}

# parameters for varius xChecks
devel        = False
convergeTest = False
untaggedFit  = False   
convergeTestDataFile = '/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting.root'
#convergeTestDataFile = [ '/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting_part%s.root'%n for n in [1,2] ]
#convergeTestDataFile = '/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting_part1.root'
convergeTestSourceDataName = 'source'
convergeTestTargetDataName = 'target'
convergeTestWeights        = 'trivialWeights'

# import data parameters obrtained from an initial fit to data
if not convergeTest: from P2VV.Utilities.MCReweighting import parValuesNoKKBinsWideKKWindow as dataParameters,\
                                                              parValuesMc2011Gen as monteCarloParamters
else: 
    from P2VV.Utilities.MCReweighting import parValuesMc2012Fit as dataParameters, parValuesMc2012Fit as monteCarloParamters
    mcTuplePath = convergeTestDataFile
    mcTupleName = convergeTestSourceDataName


###########################################################################################################################
## Begin iterative procedure  ##
################################
# import stuff  initialize objects.
from P2VV.Utilities.MCReweighting import MatchPhysics, MatchWeightedDistributions, \
                                        compareDistributions, BuildBs2JpsiKK2011sFit                                      
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.Utilities.Plotting import plot
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from P2VV import RooFitDecorators
from ROOT import RooArgSet, TFile, TCanvas
from math import pi, sqrt

# define a workspace
worksp = RooObject( workspace = 'iterativeProcedure' ).ws()

#  build data pdf and prepare the sFit ( This pdf is not multiplied by the angular acceptance ).
Bs2JpsiKK2011sFit = BuildBs2JpsiKK2011sFit( dataSetPath          = sDataPath    if not convergeTest else convergeTestDataFile ,
                                            dataSetName          = sDataName    if not convergeTest else convergeTestTargetDataName ,
                                            weightsName          = sWeightsName if not convergeTest else convergeTestWeights  ,
                                            timeEffHistFile      = timeEffPath, 
                                            timeEffHistUBName    = timeEffUBName,
                                            timeEffHistExclBName = timeEffExclBName, 
                                            KKmassBins           = None,  # '4KKMassBins',
                                            doUntaggedFit        = True if untaggedFit else False
                                            )

# initialise physics matching class and build MC pdf
PhysicsReweight = MatchPhysics( mcTuplePath, 
                                nTupleName          = mcTupleName,
                                mcParameters        = monteCarloParamters, 
                                multiplyWithTimeAcc = True ,
                                timeEffFiles        = timeEffPath,
                                timeEffHistNames    = dict( UB = timeEffUBName, ExclB = timeEffExclBName )
                                )
PhysicsReweight.buildMonteCarloPdf( dataPdfBuilder=Bs2JpsiKK2011sFit.getPdfBuilderObject() )

# get pdfs.
dataPdf = Bs2JpsiKK2011sFit.getPdf()
mcPdf   = PhysicsReweight.getPdf()

# get observables
angles     = Bs2JpsiKK2011sFit.getObservables('angles') 
time       = Bs2JpsiKK2011sFit.getObservables('time') # reconstructed decay time 
KKMassCat  = Bs2JpsiKK2011sFit.getPdfBuilderObject()['observables']['KKMassCat']
condObsSet = Bs2JpsiKK2011sFit.getPdf().ConditionalObservables().union( set([KKMassCat]) )

# get datasets.
sWeightedData = Bs2JpsiKK2011sFit.getDataSet()
projDataSet   = sWeightedData.reduce( RooArgSet(condObsSet) )
assert False


from ROOT import RooFit
fr = worksp['truetime'].frame(RooFit.Range(0,1))
mcPdf.plotOn(fr)
fr.Draw()

# initialise kinematic reweighting classs
KinematicReweight = MatchWeightedDistributions( inTree         = None,          # Source: Distribution to be altered 
                                                outTree        = sWeightedData, # Target: Distribution to be matched with
                                                whichVars      = ['Kminus_P'],  # Select which variables enter the transformation
                                                inWeightName   = physWeightName,
                                                outWeightName  = sWeightsName if not convergeTest else convergeTestWeights,
                                                observables    = PhysicsReweight.getMcObsSet(),   # true/reco decay time + angles
                                                spectatorVars  = PhysicsReweight.getNtupleVars(), # vraiables copied to the new dataset
                                                nBins          = 1000          # controls the preceision of the reweighting 
                                                ) 
assert False
# perform initial fit on data with the nominal angular acceptance
if initialFitOnData:
    Bs2JpsiKK2011sFit.doFit( angAccFile=angEffMomentsFile )
    if makePlots and plotAfterFitting:
            canvs['nomFit'] = TCanvas( 'initFit', 'intFit' )
            canvs['nomFit'].Divide(2,2)
            for can, obs, Logy in zip( [ canvs['nomFit'].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
                plot( can, obs,Bs2JpsiKK2011sFit.getDataSet(), Bs2JpsiKK2011sFit.getPdf(), plotResidHist=True, logy=Logy,  
                      pdfOpts=dict( ProjWData=projDataSet ) )
            canvs['nomFit'].Print('sFit_nom.pdf')
    # Update the data parameter values with the ones obtained from the fit 
    Bs2JpsiKK2011sFit.updateDataParameters( dataParameters )

# start looping.
for iterNumb in range( 1, NumbOfIterations + 1 ):
    print 'P2VV - INFO: Begin iteratitive procedure. Iteration step: ' +str(iterNumb) + '.'  

################################################################################################################################
## Match mc physics to sData and reweight Kaon momenta ##
#########################################################
    # calculate physics matcing weights.
    PhysicsReweight.calculateWeights( iterNumb, dataParameters )
    PhysicsReweight.writeWeights( weightsName=physWeightName) # mc dataset is internally updated with the weights.
    
    # reweight track momenta
    if devel:pass
    else:
        KinematicReweight.reweight( iterNumb, PhysicsReweight.getDataset() )
        reweightedData = KinematicReweight.getDataSet()
    
    if makePlots: # super impose data after each reweighting step
        
        if devel:
            mcData =  PhysicsReweight.getDataset() 
            reweightedData = TFile.Open('/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/temp_reweightedData_1stIter.root').Get('MomRewMC_1_Iter')
        

        compPlots = compareDistributions( mcData          = PhysicsReweight.getDataset(),
                                          momRewData      = reweightedData,
                                          sData           = sWeightedData, 
                                          itNumb          = iterNumb,
                                          physWeightsName = physWeightName,
                                          nullTest        = True if convergeTest else False 
                                          )
        for canv, name in zip( compPlots, 
                               ['anglesTime_%s.pdf'%iterNumb, 'KaonMomenta_%s.pdf'%iterNumb, 'muonMomenta_%s.pdf'%iterNumb, \
                                                         'assymKaonmomenta_%s.pdf'%iterNumb, 'assymMuonMomenta_%s.pdf'%iterNumb]
                               ): canv.Print(name)

##################################################################################################################################
## Compute angular efficiency moments for the new reweighted MC sample.##
#########################################################################
    nominalEffMoms = 'hel_UB_UT_trueTime_BkgCat050_KK30' # efficeincy moments output file 
    effMomentsFile = nominalEffMoms + '_Phys_%s_Iteration'%iterNumb 
    effWeightsFile = nominalEffMoms + '_weights_%s_Iteration'%iterNumb
    
    angleFuncs  = PhysicsReweight.getAngleFunctions()    # grab angular functions from mc pdf
    PhysicsReweight.setDataFitParameters(dataParameters) # set data pars to pdf (reweighted data has the data physics.)

    # build and write  effciency moments.
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = mcPdf, IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func ))
    
    scaleFactor = 1 / 16. / sqrt(pi)
    physMoments.initCovariances()
    physMoments.compute(reweightedData) 
    physMoments.write( effMomentsFile , Scale=scaleFactor )
    physMoments.convertEffWeightsToMoments( OutputFilePath = effWeightsFile , Scale = scaleFactor )
    
#####################################################################################################################################
## Perform sFit on data using the new angular acceptance.##
#########################################################################
    # multiply pdf with the new acceptance and do fit.
    Bs2JpsiKK2011sFit.doFit( itNum=iterNumb, angAccFile=effWeightsFile )
    
    if makePlots and plotAfterFitting: # plot
        canvs['%siter_sFit'%iterNumb] = TCanvas( 'sFit, %s AngAccCorr'%iterNumb, 'sFit, %s AngAccCorr'%iterNumb )
        canvs['%siter_sFit'%iterNumb].Divide(2,2)
        for can, obs, Logy in zip( [ canvs['%siter_sFit'%iterNumb].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
            plot( can, obs, Bs2JpsiKK2011sFit.getDataSet(), dataPdf_AngEff, plotResidHist=True, logy=Logy, 
                  pdfOpts=dict( ProjWData=projDataSet ) )
        canvs['%siter_sFit'%iterNumb].Print( 'sFit_%siter.pdf'%iterNumb )

    # update data physics parameters dictionary
    Bs2JpsiKK2011sFit.updateDataParameters( dataParameters, itNum=iterNumb ) 
   

    












# speed up development
# reweightedData = TFile.Open('/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/temp_reweightedData_1stIter.root').Get('MomRewMC_1_Iter') # Speed up
# KinematicReweight._recalculatedData = reweightedData

# Question
# C_SP factor in MC pdf ???

# Improvement ideas:
# Unify the two classes into 1.
# plot after mimicing

