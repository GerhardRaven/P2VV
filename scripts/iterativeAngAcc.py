

############################################################################################################
## Specify paths and paramters ##
################################

# mc input File    
mcTuplePath = '/project/bfys/vsyropou/data/P2VVDataSets2011_MC_forReweighting.root'

# sData input file
sDataPath    = '/project/bfys/vsyropou/data/P2VVDataSets2011Reco12_wideKKMass_noKKMassBins_2TagCats_forReweighting.root'
sDataName    = 'JpsiKK_sigSWeight' # 'DecayTree' # 'JpsiKK_sigSWeight', # JpsiKK
sWeightsName = 'N_sigMass_sw'      # 'sWeight'   # 'N_sigMass_sw'       # 'weightVar'
      
# time acceptance
timeEffPath       = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
timeEffUBName     = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
timeEffExclBName  = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
angEffMomentsFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

# parameters obtained from the initial sFit to data, no KKMass bins, wide KK mass window. 
dataParameters = dict(  AperpMag2        = 0.24663
                       ,AperpPhase       = 3.032
                       ,A0Mag2           = 0.52352
                       ,A0Phase          = 0 # cosntrained
                       ,AparPhase        = 3.2176
                       ,f_S_bin0         = 0.046222
                       ,ASOddPhase_bin0  = -0.066683
                       ,dM               = 17.677
                       ,dGamma           = 0.1023
                       ,Gamma            = 0.67283
                       ,phiCP            = 0.085636
                       ,lambdaCP         = 0.92737
                       )

# varius flags / names
NumbOfIterations = 5 # number of iterations.
physWeightName   = 'weightPhys'
makePlots        = True
plotAfterFitting = False
doBaselineFit    = False
canvs            = {} if makePlots else '' # keep reference to the plots

# parameters for varius xChecks
doNothingTest  = False
sDataXcheck    =  '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/xChecks/P2VVDataSets2011_MC_forReweighting_RenameTrueTime2Time_TrivialWeights.root' 
sWeightsNameXcheck = 'trivialWeights'
devel = True
if doNothingTest: 
    dataParameters = dict( AperpMag2        = 0.16
                           ,AperpPhase       = -0.17
                           ,A0Mag2           =  0.6
                           ,A0Phase          = 0 # cosntrained
                           ,AparPhase        = 2.50
                           ,f_S_bin0         = 0
                           ,ASOddPhase_bin0  = 0
                           ,dM               = 17.8
                           ,dGamma           = 0.06
                           ,Gamma            = 0.679
                           ,phiCP            = -0.04
                           ,lambdaCP         = 1.
                           )

###########################################################################################################################
## Begin iterative procedure ##
################################
# initialize objects.
from P2VV.Utilities.MCReweighting import MatchMCphysics2Data, MatchWeightedDistributions, \
                                        compareDistributions, BuildBs2JpsiKK2011sFit
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.Utilities.Plotting import plot
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from ROOT import RooArgSet, TFile, TCanvas
from math import pi, sqrt

# define a workspace
worksp = RooObject( workspace = 'iterativeAngularAcceptance' ).ws()

#  build data pdf and prepare the sFit ( This pdf is not multiplied by the angular acceptance ).
Bs2JpsiKK2011sFit = BuildBs2JpsiKK2011sFit( dataSetPath          = sDataPath,
                                            dataSetName          = sDataName,
                                            weightsName          = sWeightsName,
                                            timeEffHistFile      = timeEffPath, 
                                            timeEffHistUBName    = timeEffUBName,
                                            timeEffHistExclBName = timeEffExclBName, 
                                            KKmassBins           = None  # '4KKMassBins',
                                            )

# initialise physics matching class and build MC pdf
MatchPhysics = MatchMCphysics2Data( mcTuplePath )
MatchPhysics.buildMonteCarloPdf( dataPdfBuilder=Bs2JpsiKK2011sFit.getPdfBuilderObject() )

# get pdfs.
dataPdf = Bs2JpsiKK2011sFit.getPdf()
mcPdf   = MatchPhysics.getPdf()

# get observables
angles     = Bs2JpsiKK2011sFit.getObservables('angles') 
time       = Bs2JpsiKK2011sFit.getObservables('time') # reconstructed decay time 
KKMassCat  = Bs2JpsiKK2011sFit.getPdfBuilderObject()['observables']['KKMassCat']
condObsSet = Bs2JpsiKK2011sFit.getPdf().ConditionalObservables().union( set([KKMassCat]) )

# get datasets.
sWeightedData = Bs2JpsiKK2011sFit.getDataSet() if not doNothingTest else TFile.Open(sDataXcheck).Get('MC4xCheck')
projDataSet = sWeightedData.reduce( RooArgSet(condObsSet) ) if not doNothingTest else None

# initialise kinematic reweighting classs
matchMC2Data = MatchWeightedDistributions( inTree         = None,          # distribution to be matched with.
                                           outTree        = sWeightedData, # distribution to be matched with.
                                           whichVars      = ['Kminus_P'],
                                           inWeightName   = physWeightName,
                                           outWeightName  = sWeightsName if not doNothingTest else sWeightsNameXcheck, 
                                           observables    = MatchPhysics.getMcObsSet(),   # true/reco decay time + angles
                                           spectatorVars  = MatchPhysics.getNtupleVars(), # vraiables copied to the new dataset
                                           nBins          = 1000          # controls the preceision of the reweighting 
                                             ) 


# start looping.
for iterNumb in range( 1, NumbOfIterations + 1 ):
    print 'P2VV - INFO: Begin iteratitive procedure. Iteration step: ' +str(iterNumb) + '.'  

################################################################################################################################
## Match mc physics to sData and reweight Kaon momenta ##
#########################################################
    # calculate physics matcing weights.
    MatchPhysics.calculateWeights( iterNumb, dataParameters )
    MatchPhysics.writeWeights( weightsName=physWeightName) # mc dataset is internally updated with the weights.
    
    # reweight track momenta
    if devel:pass
    else:
        matchMC2Data.reweight( iterNumb, MatchPhysics.getDataset() )
        reweightedData = matchMC2Data.getDataSet()
    
    if makePlots: # super impose data after each reweighting step
        
        if devel:
            mcData =  MatchPhysics.getDataset() 
            reweightedData = TFile.Open('../temp_reweightedData_1stIter.root').Get('MomRewMC_1_Iter')
        

        compPlots = compareDistributions( mcData          = MatchPhysics.getDataset(),
                                          momRewData      = reweightedData,
                                          sData           = sWeightedData, 
                                          itNumb          = iterNumb,
                                          physWeightsName = physWeightName
                                          )
        for canv, name in zip( compPlots, ['anglesTime_%s.pdf'%iterNumb, 'KaonMomenta_%s.pdf'%iterNumb, 'muonMomenta_%s.pdf'%iterNumb, \
                                                                    'assymKaonmomenta_%s.pdf'%iterNumb, 'assymMuonMomenta_%s.pdf'%iterNumb]
                               ): canv.Print(name)
    if devel: assert False
##################################################################################################################################
## Compute angular efficiency moments for the new reweighted MC sample.##
#########################################################################
    nominalEffMoms = 'hel_UB_UT_trueTime_BkgCat050_KK30' # efficeincy moments output file 
    effMomentsFile = nominalEffMoms + '_Phys_%s_Iteration'%iterNumb 
    effWeightsFile = nominalEffMoms + '_weights_%s_Iteration'%iterNumb
    
    angleFuncs  = MatchPhysics.getAngleFunctions()    # grab angular functions from mc pdf
    MatchPhysics.setDataFitParameters(dataParameters) # set data pars to pdf (reweighted data have the data physics.)

    # build and write  effciency moments.
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = mcPdf, IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func ))
    
    scaleFactor = 1 / 16. / sqrt(pi)
    physMoments.initCovariances()
    physMoments.compute(reweightedData) 
    physMoments.write( effMomentsFile , Scale=scaleFactor )
    physMoments.convertPhysMomsToEffWeights( effWeightsFile , Scale=scaleFactor )
    
#####################################################################################################################################
## Perform sFit on data using the new angular acceptance.##
#########################################################################
    
    # perform the sFit with the nominal ang. acceptance
    if doBaselineFit:Bs2JpsiKK2011sFit.doFit( iterNumb, angEffMomentsFile )
    
    if makePlots and plotAfterFitting:
            canvs['nomFit'] = TCanvas( 'sFit, no AngAccCorr', 'sFit, no AngAccCorr' )
            canvs['nomFit'].Divide(2,2)
            for can, obs, Logy in zip( [ canvs['nomFit'].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
                plot( can, obs,Bs2JpsiKK2011sFit.getDataSet(), dataPdf_StandardAngEff, plotResidHist=True, logy=Logy,  
                      pdfOpts=dict( ProjWData=projDataSet ) )
            canvs['nomFit'].Print('sFit_nom.pdf')

    # multiply pdf with the new acceptance and do fit.
    Bs2JpsiKK2011sFit.doFit( iterNumb, effWeightsFile )

    if makePlots and plotAfterFitting: # plot
        canvs['%siter_sFit'%iterNumb] = TCanvas( 'sFit, %s AngAccCorr'%iterNumb, 'sFit, %s AngAccCorr'%iterNumb )
        canvs['%siter_sFit'%iterNumb].Divide(2,2)
        for can, obs, Logy in zip( [ canvs['%siter_sFit'%iterNumb].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
            plot( can, obs, Bs2JpsiKK2011sFit.getDataSet(), dataPdf_AngEff, plotResidHist=True, logy=Logy, 
                  pdfOpts=dict( ProjWData=projDataSet ) )
        canvs['%siter_sFit'%iterNumb].Print( 'sFit_%siter.pdf'%iterNumb )


    Bs2JpsiKK2011sFit.updateDataParameters( iterNumb, dataParameters ) # update data physics parameters dictionary














# speed up development
# reweightedData = TFile.Open('/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/temp_reweightedData_1stIter.root').Get('MomRewMC_1_Iter') # Speed up
# matchMC2Data._recalculatedData = reweightedData

# Question
# C_SP factor in MC pdf ???

# Improvement ideas:
# Unify the two classes into 1.
# plot after mimicing

