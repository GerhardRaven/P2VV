
############################################################################################################
## Specify paths and paramters ##
################################
globalOutputFolder = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/'

# mc input File    
mcTuplePath = '/project/bfys/vsyropou/data/P2VVDataSets2011_MC_forReweighting.root'
# sData input file
sDataPath    = '/project/bfys/vsyropou/data/P2VVDataSets2011Reco12_wideKKMass_noKKMassBins_2TagCats_forReweighting.root'
sDataName    = 'JpsiKK'        # 'DecayTree' # 'JpsiKK_sigSWeight', # JpsiKK
sWeightsName = 'N_sigMass_sw'  # 'sWeight'   # 'N_sigMass_sw'       # 'weightVar'
      
# time acceptance
timeEffPath       = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
timeEffUBName     = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
timeEffExclBName  = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
angEffMomentsFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

# parameters obtained from the initial sFit to data, no KK mass binning, wide KK mass window. 
dataParameters = dict(  AperpMag2   = 0.246729
                       ,AperpPhase  = 3.0270
                       ,A0Mag2      =  0.52348
                       ,A0Phase     = 0 # cosntrained
                       ,AparPhase   = 3.2133
                       ,f_S         = 0.046284
                       ,ASOddPhase  = -0.06552
                       ,dM          = 17.676
                       ,dGamma      = 0.10220
                       ,Gamma       = 0.67283
                       ,phiCP       = 0.08455
                       ,lambdaCP    = 0.92747
                       )

# varius flags / names
NumbOfIterations = 2 # desired number of iterations.
makePlots        = True
physWeightName   = 'weightPhys'
doBaselineFit    = False




############################################################################################################
## Begin iterative procedure ##
################################
# initialize objects.
from P2VV.Utilities.MCReweighting import MatchMCphysics2Data, MatchWeightedDistributions, \
                                         CompareWeightedDistributions, BuildBs2JpsiKK2011sFit
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.Utilities.Plotting import plot
from P2VV.RooFitWrappers import RooObject, RealEffMoment
#from P2VV.Load import LHCbStyle
from ROOT import RooArgSet, TFile, TCanvas
from math import pi, sqrt

worksp = RooObject( workspace = 'iterativeAngularAcceptance' ).ws()

#  build data pdf for sFiting ( This pdf is not multiplied by the angular acceptance ).
Bs2JpsiKK2011sFit = BuildBs2JpsiKK2011sFit( dataSetPath          = sDataPath,
                                            dataSetName          = sDataName,
                                            weightsName          = sWeightsName,
                                            timeEffHistFile      = timeEffPath, 
                                            timeEffHistUBName    = timeEffUBName,
                                            timeEffHistExclBName = timeEffExclBName, 
                                            KKmassBins           = None  # '4KKMassBins',
                                            )

# build MC pdf.
MatchPhysics = MatchMCphysics2Data( mcTuplePath )
MatchPhysics.buildMonteCarloPdf( dataPdfBuilder=Bs2JpsiKK2011sFit.getPdfBuilderObject() )

# get pdfs.
dataPdf = Bs2JpsiKK2011sFit.getPdf()
mcPdf   = MatchPhysics.getPdf()

# get datasets.
mcData        = MatchPhysics.getInitialMCafterSel() 
sWeightedData = Bs2JpsiKK2011sFit.getDataSet()

# get observables
angles     = Bs2JpsiKK2011sFit.getObservables('angles') 
time       = Bs2JpsiKK2011sFit.getObservables('time') # reconstructed decay time 
KKMassCat  = Bs2JpsiKK2011sFit.getPdfBuilderObject()['observables']['KKMassCat']
mcObsSet   = MatchPhysics.getMcObsSet()               # true/reco decay time + angles 
condObsSet = Bs2JpsiKK2011sFit.getPdf().ConditionalObservables().union( set([KKMassCat]) )

# projection dataset for the simultaneous pdf.
projDataSet = sWeightedData.reduce( RooArgSet(condObsSet) )

# prepare canvases container.
canvs = {}



# THIS IS TEMP: Speed up development.
#reweightedData = TFile.Open('/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/temp_reweightedData_1stIter.root').Get('MomRewMC_1_Iter') # Speed up



# start looping.
for iterNumb in range( 1, NumbOfIterations + 1 ):
    print 'P2VV - INFO: Iteration number ' +str(iterNumb) + '.'  
############################################################################################################
## Match Mc physics to sData. ##
################################
    # calculate and write physics matcing weights
    MatchPhysics.calculateWeights( iterNumb, dataParameters, mcData )
    MatchPhysics.combineWeights() # weights from all iterations are combiend to one weights set.
    MatchPhysics.writeWeights( weightsName=physWeightName) # mcData is internally updated with the weights.
    
    if makePlots: # check the effect of physics reweighting. 
        MatchPhysics.compareAngles()        # save by default to anglesComparision_?.pdf'
        MatchPhysics.compareTrackMomenta()  # save by default to KaonComparision_?.pdf' and muonComparision_?.pdf
   
 ############################################################################################################
 ## Reweight Kaon momenta of the previously reweighted MC to match the Kaon momenta of sData.##
 ##############################################################################################
    reweightArgs = dict( inTree         = mcData,        # distribution to be modified. 
                         outTree        = sWeightedData, # distribution to be matched with.
                         whichVars      = ['Kminus_P'],
                         inWeightName   = physWeightName,
                         outWeightName  = sWeightsName,
                         observables    = mcObsSet,
                         spectatorVars  = MatchPhysics.getNtupleVars(), # vraiables copied to the new dataset.
                         nBins          = 1000,          # controls the preceision of the reweighting 
                         itNum          = iterNumb
                         )
    matchMC2Data = MatchWeightedDistributions( **reweightArgs )
    matchMC2Data.reweight()
    reweightedData = matchMC2Data.getDataSet()
    
    if makePlots:
        #matchMC2Data._recalculatedData = reweightedData
        matchMC2Data.compareAngles()
        matchMC2Data.compareTrackMomenta()
    #assert False  
##############################################################################################################
## Compute angular efficiency moments for the new reweighted MC sample.##
#########################################################################
    momentsFile   = globalOutputFolder + 'hel_UB_UT_trueTime_BkgCat050_KK30' # efficeincy moments output file. 
    mc_angleFuncs = MatchPhysics.getAngleFunctions() # grab angular functions from mc pdf.
    
    MatchPhysics.setDataFitParameters(dataParameters) # set data pars to pdf (reweighted data have the data physics.)

    # build and write  effciency moments.
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = mcPdf, IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in mc_angleFuncs.functions.itervalues() for func in complexFunc if func ))
    
    scaleFactor = 1 / 16. / sqrt(pi)
    physMoments.initCovariances()
    physMoments.compute(reweightedData) 
    physMoments.write( momentsFile + '_Phys_%s_Iteration'%iterNumb , Scale=scaleFactor )
    physMoments.convertPhysMomsToEffWeights( momentsFile + '_weights_%s_Iteration'%iterNumb , Scale=scaleFactor )

############################################################################################################
## Perform sFit on data using the new angular acceptance.##
#########################################################################
    if doBaselineFit: # perform the sFit with the nominal ang. acceptance and plot.
        dataPdf_StandardAngEff = Bs2JpsiKK2011sFit.multiplyPdfWithAcc( angEffMomentsFile )
        Bs2JpsiKK2011sFit.doFit( iterNumb, dataPdf_StandardAngEff )
        if makePlots:
            canvs['nomFit'] = TCanvas( 'sFit, no AngAccCorr', 'sFit, no AngAccCorr' )
            canvs['nomFit'].Divide(2,2)
            for can, obs, Logy in zip( [ canvs['nomFit'].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
                plot( can, obs,Bs2JpsiKK2011sFit.getDataSet(), dataPdf_StandardAngEff, plotResidHist=True, logy=Logy,  
                      pdfOpts=dict( ProjWData=projDataSet ) )
            canvs['nomFit'].Print('sFit_nom.pdf')

    # multiply pdf with the new acceptance, do Fit and plot.
    angAccFile     = momentsFile + '_weights_%s_Iteration'%iterNumb
    dataPdf_AngEff = Bs2JpsiKK2011sFit.multiplyPdfWithAcc( angAccFile, iterNumb )
    Bs2JpsiKK2011sFit.doFit( iterNumb, dataPdf_AngEff )
    if makePlots:
        canvs['%siter_sFit'%iterNumb] = TCanvas( 'sFit, %s AngAccCorr'%iterNumb, 'sFit, %s AngAccCorr'%iterNumb )
        canvs['%siter_sFit'%iterNumb].Divide(2,2)
        for can, obs, Logy in zip( [ canvs['%siter_sFit'%iterNumb].cd(i) for i in [1,2,3,4]],  angles + time, 3*[False] + [True] ):
            plot( can, obs, Bs2JpsiKK2011sFit.getDataSet(), dataPdf_AngEff, plotResidHist=True, logy=Logy, 
                  pdfOpts=dict( ProjWData=projDataSet ) )
        canvs['%siter_sFit'%iterNumb].Print( 'sFit_%siter.pdf'%iterNumb )


    Bs2JpsiKK2011sFit.setFitParameters( iterNumb, dataParameters ) # update dataParameters with the new onces from the fit.
    mcData = reweightedData                                        # set the reweighted data as input for the next iteration.










assert False
# Question
# C_SP factor in MC pdf ???


# Improvements idea:
# Modulate the building of pdfs in the classes matchMCphysics2Data and matchWeightedDistributions
# Unify the two classes into 1.
# Speed up writting in MatchWeightedDistributions:ReweightAndTransformAngles
# plot after mimicing






# compPhys = CompareWeightedDistributions( mcData, mcData, 'B_P', weight = physWeightName,
#                                          rangeX = [0,4e5], 
#                                          title  = 'B_P before and after physics reweight',
#                                          save   = 'PhysRew_BP_%s.pdf'%iterNumb )
   
