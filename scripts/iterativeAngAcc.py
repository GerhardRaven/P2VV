
############################################################################################################
## Specify paths and paramters ##
################################
globalOutputFolder = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/'

# mc input File 
mcTupleFile = 'Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130628.root'
mcTuplePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/' + mcTupleFile
# sData input file
sDataPath    = '/project/bfys/vsyropou/data/P2VVDataSets2011Reco12_noKKMassBins_2TagCats_Kmom.root'
             # '/project/bfys/vsyropou/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_sWeights.root'
             # '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_6KKMassBins_2TagCats.root' 
             # '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_noKKMassBins_2TagCats.root'
sDataName    = 'JpsiKK_sigSWeight' # 'DecayTree' # 'JpsiKK_sigSWeight'
sWeightsName = 'N_sigMass_sw'      # 'sWeight'   # 'N_sigMass_sw'      # 'weightVar'
      
# time acceptance
timeEffPath       = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
timeEffUBName     = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
timeEffExclBName  = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
# angEffMomentsFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

# parameters obtained from the initial sFit to data, no KK mass binning, wide KK mass window. 
dataParameters = dict(  AperpMag2   =  0.246729
                       ,AperpPhase  =  3.08416
                       ,A0Mag2      =  0.523498
                       ,A0Phase     =  0 # cosntrained
                       ,AparPhase   =  3.2097
                       ,dM          = 17.688
                       ,dGamma      = 0.102546
                       ,Gamma       = 0.672938
                       ,phiCP       = 0.0899423
                       ,lambdaCP    = 0.901019
                       )
# varius flags / names
NumbOfIterations = 3 # desired number of iterations.
makePlots = False
physWeightName = 'weightPhys'
############################################################################################################




############################################################################################################
## Begin iterative procedure ##
################################
from P2VV.Utilities.MCReweighting import MatchMCphysics2Data, MatchWeightedDistributions, \
                                         CompareWeightedDistributions, BuildBs2JpsiKK2011sFit
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.Utilities.Plotting import plot
from P2VV.RooFitWrappers import RooObject, RealEffMoment
from ROOT import TFile, TCanvas
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
MatchPhysics = MatchMCphysics2Data( mcTuplePath, modelSwave=False )
MatchPhysics.buildMonteCarloPdf( dataPdfBuilder=Bs2JpsiKK2011sFit.getPdfBuilderObject() )

# get pdfs.
dataPdf = Bs2JpsiKK2011sFit.getPdf()
mcPdf   = MatchPhysics.getPdf()

# get datasets
mcData       = MatchPhysics.getInitialMCafterSel() 
sWeightedData = Bs2JpsiKK2011sFit.getDataSet()

# get observables
angles   = Bs2JpsiKK2011sFit.getObservables('angles') 
time     = Bs2JpsiKK2011sFit.getObservables('time') # reconstructed decay time 
mcObsSet = MatchPhysics.getMcObsSet()               # true/reco decay time + angles 

# THIS IS TEMP: Speed up development.
#mcData = TFile.Open('/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/temp_mcDataAfterSel.root').Get('DecayTree') # speed up 
#reweightedData = TFile.Open('/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/temp_reweightedData_1stIter.root').Get('MomRewMC_1_Iter') # Speed up

# start looping.
for iterNumb in range(1,NumbOfIterations):
    print 'P2VV - INFO: Iteration number ' +str(iterNumb) + '.'  
############################################################################################################
## Match Mc physics to sData. ##
################################

    # calculate and write physics matcing weights
    MatchPhysics.calculateWeights( iterNumb, dataParameters, mcData )
    MatchPhysics.combineWeights() # weights from all iterations are combiend to one weights set.
    MatchPhysics.writeWeights( weightsName=physWeightName) # mcData is internally updated with a the weights column
     
    if makePlots: # check the effect of physics reweighting. 
        compPhys = CompareWeightedDistributions( mcData, mcData, 'B_P', weight = MatchPhysics.getWeightName(), 
                                                                          save = 'PhysRew_%s.pdf'%iterNumb )
   
 ############################################################################################################
 ## Reweight Kaon momenta of the previously reweighted MC to match the Kaon momenta of sData.##
 ##############################################################################################
    reweightArgs = dict( inTree         = mcData, # distribution to be modified 
                         outTree        = sWeightedData, # distribution to be matched with
                         whichVars      = ['Kminus_P'],
                         inWeightName   = physWeightName,
                         outWeightName  = sWeightsName,
                         observables    = mcObsSet,
                         spectatorVars  = MatchPhysics.getNtupleVars(), # vraiables copied to the new dataset.
                         nBins          = 500,  # controls the preceision of the reweighting 
                         itNum          = iterNumb
                          )
    matchMC2Data = MatchWeightedDistributions( **reweightArgs )
    matchMC2Data.reweight()
    reweightedData = matchMC2Data.getDataSet()
    
    if makePlots:  # check if the momentum reweighting worked. # check B_P between sWweighted data and reweighted mc.
        compRew1 = CompareWeightedDistributions( sWeightedData, reweightedData, 'Kminus_P', weight = sWeightsName, 
                                                                        save = 'MomRewCheck_%s.pdf'%iterNumb )
         # Frist add B_P variable to sData and then switch this on.
         # compRew2 = CompareWeightedDistributions( sWeightedData, reweightedData, 'B_P', weight = sWeightsName, 
         #                                                                save = 'MomRew_B_P_%s.pdf'%iterNumb )

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

# ############################################################################################################
# ## Perform sFit on data using the new angular acceptance.##
# #########################################################################
    # miltiply pdf with the new acceptance.
    angAccFile     = momentsFile + '_weights_%s_Iteration'%iterNumb
    dataPdf_AngEff = Bs2JpsiKK2011sFit.multiplyPdfWithAcc( iterNumb, angAccFile )

    # do Fit and update dataParameters with the new onces from the fit
    Bs2JpsiKK2011sFit.doFit( iterNumb, dataPdf_AngEff )
    Bs2JpsiKK2011sFit.setFitParameters( iterNumb, dataParameters )
    
    if makePlots: # sFit plots WITH the new acceptance.
        ObjNam = lambda obj: obj.GetName() +'iterNumb%s'%iterNumb # plots nameing
        for can, obs in zip( [ TCanvas(ObjNam(o),ObjNam(o)) for o in angles+[time] ], angles+[time] ):
            plot(can, obs, data=sWeightedData, pdf=dataPdf_AngEff.getPdf('bin0')  )
            can.Print( ObjNam(obs) + '.pdf' )

    mcData = reweightedData 

# super impose all the pdfs
Bs2JpsiKK2011sFit.plotPdfs()










# some xChecks
# Bs2JpsiKK2011sFit._PDFS['%s_iteration'%iterNumb] = dataPdf_AngEff.getPdf('bin0')    
# Bs2JpsiKK2011sFit.acceptancePlots()



# Improvements idea:
#  Modulate the building of pdfs in the classes matchMCphysics2Data and matchWeightedDistributions
# Unify the two classes into 1.
# Speed up writting in MatchWeightedDistributions:ReweightAndTransformAngles

