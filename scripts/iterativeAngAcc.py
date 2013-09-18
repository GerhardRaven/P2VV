




############################################################################################################
## Specify paths and paramters ##
################################
globalOutputFolder = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/'

# Mc input File 
mcTupleFile = 'Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130628.root'
mcTuplePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/' + mcTupleFile
# sData input file
sDataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_noKKMassBins_2TagCats.root'
             # '/project/bfys/vsyropou/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_sWeights.root'
             # '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_6KKMassBins_2TagCats.root' 
             # '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_noKKMassBins_2TagCats.root'
sDataName    = 'JpsiKK_sigSWeight' # 'DecayTree' # 'JpsiKK_sigSWeight'
sWeightsName = 'N_sigMass_sw' # 'sWeight'   # 'N_sigMass_sw' # 'weightVar'
sDataType    =  'RooDataSet' # 'TTree'     # 'RooDataSet'
      
# Time acceptance
timeEffPath       = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
timeEffUBName     = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
timeEffExclBName  = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
#angEffMomentsFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

# Parameters obtained from the initial sFit to data, on KK mass binning, wide KK mass window. 
initDataParameters = dict(ParNamePrefix   = 'data_initVals'
                       ,AperpMag2   =  0.246729
                       ,AperpPhase  =  3.08416
                       ,A0Mag2      =  0.523498
                       ,A0Phase     =  0 # cosntrained
                       ,AparPhase   =  3.2097
                       ,ASOddPhase  =-0.0988961
                       ,f_S         = 0.0471
                       ,C_SP        = 0.326
                       ,dM          = 17.688
                       ,dGamma      = 0.102546
                       ,Gamma       = 0.672938
                       ,phiCP       = 0.0899423
                       ,lambdaCPSq  = 0.901019**2
                       )
# varius flags / parameters
makePlots = False
physWeightName = 'weightPhys'

############################################################################################################
## Begin iterative procedure ##
################################
from P2VV.Utilities.MCReweighting import MatchMCphysics2Data, MatchWeightedDistributions, \
                                         CompareWeightedDistributions, BuildBs2JpsiKK2011sFit
from P2VV.RooFitWrappers import RooObject
from ROOT import TFile

worksp = RooObject( workspace = 'iterativeAngularAcceptance' ).ws()


#  Build data pdf for sFiting ( This pdf is not multiplied by the angular acceptance )
Bs2JpsiKK2011sFit = BuildBs2JpsiKK2011sFit( dataSetPath          = sDataPath,
                                            dataSetName          = sDataName,
                                            timeEffHistFile      = timeEffPath, 
                                            timeEffHistUBName    = timeEffUBName,
                                            timeEffHistExclBName = timeEffExclBName, 
                                            KKmassBins           = None  , # '4KKMassBins',
                                            parFileIn            = None  , 
                                            parFileOut           = None  , 
                                            )
# Build MC pdf
matchPhysics = MatchMCphysics2Data( mcTuplePath )
matchPhysics.buildMonteCarloPdf( dataPdfBuilder=Bs2JpsiKK2011sFit.getPdfBuilderObject() )

mcPdf = matchPhysics.getPdf() # Monte Carlo pdf
dataPdf = Bs2JpsiKK2011sFit.getPdf() 


assert False

for iterNumb in range(1,3):
    print 'P2VV - INFO: Iteration number ' +str(iterNumb) + '.'  

############################################################################################################
## Match Mc physics to sData. ##
################################
    if iterNumb==1: 
        dataParameters = initDataParameters # Specify pdf paramters from the sFit at the end of the loop
        mcData = matchPhysics.getInitialMCafterSel()
    else: 
       mcData =  reweightedData

    # Calculate and write physics matcing weights
    matchPhysics.calculateWeights( iterNumb, dataParameters, mcData )
    matchPhysics.combineWeights() # Weights from all iterations are combiend.
    matchPhysics.writeWeights( weightsName=physWeightName) # RooDataSet object.
    physReweightOutput = mcData.buildTree()
     
    # Check the effect of physics reweighting.
    if makePlots: 
        c,a = compareWeightedDistributions( physReweightOutput, physReweightOutput, 'B_P', \
                                                weight    = matchPhysics.getWeightName(),  \
                                                assymPlot = True)
   
 ############################################################################################################
 ## Reweight Kaon momenta of the previously reweighted MC to match the Kaon momenta of sData.##
 ##############################################################################################
    reweightArgs = dict( inTree         = physReweightOutput, # Distribution to be modified 
                         outTree        = TFile.Open(sDataPath).Get(sDataName), # Distribution to be matched with
                         whichVars      = ['Kminus_P'],
                         inWeightName   = matchPhysics.getWeightName(),
                         outWeightName  = sWeightsName,
                         PhysWeightName = physWeightName,
                         nBins          = 500,  # SET THIS TO 1000
                         itNum          = iterNumb,
                          )
    matchMC2Data = matchWeightedDistributions( **reweightArgs )
    matchMC2Data.mimicWeights()
    matchMC2Data.reweightMC( matchPhysics.getObservables(), copyVars= matchPhysics.getNtupleVars() )
    reweightedData = matchMC2Data.getDataSet()
    
    if makePlots: 
        # Check if the reweighting worked.
        comp = compareWeightedDistributions( TFile.Open(sDataPath).Get(sDataName), \
                                                  reweightedData.buildTree(), 'Kminus_P', weight = sWeightsName)
        #Check B_P  between sWeighted data and mc.
 




    assert False



    reweightedData = TFile.Open('/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/MomRewMC_1_Iter.root').Get('MomRewMC_1_Iter')
################################################################################################################
## Compute angular efficiency moments for the new reweighted MC sample.##
#########################################################################
    # Specify efficeincy moments output file.
    effWieghtsOutputFile = 'hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights_%s_Iteration'%iterNumb 
    
    # Grab some stuff from the MC pdf.
    angleFuncs = matchPhysics.getAngleFunctions()
    angles     = matchPhysics.getObservables()[1:4]
    
    from P2VV.Utilities.DataMoments import RealMomentsBuilder
    from P2VV.RooFitWrappers import RealEffMoment
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = mcPdf, IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func ))
    physMoments.initCovariances()

    # Compute moments from data set and write them in a file
       #data = matchPhysics.getInitialMCafterSel() # For xcheck you can use the initial MC dataset (after selection).
    physMoments.compute(reweightedData)
       #PhysMoments.Print(  Scale = 1 / 16. / sqrt(pi) ) # Optionaly print the intermediate result.
    # Convert efficinecy moments to efficiency weights and write them in a file.
    scaleFactor = 0.12500024417259972 # CHECK THISN ot so sure how and why this factor is necessary, but it gives the correct result.
    physMoments.convertPhysicsMomentsToBasisMoments(globalOutputFolder + effWieghtsOutputFile,scale=scaleFactor)
    
# ############################################################################################################
# ## Perform sFit on data using the new angular acceptance.##
# #########################################################################
   
    # #Read efficiency moments file
    from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles     as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = observables['cpsi'], ctheta = observables['ctheta'], phi = observables['phi'] )
    
    # #Multiply pure pdf with angular acceptance
    from P2VV.Utilities.DataMoments import angularMomentIndices
    moments = RealMomentsBuilder()
    moments.appendPYList( angleFuncs.angles, angularMomentIndices(multiplyByAngEff,angleFuncs ) )
    moments.read(globalOutputFolder + effWieghtsOutputFile)
    moments.Print()
    dataPdfAngAcc = moments * dataPdf


    # do Fit
    Bs2JpsiKK2011sFit.doFit( iterNumb )

    # grab the fit result
    sFitResult = Bs2JpsiKK2011sFit.getFitResult( iterNumb )   

    sFitParameters = dict( prefix   = 'data_%s_Iter'%iterNumb
                            ,AperpMag2   = sFitResult.find('AperpMag2')
                            ,AperpPhase  = sFitResult.find('AperpPhase')
                            ,A0Mag2      = sFitResult.find('A0Mag2')
                            ,A0Phase     =  0 # cosntrained
                            ,AparPhase   = sFitResult.find('AparPhase')
                            ,ASOddPhase  = sFitResult.find('ASOddPhase')
                            ,f_S         = sFitResult.find('f_S')
                            ,C_SP        = sFitResult.find('C_SP')
                            ,dM          = sFitResult.find('')
                            ,dGamma      = sFitResult.find('')
                            ,Gamma       = sFitResult.find('')
                            ,phiCP       = sFitResult.find('')
                            ,lambdaCPSq  = sFitResult.find('')
                            )





# Improvements idea:
#  Modulate the building of pdfs in the classes matchMCphysics2Data and matchWeightedDistributions
# Unify the two classes into 1.


# Temp shortcuts
# dataPdf.getPdf(worksp.cat('KKMassCat').lookupType(0).GetName()).Print()


# >>>dataPdf.getPdf(worksp.cat('KKMassCat').lookupType(0).GetName()).Print()
# RooBTagDecay::data_sig_t_angles_bin0[ time=time iTag0=iTagOS iTag1=iTagSS fTag=NULL tagCat0=tagCatP2VVOS tagCat1=tagCatP2VVSS tau=data_MeanLifetime dGamma=data_dGamma dm=data_dM dilutions0=(data_tagDilution0_0,data_tagDilutionOS1) dilutions1=(data_tagDilution1_0,data_tagDilutionSS1) ADilWTags0=(data_ADilWTag0_0,data_ADilWTagOS1) ADilWTags1=(data_ADilWTag1_0,data_ADilWTagSS1) ANorm=NULL avgCEvenSum=data_avgCEvenSum avgCOddSum=data_avgCOddSum coshCoef=data_coshCoef_bin0 sinhCoef=data_sinhCoef_bin0 cosCoef=data_cosCoef_bin0 sinCoef=data_sinCoef_bin0 createdVars=() tagCatCoefs0=(data_sig_t_angles_tagCatCoef0,data_tagCatCoef0-1) tagCatCoefs1=(data_tagCatCoef1-0,data_tagCatCoef1-1) avgCEvens0=(data_sig_t_angles_avgCEven0,data_avgCEvenSSTagged) avgCOdds0=(data_sig_t_angles_avgCOdd0,data_avgCOddSSTagged) avgCEvens1=(data_avgCEvenOSTagged,data_avgCEvenTagged) avgCOdds1=(data_avgCOddOSTagged,data_avgCOddTagged) ] = 0.402156
# >>> mcPdf.Print()
# RooBTagDecay::_sig_t_angles_tagCat_iTag[ time=truetime iTag0=NULL iTag1=iTag fTag=NULL tagCat0=NULL tagCat1=NULL tau=mc_MeanLifetime dGamma=mc_dGamma dm=mc_dM dilutions0=() dilutions1=(mc_one) ADilWTags0=() ADilWTags1=(mc_zero) ANorm=NULL avgCEvenSum=one avgCOddSum=zero coshCoef=mc_coshCoef sinhCoef=mc_sinhCoef cosCoef=mc_cosCoef sinCoef=mc_sinCoef createdVars=() avgCEvens0=(one) avgCOdds0=(zero) ] = 9.66352
