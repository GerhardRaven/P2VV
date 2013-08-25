

# What this script does:
# It takes as imput a MC sample and an sWeighted dataset.
# And does the folowing:
# (1) Constructs a MC pdf                                                                         ( DONE )
# (2) Reweights the MC sample with a weight w = pdf(t,Omega;pars_sData) / pdf(t,Omega;pars_MC)    ( DONE )
#  in order to match the physics between data and MC                                              ( DONE )
# (3) Writes the output reweighted MC sample in 'physReweightOutputFile'                          ( DONE )
# (4) Takes the physics reweighted MC sample and reweights it again to mach the K momentum        ( DONE )
# (5) Takes the physics + K momentum matched MC sample and calculates the angular acceptance      ( DONE ) 
# (6) Performs the sFit on data using the new accptance.                                          ( Under Devel )
# At this stage a loop takes you to step (1)                                                     
#----------------------------------------------------------------------------------------------------------------------------



#Import stuff
# MUST SetUp Urania FIRST and then P2VV
import sys
import os    # exit
sys.path.append(os.environ["BS2MUMUROOT"] +"/python/Bs2MuMu/")

# Global output folder and default tree name
# Modify this according to where you want to put your output datasets.
globalOutputFolder = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/'
# Default tree name
treeName = 'DecayTree'

# Specify input data
# Mc input File 
mcTupleFile = 'Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130628.root'
mcTuplePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/' + mcTupleFile
# sData input file
sDataFile   = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_sWeights.root'


# Parameters obtained from the initial sFit to data, on KK mass binning, wide KK mass window. 
initDataParameters = dict( prefix   = 'data_initVals'
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





####################################### Begin iterative procedure ##########################################
############################################################################################################
from P2VV.GeneralUtils import matchMCphysics2Data, matchWeightedDistributions, compareWeightedDistributions
from ROOT import TFile

for iterNumb in range(1,3):
    print 'P2VV - INFO: Iteration number ' +str(iterNumb) + '.'  

    # Name sufix of comparision plots
    plotsSufix = str(iterNumb) + '_Iteration' + '.pdf'

############################################################################################################
## Match Mc physics to sData. ##
################################
    # Input/Output file names 
    physReweightInputFile  = mcTuplePath if iterNumb==1 else momReweightOutputFile + '.root'
    physReweightOutputFile = globalOutputFolder + mcTupleFile.partition('.')[0] + '_physWeights'\
                                                + '_{0}'.format(iterNumb) + '_Iteration.root'
    print 'P2VV - INFO: Physics matching input file: ',  physReweightInputFile
    print 'P2VV - INFO: Physics matching output file: ', physReweightOutputFile
    
    # Specify weight names for each iteration and create a string with all the weights.
    WeightName = 'weightPhys_' + str(iterNumb) + '_iter'
    if iterNumb==1: allWeights = WeightName  
    else:           allWeights += '*' + WeightName

    # Specify pdf paramters from the sFit
    if iterNumb==1:dataParameters = initDataParameters
    else:
        # This is reminder to grab the sFit result.
        print 'Impliment dataParameters = getDataParametersFromRooFitResult()'
        dataParameters = initDataParameters

    # Build angular pdf, calculate and write physics matcing weights
    matchPhysics = matchMCphysics2Data( physReweightInputFile )
    matchPhysics.buildMonteCarloPdf()
    matchPhysics.calculateWeights(dataParameters)
    matchPhysics.writeWeightsToFile( physReweightOutputFile, weightsName=WeightName ) 
    
     # Comment in to check the effect of physics reweighting.
     #t_mc = TFile.Open(physReweightOutputFile).wGet(treeName)
     #c,a = compareWeightedDistributions(t_mc, t_mc, 'B_P', weight=allWeights,\
     #                                                      assymPlot=True,   \
     #                                                      Save=[True,'PhysRew_' + plotsSufix] )

 ############################################################################################################
 ## Reweight Kaon momenta of the previously reweighted MC to match the Kaon momenta of sData.##
 ##############################################################################################
 # Specidy input/output files
    momReweightInputFile   = physReweightOutputFile 
    momReweightOutputFile  = globalOutputFolder + mcTupleFile.partition('.')[0] + '_physWeights_momReWeight' \
                                                 + '_%s'%iterNumb + '_Iteration'

    print 'P2VV - INFO: Kinematic matching input file: ',  momReweightInputFile
    print 'P2VV - INFO: Kinematic matching output file: ', momReweightOutputFile
    
    reweightArgs = dict( sDInfo     = dict(path=sDataFile,            name=treeName, weight='sWeight'    )
                          ,mcInfo    = dict(path=momReweightInputFile, name=treeName, weight=allWeights   )
                          ,whichVars = 'KaonMomenta'
                          ,nBins     = 500  ### WATCH OUT IN THE END YOU WILL SET THIS TO 1000
                          ,itNum     = iterNumb
                          )
     matchMC2Data = matchWeightedDistributions( momReweightOutputFile, **reweightArgs )
     matchMC2Data.mimicWeights()
     matchMC2Data.reweightMC(obsSet=matchPhysics.getObservables())
   
     # Check if K_minus before and after reweighting DO match.
     t_mc = TFile.Open(momReweightOutputFile + '.root').Get('T')
     t_sD = TFile.Open(sDataFile).Get(treeName)
     c,a = compareWeightedDistributions(t_mc, t_sD, 'Kminus_P_mod%s'%iterNumb, sVar='Kminus_P',  \
                                            weight=allWeights, \
                                            sWeight='sWeight', \
                                            rangeX=[0,10e4],   \
                                            assymPlot=True,    \
                                            Save=[True,'KmomRew_' + plotsSufix] 
                                        )
  
############################################################################################################
## Compute angular efficiency moments for the new reweighted MC sample.##
#########################################################################
    # Specify efficeincy moments output file.
    effWieghtsOutputFile = 'hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights_%s_Iteration'%iterNumb 
    
    # Grab some stuff from the MC pdf.
    pdf        = matchPhysics.getPdf()
    angleFuncs = matchPhysics.getAngleFunctions()
    angles     = matchPhysics.getObservables()[1:4]
    data       = TFile.Open(momReweightOutputFile + '_RDS.root').Get('MomRewMC_%s_Iter'%iterNumb)
    
    from P2VV.GeneralUtils import RealMomentsBuilder
    from P2VV.RooFitWrappers import RealEffMoment
    physMoments = RealMomentsBuilder( Moments = ( RealEffMoment( Name = func.GetName(), BasisFunc = func,
                                                                 Norm = 1., PDF = pdf, IntSet = [ ], NormSet = angles )\
                                                      for complexFunc in angleFuncs.functions.itervalues() for func in complexFunc if func ))
    physMoments.initCovariances()

    # Compute moments from data set and write them in a file
       #data = matchPhysics.getInitialMCafterSel() # For xcheck you can use the initial MC dataset (after selection).
    physMoments.compute(data)
       #PhysMoments.Print(  Scale = 1 / 16. / sqrt(pi) ) # Optionaly print the intermediate result.
    # Convert efficinecy moments to efficiency weights and write them in a file.
    scaleFactor = 0.12500024417259972 # Not so sure how and why this factor is necessary, but it gives the correct result.
    physMoments.convertPhysicsMomentsToBasisMoments(globalOutputFolder + effWieghtsOutputFile,scale=scaleFactor)
    

############################################################################################################
## Perform sFit on data using the new angular acceptance.##
#########################################################################
# set script parameters




assert False






