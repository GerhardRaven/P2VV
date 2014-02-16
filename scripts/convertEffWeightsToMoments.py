# script parameters
weightsFile = 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm'
momentsFile = 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
funcNames = [ 'Re_ang_A0_A0', 'Re_ang_Apar_Apar', 'Re_ang_Aperp_Aperp'
             , 'Im_ang_Apar_Aperp', 'Re_ang_A0_Apar', 'Im_ang_A0_Aperp'
             , 'Re_ang_AS_AS', 'Re_ang_Apar_AS', 'Im_ang_Aperp_AS', 'Re_ang_A0_AS' ]

# read weights
from P2VV.Utilities.DataMoments import readMoments, printMoments
weights = { }
weightCorrs = { }
readMoments( weightsFile, BasisFuncNames = funcNames, Moments = weights, Correlations = weightCorrs )
printMoments( BasisFuncNames = funcNames, Moments = weights, Correlations = weightCorrs )

# convert weights to moments
from P2VV.Utilities.DataMoments import convertEffWeightsToMoments
convertEffWeightsToMoments( Weights = weights, OutputFilePath = momentsFile )
