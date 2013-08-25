# script parameters
momentsFile1 = 'MC11_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm'
momentsFile2 = 'MC12_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm'
funcNames = [ 'Re_ang_A0_A0', 'Re_ang_Apar_Apar', 'Re_ang_Aperp_Aperp'
             , 'Im_ang_Apar_Aperp', 'Re_ang_A0_Apar', 'Im_ang_A0_Aperp'
             , 'Re_ang_AS_AS', 'Re_ang_Apar_AS', 'Im_ang_Aperp_AS', 'Re_ang_A0_AS' ]
negativeFuncs = [ 'Im_ang_A0_Aperp', 'Im_ang_Apar_Aperp' ]

# moment containers
moments1 = { }
moments2 = { }
correlations1 = { }
correlations2 = { }

# read moments
print 'input moments 1:'
from P2VV.Utilities.DataMoments import readMoments, printMoments
readMoments( momentsFile1, BasisFuncNames = funcNames, Moments = moments1, Correlations = correlations1 )
printMoments( BasisFuncNames = funcNames, Moments = moments1, Correlations = correlations1 )

print 'input moments 2:'
from P2VV.Utilities.DataMoments import readMoments, printMoments
readMoments( momentsFile2, BasisFuncNames = funcNames, Moments = moments2, Correlations = correlations2 )
printMoments( BasisFuncNames = funcNames, Moments = moments2, Correlations = correlations2 )

from math import sqrt
print 'differences:'
for name in funcNames :
    val1 = moments1[name][0]
    val2 = moments2[name][0]
    if name in negativeFuncs :
        val1 = -val1
        val2 = -val2
    diffVal = val2 - val1
    diffErr = sqrt( moments1[name][1]**2 + moments2[name][1]**2 )
    diffSig = ( diffVal / diffErr ) if diffErr != 0. else 0. if abs(diffVal) < 1.e-10 else float('inf')
    print '  {0:20s}  {1:+7.4f} +/- {2:+7.4f}    {3:+7.4f} +/- {4:+7.4f}    {5:+7.4f} +/- {6:+7.4f} ({7:+.1f} sigma)'\
        .format( name, val1, moments1[name][1], val2, moments2[name][1], diffVal, diffErr, diffSig )
