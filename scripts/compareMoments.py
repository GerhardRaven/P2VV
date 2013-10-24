# script parameters
momentsFile1 = 'Sim08_2011_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm'
momentsFile2 = 'Sim08_2012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm'
funcNames = [ 'Re_ang_A0_A0', 'Re_ang_Apar_Apar', 'Re_ang_Aperp_Aperp'
             , 'Im_ang_Apar_Aperp', 'Re_ang_A0_Apar', 'Im_ang_A0_Aperp'
             , 'Re_ang_AS_AS', 'Re_ang_Apar_AS', 'Im_ang_Aperp_AS', 'Re_ang_A0_AS' ]

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

from ROOT import TMatrixD
diffInd   = 0
diffNames = [ name for name in funcNames if moments1[name][1] != 0. ]
diffs     = TMatrixD( len(diffNames), 1 )
diffCovs  = TMatrixD( len(diffNames), len(diffNames) )

from math import sqrt
print 'differences:'
for name in funcNames :
    val1 = moments1[name][0]
    val2 = moments2[name][0]
    diffVal = val2 - val1
    diffErr = sqrt( moments1[name][1]**2 + moments2[name][1]**2 )
    diffSig = ( diffVal / diffErr ) if diffErr != 0. else 0. if abs(diffVal) < 1.e-10 else float('inf')
    print '  {0:20s}  {1:+7.4f} +/- {2:+7.4f}    {3:+7.4f} +/- {4:+7.4f}    {5:+7.4f} +/- {6:+7.4f} ({7:+.1f} sigma)'\
        .format( name, val1, moments1[name][1], val2, moments2[name][1], diffVal, diffErr, diffSig )

    if name in diffNames :
        diffs[diffInd][0] = diffVal
        for it, name2 in enumerate(diffNames) :
            diffCovs[diffInd][it] = moments1[name][1] * moments1[name2][1] * correlations1[name][name2]\
                                    + moments2[name][1] * moments2[name2][1] * correlations2[name][name2]
        diffInd += 1

diffsTrans = TMatrixD( 1, len(diffNames) )
diffsTrans.Transpose(diffs)
diffCovsInv = TMatrixD(diffCovs).Invert()
chiSq = diffsTrans * diffCovsInv * diffs

from ROOT import TMath
print 'differences chi^2 / #dof = %.1f / %d = %.1f (%.1f%%)'\
      % ( chiSq[0][0], len(diffNames), chiSq[0][0] / float( len(diffNames) ), TMath.Prob( chiSq[0][0], len(diffNames) ) * 100. )
