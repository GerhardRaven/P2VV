# script parameters
momentsFilesIn = [ 'Sim08_2011_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm', 'Sim08_2012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm' ]
momentsFileOut = 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm'
funcNames = [ 'Re_ang_A0_A0', 'Re_ang_Apar_Apar', 'Re_ang_Aperp_Aperp'
             , 'Im_ang_Apar_Aperp', 'Re_ang_A0_Apar', 'Im_ang_A0_Aperp'
             , 'Re_ang_AS_AS', 'Re_ang_Apar_AS', 'Im_ang_Aperp_AS', 'Re_ang_A0_AS' ]

# moment containers
moments = [ ]
correlations = [ ]

# read moments
from P2VV.Utilities.DataMoments import readMoments, printMoments
for momsFile in momentsFilesIn :
    moments.append( { } )
    correlations.append( { } )
    readMoments( momsFile, BasisFuncNames = funcNames, Moments = moments[-1], Correlations = correlations[-1] )
    printMoments( BasisFuncNames = funcNames, Moments = moments[-1], Correlations = correlations[-1] )

# combine moments
from ROOT import TMatrixD
combNames = [ name for name in funcNames if moments[0][name][1] > 0. ]
momVec    = TMatrixD( len(combNames), 1              )
covMat    = TMatrixD( len(combNames), len(combNames) )
for moms, corrs in zip( moments, correlations ) :
    mom = TMatrixD( len(combNames), 1              )
    cov = TMatrixD( len(combNames), len(combNames) )
    for it1, name1 in enumerate(combNames) :
        mom[it1][0] = moms[name1][0]
        for it2, name2 in enumerate(combNames) :
            cov[it1][it2] = moms[name1][1] * moms[name2][1] * corrs[name1][name2]
    cov.Invert()
    momVec += cov * mom
    covMat += cov
covMat.Invert()
momVec = covMat * momVec

momsNew  = { }
corrsNew = { }
for name in funcNames :
    if moments[0][name][1] > 0. : continue
    momsNew[name]  = moments[0][name]
    corrsNew[name] = correlations[0][name]

from math import sqrt
for it1, name1 in enumerate(combNames) :
    val = momVec[it1][0]
    err = sqrt( covMat[it1][it1] )
    momsNew[name1]  = ( val, err, abs( val / err ) )
    corrsNew[name1] = { }
    for name in funcNames :
        if moments[0][name][1] > 0. : continue
        corrsNew[name1][name] = correlations[0][name1][name]
    for it2, name2 in enumerate(combNames) :
        corrsNew[name1][name2] = covMat[it1][it2] / err / sqrt( covMat[it2][it2] )

# write combined moments to file
from P2VV.Utilities.DataMoments import writeMoments
writeMoments( momentsFileOut, BasisFuncNames = funcNames, Moments = momsNew, Correlations = corrsNew )
printMoments( BasisFuncNames = funcNames, Moments = momsNew, Correlations = corrsNew )
