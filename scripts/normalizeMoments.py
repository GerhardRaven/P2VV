# script parameters
momentsFileOrig = 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Phys'
momentsFileNorm = 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm'
normMoment = 'Re_ang_A0_A0'

# moment containers
funcNames = [ ]
moments = { }
correlations = { }

# read moments
print 'input moments:'
from P2VV.Utilities.DataMoments import readMoments, printMoments
readMoments( momentsFileOrig, BasisFuncNames = funcNames, Moments = moments, Correlations = correlations, ProcessAll = True )
printMoments( BasisFuncNames = funcNames, Moments = moments, Correlations = correlations )

# calculate moments with new normalization
from math import sqrt
valNorm = moments[normMoment][0]
errNorm = moments[normMoment][1]
momsNew  = { }
corrsNew = { }
# moments value and error
for name in funcNames :
    valNew   = moments[name][0] / valNorm
    err      = moments[name][1]
    corrNorm = correlations[name][normMoment]
    varNew   = ( err**2 + valNew**2 * errNorm**2 - 2. * valNew * err * errNorm * corrNorm ) / valNorm**2
    momsNew[name] = ( valNew, sqrt(varNew), ( sqrt( valNew**2 / varNew ) ) if varNew != 0. else 0. if valNew == 0. else float('inf') )
    corrsNew[name] = { }

# moments correlations
for name0 in funcNames :
    for name1 in funcNames :
        valsNew   = ( momsNew[name0][0], momsNew[name1][0] )
        errsNew   = ( momsNew[name0][1], momsNew[name1][1] )
        errs      = ( moments[name0][1], moments[name1][1] )
        corr01    = correlations[name0][name1]
        corr0Norm = correlations[name0][normMoment]
        corr1Norm = correlations[name1][normMoment]
        corr01New = ( errs[0] * errs[1] * corr01 + valsNew[0] * valsNew[1] * errNorm**2\
                     - valsNew[1] * errs[0] * errNorm * corr0Norm - valsNew[0] * errs[1] * errNorm * corr1Norm ) / valNorm**2
        if abs( errsNew[0] ) > 1.e-10 and abs( errsNew[1] ) > 1.e-10 :
            corr01New /= errsNew[0] * errsNew[1]
        elif abs( corr01New ) < 1.e-10 and name0 == normMoment and name1 == normMoment :
            corr01New = 1.
        elif abs( corr01New ) < 1.e-10 and ( name0 == normMoment or name1 == normMoment ) :
            corr01New = 0.
        else :
            corr01New = -999.
        corrsNew[name0][name1] = corr01New
        corrsNew[name1][name0] = corr01New

print 'output moments:'
from P2VV.Utilities.DataMoments import writeMoments
writeMoments( momentsFileNorm, BasisFuncNames = funcNames, Moments = momsNew, Correlations = corrsNew )
printMoments( BasisFuncNames = funcNames, Moments = momsNew, Correlations = corrsNew )
