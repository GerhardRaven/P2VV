from ROOT import TMatrixT
from math import sqrt

def propagateScaleFactor(result, suffix = ''):
    fpf = result.floatParsFinal()
    cov = result.covarianceMatrix()
    C = TMatrixT('double')(3, 3)
    J = TMatrixT('double')(1, 3)
    
    frac = fpf.find('timeResFrac2' + suffix)
    if frac:
        frac = frac.getVal()
        names = [n + suffix for n in ['timeResFrac2', 'timeResSigmaSF_1', 'timeResSigmaSF_2']]
        indices = [fpf.index(n) for n in names]
    else:
        names = [n for n in ['timeResFrac2', 'timeResSigmaSF_1' + suffix, 'timeResSigmaSF_2' + suffix]]
        indices = [fpf.index(n) for n in names]
        frac = fpf.find('timeResFrac2').getVal()
    sf1 = fpf.find('timeResSigmaSF_1' + suffix).getVal()
    sf2 = fpf.find('timeResSigmaSF_2' + suffix).getVal()
        
    # Make our own small covariance matrix
    for i, k in enumerate(indices):
        for j, l in enumerate(indices):
            C[i][j] = cov[k][l]
        
    # Jacobian for calculation of sf
    J[0][0] = sf1 + sf2
    J[0][1] = 1 - frac
    J[0][2] = frac
        
    # Calculate J * C * J^T
    JT = J.Clone().T()
    tmp = TMatrixT('double')(3, 1)
    tmp.Mult(C, JT)
    r = TMatrixT('double')(1, 1)
    r.Mult(J, tmp)
        
    sf = (1 - frac) * sf1 + frac * sf2
    sf_e = sqrt(r[0][0])
    return (sf, sf_e)
