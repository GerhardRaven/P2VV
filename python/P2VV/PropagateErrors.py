from ROOT import TMatrixT
from math import sqrt, exp

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

def dilutionError(result, derivs):
    """
    propagate parameter errors to the dilution
    result = fit result
    derivs = [(name, deriv), ...]
    """
    n = len(derivs)
    fpf = result.floatParsFinal()
    cov = result.covarianceMatrix()
    C = TMatrixT('double')(n, n)
    J = TMatrixT('double')(1, n)

    indices = [fpf.index(name) for e[0] in derivs if fpf.find(e[0])]
    # Make our own small covariance matrix
    for i, k in enumerate(indices):
        for j, l in enumerate(indices):
            C[i][j] = cov[k][l]
        
    # Jacobian for calculation of sf
    args = [a[0] for a in derivs]
    for i in range(n):
        J[0][i] = derivs[i][1](*[fpf.find(a).getVal() for a in args])
        
    # Calculate J * C * J^T
    JT = J.Clone().T()
    tmp = TMatrixT('double')(3, 1)
    tmp.Mult(C, JT)
    r = TMatrixT('double')(1, 1)
    r.Mult(J, tmp)
        
    return sqrt(r[0][0])

class Parameter(object):
    def __init__(self, name, value, error):
        self.__n = name
        self.__v = value
        self.__e = error
    
    def name(self):
        return self.__n
    
    def value(self):
        return self.__v
    
    def error(self):
        return self.__e

    def setValue(self, v):
        self.__v = v

    def setError(self, e):
        self.__e = e
    
class ErrorSFC(object):
    '''
    Class to calculate errors on D^2 for a fit of a Double Gauss to the prompt sample
    '''
    def __init__(self, dms, sfc, frac, sf2, cv):
        self.__dms = dms
        self.__sfc = sfc
        self.__frac = frac
        self.__sf2 = sf2
        
        # Covariance matrix from the fit, add a column and row for dms
        from ROOT import TMatrixT
        self.__C = TMatrixT('double')(4, 4)
        self.__C[0][0] = dms.error() ** 2
        for i in range(3):
            for j in range(3):
                self.__C[i + 1][j + 1] = cv[i][j]
            self.__C[0][i + 1] = 0.
            self.__C[i + 1][0] = 0.
        
        # Derivatives of D^2 wrt dms, sfc, frac and sf2
        import dilution
        self.__d = [dilution.dDc2_ddms_c, dilution.dDc2_dsfc_c, dilution.dDc2_df_c, dilution.dDc2_dsf2_c]
        # Jacobian matrix to propagate from dms, sfc, frac and sf2 do D^2
        self.__J = TMatrixT('double')(1, 4)
        self.__JT = TMatrixT('double')(4, 1)
        self.__sw = 0
        
    def reset(self):
        self.__J = TMatrixT('double')(1, 4)
        self.__JT = TMatrixT('double')(4, 1)
        self.__sw = 0        
    
    def __call__(self, w, st):
        # Calculate derivatives for this event and add to the sum stored in the
        # Jacobian matrix (since the covariance does not depend on sigma_t).
        for i in range(4):
            d = self.__d[i](st, self.__dms.value(), self.__sfc.value(),
                            self.__frac.value(), self.__sf2.value())
            self.__J[0][i] += w * d
            self.__JT[i][0] += w * d
        self.__sw += w
    
    def error(self):
        # Calculate error on D^2
        tmp = TMatrixT('double')(4, 1)
        tmp.Mult(self.__C, self.__JT)
        r = TMatrixT('double')(1, 1)
        r.Mult(self.__J, tmp)
        return 1 / self.__sw * sqrt(r[0][0])

class ErrorCDG(object):
    def __init__(self, st_mean, dms, sfc_offset, sfc_slope,
                 frac, sf2_offset, sf2_slope, cv):
        self.__sfc_offset = sfc_offset
        self.__sfc_slope = sfc_slope
        self.__sf2_offset = sf2_offset
        self.__sf2_slope = sf2_slope
        self.__dms = dms
        self.__st_mean = st_mean
        self.__sfc = 0
        self.__frac = frac.value()
        self.__sf2 = 0
        # Covariance matrix, add row and column for dms
        from ROOT import TMatrixT
        self.__Cc = TMatrixT('double')(6, 6)
        self.__Cc[0][0] = dms.error() ** 2
        for i in range(5):
            for j in range(5):
                self.__Cc[i + 1][j + 1] = cv[i][j]
            self.__Cc[0][i + 1] = 0.
            self.__Cc[i + 1][0] = 0.

        # Derivatives of D^2 wrt dms, sfc, frac and sf2
        import dilution
        self.__d = [dilution.dDc2_ddms_c, dilution.dDc2_dsfc_c, dilution.dDc2_df_c, dilution.dDc2_dsf2_c]
        self.reset()
        
    def reset(self):
        # Jacobian matrix
        self.__J = TMatrixT('double')(1, 4)
        self.__JT = TMatrixT('double')(4, 1)
        # Temporary matrices to perform multiplication and contrain D^2 error^2
        self.__tmp = TMatrixT('double')(4, 1)
        self.__r = TMatrixT('double')(1, 1)
        # Matrices to propagate from calibration parameters to sfc, frac, sf2
        self.__Jc = TMatrixT('double')(4, 6)
        self.__tmp_c = TMatrixT('double')(6, 4)
        # Prefill constant terms of Jacobian matrix
        self.__Jc[0][0] = 1 #ddms_ddms
        self.__Jc[1][1] = 1 #dsfc_do1
        self.__Jc[2][3] = 1 #dfrac_dfrac
        self.__Jc[3][4] = 1 #dsf2_do2
        # Transpose
        self.__JcT = TMatrixT('double')(TMatrixT('double').kTransposed, self.__Jc)
        self.__C = TMatrixT('double')(4, 4)
        # "Final" variables.
        self.__D2_err2 = 0
        self.__sw = 0
        
    def __call__(self, w, st):
        # Fill the Jacobian to propagate from calibration parameters to sfc, frac and sf2
        self.__Jc[1][2] = st
        self.__Jc[3][5] = st
        self.__JcT.Transpose(self.__Jc)
        self.__Jc.Print()
        # Calculate propagated covariance matrix
        self.__tmp_c.Mult(self.__Cc, self.__JcT)
        self.__C.Mult(self.__Jc, self.__tmp_c)
        self.__C.Print()
        # Calculate values of sfc and sf2
        self.__sfc = self.__sfc_offset.value() + self.__sfc_slope.value() * (st - self.__st_mean)
        self.__sf2 = self.__sf2_offset.value() + self.__sf2_slope.value() * (st - self.__st_mean)
        # Calculate Jacobian of dms, sfc, frac, sf2 to D^2
        for i in range(4):
            d = self.__d[i](st, self.__dms.value(), self.__sfc,
                            self.__frac, self.__sf2)
            self.__J[0][i] = w * d
            self.__JT[i][0] = w * d
        # Propagate to D^2
        self.__J.Print()
        self.__tmp.Mult(self.__C, self.__JT)
        self.__r.Mult(self.__J, self.__tmp)
        # Add contribution of this event to the sum
        print self.__r[0][0]
        self.__D2_err2 += self.__r[0][0]
        self.__sw += w
        
    def error(self):
        return 1 / self.__sw * sqrt(self.__D2_err2)

class ErrorSG(object):
    def __init__(self, dms, sf):
        self.__dms = dms
        self.__sf = sf
        self.__edms = 0
        self.__esf = 0
        self.__sw = 0
        
    def reset(self):
        self.__edms = 0
        self.__esf = 0
        self.__sw = 0
    
    def __call__(self, w, st):
        self.__edms += w * self.__ddms(st)
        self.__esf += w * self.__dsf(st)
        self.__sw += w
        
    def __ddms(self, st):
        sf2 = self.__sf.value() ** 2
        st2 = st ** 2
        dms2 = self.__dms.value() ** 2
        return - st2 * self.__dms.value() * sf2 * exp(- dms2 * sf2 * st2  / 2.)
    
    def __dsf(self, st):
        sf2 = self.__sf.value() ** 2
        st2 = st ** 2
        dms2 = self.__dms.value() ** 2
        return - st2 * dms2 * self.__sf.value() * exp(- dms2 * sf2 * st2  / 2.)
    
    def error(self):
        dmse2 = self.__dms.error() ** 2
        sfe2 = self.__sf.error() ** 2
        return 1 / self.__sw * sqrt(self.__edms ** 2 * dmse2 + self.__esf ** 2 * sfe2)
