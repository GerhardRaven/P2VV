# P_syst = P_MCS * P_MCP^-1 * P_D
import ROOT
from ROOT import Math

Matrix = ROOT.TMatrixT('double')
Vector = ROOT.TVectorT('double')
from array import array
from itertools import chain, product

def get_width_parameters(signal, pmc, data, cov):
    __m = lambda v: Matrix(len(v), 1, array('d', v))
    m = Matrix(cov)
    ev = Vector(5)
    DI = m.EigenVectors(ev)
    D = Matrix(Matrix.kInverted, DI)
    s = __m(signal)
    p = __m(pmc)

    pd = Matrix(D, Matrix.kMult, p)
    sd = Matrix(D, Matrix.kMult, s)

    T = Matrix(5, 5)
    for i, j in product(range(len(pmc)), range(len(pmc))):
        if i == j:
            T[i][i] = sd[i][0] / pd[i][0]
        else:
            T[i][j] = 0.
    tmp = Matrix(T, Matrix.kMult, D)
    tmp = Matrix(DI, Matrix.kMult, tmp)
    return Matrix(tmp, Matrix.kMult, __m(data))

def get_offset_parameters(signal, pmc, data):
    n = len(signal)
    pars = [0. for i in range(n)]
    for i in range(n):
        pars[i] = data[i] + (signal[i] - pmc[i])
    return pars

results = {'2011' : {'PMC'  : {'datatype' : 'MC2011_Sim08a_incl_Jpsi',
                               'result'   : '1bin_15500.00fs_simple/2027465761870101697/results/time_result_double_RMS_Gauss_quadratic_no_offset',
                               'st_mean'  : 0.03599},
                     'MC'   : {'datatype'  : 'MC2011_Sim08a',
                               'result'   : '1bin_19000.00fs_simple/m4074732224475561764/results/time_result_double_RMS_Rest_quadratic_no_offset',
                               'st_mean'  : 0.03318},
                     'data' : {'datatype' : '2011_Reco14',
                               'result'   : '1bin_15500.00fs_simple/2027465761870101697/results/time_result_double_RMS_Gauss_quadratic_no_offset',
                               'st_mean'  : 0.03502}
                               },
           '2012' : {'PMC'  : {'datatype' : 'MC2012_incl_Jpsi',
                               'result'   : '1bin_15500.00fs_simple/2027465761870101697/results/time_result_double_RMS_Gauss_quadratic_no_offset',
                               'st_mean'  : 0.03674},
                     'MC'   : {'datatype' : 'MC2012',
                               'result'   : '1bin_19000.00fs_simple/m4074732224475561764/results/time_result_double_RMS_Rest_quadratic_no_offset',
                               'st_mean'  : 0.03355},
                     'data' : {'datatype' : '2012',
                               'result'   : '1bin_15500.00fs_simple/2027465761870101697/results/time_result_double_RMS_Gauss_quadratic_no_offset',
                               'st_mean'  : 0.0349}
                               }
          }

def transform_width_pars(pars, xb, xbp):
    pars[0] = pars[0] + pars[1] * (xb - xbp)
    pars[3] = pars[3] + pars[4] * (xb - xbp)
    return pars

def transform_offset_pars(par_list, xb, xbp):
    mu = par_list.find('timeResMu')
    pars = [0., 0., 0.]
    if mu:
        pars[0] = mu.getVal()
    else:
        o = par_list.find('timeResMu_offset').getVal()
        s = par_list.find('timeResMu_slope').getVal()
        q = par_list.find('timeResMu_quad').getVal()
        d = xb - xbp
        pars[0] = o + s * d - q * d * (xbp + 3 * xb)
        pars[1] = s - 2 * q * d
        pars[2] = q
    return pars

def transform_cov(parameters, result, xb, xbp):
    from ROOT import RooArgList
    args = RooArgList(*tuple(result.floatParsFinal().find(p) for p in parameters))
    cov = result.reducedCovarianceMatrix(args)
    n = len(parameters)
    J = Matrix(n, n)
    for i in range(n):
        J[i][i] = 1.
    J[0][1] = xb - xbp
    J[3][4] = xb - xbp
    JT = Matrix(Matrix.kTransposed, J)
    tmp = Matrix(cov, Matrix.kMult, JT)
    return Matrix(J, Matrix.kMult, tmp)
    
parameters = ['sf_mean_slope', 'sf_mean_quad', 'timeResFrac2', 'sf_sigma_slope', 'sf_sigma_quad']
from P2VV.Utilities.Resolution import input_data
for year, year_info in results.iteritems():
    print '\n' + year
    rs = {}
    def __get_result(dt, res):
        from P2VV.CacheUtils import CacheFiles
        cache_files = CacheFiles(*input_data[dt]['cache'].rsplit('/', 1))
        return cache_files.getFromCache(res)[0]
    rs = dict((t, __get_result(info['datatype'], info['result'])) for t, info in year_info.iteritems())
    xb = year_info['data']['st_mean']
    def __gwp(t):
        xbp = year_info[t]['st_mean']
        return transform_width_pars([rs[t].floatParsFinal().find(p).getVal() for p in parameters], xb, xbp)

    tcov = transform_cov(parameters, rs['PMC'], xb, year_info['PMC']['st_mean'])
    wp = get_width_parameters(__gwp('MC'), __gwp('PMC'), __gwp('data'), tcov)
    for pn, i in zip(parameters, range(5)):
        print '{0:<20s}: {1:>6.4f}'.format(pn, wp[i][0])

    def __gop(t):
        xbp = year_info[t]['st_mean']
        return transform_offset_pars(rs[t].floatParsFinal(), xb, xbp)
    op = get_offset_parameters(__gop('MC'), __gop('PMC'), __gop('data'))
    for pn, p in zip(['timeResMu_offset', 'timeResMu_slope', 'timeResMu_quad'], op):
        print '{0:<20s}: {1:>6.4f}'.format(pn, p)

