from ROOT import TF1
from ROOT import TGraph
from ROOT import TGraph2D
from ROOT import RooArgSet
from ROOT import TMatrixTSym
from array import array
from itertools import product
from math import isnan

__points = {}
__graphs = {}
__results = {}

def nll_values(nll, nll_min, pars, N):
    zeros = [(p.getVal(), (2 * p.getError() / N)) for p in pars]
    key = (frozenset(pars), N)
    if key in __points:
        return __points[key]
    
    points = []
    for idx in product(*[range(N) for i in range(len(pars))]):
        vals = []
        for j, (i, x) in enumerate(zip(idx, pars)):
            x0, dx = zeros[j]
            v = x0 + (i - (N - 1) / 2) * dx
            x.setVal(v)
            vals.append(v)
        nll_val = nll.getVal()
        if isnan(nll_val):
            continue
        points.append(tuple(vals + [nll.getVal() - nll_min]))
    
    for i, x in enumerate(pars):
        x.setVal(zeros[i][0])
    __points[key] = points
    return points
    
def fit_parabola(nll, nll_min, x, N):
    print 'Calculating diagonal term for {0}.'.format(x.GetName())
    points = nll_values(nll, nll_min, [x], N)
    x0 = x.getVal()
    g_x = TGraph(N, array('d', [p[0] for p in points]), array('d', [p[1] for p in points]))
    __graphs[(frozenset((x, x)), N)] = g_x
    f_x = TF1('f_{0}'.format(x.GetName()), "[3] * (x - {0})^3 + [2] * (x - {0})^2 + [1] * (x - {0}) + [0]". format(x0))
    for i in range(3):
        result_x = g_x.Fit(f_x, 'QS+')
        par = result_x.Parameter(2)
        if result_x.Status() == 0 and par > 0:
            break
    __results[(frozenset([x]), N)] = result_x
    return tuple(result_x.Parameter(i) for i in range(1, 4))

def fit_cross(nll, nll_min, x, y, diagonals, N):
    print 'Calculating cross term for {0} and {1}'.format(x.GetName(), y.GetName())
    x0 = x.getVal()
    y0 = y.getVal()
    points = nll_values(nll, nll_min, [x, y], N)
    
    g_xy = TGraph2D(len(points))
    g_xy.SetName('g_{0}_{1}'.format(x.GetName(), y.GetName()))
    for i, p in enumerate(points):
       g_xy.SetPoint(i, p[0], p[1], p[2])
    
    __graphs[frozenset((x, y))] = g_xy
    
    from ROOT import TF2
    paraboloid = TF2("f_xy",
                     ("[0] + " + \
                      "[1] * (x - {0}) + " + \
                      "[2] * (x - {0})^2 + " + \
                      "[3] * (x - {0})^3 + " + \
                      "[4] * (y - {1}) + " + \
                      "[5] * (y - {1})^2 + " + \
                      "[6] * (y - {1})^3 + " + \
                      "[7] * (x - {0}) * (y - {1}) + " + \
                      "[8] * (x - {0})^2 * (y - {1}) + " + \
                      "[9] * (y - {1})^2 * (x - {0})").format(x0, y0),
                      points[0][0], points[-1][0], points[0][1], points[-1][1])
    for i, v in enumerate(diagonals[x.GetName()]):
        paraboloid.FixParameter(i + 1, v)
    
    for i, v in enumerate(diagonals[y.GetName()]):
        paraboloid.FixParameter(i + len(diagonals[x.GetName()]) + 1, v)

    for i in range(3):
        result_xy = g_xy.Fit(paraboloid, 'QS')
        par = result_xy.Parameter(7)
        if result_x.Status() == 0 and par > 0:
            break
    
    __results[(frozenset([x, y]), N)] = result_xy.Get().Clone()
    return par
    
def hessian(nll, parameters, N = 9):
    values = dict([(p.GetName(), (p.getVal(), p.getError())) for p in parameters])
    nll_min = nll.getVal()
    __parameters = RooArgSet()
    for p in parameters:
        __parameters.add(p)
    nll_params = nll.getObservables(__parameters)
    
    H = TMatrixTSym('double')(parameters.getSize())
    diagonals = {}
    for i, x in enumerate(parameters):
        p = nll_params.find(x.GetName())
        r = fit_parabola(nll, nll_min, p, N)
        diagonals[p.GetName()] = r
        H[i][i] = 2 * r[1]
    
    for i, x in enumerate(parameters):
        x = nll_params.find(x.GetName())
        for j in range(i + 1, parameters.getSize()):
            y = nll_params.find(parameters.at(j).GetName())
            H[i][j] = H[j][i] = fit_cross(nll, nll_min, x, y, diagonals, N)
    
    return H

def graphs():
    return __graphs

def results():
    return __results


for i, p in enumerate(time_result.floatParsFinal()):
    p = nll_params.find(p)
    k = (frozenset([p]), 9)
    g = Hessian.graphs()[k]
    canvas.cd(i + 1)
    g.Draw("A*")
