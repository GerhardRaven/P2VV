#!/usr/bin/env python
from array import array
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog data_type')
(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.usage
    sys.exit(-2)
elif args[0] not in ['2011', '2012', 'MC11a']:
    print parser.usage
    sys.exit(-2)

from ROOT import *
from collections import defaultdict
from P2VV.Load import LHCbStyle
from P2VV.RooFitDecorators import *

if args[0] == 'MC11a':
    input_file = 'Bs2JpsiPhi_MC11a_Prescaled_st.root'
elif args[0] == '2011':
    input_file = 'Bs2JpsiPhi_2011_Prescaled.root'
elif args[0] == '2012':
    input_file = 'Bs2JpsiPhi_2012_Prescaled.root'
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
directory = os.path.join(prefix, 'p2vv/data')

from P2VV.CacheUtils import CacheFiles
cfs = CacheFiles(directory, input_file)
cache_files = cfs.getCacheFiles()

sdatas = defaultdict(dict)
results = defaultdict(dict)

def add_keys(input_file, class_name, keys = None, path = None):
    if keys == None or type(keys) != dict:
        keys = {}
    if path == None or path == '/':
        d = input_file
        path = ''
    else:
        d = input_file.Get(path)
    for key in d.GetListOfKeys():
        if path:
            new_path = path + '/' + key.GetName()
        else:
            new_path = key.GetName()
        if new_path not in keys and key.GetClassName() == class_name:
            keys[new_path] = key.ReadObj()
        if key.GetClassName() == 'TDirectoryFile':
            add_keys(input_file, class_name, keys, new_path)
    return keys

dirs = {}
for f in cache_files:
    add_keys(f, 'TDirectoryFile', keys = dirs)

interesting = {}
for key, d in dirs.iteritems():
    if len(key.split('/')) != 2:
        continue
    if key.find('simul') == -1:
        continue
    cut = d.Get('cut')
    if not cut:
        continue
    interesting[key] = d

titles = {}
for k, d in interesting.items():
    cut = d.Get('cut')
    cut = str(cut)
    cuts = [c.strip() for c in cut.split('&&')]
    dc = []
    for c in cuts:
        if c.startswith('sel'):
            break
        else:
            dc.append(c)
    if dc:
        dc = ' && '.join(dc)
    else:
        dc = 'no extra cut'
    titles[k] = dc
    print k, dc

if args[0] == '2011':
    ## good = {'1243785060103642893' : 4, 'm2334064025374600976' : 3,
    ##         '1626518906014697943' : 2, 'm3832912631969227654' : 1,
    ##         '4086600821164745518' : 6, 'm6573713017788044320' : 5}
    good = {'m934737057402830078' : 1}
    sig_name = 'psi_ll'
elif args[0] == '2012':
    good = {'m934737057402830078' : 1}
    sig_name = 'psi_ll'
elif args[0] == 'MC11a':
    good = {'389085267962218368' : 4, 'm3019457528953402347' : 3,
            'm7780668933605436626' : 2, 'm8376372569899625413' : 1,
            'm1545059518337894505' : 6, 'm8342219958663192955' : 5}
    sig_name = 'signal'

PDFs = defaultdict(dict)
for k, cache_dir in filter(lambda k: k[0].split('/')[0] in ['9bins_14.10fs_simul'], interesting.iteritems()):
    hd = k.split('/')[-1]
    try:
        index = good[hd]
    except KeyError:
        continue
    sdata_dir = cache_dir.Get('sdata')
    for e in sdata_dir.GetListOfKeys():
        if e.GetClassName() == 'RooDataSet':
            sdatas[index][e.GetName()] = e.ReadObj()
    rd = cache_dir.Get('results')
    for e in rd.GetListOfKeys():
        if e.GetClassName() == 'RooFitResult':
            results[k][e.GetName()] = e.ReadObj()
    ## pdf_dir = cache_dir.Get('PDFs')
    ## for e in pdf_dir.GetListOfKeys():
    ##     PDFs[index][e.GetName()] = e.ReadObj()

from ROOT import kGray
from ROOT import TH1D
from ROOT import TMatrixT
from math import sqrt
from P2VV.PropagateErrors import propagateScaleFactor

__canvases = []
__histos = []
__fit_funcs = []

fit_type = 'double_Comb_Gauss'

__fit_results = defaultdict(list)
from array import array

def fr_latex(frs):
    frs = sorted(frs, key = lambda fr: fr.NPar())
    names = ' & '.join(['result ' + fr.GetName().split('_', 1)[1] for fr in frs])
    s = "\\begin{tabular}{|l|%s|}\n\hline\nparameter & %s \\\\ \n\hline\hline\n" % ('|'.join(['r' for i in range(len(frs))]), names)
    pars = dict((i, [frs[i].ParName(j) for j in range(frs[i].NPar())]) for i in range(len(frs)))
    max_pars = 0
    max_index = None
    for i, p in pars.iteritems():
        if len(p) > max_pars:
            max_pars = len(p)
            max_index = i
    assert(i != None)
    max_pars = pars[max_index]

    for i, par in enumerate(max_pars):
        s += "{0:<5}".format(par)
        for j, fr in enumerate(frs):
            if par in pars[j]:
                s += "& ${0:>10.2e} \pm {1:>10.2e}$".format(fr.Value(i), fr.ParError(i))
            else:
                s += "& {0:>20}".format('-')
        s += '\\\\ \n'

    s += "{0:<5} ".format('$\chi^2$/\#DoF')
    for fr in frs:
        s += '& ${0:>10.3}$'.format(fr.Chi2() / fr.Ndf())
    s += '\\\\ \n'    
    s += "\hline\n\end{tabular}\n"
    return s

def draw_res_graph(res_graph, hist_events):
    res_max = TMath.MaxElement(res_graph.GetN(), res_graph.GetY()) * 1.10
    scale = res_max / hist_events.GetMaximum()
    hist_events = hist_events.Clone()
    __histos.append(hist_events)
    hist_events.Scale(scale)
    hist_events.GetYaxis().SetRangeUser(0, res_max * 1.10)
    res_graph.GetYaxis().SetRangeUser(0, res_max * 1.10)
    hist_events.Draw('hist') 
    from ROOT import kGray
    hist_events.SetFillColor(kGray + 1)
    res_graph.Draw('P')
    return hist_events
    
for key, fit_results in sorted(results.items(), key = lambda e: good[e[0].split('/')[-1]]):
    index = good[key.split('/')[-1]]
    full_sdata = sdatas[index]['sig_sdata']
    st = full_sdata.get().find('sigmat')
    st_cat = full_sdata.get().find('sigmat_cat')
    st_binning = st.getBinning('st_binning')
    split_bounds = array('d', [st_binning.binLow(0)] + [st_binning.binHigh(k) for k in range(st_binning.numBins())])
    
    name = 'canvas_%s' % index
    canvas = TCanvas(name, titles[key], 1000, 500)
    canvas.Divide(2, 1)
    __canvases.append(canvas)
    name = 'hist_events_%s' % index
    hist_events = TH1D(name, name, len(split_bounds) - 1, array('d', [v for v in split_bounds]))
    __histos.append(hist_events)
    mass_result = mass_fpf = fit_results['sWeight_mass_result']
    mass_fpf = mass_result.floatParsFinal()
    time_result = fit_results['time_result_%s' % fit_type]
    time_fpf = time_result.floatParsFinal()
    
    res_x = array('d')
    comb = array('d')
    comb_e = array('d')
    sf2s = array('d')
    sf2_es = array('d')
    total = full_sdata.sumEntries()
    for b, ct in enumerate(st_cat):
        d = split_bounds[b + 1] - split_bounds[b]
        bin_name = '_'.join(('N', sig_name, ct.GetName()))
        events = mass_fpf.find(bin_name)
        hist_events.SetBinContent(b + 1, events.getVal() / d)
        hist_events.SetBinError(b + 1, events.getError() / d)
        
        if fit_type.startswith('double_Comb'):
            sf_comb = time_fpf.find('timeResComb_%s' % ct.GetName())
            sf, sf_e = sf_comb.getVal(), sf_comb.getError()
            tmp = time_fpf.find('timeResSigmaSF_2_%s' % ct.GetName())
            sf2, sf2_e = tmp.getVal(), tmp.getError()
        elif fit_type == 'double':
            from P2VV.PropagateErrors import propagateScaleFactor
            sf, sf_e = propagateScaleFactor(time_result, '_' + ct.GetName())
        elif fit_type == 'single':
            sf_var = time_fpf.find('sigmaSF_%s' % ct.GetName())
            sf, sf_e = sf_var.getVal(), sf_var.getError()
        
        range_cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName())
        mean = full_sdata.mean(st, range_cut)
        res_x.append(mean)
        comb.append(sf)
        comb_e.append(sf_e)
        sf2s.append(sf2)
        sf2_es.append(sf2_e)
    
    res_ex = array('d', [0 for i in range(len(res_x))])
    res_graph = TGraphErrors(len(res_x), res_x, comb, res_ex, comb_e)
    res_graph.SetName('res_graph_%d' % index)
    sf2_graph = TGraphErrors(len(res_x), res_x, sf2s, res_ex, sf2_es)
    sf2_graph.SetName('sf2_graph_%d' % index)
    __histos.extend([res_graph, sf2_graph])
    
    from ROOT import TF1
    fit_funcs = {'pol1' : ('pol1', 'S0+'),
                 'pol2' : ('pol2', 'S0+'),
                 'pol1_mean_param' : ('[0] + [1] + [2] * (x - [0])', 'S+'),
                 'pol2_no_offset' : ('x ++ x * x', 'S0+'),
                 'pol2_mean_param' : ('[0] + [1] + [2] * (x - [0]) + [3] * (x - [0])^2', 'S0+')}
    print titles[key]
    st_mean = full_sdata.mean(st)
    for g in (res_graph, sf2_graph):
        frs = []
        for i, (name, (func, opts)) in enumerate(fit_funcs.iteritems()):
            fit_func = TF1(name, func, split_bounds[0], split_bounds[-1])
            if name.endswith('mean_param'):
                fit_func.FixParameter(0, st_mean)
            print name
            fit_result = g.Fit(fit_func, opts, "L")
            fit_result.SetName('result_' + name)
            print 'Chi2 / nDoF = %5.3f\n' % (fit_result.Chi2() / fit_result.Ndf())
            __fit_results[g.GetName().rsplit('_', 1)[0]].append(fit_result)
            frs.append(fit_result.Get())
        print fr_latex(frs)
    
    print ''
    
    canvas.cd(1)
    sf1_hist = draw_res_graph(res_graph, hist_events)
    sf1_hist.GetXaxis().SetTitle('estimated decay time resolution [ps]')
    sf1_hist.GetYaxis().SetTitle('combined scale factor')
    
    canvas.cd(2)
    sf2_hist = draw_res_graph(sf2_graph, hist_events)
    sf2_hist.GetXaxis().SetTitle('estimated decay time resolution [ps]')
    sf2_hist.GetYaxis().SetTitle('2nd scale factor')

## def chunks(l, n):
##     """ Yield successive n-sized chunks from l.
##     """
##     for i in range(0, len(l), n):
##         yield l[i:i+n]

## for frs in chunks(__fit_results, 2):
##     print fr_latex(frs)
##     print ''

## graphs = [g for g in __histos if g.GetName().find('graph') != -1]

## func = 'pol2'
## canvas = TCanvas('canvas', 'canvas', 500, 500)
## colors = [kRed, kGreen, kBlue, kOrange, kBlack]
## __extra_results = []
## __funcs = []
## for i, g in enumerate(graphs):
##     g.SetMarkerColor(colors[i])
##     g.SetLineColor(colors[i])
##     fit_func = TF1('fit_func_%s_%d' % (func, i) , func, split_bounds[0], split_bounds[-1])
##     __funcs.append(fit_func)
##     fit_func.SetLineColor(colors[i])
##     fit_result = g.Fit(fit_func, "S+", "L")
##     __extra_results.append(fit_result)
##     if i == 0:
##         g.Draw('AP')
##     else:
##         g.Draw('P')

## Make datasets for simultaneous fit
import ROOT
         
######################################################################

class SimFunction( ROOT.TPyMultiGenFunction ):
    def __init__(self, x, c, ce, s2, s2e, sim = True):
        ROOT.TPyMultiGenFunction.__init__(self, self)
        self.__x = x
        self.__c = c
        self.__ce = ce
        self.__s2 = s2
        self.__s2e = s2e
        self.__sim = sim
        
        self.__mean = sum(x) / float(len(x))
        
    def NDim(self):
        return 3 if self.__sim else 4
    
    def __fun(self, x, par):
        return par[0] + par[1] * (x - self.__mean)
    
    def DoEval(self, x):
        # calculate chisquare
        chisq, delta = 0., 0.
        if self.__sim:
            par_c = [x[0], x[1]]
            par_s2 = [x[0], x[2]]
        else:
            par_c = [x[0], x[1]]
            par_s2 = [x[2], x[3]]
        for z, ze, par in [(self.__c, self.__ce, par_c),
                           (self.__s2, self.__s2e, par_s2)]:
            for i in range(len(self.__x)):
                delta  = (z[i]-self.__fun(self.__x[i], par)) / ze[i]
                chisq += delta*delta
        
        return chisq

def sim_fit(x, c, ce, s2, s2e, sim = True):
    fitter = ROOT.Fit.Fitter()
    fitter.Config().MinimizerOptions().SetMinimizerType("Minuit2")
    fitter.Config().MinimizerOptions().SetMinimizerAlgorithm("Migrad")
    fitter.Config().MinimizerOptions().SetPrintLevel(2) 
    
    fcn = SimFunction(x, c, ce, s2, s2e, sim)
    
    start = array('d', (0, -3, -3))
    fitter.FitFCN(fcn, start)
    
    return fitter.Result()

## sim_result = sim_fit(res_x, comb, comb_e, sf2s, sf2_es, True)
## sep_result = sim_fit(res_x, comb, comb_e, sf2s, sf2_es, False)

from ROOT import TFile
sig_file = TFile.Open("/bfys/raaij/p2vv/data/P2VVDataSets_4KKMassBins_noTagCats.root")
sig_data = sig_file.Get("JpsiKK_splotdata_weighted_sigMass")
sig_st = sig_data.get().find("sigmat")
st_data = []
for i in range(sig_data.numEntries()):
    r = sig_data.get(i)
    st_data.append((sig_st.getVal(), sig_data.weight()))

from operator import itemgetter
st_data = sorted(st_data, key = itemgetter(0))
binning = st.getBinning('st_binning')

bins = []
for i in range(st_binning.numBins()):
    bins.append(st_binning.binLow(i))
bins.append(st_binning.binHigh(st_binning.numBins() - 1))

binned_data = []
it = iter(bins)
boundary = it.next()
bin_list = []
for v in st_data:
    try:
        if not bin_list and v[0] < boundary:
            continue
        if v[0] > boundary:
            if bin_list:
                binned_data.append(bin_list)
                bin_list = []
            boundary = it.next()
        bin_list.append(v)
    except StopIteration:
        break

def calib_sf1(pars, st):
    sfc_pars = pars[0]
    sf2_pars = pars[1]
    frac = pars[2]
    sfc = sum(sfc_pars[i + 1] * pow(st - sfc_pars[0], i) for i in range(len(sfc_pars) - 1))
    sf2 = calib_sf2(sf2_pars, st) / st
    sf1 = (sfc - frac * sf2) / (1 - frac)
    return sf1 * st

def calib_sf2(pars, st):
    return sum(pars[i + 1] * pow(st - pars[0], i) for i in range(len(pars) - 1)) * st

def sf1(p, st):
    return (p[0].value() - p[1].value() * p[2].value()) / (1 - p[1].value()) * st

from P2VV.PropagateErrors import Parameter, ErrorSFC, ErrorSG, ErrorCDG

def make_parameter(result, name):
    p = result.floatParsFinal().find(name)
    return Parameter(name, p.getVal(), p.getError())

dms = Parameter('dms', 17.768,  0.024)

from P2VV.Dilution import dilution
# Dilution of double Gauss with scalefactors calibrated.
cdg_rd = dirs['9bins_14.10fs_simul/m934737057402830078/results']
result_cdg = rd.Get("time_result_double_Comb_Gauss_linear")
sfc_offset = make_parameter(result_cdg, "sfc_offset")
sfc_slope = make_parameter(result_cdg, "sfc_slope")
cdg_frac = make_parameter(result_cdg, "timeResFrac2")
sf2_offset = make_parameter(result_cdg, "sf2_offset")
sf2_slope = make_parameter(result_cdg, "sf2_slope")
cdg_cv = result_cdg.reducedCovarianceMatrix(RooArgList(*[result_cdg.floatParsFinal().find(p.name()) for p in [sfc_offset, sfc_slope, cdg_frac, sf2_offset, sf2_slope]]))
error_cdg = ErrorCDG(st_mean, dms, sfc_offset, sfc_slope, cdg_frac, sf2_offset, sf2_slope, cdg_cv)
calib_dilutions = []
for bin_data in binned_data:
    d = dilution(bin_data, [([(st_mean, sfc_offset.value(), sfc_slope.value()),
                              (st_mean, sf2_offset.value(), sf2_slope.value()), cdg_frac.value()],
                              (1 - cdg_frac.value())), ((st_mean, sf2_offset.value(), sf2_slope.value()), cdg_frac.value())],
                              (calib_sf1, calib_sf2), error_cdg)
    calib_dilutions.append(d)

# Dilution of single Gauss fit.
sg_dilutions = []
sf_sg = Parameter('sf_sg', 1.45, 0.06)
sf_sge = ErrorSG(dms, sf_sg)
for bin_data in binned_data:
    d = dilution(bin_data, [((sf_sg,), 1)], (lambda p, st: p[0].value() * st,), sf_sge)
    sg_dilutions.append(d)

# Dilution of Double Gauss fit
sfc_dir = result_sfc = dirs['1bin_9500.00fs_simple/m934737057402830078/results']
result_sfc = sfc_dir.Get("time_result_double_Comb_Gauss")
sfc = make_parameter(result_sfc, "timeResComb")
dg_frac = make_parameter(result_sfc, "timeResFrac2")
sf2 = make_parameter(result_sfc, "timeResSigmaSF_2")
matrix = result_sfc.reducedCovarianceMatrix(RooArgList(*[result_sfc.floatParsFinal().find(p.name()) for p in [sfc, dg_frac, sf2]]))

dg_dilutions = []
error_dg = ErrorSFC(dms, sfc, dg_frac, sf2, matrix)
for bin_data in binned_data:
    d = dilution(bin_data, [((sfc, dg_frac, sf2), (1 - dg_frac.value())), ((0, sf2), dg_frac.value())], (sf1, lambda p, st: p[1].value() * st), error_dg)
    dg_dilutions.append(d)

# Plot Dilutions
means = array('d')
for bin_data in binned_data:
    means.append(sum(e[0] * e[1] for e in bin_data) / sum(e[1] for e in bin_data))
    
graphs = []
canvas = TCanvas('dilution_canvas', 'dilution_canvas', 550, 500)
canvas.SetLeftMargin(0.18)
colors = [kGreen, kBlue, kBlack]
first = True
for color, (ds, name) in zip(colors, [(calib_dilutions, 'calibrated'), (dg_dilutions, 'double'),
                                (sg_dilutions, 'single')]):
    graph = TGraphErrors(len(means), means, array('d', [d[0] for d in ds]),
                         array('d', len(ds) * [0]), array('d', [d[1] for d in ds]))
    graph.SetName(name)
    if first:
        graph.Draw("AP")
        graph.GetXaxis().SetTitle('estimated decay time error [ps]')
        graph.GetYaxis().SetTitle('dilution')
        graph.GetYaxis().SetTitleOffset(1.2)
        first = False
    else:
        graph.Draw("P, same")
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graphs.append(graph)
