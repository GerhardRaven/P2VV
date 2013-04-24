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
elif args[0] not in ['2011', 'MC11a']:
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
for k, d in sorted(interesting.items(), key = lambda e: int(e[0].split('bins')[0])):
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
elif args[0] == 'MC11a':
    good = {'389085267962218368' : 4, 'm3019457528953402347' : 3,
            'm7780668933605436626' : 2, 'm8376372569899625413' : 1,
            'm1545059518337894505' : 6, 'm8342219958663192955' : 5}
    sig_name = 'signal'

PDFs = defaultdict(dict)
for k, cache_dir in filter(lambda k: k[0].split('bins')[0] in ['9', '10'], interesting.iteritems()):
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

fit_type = 'double_Comb'

__fit_results = []
from array import array

def fr_latex(frs):
    frs = sorted(frs, key = lambda fr: fr.NPar())
    names = ' & '.join(['result ' + fr.GetName().rsplit('_', 1)[1] for fr in frs])
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

for key, fit_results in sorted(results.items(), key = lambda e: good[e[0].split('/')[-1]]):
    index = good[key.split('/')[-1]]
    full_sdata = sdatas[index]['sig_sdata']
    st = full_sdata.get().find('sigmat')
    st_cat = full_sdata.get().find('sigmat_cat')
    st_binning = st.getBinning('st_binning')
    split_bounds = array('d', [st_binning.binLow(0)] + [st_binning.binHigh(k) for k in range(st_binning.numBins())])
    
    name = 'canvas_%s' % index
    canvas = TCanvas(name, titles[key], 500, 500)
    __canvases.append(canvas)
    name = 'hist_events_%s' % index
    hist_events = TH1D(name, name, len(split_bounds) - 1, array('d', [v for v in split_bounds]))
    __histos.append(hist_events)
    mass_result = mass_fpf = fit_results['sWeight_mass_result']
    mass_fpf = mass_result.floatParsFinal()
    time_result = fit_results['time_result_%s' % fit_type.split('_')[0]]
    time_fpf = time_result.floatParsFinal()
    
    res_x = array('d')
    res = array('d')
    res_e = array('d')
    total = full_sdata.sumEntries()
    for b, ct in enumerate(st_cat):
        d = split_bounds[b + 1] - split_bounds[b]
        bin_name = '_'.join(('N', sig_name, ct.GetName()))
        events = mass_fpf.find(bin_name)
        hist_events.SetBinContent(b + 1, events.getVal() / d)
        hist_events.SetBinError(b + 1, events.getError() / d)
        
        if fit_type == 'double':
            from P2VV.PropagateErrors import propagateScaleFactor
            sf, sf_e = propagateScaleFactor(time_result, '_' + ct.GetName())
        elif fit_type == 'double_Comb':
            sf_comb = time_fpf.find('timeResComb_%s' % ct.GetName())
            sf, sf_e = sf_comb.getVal(), sf_comb.getError()
        elif fit_type == 'single':
            sf_var = time_fpf.find('sigmaSF_%s' % ct.GetName())
            sf, sf_e = sf_var.getVal(), sf_var.getError()
        
        range_cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName())
        mean = full_sdata.mean(st, range_cut)
        mean *= total / full_sdata.sumEntries(range_cut)
        res_x.append(mean)
        res.append(mean * sf)
        res_e.append(mean * sf_e)
    
    res_ex = array('d', [0 for i in range(len(res_x))])
    res_graph = TGraphErrors(len(res_x), res_x, res, res_ex, res_e)
    res_graph.SetName('res_graph_%d' % index)
    __histos.append(res_graph)
    scale = 0.1 / hist_events.GetMaximum()
    hist_events.Scale(scale)
    hist_events.GetYaxis().SetRangeUser(0, 0.11)
    
    from ROOT import TF1
    fit_funcs = {'pol1' : ('pol1', 'S0+'), 'pol2' : ('pol2', 'S+'),
                 'x ++ x * x' : ('pol2_no_offset', 'S0+'),
                 '[0] + [1] + [2] * (x - [0]) + [3] * (x - [0])^2' : ('pol2_mean_param', 'S0+')}
    print titles[key]
    frs = []
    for i, (func, opts) in enumerate(fit_funcs.iteritems()):
        fit_func = TF1('fit_func_%s' % opts[0], func, split_bins[0], split_bins[-1])
        if opts[0] == 'pol2_mean_param':
            fit_func.FixParameter(0, full_sdata.mean(st))
        fit_result = res_graph.Fit(fit_func, opts[1], "L")
        print 'Chi2 / nDoF = %5.3f' % (fit_result.Chi2() / fit_result.Ndf())
        __fit_results.append(fit_result)
        frs.append(fit_result.Get())
    print fr_latex(frs)
    
    print ''
    res_graph.GetYaxis().SetRangeUser(0, 0.11)
    hist_events.Draw('hist') 
    hist_events.GetXaxis().SetTitle('estimated decay time error [ps]')
    hist_events.GetYaxis().SetTitle('decay time resulution [ps]')
    from ROOT import kGray
    hist_events.SetFillColor(kGray + 1)
    res_graph.Draw('P')


## def chunks(l, n):
##     """ Yield successive n-sized chunks from l.
##     """
##     for i in range(0, len(l), n):
##         yield l[i:i+n]

## for frs in chunks(__fit_results, 2):
##     print fr_latex(frs)
##     print ''

graphs = [g for g in __histos if g.GetName().find('graph') != -1]

func = 'pol2'
canvas = TCanvas('canvas', 'canvas', 500, 500)
colors = [kRed, kGreen, kBlue, kOrange, kBlack]
__extra_results = []
__funcs = []
for i, g in enumerate(graphs):
    g.SetMarkerColor(colors[i])
    g.SetLineColor(colors[i])
    fit_func = TF1('fit_func_%s_%d' % (func, i) , func, split_bounds[0], split_bounds[-1])
    __funcs.append(fit_func)
    fit_func.SetLineColor(colors[i])
    fit_result = g.Fit(fit_func, "S+", "L")
    __extra_results.append(fit_result)
    if i == 0:
        g.Draw('AP')
    else:
        g.Draw('P')
