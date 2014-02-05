#!/usr/bin/env python
from array import array
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog data_type top_dir sub_dir')
(options, args) = parser.parse_args()

from P2VV.Utilities.Resolution import input_data, prefix

if len(args) not in [1, 3]:
    print parser.print_usage()
    sys.exit(-2)
elif args[0] not in input_data.keys():
    print parser.print_usage()
    print "Possible samples are: %s" % ' '.join(input_data.keys())
    sys.exit(-2)

from ROOT import *
from collections import defaultdict
from P2VV.Load import LHCbStyle
from P2VV.RooFitDecorators import *

directory = os.path.join(prefix, 'p2vv/data')

from P2VV.CacheUtils import CacheFiles
cfs = CacheFiles(*input_data[args[0]]['cache'].rsplit('/', 1))
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
    
if len(args) == 1:
    print 'possible top level directories are:'
    for k in interesting.iterkeys():
        sub_dirs = filter(lambda x: len(x.split('/')) == 3 and x.startswith(k) and x.endswith('bins'), dirs.keys())
        print k + ':', ' '.join([s.split('/')[-1] for s in sub_dirs])
    sys.exit(0)


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

top_dir = args[1]
sub_dir = args[2]
if args[0] in ['MC2011_Sim08a', 'MC2012']:
    single_dir = '1bin_19000.00fs_simple'
else:
    single_dir = '1bin_15500.00fs_simple' 

interesting.update(dict([(k, dirs[k]) for k in set([k.replace(top_dir, single_dir)
                                                    for k in interesting.iterkeys()]) if k in dirs]))

if args[0] == '2011':
    ## good = {'1243785060103642893' : 4, 'm2334064025374600976' : 3,
    ##         '1626518906014697943' : 2, 'm3832912631969227654' : 1,
    ##         '4086600821164745518' : 6, 'm6573713017788044320' : 5}
    good = {'2027465761870101697' : 1}
    sig_name = 'psi_ll'
elif args[0] == '2012':
    good = {'2027465761870101697' : 1}
    ## 'm934737057402830078'
    sig_name = 'psi_ll'
elif args[0] == '2011_Reco14':
    good = {'2027465761870101697' : 1}
    sig_name = 'psi_ll'
elif args[0] == 'MC11a':
    good = {'389085267962218368' : 4, 'm3019457528953402347' : 3,
            'm7780668933605436626' : 2, 'm8376372569899625413' : 1,
            'm1545059518337894505' : 6, 'm8342219958663192955' : 5}
    sig_name = 'signal'
elif args[0] == 'MC2012':
    good = {'m8673867318518908571' : 1}
    sig_name = 'signal'
elif args[0] == 'MC2011_Sim08a':
    good = {'m8673867318518908571' : 1}
    sig_name = 'signal'
elif args[0] == 'MC2011_Sim08a_incl_Jpsi':
    good = {'2027465761870101697' : 1}
    sig_name = 'psi_ll'
elif args[0] == 'MC2012_incl_Jpsi':
    good = {'2027465761870101697' : 1}
    sig_name = 'psi_ll'

from itertools import chain
PDFs = defaultdict(dict)
for k, cache_dir in filter(lambda k: k[0].split('/')[0] in [top_dir, single_dir], interesting.iteritems()):
    hd = k.split('/')[-1]
    try:
        index = good[hd]
    except KeyError:
        continue

    def add_sdata(d, name):
        for e in d.GetListOfKeys():
            if e.GetClassName() == 'RooDataSet':
                sdatas[index]['/'.join((name, e.GetName())) if name else e.GetName()] = e.ReadObj()

    def add_results(d):
        for e in d.GetListOfKeys():
            if e.GetClassName() == 'RooFitResult':
                results[k][e.GetName()] = e.ReadObj()

    add_sdata(cache_dir.Get('sdata'), '')
    add_results(cache_dir.Get('results'))
    sub_cache_dir = cache_dir.Get(sub_dir)
    if sub_cache_dir:
        add_sdata(sub_cache_dir.Get('sdata'), sub_dir)
        add_results(sub_cache_dir.Get('results'))

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

fit_type = 'double_RMS_Gauss'
if args[0] in ['MC2012', 'MC2011_Sim08a']:
    fit_type = 'double_RMS_Rest'

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
    from ROOT import kGray
    hist_events.SetFillColor(kGray + 1)
    hist_events.Draw('hist') 
    res_graph.Draw('P')
    hist_events.GetXaxis().SetTitle('estimated decay time resolution [ps]')
    return hist_events

if top_dir.split('_')[0].startswith('momentum'):
    split_cat_name = 'B_P_cat'
    binning_name = 'momentum_binning'
    obs_name = 'B_P'
else:
    split_cat_name = 'sigmat_cat'
    binning_name = 'st_binning'
    obs_name = 'sigmat'
    
ffs = defaultdict(dict)
for key, fit_results in sorted(results.items(), key = lambda e: good[e[0].split('/')[-1]]):
    if top_dir not in key:
        continue
    index = good[key.split('/')[-1]]
    full_sdata = sdatas[index][sub_dir + '/sig_sdata']
    observable = full_sdata.get().find(obs_name)
    split_cat = full_sdata.get().find(split_cat_name)
    binning = observable.getBinning(binning_name)
    split_bounds = array('d', [binning.binLow(0)] + [binning.binHigh(k) for k in range(binning.numBins())])
    
    name = 'canvas_%s' % index
    canvas = TCanvas(name, titles[key], 1200, 400)
    
    canvas.Divide(2, 1)
    __canvases.append(canvas)
    name = 'hist_events_%s' % index
    hist_events = TH1D(name, name, len(split_bounds) - 1, array('d', [v for v in split_bounds]))
    __histos.append(hist_events)
    mass_result = mass_fpf = fit_results['sWeight_mass_result']
    mass_fpf = mass_result.floatParsFinal()
    for suffix in (fit_type, fit_type + '_pee'):
        if 'time_result_' + suffix in fit_results:
            time_result = fit_results['time_result_' + suffix]
            break
    else:
        continue
    time_fpf = time_result.floatParsFinal()

    split_mean = time_fpf.find('timeResMu_st_bin_0')

    res_x = array('d')
    comb = array('d')
    comb_e = array('d')
    sfos = array('d')
    means = array('d')
    mean_es = array('d')
    sfo_es = array('d')
    total = full_sdata.sumEntries()
    for b, ct in enumerate(split_cat):
        d = split_bounds[b + 1] - split_bounds[b]
        bin_name = '_'.join(('N', sig_name, ct.GetName()))
        events = mass_fpf.find(bin_name)
        hist_events.SetBinContent(b + 1, events.getVal() / d)
        hist_events.SetBinError(b + 1, events.getError() / d)
        
        if fit_type.startswith('double_Comb'):
            sf_comb = time_fpf.find('timeResComb_%s' % ct.GetName())
            sf, sf_e = sf_comb.getVal(), sf_comb.getError()
            tmp = time_fpf.find('timeResSigmaSF_2_%s' % ct.GetName())
            sfo, sfo_e = tmp.getVal(), tmp.getError()
        elif fit_type.startswith('double_RMS'):
            sf_av = time_fpf.find('timeResSFMean_%s' % ct.GetName())
            sf_sigma = time_fpf.find('timeResSFSigma_%s' % ct.GetName())
            sf, sf_e = sf_av.getVal(), sf_av.getError()
            sfo, sfo_e = sf_sigma.getVal(), sf_sigma.getError()
        elif fit_type == 'double':
            from P2VV.PropagateErrors import propagateScaleFactor
            sf, sf_e = propagateScaleFactor(time_result, '_' + ct.GetName())
        elif fit_type == 'single':
            sf_var = time_fpf.find('sigmaSF_%s' % ct.GetName())
            sf, sf_e = sf_var.getVal(), sf_var.getError()
        
        if split_mean:
            mean_var = time_fpf.find('timeResMu_%s' % ct.GetName())
            means.append(mean_var.getVal())
            mean_es.append(mean_var.getError())

        range_cut = '{0} == {0}::{1}'.format(split_cat.GetName(), ct.GetName())
        mean = full_sdata.mean(observable, range_cut)
        res_x.append(mean)
        comb.append(sf)
        comb_e.append(sf_e)
        sfos.append(sfo)
        sfo_es.append(sfo_e)
    
    res_ex = array('d', [0 for i in range(len(res_x))])
    res_graph = TGraphErrors(len(res_x), res_x, comb, res_ex, comb_e)
    res_graph.SetName('res_graph_%d' % index)
    sfo_graph = TGraphErrors(len(res_x), res_x, sfos, res_ex, sfo_es)
    sfo_graph.SetName('sfo_graph_%d' % index)
    __histos.extend([res_graph, sfo_graph])

    if split_mean:
        mean_graph = TGraphErrors(len(res_x), res_x, means, res_ex, mean_es)
        mean_graph.SetName('mean_graph_%d' % index)
        __histos.append(mean_graph)

    from ROOT import TF1
    fit_funcs = {'pol0' : ('pol0', 'S0+'),
                 'pol1' : ('pol1', 'S0+'),
                 'pol2' : ('pol2', 'S0+'),
                 'pol1_mean_param' : ('[1] + [2] * (x - [0])', 'S0+'),
                 'pol2_mean_param_no_offset' : ('([1] + [2] * (x - [0])) * x', 'S0+'),
                 'pol2_mean_param' : ('[1] + [2] * (x - [0]) + [3] * (x - [0])^2', 'S0+')}
    print titles[key]
    mean = full_sdata.mean(observable)
    graphs = [res_graph, sfo_graph]
    if split_mean:
        graphs.append(mean_graph)
    for g in graphs:
        print ('-' * 15) + g.GetName() + ('-' * 15)
        frs = []
        for i, (name, (func, opts)) in enumerate(fit_funcs.iteritems()):
            fit_func = TF1(name, func, split_bounds[0], split_bounds[-1])
            if 'mean_param' in name:
                fit_func.FixParameter(0, mean)
            print name
            fit_result = g.Fit(fit_func, opts, "L")
            fit_result.SetName('result_' + name)
            print 'Chi2 / nDoF = %5.3f\n' % (fit_result.Chi2() / fit_result.Ndf())
            __fit_results[g.GetName().rsplit('_', 1)[0]].append(fit_result)
            frs.append(fit_result.Get())
            ffs[g.GetName()][name] = fit_func
        print fr_latex(frs)
    
    print ''

    def draw_calib_graph(prefix, calib_type, formula, pars, color):
        if key.replace(top_dir, single_dir) not in results:
            return
        ## Try extra _pee suffix for backwards compatibility
        for rn in [(calib_type,), (calib_type, 'pee')]:
            time_name = time_result.GetName().replace('_pee', '')
            calib_result = results[key.replace(top_dir, single_dir)].get('_'.join((time_name,) + rn), None)
            if calib_result:
                break
        else:
            return
        calib_fpf = dict([(p.GetName(), [p.getVal(), p.getError()]) for p in calib_result.floatParsFinal()])
        calib_func = TF1(prefix + '_simul', formula,
                         (split_bounds[0] + split_bounds[1]) / 2.,
                         (split_bounds[-2] + split_bounds[-1]) / 2.)
        for i, par_name in enumerate(pars):
            par = calib_fpf[prefix + '_' + par_name]
            calib_func.SetParameter(i + 1, par[0])
            calib_func.SetParError(i, par[1])
        calib_func.SetParameter(0, mean)
        calib_func.SetLineColor(color)
        calib_func.Draw('same')
        __histos.append(calib_func)
    
    canvas.cd(1)
    sf1_hist = draw_res_graph(res_graph, hist_events)
    sf1_hist.GetYaxis().SetTitle('#bar{sf}')
    sf1_hist.GetYaxis().SetTitleOffset(1.05)
    
    canvas.cd(2)
    sfo_hist = draw_res_graph(sfo_graph, hist_events)
    sfo_hist.GetYaxis().SetTitle('sf_{#sigma}')
    sfo_hist.GetYaxis().SetTitleOffset(1.05)

    from itertools import product
    par_defs = [('linear', '[1] + [2] * (x - [0])', ('offset', 'slope'), kBlue),
                ('quadratic', '[1] + [2] * (x - [0]) + [3] * (x - [0])^2', ('offset', 'slope', 'quad'), kGreen),
                ('quadratic_no_offset', '([1] + [2] * (x - [0])) * x', ('slope', 'quad'), kOrange)]
    for (pad, prefix), args in product(((1, 'sf_mean'), (2, 'sf_sigma')), par_defs):
        canvas.cd(pad)
        draw_calib_graph(*tuple((prefix, ) + args))

    if split_mean:
        name = 'mean_canvas_%s' % index
        mean_canvas = TCanvas(name, titles[key], 600, 400)
        mean_canvas.SetLeftMargin(0.15)
        __canvases.append(mean_canvas)
        mean_graph.Draw("AP")
        mean_graph.GetXaxis().SetTitle("estimated decay time error [ps]")
        mean_graph.GetYaxis().SetTitle("#mu [ps]  ")
        mean_graph.GetYaxis().SetTitleOffset(1.07)
        def __dg(n, c):
            g = ffs[mean_graph.GetName()][n]
            g.SetLineColor(c)
            g.SetRange((split_bounds[0] + split_bounds[1]) / 2.,
                       (split_bounds[-2] + split_bounds[-1]) / 2.)
            g.Draw('same')
        for n, c in [('pol1_mean_param', kBlue), ('pol2_mean_param', kGreen), ('pol0', 7)]:
            __dg(n, c)
        mean_canvas.Update()
    
    canvas.Update()

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
