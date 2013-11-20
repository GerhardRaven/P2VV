#!/usr/bin/env python
from array import array
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog data_type')
(options, args) = parser.parse_args()

from P2VV.Utilities.Resolution import input_data, prefix

if len(args) != 1:
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
elif args[0] == 'MC2012':
    good = {'m4074732224475561764' : 1}
    sig_name = 'signal'
elif args[0] == 'MC2011_Sim08a':
    good = {'m4074732224475561764' : 1}
    sig_name = 'signal'
elif args[0] == 'MC2011_Sim08a_incl_Jpsi':
    good = {'2027465761870101697' : 1}
    sig_name = 'psi_ll'

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

fit_type = 'double_RMS_Gauss'

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
    
ffs = defaultdict(dict)
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

    split_mean = time_fpf.find('timeResMu_st_bin_0')

    res_x = array('d')
    comb = array('d')
    comb_e = array('d')
    sfos = array('d')
    means = array('d')
    mean_es = array('d')
    sfo_es = array('d')
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
            sfo, sfo_e = tmp.getVal(), tmp.getError()
        elif fit_type.startswith('double_RMS'):
            sf_av = time_fpf.find('timeResSFMean_%s' % ct.GetName())
            sf, sf_e = sf_av.getVal(), sf_av.getError()
            sf_sigma = time_fpf.find('timeResSFSigma_%s' % ct.GetName())
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

        range_cut = '{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName())
        mean = full_sdata.mean(st, range_cut)
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
    fit_funcs = {'pol1' : ('pol1', 'S0+'),
                 'pol2' : ('pol2', 'S0+'),
                 'pol1_mean_param' : ('[1] + [2] * ((x - [0]) / [3])', 'S0+'),
                 'pol2_no_offset' : ('x ++ x * x', 'S0+'),
                 'pol2_mean_param' : ('[1] + [2] * ((x - [0]) / [3]) + [4] * ((x - [0]) / [3])^2', 'S0+')}
    print titles[key]
    st_mean = full_sdata.mean(st)
    graphs = [res_graph, sfo_graph]
    if split_mean:
        graphs.append(mean_graph)
    for g in graphs:
        frs = []
        for i, (name, (func, opts)) in enumerate(fit_funcs.iteritems()):
            fit_func = TF1(name, func, split_bounds[0], split_bounds[-1])
            if name.endswith('mean_param'):
                fit_func.FixParameter(0, st_mean)
                fit_func.FixParameter(3, res_x[-1] - res_x[0])
            print name
            fit_result = g.Fit(fit_func, opts, "L")
            fit_result.SetName('result_' + name)
            print 'Chi2 / nDoF = %5.3f\n' % (fit_result.Chi2() / fit_result.Ndf())
            __fit_results[g.GetName().rsplit('_', 1)[0]].append(fit_result)
            frs.append(fit_result.Get())
            ffs[g.GetName()][name] = fit_func
        print fr_latex(frs)
    
    print ''

    def draw_calib_graph(prefix):
        calib_result = fit_results[time_result.GetName() + '_linear']
        if not calib_result:
            return
        calib_fpf = dict([(p.GetName(), [p.getVal(), p.getError()]) for p in calib_result.floatParsFinal()])
        offset = calib_fpf['_'.join((prefix, 'offset'))]
        slope = calib_fpf['_'.join((prefix, 'slope'))]
        slope[0] /= (split_bounds[-1] - split_bounds[0])
        calib_func = TF1(prefix + '_simul', '[0] + [1] * (x - [2])',
                         (split_bounds[0] + split_bounds[1]) / 2.,
                         (split_bounds[-2] + split_bounds[-1]) / 2.)
        for i, (v, e) in [(0, offset), (1, slope), (2, (st_mean, 0))]:
            calib_func.SetParameter(i, v)
            calib_func.SetParError(i, e)
        from ROOT import kBlue
        calib_func.SetLineColor(kBlue)
        calib_func.Draw('same')
        __histos.append(calib_func)
    
    canvas.cd(1)
    sf1_hist = draw_res_graph(res_graph, hist_events)
    sf1_hist.GetYaxis().SetTitle('#bar{sf}')
    sf1_hist.GetYaxis().SetTitleOffset(1.05)
    draw_calib_graph('sf_mean')
    
    canvas.cd(2)
    sfo_hist = draw_res_graph(sfo_graph, hist_events)
    sfo_hist.GetYaxis().SetTitle('sf_{#sigma}')
    sfo_hist.GetYaxis().SetTitleOffset(1.05)
    draw_calib_graph('sf_sigma')

    if split_mean:
        name = 'mean_canvas_%s' % index
        mean_canvas = TCanvas(name, titles[key], 600, 400)
        mean_canvas.SetLeftMargin(0.15)
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
        for n, c in [('pol1_mean_param', kGreen), ('pol2_mean_param', kBlue)]:
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
