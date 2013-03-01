import os
from array import array
from ROOT import *
from collections import defaultdict
from P2VV.Load import LHCbStyle
from P2VV.RooFitDecorators import *

f = TFile.Open('/stuff/PhD/p2vv/data/Bs2JpsiPhi_2011_Prescaled_st_bins.root')

sdatas = defaultdict(dict)
results = defaultdict(dict)

def add_keys(class_name, keys = None, path = None):
    if keys == None or type(keys) != dict:
        keys = {}
    if path == None or path == '/':
        d = gFile
        path = ''
    else:
        d = gFile.Get(path)
    for key in d.GetListOfKeys():
        if path:
            new_path = path + '/' + key.GetName()
        else:
            new_path = key.GetName()
        if new_path not in keys and key.GetClassName() == class_name:
            keys[new_path] = key.ReadObj()
        if key.GetClassName() == 'TDirectoryFile':
            add_keys(class_name, keys, new_path)
    return keys

dirs = add_keys('TDirectoryFile')
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

for k, d in sorted(interesting.items(), key = lambda e: int(e[0].split('bins')[0])):
    cut = d.Get('cut')
    print k, cut.String().Data().split('&&')[0]

good = {'1243785060103642893' : 4, 'm2334064025374600976' : 3,
        '1626518906014697943' : 2, 'm3832912631969227654' : 1}
fd = {}
PDFs = defaultdict(dict)
for k, cache_dir in filter(lambda k: k[0].split('bins')[0] in ['9', '10'], interesting.iteritems()):
    hd = k.split('/')[-1]
    try:
        index = good[hd]
    except KeyError:
        continue
    fd[index] = cache_dir.GetDirectory('..')
    sdata_dir = cache_dir.Get('sdata')
    for e in sdata_dir.GetListOfKeys():
        if e.GetClassName() == 'RooDataSet':
            sdatas[index][e.GetName()] = e.ReadObj()
    rd = cache_dir.Get('results')
    for e in rd.GetListOfKeys():
        if e.GetClassName() == 'RooFitResult':
            results[index][e.GetName()] = e.ReadObj()
    pdf_dir = cache_dir.Get('PDFs')
    for e in pdf_dir.GetListOfKeys():
        PDFs[index][e.GetName()] = e.ReadObj()

from ROOT import kGray
from ROOT import TH1D
from ROOT import TMatrixT
from math import sqrt
from P2VV.PropagateErrors import propagateScaleFactor

__canvases = []
__histos = []
__fit_funcs = []

fit_type = 'single'

__fit_results = []
from array import array



for i, fit_results in results.iteritems():
    full_sdata = sdatas[i]['sig_sdata']
    st = full_sdata.get().find('sigmat')
    st_cat = full_sdata.get().find('sigmat_cat')
    st_binning = st.getBinning('st_binning')
    split_bounds = array('d', [st_binning.binLow(0)] + [st_binning.binHigh(k) for k in range(st_binning.numBins())])
    
    name = 'canvas_%d' % i
    canvas = TCanvas(name, name, 500, 500)
    __canvases.append(canvas)
    name = 'hist_events_%d' % i
    hist_events = TH1D(name, name, len(split_bounds) - 1, array('d', [1000 * v for v in split_bounds]))
    __histos.append(hist_events)
    
    mass_fpf = fit_results['sWeight_mass_result'].floatParsFinal()
    time_fpf = fit_results['time_result_%s' % fit_type.split('_')[0]].floatParsFinal()
    
    res_x = array('d')
    res = array('d')
    res_e = array('d')
    for index, ct in enumerate(st_cat):
        d = split_bounds[index + 1] - split_bounds[index]
        bin_name = '_'.join(('N_psi_ll', ct.GetName()))
        events = mass_fpf.find(bin_name)
        hist_events.SetBinContent(index + 1, events.getVal() / d)
        hist_events.SetBinError(index + 1, events.getError() / d)
        
        if fit_type == 'double':
            from PropagateErrors import propagateScaleFactor
            sf, sf_e = propagateScaleFactor(time_result, '_' + ct.GetName())
        elif fit_type == 'double_Comb':
            sf_comb = time_fpf.find('timeResComb_%s' % ct.GetName())
            sf, sf_e = sf_comb.getVal(), sf_comb.getError()
        elif fit_type == 'single':
            sf_var = time_fpf.find('sigmaSF_%s' % ct.GetName())
            sf, sf_e = sf_var.getVal(), sf_var.getError()
        
        mean = sdatas[i][ct.GetName()].mean(st)
        res_x.append(1000 * mean)
        res.append(1000 * mean * sf)
        res_e.append(1000 * mean * sf_e)
    
    res_ex = array('d', [0 for i in range(len(res_x))])
    res_graph = TGraphErrors(len(res_x), res_x, res, res_ex, res_e)
    __histos.append(res_graph)
    scale = 100 / hist_events.GetMaximum()
    hist_events.Scale(scale)
    hist_events.GetYaxis().SetRangeUser(0, 110)
    
    from ROOT import TF1
    fit_funcs = {'pol1' : 'S0+', 'pol2' : 'S+'}
    for i, (func, opts) in enumerate(fit_funcs.iteritems()):
        fit_func = TF1('fit_func_%s' % func , func, split_bounds[0], split_bounds[-1])
        fit_result = res_graph.Fit(fit_func, opts, "L")
        __fit_results.append(fit_result)
        ## fr = fit_result.Get()
        ## fr.SetName('fit_result_%s_%s' % (func, args[1]))
        ## results.append(fr)
    
    res_graph.GetYaxis().SetRangeUser(0, 110)
    hist_events.Draw('hist')
    from ROOT import kGray
    hist_events.SetFillColor(kGray + 1)
    res_graph.Draw('P')
    res_graph.GetXaxis().SetTitle('estimated decay time error [fs]')
    res_graph.GetYaxis().SetTitle('decay time resulution [fs]')

from ROOT import RooBinning
plot_bounds = array('d', [-1.5 + i * 0.1 for i in range(12)] + [-0.3 + i * 0.01 for i in range(60)] + [0.3 + i * 0.1 for i in range(57)] + [6 + i * 0.4 for i in range(6)])
zoom_bounds = array('d', [-0.2 + i * 0.005 for i in range(81)])

plot_binning = RooBinning(len(plot_bounds) - 1, plot_bounds)
plot_binning.SetName('full')
zoom_binning = RooBinning(len(zoom_bounds) - 1, zoom_bounds)
zoom_binning.SetName('zoom')

from ROOT import kDashed, kRed, kGreen, kBlue, kBlack, kOrange
from P2VV.GeneralUtils import plot

print 'plotting'


binnings = [plot_binning, zoom_binning]
plotLog = [True, False]

plots = defaultdict(list)
assert(False)

for i, fit_results in results.iteritems():
    full_sdata = sdatas[i]['sig_sdata']
    st = full_sdata.get().find('sigmat')
    st_cat = full_sdata.get().find('sigmat_cat')
    st_binning = st.getBinning('st_binning')
    split_bounds = array('d', [st_binning.binLow(0)] + [st_binning.binHigh(k) for k in range(st_binning.numBins())])
    time = full_sdata.get().find('time')
    
    time_result = None
    for r, result in fit_results.iteritems():
        if r == 'time_result_' + fit_type:
            time_result = result
            break
    assert(time_result)
    
    time_pdf = PDFs[i]['time_pdf']
    parameters = time_pdf.getParameters(full_sdata)
    
    pdf_set = set([par.GetName() for par in parameters if not par.isConstant()])
    result_set = set([par.GetName() for par in time_result.floatParsFinal()])
    assert(pdf_set == result_set)
    
    for par in time_result.floatParsFinal():
        pdf_par = parameters.find(par.GetName())
        if not pdf_par:
            print par.GetName()
            continue
        pdf_par.setVal(par.getVal())
        pdf_par.setError(par.getError())
    
    time_pdf.recursiveRedirectServers(RooArgSet(st, st_cat))
    
    for l, (bins, pl) in enumerate(zip(binnings, plotLog)):
        if fd[i].GetName().find('simul') != -1:
            projSet = RooArgSet(st, st_cat)
            r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
            for ct in st_cat:
                time_pdf.indexCat().setLabel(ct.GetName())
                name = 'time_canvas_%s_%d' % (ct.GetName(), l)
                canvas = TCanvas(name, name, 600, 400)
                __canvases.append(canvas)
                p = canvas.cd(1)
                pd = sdatas[i][ct.GetName()]
                pdfOpts  = dict(ProjWData = (projSet, full_sdata, True))
                ps = plot(p, time, pdf = time_pdf, data = pd
                          , frameOpts = dict(Range = r, Title = "")
                          , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                          , pdfOpts  = dict(LineWidth = 4, Slice = (st_cat, ct.GetName()), **pdfOpts)
                          , logy = pl
                          , plotResidHist = False
                          ## , components = { 'wpv_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                          ##                  , 'prompt_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                          ##                  , 'sig_*'     : dict( LineColor = kOrange,  LineStyle = kDashed )
                          ##                  }
                          )
                for frame in ps:
                    plot_name = '_'.join((time.GetName(), bins.GetName(), ct.GetName(), frame.GetName()))
                    frame.SetName(plot_name)
                plots[i].append(ps)
        else:
            canvas = TCanvas('time_canvas_%d' % l, 'time_canvas_%d' % l, 600, 400)
            __canvases.append(canvas)
            p = canvas.cd(1)
            r = (bins.binLow(0), bins.binHigh(bins.numBins() - 1))
            projSet = RooArgSet(st)
            pdfOpts  = dict(ProjWData = (projSet, full_sdata, True))
            ps = plot(p, t, pdf = time_pdf, data = pd
                      , frameOpts = dict(Range = r, Title = "")
                      , dataOpts = dict(MarkerSize = 0.8, Binning = bins, MarkerColor = kBlack)
                      , pdfOpts  = dict(LineWidth = 4, **pdfOpts)
                      , logy = pl
                      , plotResidHist = False)
            for frame in ps:
                plot_name = '_'.join((t.GetName(), bins.GetName(), frame.GetName()))
                frame.SetName(plot_name)
            plots[i].append(ps)
