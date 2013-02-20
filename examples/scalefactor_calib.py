import os
from array import array
from ROOT import *
from collections import defaultdict
from P2VVLoad import LHCbStyle

f = TFile.Open('/stuff/PhD/p2vv/data/Bs2JpsiPhi_2011_Prescaled_st_bins.root')

order = {1 : '8431586533064055611', 2 : 'm5253414056674419702' , 3 : '8922847604268075573', 4 : 'm2133574294912993512'}

sdatas = defaultdict(list)
datas = defaultdict(list)
results = {}

st_bins = array('d', [0.01 + i * 0.01 for i in range(7)])
directory = '%sbins_%4.2ffs' % (len(st_bins) - 1, (1000 * (st_bins[1] - st_bins[0])))
for i, cd in order.iteritems():
    cache_dir = f.Get(os.path.join(directory, cd))
    for j in range(len(st_bins) - 1):
        sds_name = 'sdata/DecayTree_%02d_weighted_psi_background' % j
        sdatas[i].append(cache_dir.Get(sds_name))
        ds_name = 'data/DecayTree_%02d' % j
        datas[i].append(cache_dir.Get(ds_name))
    results[i] = defaultdict(list)
    for t in ['time', 'mass']:
        rs = set()
        rd = cache_dir.Get('%s_results' % t)
        for e in rd.GetListOfKeys():
            if e.GetClassName() == 'RooFitResult':
                rs.add(os.path.join(rd.GetName(), e.GetName()))
        for e in rs:
            results[i][t].append(cache_dir.Get(e))
        results[i][t] = sorted(results[i][t], key = lambda r: int(r.GetName().split('_', 1)[0]))


from ROOT import kGray
from ROOT import TH1D
from ROOT import TMatrixT
from math import sqrt
from PropagateErrors import propagateScaleFactor

__canvases = []
__histos = []
__fit_funcs = []

for i, res in results.iteritems():
    name = 'canvas_%d' % i
    canvas = TCanvas(name, name, 500, 500)
    __canvases.append(canvas)
    name = 'hist_res_%d' % i
    hist_res = TH1D(name, name, len(st_bins) - 1, array('d', [1000 * v for v in st_bins]))
    __histos.append(hist_res)
    name = 'hist_events_%d' % i
    hist_events = TH1D(name, name, len(st_bins) - 1, array('d', [1000 * v for v in st_bins]))
    __histos.append(hist_events)

    for index, rs in enumerate(res['time']):
        # Fill events histo
        events = res['mass'][index].floatParsFinal().find('N_psi_background')
        hist_events.SetBinContent(index + 1, events.getVal())
        hist_events.SetBinError(index + 1, events.getError())
    
        if rs.status() != 0:
            hist_res.SetBinContent(index + 1, 0)
            hist_res.SetBinError(index + 1, 0)
            continue

        fpf = rs.floatParsFinal()
        if fpf.find('timeResSigma2'):
            sf, sf_e = propagateScaleFactor(rs)
        else:
            sf = fpf.find('sigmaSF').getVal()
            sf_e = fpf.find('sigmaSF').getError()

        st = sdatas[i][index].get().find('sigmat')
        mean = sdatas[i][index].mean(st)
        r = mean * sf
        r_e = mean * sf_e
        
        hist_res.SetBinContent(index + 1, 1000 * r)
        hist_res.SetBinError(index + 1, 1000 * r_e)

    scale = 100 / hist_events.GetMaximum()
    hist_events.Scale(scale)
    hist_res.Draw('pe')
    hist_events.Draw('hist, same')
    hist_events.SetFillColor(kGray + 1)
    hist_res.Draw('pe, same')
    hist_res.GetXaxis().SetTitle('estimated decay time error [fs]')
    hist_res.GetYaxis().SetTitle('decay time resulution [fs]')
    hist_res.GetYaxis().SetRangeUser(0, 120)

    name = 'fit_func_%d' % i
    fit_func = TF1(name, "pol1", split_bins[0], split_bins[-1])
    __fit_funcs.append(fit_func)
    fr = hist_res.Fit(fit_func, "S")
    results[i]['resolution'].append(fr)

def add_keys(keys, path):
    print path
    if not path:
        d = gFile
    else:
        d = gFile.Get(path)
    for key in d.GetListOfKeys():
        if path:
            new_path = path + '/' + key.GetName()
        else:
            new_path = key.GetName()
        if key.GetClassName() == 'TDirectoryFile':
            add_keys(keys, new_path)
        elif key.GetClassName() == 'RooFitResult':
            keys.append(new_path)
