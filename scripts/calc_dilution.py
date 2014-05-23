#!/usr/bin/env python
from array import array
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog data_type')
(options, args) = parser.parse_args()

prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/project/bfys/raaij'
input_data = {'2011' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'),
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled_cache.root'),
                        'result_dir' : '1bin_15500.00fs_simple/2027465761870101697/results',
                        'dataset' : 'JpsiKK_sigSWeight'},
              '2011_Reco14' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'),
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Reco14_Prescaled_cache.root'),
                               'result_dir' : '1bin_15500.00fs_simple/2027465761870101697/results',
                               'dataset' : 'JpsiKK_sigSWeight'},
              '2012' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'),
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_Prescaled_cache.root'),
                               'result_dir' : '1bin_15500.00fs_simple/2027465761870101697/results',
                               'dataset' : 'JpsiKK_sigSWeight'}              
              }

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)
elif args[0] not in input_data:
    print parser.print_usage()
    sys.exit(-2)

input_data = input_data[args[0]]

from ROOT import TFile
sig_file = TFile.Open(input_data['data'])
sig_data = sig_file.Get(input_data['dataset'])
sig_st = sig_data.get().find("sigmat")
st_mean = sig_data.mean(sig_st)

st_data = []
for i in range(sig_data.numEntries()):
    r = sig_data.get(i)
    st_data.append((sig_st.getVal(), sig_data.weight()))

from operator import itemgetter
st_data = sorted(st_data, key = itemgetter(0))

from P2VV.Utilities.Resolution import SplitSigmat
split = SplitSigmat(args[0], sig_st)
bins = split.binning(sig_st)

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
means = [(sum((e[0] * e[1]) for e in bin_data)) / sum((e[1] for e in bin_data)) for bin_data in binned_data]
means = array('d', means)

# Dilution of double Gauss with scalefactors calibrated.
cache_file = TFile.Open(input_data['cache'])
simple_result_dir = cache_file.Get(input_data['result_dir'])

from P2VV.Dilution import DilutionCSFC, DilutionCSFS, DilutionSG, DilutionSFS, DilutionSFC
dilutions = {}
for calc_type, cargs, name in [(DilutionCSFS, (st_mean, simple_result_dir.Get('time_result_double_RMS_Gauss_quadratic_no_offset')), 'double_sig_calibrated'),
                               (DilutionSFS, (simple_result_dir.Get('time_result_double_RMS_Gauss_linear_no_offset'),), 'double_sig'),
                               (DilutionSG, (), 'single')]:                      
    if cargs and not cargs[-1]:
        continue
    calc = calc_type(*cargs)
    dilutions[name] = [calc.dilution(bin_data) for bin_data in binned_data]
    print 'dilution for %s = %5.4f +- %5.4f' % (tuple([name]) + calc.dilution(st_data))

graphs = []
from ROOT import TCanvas, TGraphErrors, TLegend
canvas = TCanvas('dilution_canvas', 'dilution_canvas', 550, 500)
canvas.SetLeftMargin(0.18)
from ROOT import kGreen, kBlue, kBlack, kOrange, kWhite
colors = [kGreen, kBlue, kBlack, kOrange]
first = True

legend = TLegend(0.7, 0.7, 0.89, 0.82)
legend.SetFillColor(kWhite)
legend.SetBorderSize(0)
from itertools import izip
for color, (name, ds) in izip(colors, dilutions.iteritems()):
    graph = TGraphErrors(len(means), means, array('d', [d[0] for d in ds]),
                         array('d', len(ds) * [0]), array('d', [d[1] for d in ds]))
    graph.SetName(name)
    legend.AddEntry(graph, name, 'lep')
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

legend.Draw()
