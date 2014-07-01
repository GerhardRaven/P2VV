#!/usr/bin/env python
from array import array
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = '%prog data_type')
(options, args) = parser.parse_args()

prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/project/bfys/raaij'
input_data = {'2011' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140427.root'),
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled_cache.root'),
                        'result_dir' : '1bin_15500.00fs_simple/2027465761870101697/results',
                        'dataset' : 'JpsiKK_sigSWeight'},
              '2011_Reco14' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140427.root'),
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Reco14_Prescaled_cache.root'),
                               'result_dir' : '1bin_15500.00fs_simple/2027465761870101697/results',
                               'dataset' : 'JpsiKK_sigSWeight'},
              '2012' : {'data' : os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140427.root'),
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
sig_data = sig_data.reduce('runPeriod == runPeriod::p%s' % args[0][:4])

sig_st = sig_data.get().find("sigmat")
sig_p = sig_data.get().find("B_P")
st_mean = sig_data.mean(sig_st)

datas = dict([(o, []) for o in (sig_st, sig_p)])
for i in range(sig_data.numEntries()):
    r = sig_data.get(i)
    for o, l in datas.iteritems():
        l.append((sig_st.getVal(), o.getVal(), sig_data.weight()))
    
from operator import itemgetter
for o, l in datas.iteritems():
    datas[o] = sorted(l, key = itemgetter(1))

from P2VV.Utilities.Resolution import SplitSigmat

# Dilution of double Gauss with scalefactors calibrated.
cache_file = TFile.Open(input_data['cache'])
simple_result_dir = cache_file.Get(input_data['result_dir'])

from P2VV.Dilution import DilutionCSFC, DilutionCSFS, DilutionSG, DilutionSFS, DilutionSFC
dilutions = {}
calculators = dict(((name, ct(*ca)) for ct, ca, name in [(DilutionCSFS, (st_mean, simple_result_dir.Get('time_result_double_RMS_Gauss_quadratic_no_offset')), 'double_sig_calibrated'),
                                                         (DilutionSFS, (simple_result_dir.Get('time_result_double_RMS_Gauss_linear_no_offset'),), 'double_sig'),
                                                         (DilutionSG, (), 'single')] if ((ca and ca[-1]) or not ca)))
for name, calc in calculators.iteritems():
    print 'dilution for %s = %5.4f +- %5.4f' % (tuple([name]) + calc.dilution(((v[0], v[-1]) for v in datas[sig_st])))

n_bins = 12
from array import array
binnings = {sig_p : [0 + i * 5e5 / float(n_bins) for i in range(n_bins + 1)],
            sig_st : [0.007 + i * (0.08 - 0.007) / float(n_bins) for i in range(n_bins + 1)] }

def bin_data(binning, data):
    it = iter(binning)
    boundary = it.next()
    bin_list = []
    binned_data = []
    for v in data:
        try:
            if v[1] > boundary:
                if bin_list:
                    binned_data.append(bin_list)
                    bin_list = []
                boundary = it.next()
            bin_list.append(v)
        except StopIteration:
            break
    return binned_data

binned_datas = dict([(o, bin_data(binnings[o], l)) for o, l in datas.iteritems()])

means = dict([(o.GetName(), []) for o in datas.iterkeys()])
for o, l in binned_datas.iteritems():
    for bd in l: 
        n = sum((e[1] * e[-1]) for e in bd)
        d = sum(e[-1] for e in bd)
        means[o.GetName()].append(n / d)
entries = dict([(o.GetName(), [sum(e[-1] for e in bd) for bd in bds]) for o, bds in binned_datas.iteritems()])
        
calculator = calculators['double_sig_calibrated']
assert(calculator)
dilutions = dict([(o.GetName(), []) for o in datas.iterkeys()])
for obs, binned_data in binned_datas.iteritems():
    for bd in binned_data:
        dilutions[obs.GetName()].append(calculator.dilution(((v[0], v[-1]) for v in bd)))


import shelve 
dbase = shelve.open('dilutions.db') 
dbase[args[0]] = {'means' : means, 'dilutions' : dilutions, 'entries' : entries}
dbase.close()
