#!/usr/bin/env python
import gc
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog <2011|2012|combined>')

(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)

years = args[0]
if not years in ['2011_Reco14', '2012', 'combined']:
    print parser.print_usage()
    sys.exit(-2)
    
obj = RooObject( workspace = 'w')
w = obj.ws()

from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle

t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-1.5, 14))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

# signal data
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
input_file = os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats.root')
from ROOT import TFile
input_file = TFile(input_file)
ds_name = 'JpsiKK_sigSWeight'
data = input_file.Get(ds_name)
if years == '2011_Reco14':
    data = data.reduce('runPeriod == runPeriod::p2011')
elif years == '2012':
    data = data.reduce('runPeriod == runPeriod::p2012')

# calibration data
from P2VV.Utilities.Resolution import input_data
from P2VV.CacheUtils import Cache
cut = 'sel == 1 && hlt1_unbiased == 1 && hlt2_unbiased == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
cut += ' && sel_cleantail == 1'
hash_str = cut
hd = ('%d' % hash(hash_str)).replace('-', 'm')
directory = '1bin_%4.2ffs_simple/%s' % (1000 * (t.getMax() - t.getMin()), hd)
cache = Cache(input_data[years]['cache'].rsplit('/', 1), directory)
data, sdatas = cache.read_data()
if component == 'signal':
    data = single_bin_sig_sdata = sdatas['sig_sdata']
else:
    data = single_bin_bkg_sdata = sdatas['bkg_sdata']
