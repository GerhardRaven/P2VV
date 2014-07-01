#!/usr/bin/env python
import gc
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog <resolution|physics> <2011|2012|combined> <signal|background>')

(options, args) = parser.parse_args()

if len(args) != 3:
    print parser.print_usage()
    sys.exit(-2)

sample_type = args[0]
if not sample_type in ['resolution', 'physics']:
    print parser.print_usage()
    sys.exit(-2)

years = args[1]
if not years in ['2011_Reco14', '2012', 'combined']:
    print parser.print_usage()
    sys.exit(-2)

component = args[2]
if not component in ['background', 'signal']:
    print parser.print_usage()
    sys.exit(-2)
    
if years == 'combined' and type == 'resolution':
    print 'Combined sample for resolution does not make sense'
    sys.exit(2)

from P2VV.RooFitWrappers import *
from P2VV.Load import LHCbStyle

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-1.5, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
mphi = RealVar('mdau2', Title = 'phi mass', Unit = 'MeV', Observable = True, MinMax = (990, 1050))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))
zerr = RealVar('B_s0_bpv_zerr', Title = 'Best PV Z error', Unit = 'mm', Observable = True, MinMax = (0, 1))
from math import pi
cpsi = RealVar('helcosthetaK', Title = 'cpsi', Observable = True, MinMax = (-1, 1))
ctheta = RealVar('helcosthetaL', Title = 'ctheta', Observable = True, MinMax = (-1, 1))
phi = RealVar('helphi', Title = 'helphi', Observable = True, MinMax = (-pi, pi))
angles = [cpsi, ctheta, phi]

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories
excl_biased = Category('triggerDecisionBiasedExcl', States = {'ExclBiased' : 1, 'Unbiased' : 0})
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, mphi, st, excl_biased, selected, zerr,
               hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased] + angles

project_vars = [st, excl_biased]

if sample_type == 'physics':
    prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
    input_file = os.path.join(prefix, 'p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats.root')
    from ROOT import TFile
    input_file = TFile(input_file)
    if component == 'signal':
        ds_name = 'JpsiKK_sigSWeight'
    else:
        ds_name = 'JpsiKK_bkgSWeight'
    data = input_file.Get(ds_name)
    if years == '2011_Reco14':
        data = data.reduce('runPeriod == runPeriod::p2011')
    elif years == '2012':
        data = data.reduce('runPeriod == runPeriod::p2012')
else:
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

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True,
               Verbose = False, Optimize = 2, Minimizer = 'Minuit2')

from P2VV.Parameterizations.SigmatPDFs import DoubleLogNormal
dln = DoubleLogNormal(st)
ln = dln.pdf()

# Fit
result = ln.fitTo(data, **fitOpts)

# Plot
from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 600, 400)
p = canvas.cd(1)
from P2VV.Utilities.Plotting import plot
from ROOT import kBlack, kBlue
plot(p, st, pdf = ln, data = data,
     dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = 60),
     pdfOpts  = dict(LineWidth = 3),
    xTitle = 'estimated decay-time error [ps]',
    yTitle = 'Candidates / (0.002 ps)')
