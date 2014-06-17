import optparse
import sys
import os

parser = optparse.OptionParser(usage = 'usage: %prog data_type')
(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)

from P2VV.Utilities.Resolution import input_data
from P2VV.Utilities.Resolution import plot_dir

if args[0] not in input_data.keys():
    print parser.print_usage()
    print "Possible samples are: %s" % ' '.join(input_data.keys())
    sys.exit(-2)


from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle

obj = RooObject( workspace = 'w')
w = obj.ws()

st = RealVar('sigmat',Title = '#sigma_{t}', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))
tdiff = RealVar('time_diff', Unit = 'ps', Observable = True, MinMax = (-0.5, 0.5))

from P2VV.CacheUtils import Cache
directory = '1bin_19000.00fs_simple/m4074732224475561764'
refit_directory = '1bin_15500.00fs_simple/5583982615024427910'
cache = Cache(input_data[args[0]]['cache'].rsplit('/', 1), directory)
cache_files = cache.cache_files()
sig_sdata = cache_files.getFromCache(directory + '/sdata/sig_sdata')[0]
refit_sig_sdata = cache_files.getFromCache(refit_directory + '/sdata/sig_sdata')[0]

assert(sig_sdata and refit_sig_sdata)

from ROOT import RooFormulaVar
for pn, an in [('sigmat', 'sigmat_refit'), ('time_diff', 'time_diff_refit')]:
    arg = refit_sig_sdata.get().find(an)
    if not arg:
        print 'observable %s not in refit data, skipping' % an
    formula = RooFormulaVar(pn, pn, '@0', RooArgList(arg))
    refit_sig_sdata.addColumn(formula)

from ROOT import TPaveText
year_label = TPaveText(0.71, 0.72, 0.89, 0.85, "NDC")
year_label.SetFillColor(0)
year_label.AddText(args[0].split('_')[0][-4:])
year_label.SetBorderSize(0)

from ROOT import TCanvas
from ROOT import kBlue, kOrange
st_canvas = TCanvas('st_canvas', 'st_canvas', 600, 400)
st_frame = st.frame()
sig_sdata.plotOn(st_frame, MarkerColor = kBlue, LineColor = kBlue, Rescale = 1 / sig_sdata.sumEntries())
refit_sig_sdata.plotOn(st_frame, MarkerColor = kOrange - 3, LineColor = kOrange - 3, Rescale = 1 / refit_sig_sdata.sumEntries())
st_frame.SetMinimum(0.)
st_frame.SetMaximum(0.06)

st_frame.GetXaxis().SetTitle('#sigma_{t} [ps])')
st_frame.GetXaxis().SetRangeUser(0, 0.08)
st_frame.GetYaxis().SetTitle('Candidates / (0.002 ps)')
st_frame.Draw()
year_label.Draw()
st_canvas.Print(os.path.join(plot_dir, 'sigmat_refit_%s.pdf' % args[0]), EmbedFonts = True) 

td_canvas = TCanvas('td_canvas', 'td_canvas', 600, 400)
td_canvas.SetLogy()
td_frame = tdiff.frame()
sig_sdata.plotOn(td_frame, MarkerColor = kBlue, LineColor = kBlue, Rescale = 1 / sig_sdata.sumEntries())
refit_sig_sdata.plotOn(td_frame, MarkerColor = kOrange - 3, LineColor = kOrange - 3, Rescale = 1 / refit_sig_sdata.sumEntries())
td_frame.SetMinimum(0.00001)
td_frame.SetMaximum(0.2)

td_frame.GetXaxis().SetTitle('t - t_{true} [ps])')
td_frame.GetYaxis().SetTitle('Candidates / (0.01 ps)')
td_frame.Draw()
year_label.Draw()
td_canvas.Print(os.path.join(plot_dir, 'time_diff_refit_%s.pdf' % args[0]), EmbedFonts = True) 
