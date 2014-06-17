from P2VV.Load import LHCbStyle
import subprocess
import shlex
from array import array
import shelve
import os
import optparse

parser = optparse.OptionParser(usage = '%prog')
parser.add_option("--note-labels", dest = "note_labels", default = False, action = 'store_true',
                  help = 'Use ana not axis labels.')
(options, args) = parser.parse_args()

#if options.note_labels:
from ROOT import gEnv, TStdFonts
if options.note_labels:
    gEnv.SetValue("Root.TTFont.13", "FreeSerif.otf")
    TStdFonts.SetFont(13, "Times-Roman")
else:
    TStdFonts.SetFont(13, "LinBiolinum")

try:
    from bsddb3 import dbshelve as shelve
except ImportError:
    import shelve
dbase = shelve.open('dilutions.db')

means = {}
dilutions = {}
entries = {}

for year in dbase.iterkeys():
    dilutions[year] = dbase[year]['dilutions']
    means[year] = dbase[year]['means']
    entries[year] = dbase[year]['entries']
dbase.close()

n_bins = 12
binnings = {'B_P' : [0 + i * 5e2 / float(n_bins) for i in range(n_bins + 1)],
            'sigmat' : [7 + i * (80 - 7) / float(n_bins) for i in range(n_bins + 1)] }

from math import sqrt, log
def eff_res(D):
    return 1000 * sqrt( -2*log(D) ) / 17.768

def eff_res_err(D, DE):
    return 1000. / 17.768 * 1 / (2 * sqrt(-2 * log(D))) * -2 / D * DE

factors = {'B_P' : 1. / 1000., 'sigmat' : 1000.}
graphs = dict([(year, {}) for year in ['2011_Reco14', '2012']])
hists = {}
from ROOT import TCanvas, TGraphAsymmErrors, TLegend, TH1D
for year in ['2011_Reco14', '2012']:
    mns = means[year]
    dls = dilutions[year]
    for var, mn in mns.iteritems():
        if var not in hists:
            hn = 'entries_%s' % var
            hist = TH1D(hn, hn, len(binnings[var]) - 1, array('d', binnings[var]))
            for i, n in enumerate(entries[year][var]):
                hist.SetBinContent(i + 1, n)
            hist.Scale(90 / max(entries[year][var]))
            hists[var] = hist
        ds = dls[var]
        f = factors[var]
        bins = binnings[var]
        graph = TGraphAsymmErrors(len(mn), array('d', [f * m  for m in  mn]),
                                  array('d', [eff_res(d[0]) for d in ds]),
                                  array('d', [(f * m - bins[i]) for i, m in enumerate(mn)]),
                                  array('d', [(bins[i + 1] - f * m) for i, m in enumerate(mn)]),
                                  array('d', [eff_res_err(*d) for d in ds]),
                                  array('d', [eff_res_err(*d) for d in ds]))
        graph.SetName(var)
        graphs[year][var] = graph

from ROOT import gStyle
if options.note_labels:
    size = 0.053 * 1.3
else:
    size = LHCbStyle.lhcbTSize
gStyle.SetTitleSize(size,"x")
gStyle.SetTitleSize(size,"y")
gStyle.SetTitleSize(size,"z")

from ROOT import TPaveText
labels = {'B_P' : TPaveText(0.75, 0.8, 0.9, 0.9, "NDC"),
          'sigmat' : TPaveText(0.13, 0.8, 0.3, 0.9, "NDC")}
for label in labels.itervalues():
    label.SetTextSize(size * 1.5 / 1.3 if options.note_labels else size * 1.5)
    label.SetFillColor(0)
    label.AddText("LHCb")
    label.SetBorderSize(0)

if options.note_labels:
    titles = {'B_P' : '#font[12]{p} [GeV #font[12]{c}^{-1}]',
              'sigmat' : '#sigma_{t} [fs]'}
else:
    titles = {'B_P' : 'P [GeV c^{-1}]',
              'sigmat' : '#sigma_{t} [fs]'}    

from ROOT import kBlue, kOrange, kWhite, kGray
legends = dict([(var, TLegend(0.75, 0.2, 0.9, 0.36)) for var in hists.iterkeys()])
for legend in legends.itervalues():
    legend.SetShadowColor(kWhite)
    legend.SetFillColor(kWhite)
    legend.SetBorderSize(1 if options.note_labels else 0)
    legend.SetMargin(0.35)
    legend.SetTextSize(size)

from itertools import izip
colors = [kBlue, kOrange - 3]

firsts = {}
canvases = {}
for color, (year, gs) in izip(colors, graphs.iteritems()):
    for var, graph in gs.iteritems():
        if var not in canvases:
            cn = 'canvas_%s' % var
            canvas = TCanvas(cn, cn)
            canvases[var] = canvas
            canvas.SetLeftMargin(size* 2)
            canvas.SetTopMargin(size)
            canvas.SetRightMargin(size)
            canvas.SetBottomMargin(size * 2.5)
            firsts[var] = True
        else:
            canvases[var].cd()

        if firsts[var]:
            hists[var].SetFillColor(kGray)
            hists[var].SetLineColor(kGray)
            hists[var].GetXaxis().SetTitle(titles[var])
            hists[var].GetXaxis().SetTitleSize(size)
            hists[var].GetXaxis().SetLabelSize(size)
            hists[var].GetXaxis().SetTitleOffset(1.1)
            hists[var].GetYaxis().SetRangeUser(0, 100)
            hists[var].GetYaxis().SetTitleOffset(0.9)
            hists[var].GetYaxis().SetTitleSize(size)
            hists[var].GetYaxis().SetLabelSize(size)
            hists[var].GetYaxis().SetTitle('resolution [fs]')
            hists[var].Draw()
            firsts[var] = False
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.Draw("P, same")
        legends[var].AddEntry(graph, ' ' + year[:4], 'lep')

from P2VV.Utilities.Resolution import plot_dir
names = {'B_P' : 'p', 'sigmat' : 'sigmat'}
for var, canvas in canvases.iteritems():
    canvas.cd()
    if options.note_labels:
        legends[var].Draw('NDC')
    labels[var].Draw()
    pn = 'resolution_vs_%s.pdf' % names[var]
    canvas.Print(os.path.join(plot_dir, pn), EmbedFonts = options.note_labels)
