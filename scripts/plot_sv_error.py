from P2VV.RooFitWrappers import *

obj = RooObject( workspace = 'w')
w = obj.ws()

mass = RealVar('B_s0_MM', Title = 'm(J/#psi K^{+}K^{-})', Unit = 'MeV/c^{2}', Observable = True,
               Value = 5368., MinMax = ( 5200., 5550. ))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0, 0.1))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (-1.5, 8))

excl_biased = Category('hlt1_excl_biased_dec', Observable = True,
                       States = { 'excl_biased' : 1, 'unbiased' : 0 } )

from ROOT import TFile
input_file = 'vx_errors.root'
f = TFile.Open(input_file)
datas = {}
for d in ['2011', '2012']:
    datas[d] = f.Get(d + '/' + 'sig_sdata')

from P2VV.Load import LHCbStyle

from ROOT import TCanvas
__canvs = []
from ROOT import kOrange
from ROOT import kBlue

from collections import defaultdict
frames = {}

canvases = [TCanvas(n, n, 600, 400) for n in ['pv_canvas', 'sv_canvas', 'st_canvas', 'psi_canvas']]
for color, (d, data) in zip((kOrange, kBlue), datas.items()):
    for canvas, (vn, xt) in zip(canvases, [('pv_err', '#sigma_{t,PV} [ps]'), ('sv_err', '#sigma_{t,SV} [ps]'),
                   (st.GetName(), '#sigma_{t} [ps]'), ('jpsi_vx_err', '#sigma_{t,J/#psi} [ps]')]):
        sw = data.sumEntries()
        if vn not in frames:
            var = data.get().find(vn)
            frame = var.frame()
            frame.GetXaxis().SetTitle(xt)
            frame.GetYaxis().SetTitle('a.u.')
            frame.GetYaxis().SetTitleOffset(1)
            frames[vn] = canvas, frame
        else:
            frame = frames[vn][1]
        data.plotOn(frame, MarkerColor = color, Rescale = 1 / sw)

for canvas, frame in frames.itervalues():
    canvas.cd()
    graph = frame.getObject(0)
    y = graph.GetY()
    n = graph.GetN()
    m = max(y[i] for i in range(n))
    frame.SetMaximum(1.1 * m)
    frame.Draw()

