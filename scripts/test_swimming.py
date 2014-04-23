#!/usr/bin/env python
import optparse
import sys
import os

parser = optparse.OptionParser(usage = 'usage: %prog [swim|hist]')

parser.add_option("-n", "--ntoys", dest = "n_toys", default = 1000,
                  type = int, action = 'store', help = '# toys to run')
parser.add_option("-e", "--nevents", dest = "n_events", default = 10000,
                  type = int, action = 'store', help = '# events per toy')
parser.add_option("--ntps", dest = "n_tps", default = 6,
                  type = int, action = 'store', help = 'max # tps')
parser.add_option("-d", dest = "debug", default = False,
                  action = 'store_true', help = 'enable debug')

(options, args) = parser.parse_args()
if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)
elif args[0] not in ['swim', 'hist']:
    print parser.print_usage()
    sys.exit(-2)
swimming = bool(args[0] == 'swim')

from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle

obj = RooObject( workspace = 'w')
w = obj.ws()

t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (0, 14))
pvz = RealVar('pvz', Observable = True, MinMax = (-400, 400))
tps = [RealVar('tp%d' % i, Observable = True, MinMax = t.getRange(),
               Value = t.getMin() + (t.getMax() - t.getMin()) / (options.n_tps - 1) * i) for i in range(options.n_tps)]
nPV = RealVar('nPV', Value = 1, MinMax = (0, 20))

from ROOT import RooParamBinning
binnings = []
for i in range(0, len(tps), 2):
    tp_min = tps[i]
    tp_max = tps[i + 1]
    bn = tp_min.GetName() + tp_max.GetName()
    ## binning = RooParamBinning(tp_min._target_(), tp_max._target_(), 1)
    ## binning.SetName(bn)
    ## t.setBinning(binning, bn)
    t.setRange(bn, (tp_min._target_(), tp_max._target_()))
    binnings.append(bn)
    ## binnings.append(binning)

tau_true = ConstVar(Name = 'tau_true', Value = 1.5)
tau = RealVar('tau', Value = tau_true.getVal(), MinMax = (0.1, 4))

tres_mean = RealVar('tres_mean', Value = 0)
tres_sigma = RealVar('tres_sigma', Value = 0.045)
from ROOT import RooGaussModel
tres = ResolutionModel(Name = 'tres', Type = RooGaussModel, Parameters = (t, tres_mean, tres_sigma))
## from ROOT import RooTruthModel
## tres = ResolutionModel(Name = 'tres', Type = RooTruthModel, Parameters = (t,))

from ROOT import RooGaussian as Gaussian
pvz_mean = RealVar('pvz_mean', Value = 0)
pvz_sigma = RealVar('pvz_sigma', Value = 56)
pvz_pdf = Pdf(Name = 'pvz_pdf', Type = Gaussian, Parameters = (pvz, pvz_mean, pvz_sigma))

from ROOT import TFile
pv_file = TFile("/stuff/PhD/p2vv/data/pv_data.root")
pv_data = pv_file.Get("pv_data")

tp_gen = TPGen('tp_gen', TurningPoints = tps, PVZ = pvz, PVZPdf = pvz_pdf,
                nPVs = pv_data.get().find('nPVs'), Data = pv_data)
tp_gen.setDistance(2)
tp_gen.setSmearing(3)
tp_gen.setDebug(options.debug)
## tp_decay = TPDecay('tp_decay', TPGen = tp_gen, Time = t, Tau = tau, ResolutionModel = tres)

from ROOT import RooDataSet, RooPullVar
tau_result = RealVar("tau_result", MinMax = (0.1, 4))
result_data = RooDataSet("result_data", "result_data", RooArgSet(tau_result))

# Make a dummy histogram to 
fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2,
               Offset = True, Strategy = 2)

from ROOT import RooDecay
if swimming:    
    sig_t = Pdf(Name = 'sig_t', Type = RooDecay, Parameters = (t, tau, tres, 'SingleSided'))
    sig_t.setNormRange(','.join(b for b in binnings))
else:
    # Make an exponential binning
    def next_bin(prev_bin, tau, n_bins, min_bin, max_bin):
        from math import e, log
        a = e ** (- min_bin / tau) - e ** (- max_bin / tau)
        return - tau * log(e ** (- prev_bin / tau) - a / n_bins)
    from array import array
    bins = array('d', [0])
    n_bins = 20
    for i in range(n_bins - 1):
        b = next_bin(bins[-1], tau.getVal(), n_bins, 0, t.getMax())
        bins.append(b)
    bins.append(t.getMax())

    from ROOT import TH1D
    dummy_hist = TH1D("dummy", "dummy", len(bins) - 1, bins)
    for i in range(len(bins) - 1):
        dummy_hist.SetBinContent(i + 1, 0.5)

    acc_fun = BinnedFun(Name = 'acc_fun', ObsVar = t, Histogram = dummy_hist)
    acc_res = CubicSplineGaussModel(Name = 'acc_res', ResolutionModel = tres, Efficiency = acc_fun)
    sig_t = Pdf(Name = 'sig_t', Type = RooDecay, Parameters = (t, tau, acc_res, 'SingleSided'))
    
signal = Component('signal', [sig_t], Yield = (3000, 100, 1000000))

from P2VV.Utilities.Swimming import swimming_histos

def print_data(i):
    r = data.get(i)
    print [t_data.getVal()] + [tp.getVal() for tp in r if tp.GetName().startswith('tp')]

for i in range(options.n_toys):
    tp_data = tp_gen.generate(tps, options.n_events)
    ## RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Generation))
    data = sig_t.generate([t], NumEvents = options.n_events, ProtoData = tp_data)
    t_data = data.get().find(t.GetName())

    if not swimming:
        hist = swimming_histos(bins, t, sig_t, data, [signal])[signal]
        from ROOT import TCanvas
        canvas = TCanvas('canvas', 'canvas', 600, 400)
        canvas.SetLogx(True)
        hist.Draw()
        assert(False)
        coefs = acc_fun.coefficients()[acc_fun.GetName()]
        for i, c in enumerate(coefs):
            c.setVal(hist.GetBinContent(i + 1))
        
    ## RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Integration))
    result = sig_t.fitTo(data, **fitOpts)

    tau_result.setVal(tau.getVal())
    tau_result.setError(tau.getError())
    result_data.add(RooArgSet(tau_result))

    tau.setVal(tau_true.getVal())
    tau.setError(0.1)

tau_pull = RooPullVar("tau_pull", "tau_pull", tau_result._target_(), tau_true._target_())
tau_pull = result_data.addColumn(tau_pull)
tau_pull = w.put(tau_pull)

mean = RealVar('mean', Value = 0, MinMax = (-4, 4))
sigma = RealVar('sigma', Value = 1, MinMax = (0.1, 5))
from ROOT import RooGaussian
gauss = Pdf(Name = 'gauss', Type = RooGaussian, Parameters = (tau_pull, mean, sigma))
result = gauss.fitTo(result_data, **fitOpts)

from ROOT import TPostScript
ps = TPostScript("test.eps", 113)
ps.FontEmbed()
from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 600, 400)
frame = tau_pull.frame(Range = (-5, 5))
frame.GetXaxis().SetTitle('#tau pull')
frame.GetYaxis().SetTitle('M Experiments')
result_data.plotOn(frame)
gauss.plotOn(frame)
frame.Draw()
canvas.Update()
ps.Close()
