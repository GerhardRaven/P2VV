from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

import PyCintex
gbl = PyCintex.makeNamespace('')
PyCintex.loadDictionary('libswimming')

from ROOT import set_style
style = set_style()
style.SetOptTitle(1)
style.SetOptStat(0)

ws = RooObject()
ws.setWorkspace(RooWorkspace('swimming'))

t = RealVar('tau', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-1.07, 14))
m = RealVar('m', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5155, 5450))
mpsi = RealVar('mpsi', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3170))

# now build the actual signal PDF...
from ROOT import RooTruthModel as TruthModel
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay
from ROOT import RooCBShape as CrystalBall

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

mc_res = ResolutionModel('mc_res', Type = RooTruthModel, Observables = [t])
mcpdf = Pdf('mc_pdf', Type = Decay, Observables = [t], ResolutionModel = mc_res,
            Parameters = [signal_tau], Options = ['SingleSided'])

# Time resolution model
from P2VVParameterizations.TimeResolution import ResolutionModelLP2011
tres = ResolutionModelLP2011(t).Model

# Signal time pdf
sig_t = Pdf('sig_t', Type = Decay, Observables = [t], Parameters = [signal_tau],
            ResolutionModel = tres, Options = ['SingleSided'])

# B mass pdf
m_mean  = RealVar('m_mean',  Observable = False, Unit = 'MeV', Value = 5300, MinMax = (5200, 5800))
m_sigma = RealVar('m_sigma', Observable = False, Unit = 'MeV', Value = 15, MinMax = (10, 30))
sig_m = Pdf('sig_m', Type = Gaussian, Observables = (m,), Parameters = (m_mean, m_sigma ))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',  Observable = False, Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma', Observable = False, Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha', Observable = False, Unit = '', Value = 1.36, MinMax = (0.5, 3))
mpsi_n = RealVar('mpsi_n', Observable = False, Unit = '', Value = 1, MinMax = (0.1, 2))
sig_mpsi = Pdf('sig_mpsi', Type = CrystalBall, Observables = [mpsi],
               Parameters = [mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# Create signal component
signal = Component('signal')
signal.setYield(50000,0,150000)
signal[m] = sig_m
signal[mpsi] = sig_mpsi
signal[t] = sig_t

# Create combinatorical background component
comb_background = Component('comb_background')
comb_background.setYield(50000,0,150000)

m_c = RealVar( 'm_c', Observable = False, Unit = '1/MeV',
               Value = -0.0004, MinMax = (-0.1, -0.00001))
bkg_m = Pdf('bkg_m', Observables = [m], Type = Exponential, Parameters = [m_c])
comb_background[m] = bkg_m

psi_c = RealVar( 'psi_c', Observable = False, Unit = '1/MeV',
                 Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf('bkg_mpsi', Observables = [mpsi], Type = Exponential, Parameters = [psi_c])
comb_background[mpsi] = bkg_mpsi

bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1,
                  MinMax = (0.0001, 5))
comb_t = Pdf('comb_t', Type = Decay, Observables = [t], Parameters = [bkg_tau],
             ResolutionModel = tres, Options = ['SingleSided'])
comb_background[t] = comb_t

# Create psi background component
psi_background = Component('psi_background')
psi_background.setYield(50000,0,150000)
psi_background[mpsi] = sig_mpsi
psi_background[m] = bkg_m
psi_tau = RealVar('psi_tau', Observable = False, Unit = 'ps', Value = 0.5, MinMax = (0.001, 1))
psi_t = Pdf('psi_t', Type = Decay, Observables = [t], Parameters = [psi_tau],
            ResolutionModel = tres, Options = ['SingleSided'])
psi_background[t] = comb_t

# Build PDF
pdf = buildPdf((signal, comb_background, psi_background), observables = (m,mpsi), name='pdf')

# Acceptance data
from ROOT import TFile
from ROOT import RooFit
EventRange = RooFit.EventRange

file = TFile.Open('Bu2JpsiK_biased_stripping.root')
workspace = file.Get('Bu2JpsiK_workspace')
data = workspace.data('data')
data_r = data.reduce(EventRange(0, 10000))

from ROOT import RooCmdArg
NumCPU = RooCmdArg(RooFit.NumCPU(4))
pdf.fitTo(data, NumCPU)

# Observables
observables = data.get()

# sPlot
yields = {}
for component in (signal, comb_background, psi_background):
    yields[component['Name']] = ws.ws().var(component['Yield'])

from ROOT import RooStats
from ROOT import RooArgList
yield_list = RooArgList()
for y in yields.itervalues():
    yield_list.add(y)
splot = RooStats.SPlot('data_splot', 'data_splot', data,
                       pdf._target_(), yield_list)

# Mapping for Fit
from Helpers import Mapping
mapping = Mapping({m : 'm', mpsi : 'mpsi', t : 'tau'}, data)

ranges = []
for o in observables:
    if o.GetName().find('tp') != -1:
        ranges.append(o)

weights = {}
for n, o in yields.iteritems():
    name = o.GetName() + '_sw'
    weights[n] = observables.find(name)

ranges.sort()
range_names = []
for i in range(0, len(ranges), 2):
    l = ranges[i]
    r = ranges[i + 1]
    name = l.GetName() + r.GetName()
    t._target_().setRange(name, l, r)
    range_names.append(name)

## norm_range = ','.join(range_names)
## for p in [psi_t, sig_t, comb_t]:
##     p.setNormRange(norm_range)

from ROOT import TH1F
histos = {}
for n in yields.iterkeys():
    name = 'acceptance_' + n
    histos[n] = TH1F(name, name, 110, -1.07, 14);

for i in range(data.numEntries()):
    data.get(i)
    for n, histo in histos.iteritems():
        b = 1
        interval = 0
        weight = weights[n].getVal()
        while (interval < len(ranges) / 2 and b <= histo.GetNbinsX()):
            l = ranges[2 * interval].getVal()
            r = ranges[2 * interval + 1].getVal()
            c = histo.GetBinCenter(b)
            fill = False
            if c > l:
                if c > r:
                    interval += 1
                else:
                    histo.Fill(c, weight)
                    b += 1
                    fill = True
            else:
                b += 1

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 1000, 1000)
canvas.Divide(2, 2)
for i, h in enumerate(histos.values()):
    canvas.cd(i + 1)
    h.Draw()
