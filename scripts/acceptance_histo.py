from P2VV.RooFitWrappers import *

import PyCintex
gbl = PyCintex.makeNamespace('')
PyCintex.loadDictionary('libswimming')

from ROOT import set_style
style = set_style()
style.SetOptTitle(1)
style.SetOptStat(0)

ws = RooObject( workspace = 'swimming' )

t = RealVar('tau', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-1.07, 14))
m = RealVar('m', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5155, 5450))
mpsi = RealVar('mpsi', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3170))

# now build the actual signal PDF...
from ROOT import RooTruthModel as TruthModel
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay
from ROOT import RooCBShape as CrystalBall

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5, MinMax = (1., 2.5))

mc_res = ResolutionModel('mc_res', Type = TruthModel, Parameters = [t])
mcpdf = Pdf('mc_pdf', Type = Decay,  Parameters = [t,signal_tau,mc_res,'SingleSided'])

# Time resolution model
from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution
tres = LP2011_TimeResolution(time = t)['model']

# Signal time pdf
sig_t = Pdf('sig_t', Type = Decay, Parameters = [t,signal_tau, tres, 'SingleSided'])

# B mass pdf
m_mean  = RealVar('m_mean',   Unit = 'MeV', Value = 5300, MinMax = (5200, 5800))
m_sigma = RealVar('m_sigma',  Unit = 'MeV', Value = 15, MinMax = (10, 30))
sig_m = Pdf('sig_m', Type = Gaussian, Parameters = (m,m_mean, m_sigma ))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.36, MinMax = (0.5, 3))
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 1, MinMax = (0.1, 2))
sig_mpsi = Pdf('sig_mpsi', Type = CrystalBall,  Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# Create signal component
signal = Component('signal',{m: sig_m, mpsi: sig_mpsi,t:  sig_t}, Yield= (50000,0,150000) )


# Create combinatorical background component

m_c = RealVar( 'm_c', Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.00001))
bkg_m = Pdf('bkg_m',  Type = Exponential, Parameters = [m, m_c])

psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf('bkg_mpsi', Type = Exponential, Parameters = [mpsi, psi_c])

bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1, MinMax = (0.0001, 5))
comb_t = Pdf('comb_t', Type = Decay, Parameters = [t,bkg_tau, tres, 'SingleSided'])

comb_background = Component('comb_background', { m : bkg_m, mpsi : bkg_mpsi, t : comb_t }, Yield= (50000,0,150000))


# Create psi background component
psi_tau = RealVar('psi_tau',  Unit = 'ps', Value = 0.5, MinMax = (0.001, 1))
psi_t = Pdf('psi_t', Type = Decay, Parameters = [t,psi_tau, tres, 'SingleSided'])
psi_background = Component('psi_background', { mpsi : sig_mpsi, m : bkg_m, t : comb_t }, Yield= (50000,0,150000) )


# Build PDF
pdf = buildPdf((signal, comb_background, psi_background), Observables = (m,mpsi), Name='pdf')

# Acceptance data
from ROOT import TFile
f= TFile.Open('Bu2JpsiK_biased_stripping.root')
assert f
workspace = f.Get('Bu2JpsiK_workspace')
data = workspace.data('data')
data_r = data.reduce( EventRange = (0, 10000))
pdf.fitTo(data, NumCPU = 4)

# Observables
observables = data.get()

# sPlot
yields = dict( (c['Name'], ws.ws().var(c['Yield'])) for c in (signal, comb_background, psi_background) )

from ROOT import RooStats, RooArgList
splot = RooStats.SPlot('data_splot', 'data_splot', data, pdf._target_(), RooArgList( yields.itervalues() ))

# Mapping for Fit
from Helpers import Mapping
mapping = Mapping({m : 'm', mpsi : 'mpsi', t : 'tau'}, data)

ranges = [ 0 for o in observables if o.GetName().find('tp') != -1 ]
weights = dict( ( n, observables.find(o.GetName() + '_sw') ) for n, o in yields.iteritems() )

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
histos = dict( (n, TH1F('acceptance_' + n, 'acceptance_' + n, 110, -1.07, 14)) for n in yields.iterkeys() )

for evnt in data:
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
for i, h in zip( canvas.pads(2,2), histos.values() ) : h.Draw()
