from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

dataSample = '2011'
assert(dataSample in ['2011', '2012'])

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit = 'ps', Observable = True, MinMax = (0.3, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))
sm  = RealVar('sigmam', Title = 'B mass error', Unit = 'MeV', Observable = True, MinMax = (0, 20))
P = RealVar('B_P', Title = 'momentum', Unit = 'MeV', Observable = True, MinMax = (0, 5e5))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
mphi = RealVar('mdau2', Title = 'phi mass', Unit = 'MeV', Observable = True, MinMax = (990, 1050))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))
nPV = RealVar('nPV', Title = 'number of PVs', Observable = True, MinMax = (-0.5, 10.5))
PVZ = RealVar('PVZ', Title = 'PV Z position', Observable = True, MinMax = (-300, 300))

tps = []
for i in range(6):
    tp = RealVar('tp_%d' % i, Value = 1, Observable = True, MinMax = (t.getMin(), t.getMax()))
    tps.append(tp)


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
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, st, m, sm, mpsi, mphi, excl_biased, selected, P, nPV, PVZ,
               hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased] + angles
               
project_vars = [st, excl_biased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# B mass pdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )
background = Component('background', (bkg_m.pdf(),), Yield = (10000,100,500000) )

# Create time resolution model
from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
## sig_tres = TimeResolution(Name = 'tres', time = t, sigmat = st, PerEventError = False,
##                           BiasScaleFactor = False, Cache = False, TimeResSFOffset = False,
##                           bias = dict(Value = 0, MinMax = (-1, 1), Constant = True),
##                           sigmaSF  = dict(Value = 1.45, MinMax = (0.1, 2), Constant = True))
sig_tres = TimeResolution(Name = 'tres', time = t, PerEventError = False,
                          timeResMu = dict(Name = 'timeResMu', Value = 0, Constant = True),
                          timeResSigma = dict(Name = 'timeResSigma', Value = 0.05, Constant = True))
                          

from P2VV.GeneralUtils import readData
tree_name = 'DecayTree'
if dataSample == '2011':
    input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20130131.root'
elif dataSample == '2012':
    input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_2012_ntupleB_20121212.root'
else:
    raise RuntimeError

cut = 'runNumber > 0 && sel == 1 && sel_cleantail == 1 && (hlt1_biased == 1 || hlt1_unbiased == 1) && hlt2_biased == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
data = readData(input_file, tree_name, ntupleCuts = cut, NTuple = True, observables = observables)

# Create signal component
signal = Component('signal', (sig_m.pdf(), ), Yield = (21000,100,50000))

## Build PDF
mass_pdf = buildPdf(Components = (signal, background), Observables = (m, ), Name = 'mass_pdf')
mass_pdf.Print("t")

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Offset = True,
               Verbose = False, Optimize = 2, Minimizer = 'Minuit2')

# make sweighted dataset. TODO: use mumu mass as well...
from P2VV.GeneralUtils import SData, splot

result = mass_pdf.fitTo(data, **fitOpts)

# Plot mass pdf
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
obs = [m]
for (p,o) in zip(mass_canvas.pads(len(obs)), obs):
    from P2VV.GeneralUtils import plot
    pdfOpts  = dict()
    plot(p, o, pdf = mass_pdf, data = data
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
         , frameOpts = dict(Title = dataSample)
         , plotResidHist = True
         , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                          ##, 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                          , 'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                          }
         )

for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
sig_sdata = splot.data('signal')
bkg_sdata = splot.data('background')

from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])
from ROOT import RooBinning
PV_binning = RooBinning(len(PV_bounds) - 1, PV_bounds)
PV_binning.SetName('reweigh')
nPV.setBinning(PV_binning)
nPVs = BinningCategory(Name = 'nPVs', Observable = nPV, Binning = PV_binning,
                       Data = sig_sdata, Fundamental = True)

## from P2VV.Parameterizations.OtherPDFs import AmorosoPdf

## P_pdf = AmorosoPdf(P, Name = 'P_pdf').pdf()
## P_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)

## canvas = TCanvas('canvas', 'canvas', 1000, 500)
## canvas.Divide(2, 1)
## canvas.cd(1)

## frame = P.frame()
## sig_sdata.plotOn(frame)
## P_pdf.plotOn(frame)
## frame.Draw()

from ROOT import RooGaussian
pvz_mean = RealVar(Name = 'pvz_mean', Value = 0, MinMax = (-10, 10))
pvz_sigma = RealVar(Name = 'pvz_sigma', Value = 50, MinMax = (10, 100))
pvz_pdf = Pdf(Name = 'pvz_pdf', Type = RooGaussian, Parameters = [PVZ, pvz_mean, pvz_sigma])
pvz_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)

canvas = TCanvas('canvas', 'canvas', 500, 500)
canvas.cd(1)
frame = PVZ.frame()
sig_sdata.plotOn(frame)
pvz_pdf.plotOn(frame)
frame.Draw()

## Get some more TPs
pvz_sigma.setVal(pvz_sigma.getVal() / 2)

tau = RealVar('tau', Value = 1.5, MinMax = (0.5, 2.5))

tp_decay = TPDecay('tp_decay', Time = t, Tau = tau, TurningPoints = tps,
                   ResolutionModel = sig_tres.model(), nPVs = nPVs, Data = sig_sdata,
                   PVZ = PVZ, PVZPdf = pvz_pdf)

gen_data = tp_decay.generate(RooArgSet(t, *tps), 5000)

for i in range(0, len(tps), 2):
    t._target_().setRange('r_%d' % (i / 2 + 1), tps[i]._target_(), tps[i + 1]._target_())


from ROOT import RooDecay
decay = Pdf(Name = 'decay', Type = RooDecay, Parameters = [t, tau, sig_tres.model(), 'SingleSided'])
decay.setNormRange('r_1,r_2,r_3')

RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Integration))
I = decay.createIntegral(RooArgSet(t), 'r_1')
I.recursiveRedirectServers(gen_data.get())
print I.getVal();
