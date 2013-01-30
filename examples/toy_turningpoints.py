from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from P2VVLoad import LHCbStyle
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

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
nPV = RealVar('nPV', Title = 'number of PVs', Observable = True, MinMax = (0, 12))

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

observables = [t, st, m, sm, mpsi, mphi, excl_biased, selected, P, nPV,
               hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased] + angles
               
project_vars = [st, excl_biased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )
background = Component('background', (bkg_m.pdf(),), Yield = (10000,100,500000) )

from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
if dataSample == '2011':
    input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
elif dataSample == '2012':
    input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_2012_ntupleB_20121212.root'
else:
    raise RuntimeError

cut = 'runNumber > 0 && sel == 1 && sel_cleantail == 1 && (hlt1_biased == 1 || hlt1_unbiased == 1) && hlt2_biased == 1 &&'
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
from P2VVGeneralUtils import SData, splot

result = mass_pdf.fitTo(data, **fitOpts)

# Plot mass pdf
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
mass_canvas = TCanvas('mass_canvas', 'mass_canvas', 500, 500)
obs = [m]
for (p,o) in zip(mass_canvas.pads(len(obs)), obs):
    from P2VVGeneralUtils import plot
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

from P2VVParameterizations.OtherPDFs import AmorosoPdf

P_pdf = AmorosoPdf(P, Name = 'P_pdf').pdf()
P_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)

canvas = TCanvas('canvas', 'canvas', 1000, 500)
canvas.Divide(2, 1)
canvas.cd(1)

frame = P.frame()
sig_sdata.plotOn(frame)
P_pdf.plotOn(frame)
frame.Draw()

median  = RealVar('median',   Unit = 'ps', Value = 8, MinMax = (1, 10))
k1 = RealVar('k1',  Unit = '', Value = 1.5, MinMax = (0.00001, 10))
k2 = RealVar('k2',  Unit = '', Value = 1.5, MinMax = (0.00001, 10))
frac = RealVar('frac_ln2', Value = 0.5, MinMax = (0.01, 0.99))

ln1 = LognormalPdf('ln1', Observable = sm, Median = median, Shape = k1)
ln2 = LognormalPdf('ln2', Observable = sm, Median = median, Shape = k2)

# Do our own sum pdf to have a fraction
sm_pdf = SumPdf(Name = 'ln', PDFs = [ln1, ln2], Yields = {'ln2' : frac})
sm_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)

canvas.cd(2)
frame = sm.frame()
sig_sdata.plotOn(frame)
sm_pdf.plotOn(frame)
frame.Draw()
