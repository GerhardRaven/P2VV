import os
from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

dataSample = '2011'
assert(dataSample in ['2011', '2012'])

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
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
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, mphi, st, excl_biased, selected, zerr,
               hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased] + angles
               
project_vars = [st, excl_biased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay


# B mass pdf
from P2VV.Parameterizations.MassPDFs import DoubleGauss_Signal_Mass as MassPdf
sig_m = MassPdf(m, Name = 'sig_m', AvSigParameterisation = True)
sig_m = sig_m.pdf()

# background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_m = Pdf(Name = 'bkg_m',  Type = Exponential, Parameters = [m, psi_c])

background = Component('background', (bkg_m,), Yield= (200000,500,5000000) )
signal = Component('signal', (sig_m,), Yield= (200000,500,500000) )

from P2VV.Utilities.DataHandling import readData
tree_name = 'DecayTree'
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
if dataSample == '2011':
    ## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
    input_file = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20130722.root')
elif dataSample == '2012':
    input_file = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20130722.root')
else:
    raise RuntimeError

## cut = 'runNumber > 0 && sel == 1 && sel_cleantail == 1 && (hlt1_biased == 1 || hlt1_unbiased_dec == 1) && hlt2_biased == 1 && '
cut = 'runNumber > 0 && sel == 1 && sel_cleantail == 1 && hlt1_unbiased_dec == 1 && hlt2_unbiased == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
data = readData(input_file, tree_name, ntupleCuts = cut,
                NTuple = True, observables = observables)
data_cut = data.reduce("time > 0.3")

## Build PDF
mass_pdf = buildPdf(Components = (background, signal), Observables = (m, ), Name = 'mass_pdf')
mass_pdf.Print("t")

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True,
               Verbose = False, Optimize = 2, Minimizer = 'Minuit2')

# make sweighted dataset. TODO: use mumu mass as well...
from P2VV.Utilities.SWeights import SData

result_cut = mass_pdf.fitTo(data_cut, **fitOpts)

__canvases = []

# Plot mass pdf
def plot_mass(data, prefix = ''):
    nBins = 100
    from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
    from ROOT import TCanvas
    mass_canvas = TCanvas('%smass_canvas' % prefix, '%smass_canvas' % prefix, 600, 530)
    __canvases.append(mass_canvas)
    from P2VV.Utilities.Plotting import plot
    frames = plot(mass_canvas, m, pdf = mass_pdf, data = data
                  , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = nBins)
                  , pdfOpts  = dict(LineWidth = 2)
                  , frameOpts = dict(Title = dataSample)
                  , plotResidHist = True
                  , components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                                   , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                   , 'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                   })
    
    frames[1].GetXaxis().SetTitle('M_{J/#psi K^{+}K^{-}} [MeV/c^2]')
    frames[0].GetYaxis().SetTitle('Candidates / (%3.2f MeV/c^2)' % ((m.getMax() - m.getMin()) / float(nBins)))

plot_mass(data_cut, 'cut_')

pdf_pars = mass_pdf.getParameters(data.get())
for p in sig_m.getParameters(data.get()):
    p = pdf_pars.find(p.GetName())
    p.setConstant(True)

result = mass_pdf.fitTo(data, **fitOpts)
plot_mass(data)
