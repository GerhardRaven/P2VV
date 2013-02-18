from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
ws = RooObject(workspace = 'workspace')

from math import pi

m    = RealVar('B_Mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5155, 5450))
mpsi = RealVar('Jpsi_Mass', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3025, 3170))
observables = [m, mpsi]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay
from ROOT import RooCBShape as CrystalBall

# B mass pdf
m_mean  = RealVar('m_mean',   Unit = 'MeV', Value = 5300, MinMax = (5200, 5800))
m_sigma = RealVar('m_sigma',  Unit = 'MeV', Value = 15, MinMax = (10, 30))
sig_m   = Pdf(Name = 'sig_m', Type = Gaussian,  Parameters = (m,m_mean, m_sigma ))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.36, MinMax = (0.5, 5))
mpsi_n     = RealVar('mpsi_n',  Unit = '', Value = 1, MinMax = (0.1, 2), Constant = True)

sig_mpsi   = Pdf(Name = 'sig_mpsi', Type = Gaussian,  Parameters = (mpsi,mpsi_mean, mpsi_sigma ))
## sig_mpsi   = Pdf(Name = 'sig_mpsi', Type = CrystalBall,
##                  Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# Create signal component
signal = Component('signal', (sig_m,  sig_mpsi), Yield = (10000,100,15000) )

# Create combinatorical background component
m_c = RealVar( 'm_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.00001))
bkg_m = Pdf(Name = 'bkg_m', Type = Exponential, Parameters = [m, m_c])

psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf(Name = 'bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])
comb_background = Component('comb_background', (bkg_m,  bkg_mpsi), Yield = (5000,100,15000) )

# Create psi background component
psi_background = Component('psi_background', (sig_mpsi,  bkg_m), Yield= (5000,500,15000))


# Build PDF
pdf = buildPdf((signal, comb_background, psi_background), Observables = (m, mpsi), Name='pdf')

pdf.Print("t")

from ROOT import TFile

from P2VV.GeneralUtils import readData
tree_name = 'DSTReaderAlgo/Reader_All'
input_file = '/stuff/PhD/Efficiency/JpsiK/tuples_JpsiK_Detached4.root'
data = readData(input_file, tree_name, observables = observables)

# Fit
print 'fitting data'
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True, Minimizer = 'Minuit2', Optimize = 2)
pdf.fitTo(data, **fitOpts)

from ROOT import kDashed, kRed, kGreen
from ROOT import TCanvas

print 'plotting'
m_frame = m.frame()
data.plotOn(m_frame)
pdf.plotOn(m_frame, LineWidth = 2)
pdf.plotOn(m_frame, Components = ('sig_m'), LineStyle= kDashed, LineWidth = 2)
pdf.plotOn(m_frame, Components = ('bkg_m'), LineStyle= kDashed, LineWidth = 2, LineColor= kRed)

mpsi_frame = mpsi.frame()
data.plotOn(mpsi_frame)
pdf.plotOn(mpsi_frame, LineWidth= 2)
pdf.plotOn(mpsi_frame, Components='sig_mpsi', LineStyle=kDashed, LineWidth=2)
pdf.plotOn(mpsi_frame, Components='bkg_mpsi', LineStyle=kDashed, LineWidth=2, LineColor=kRed)

canvas = TCanvas('canvas', 'canvas', 1000, 500)
canvas.Divide(2, 1)
canvas.cd(1)
m_frame.Draw()

canvas.cd(2)
mpsi_frame.Draw()
