from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

ws = RooObject( workspace = 'swimming')

from math import pi
t = RealVar('tau', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))
m = RealVar('m', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5155, 5450))
mpsi = RealVar('mpsi', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3170))
observables = [t, m, mpsi]

# now build the actual signal PDF...
from ROOT import RooTruthModel as TruthModel
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay
from ROOT import RooCBShape as CrystalBall

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.5,
                     MinMax = (1., 2.5))

mc_res = ResolutionModel('mc_res', Type = TruthModel, Parameters = [t])
mcpdf = Pdf('mc_pdf', Type = Decay,  ResolutionModel = mc_res,
            Parameters = [t,signal_tau], Options = ['SingleSided'])

# Time resolution model
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
tres = LP2011_TimeResolution(time = t)['model']

# Signal time pdf
sig_t = Pdf('sig_t', Type = Decay,  Parameters = [t,signal_tau],
            ResolutionModel = tres, Options = ['SingleSided'])

# B mass pdf
m_mean  = RealVar('m_mean',   Unit = 'MeV', Value = 5300, MinMax = (5200, 5800))
m_sigma = RealVar('m_sigma',  Unit = 'MeV', Value = 15, MinMax = (10, 30))
sig_m = Pdf('sig_m', Type = Gaussian,  Parameters = (m,m_mean, m_sigma ))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.36, MinMax = (0.5, 3))
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 1, MinMax = (0.1, 2))
sig_mpsi = Pdf('sig_mpsi', Type = CrystalBall,
               Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# Create signal component
signal = Component('signal', (  sig_m,  sig_mpsi,  sig_t ), Yield = (10000,100,15000) )

# Create combinatorical background component

m_c = RealVar( 'm_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.00001))
bkg_m = Pdf('bkg_m', Type = Exponential, Parameters = [m, m_c])


psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf('bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])
bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1, MinMax = (0.0001, 5))
comb_t = Pdf('comb_t', Type = Decay,  Parameters = [t,bkg_tau], ResolutionModel = tres, Options = ['SingleSided'])
comb_background = Component('comb_background', ( bkg_m,  bkg_mpsi,  comb_t ), Yield = (5000,100,15000) )


# Create psi background component
psi_tau = RealVar('psi_tau',  Unit = 'ps', Value = 0.5, MinMax = (0.001, 1))
psi_t = Pdf('psi_t', Type = Decay, Parameters = [t,psi_tau], ResolutionModel = tres, Options = ['SingleSided'])
psi_background = Component('psi_background', (  comb_t,  sig_mpsi,  bkg_m ), Yield= (5000,500,15000) )


# Build PDF
pdf = buildPdf((signal, comb_background, psi_background), Observables = (m,mpsi), Name='pdf')

pdf.Print("t")

# Apply acceptance (dirty way)
from ROOT import TFile

f= TFile.Open('Bu2JpsiK.root')
assert f
workspace = f.Get('Bu2JpsiK_workspace')
data = workspace.data('data')
data = data.reduce(EventRange = (0, 4000))

from Helpers import Mapping
mapping = Mapping({m : 'm', mpsi : 'mpsi', t : 'tau'}, data)

# Fit
print 'fitting data'
pdf.fitTo(data, NumCPU = 4 , Timer=1)

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
