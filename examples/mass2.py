from RooFitWrappers import *

ws = RooObject( workspace = 'workspace')

t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.8, 14))
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5250, 5450))
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3170))
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1020-8, 1020+8))

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5368, MinMax = (5363,5372) ) )
bkg_m = Background_BMass( Name = 'bkg_m', mass = m )

# J/psi mass pdf
from P2VVParameterizations.MassPDFs import Signal_PsiMass, Background_PsiMass
sig_mpsi = Signal_PsiMass(     Name = 'sig_mpsi', mass = mpsi )
bkg_mpsi = Background_PsiMass( Name = 'bkg_mpsi', mass = mpsi )

# phi mass pdf


# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t)
tres.setConstant('.*')

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as Signal_Time, LP2011_Background_Time as Background_Time
sig_t = Signal_Time(     Name = 'sig_t', time = t, resolutionModel = tres.model(), t_sig_tau = dict( Value = 1.5, Name = 't_sig_tau', MinMax=(1.0,2.0) ) )
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = tres.model() )
cmb_t = Background_Time( Name = 'cmb_t', time = t, resolutionModel = tres.model() )

# Create components
(ntot,nsig,fpsi) = (19000, 15000, 0.4)
npsi = (ntot-nsig)*fpsi
ncmb = (ntot-nsig)*(1-fpsi)
signal         = Component('signal',         ( sig_m.pdf(),  sig_mpsi.pdf(),  sig_t.pdf() ), Yield = ( nsig, 0.8*nsig, 1.2*nsig) )
psi_background = Component('psi_background', ( bkg_m.pdf(),  sig_mpsi.pdf(),  psi_t.pdf() ), Yield = ( npsi, 0.6*npsi, 1.4*npsi) )
cmb_background = Component('cmb_background', ( bkg_m.pdf(),  bkg_mpsi.pdf(),  cmb_t.pdf() ), Yield = ( ncmb, 0.6*ncmb, 1.4*ncmb) )

# Build PDF
pdf  = buildPdf((signal, cmb_background, psi_background ), Observables = (m,mpsi,t), Name='pdf')

from P2VVGeneralUtils import readData
data = readData(  '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
                , 'DecayTree'
                , True
                , [ m, mpsi, t, mphi ]
               )

# Fit
print 'fitting data'
from P2VVGeneralUtils import numCPU
pdf.fitTo(data, NumCPU = numCPU() , Timer=1, Extended = True, Verbose = False,  Optimize=0)

from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
canvas = TCanvas('canvas', 'canvas', 1000, 500)
obs = pdf.Observables()
for (p,o) in zip( canvas.pads(len(obs)), obs ) :
    from P2VVGeneralUtils import plot
    plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                        , 'psi*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                                        , 'cmb*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                        }
                         , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack )
                         , pdfOpts  = dict( LineWidth = 2 )
                         )
