from RooFitWrappers import *

ws = RooObject( workspace = 'workspace')

t    = RealVar('time',  Title = 'decay time',    Unit='ps',    Observable = True, MinMax = (-5, 14))
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5250, 5450))
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3170))
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1020-8, 1020+8))

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5368, MinMax = (5363,5372) ) )
bkg_m = Background_BMass( Name = 'bkg_m', mass = m )
#sig_m.setConstant('.*')
#bkg_m.setConstant('.*')

# J/psi mass pdf
from P2VVParameterizations.MassPDFs import Signal_PsiMass, Background_PsiMass
sig_mpsi = Signal_PsiMass(     Name = 'sig_mpsi', mass = mpsi )
bkg_mpsi = Background_PsiMass( Name = 'bkg_mpsi', mass = mpsi )

# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t)['model']

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as Signal_Time
sig_t = Signal_Time( Name = 'sig_t', time = t, resolutionModel = tres, t_sig_tau = dict( Value = 1.5, Name = 'signal_tau', MinMax=(1,2) ) )

from ROOT import RooDecay as Decay
bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1, MinMax = (0.0001, 5))
comb_t = Pdf('comb_t', Type = Decay, Parameters = [t, bkg_tau, tres, 'SingleSided'])

# Create psi background component
psi_tau = RealVar('psi_tau',  Unit = 'ps', Value = 0.5, MinMax = (0.001, 1))
psi_t = Pdf('psi_t', Type = Decay, Parameters = [t,psi_tau,tres, 'SingleSided'])

# Create components
ntot = 107182
nsig = 32020
fpsi = 0.5
npsi = (ntot-nsig)*fpsi
ncmb = (ntot-nsig)*(1-fpsi)
signal          = Component('signal',          ( sig_m.pdf(),  sig_mpsi.pdf(),  sig_t.pdf() ), Yield = ( nsig, 0.8*nsig, 1.2*nsig) )
psi_background  = Component('psi_background',  ( bkg_m.pdf(),  sig_mpsi.pdf(),  comb_t ),      Yield=  ( npsi, 0.6*npsi, 1.4*npsi) )
comb_background = Component('comb_background', ( bkg_m.pdf(),  bkg_mpsi.pdf(),  comb_t ),      Yield = ( ncmb, 0.6*ncmb, 1.4*ncmb) )

# Build PDF
pdf = buildPdf((signal, comb_background, psi_background ), Observables = (m,mpsi), Name='pdf')

from ROOT import TFile
f = TFile.Open('/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root')
#f = TFile.Open('/tmp/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root')
tree = f.Get('DecayTree')
noNAN =  ' && '.join( '( %s==%s )' % (i.GetName(),i.GetName()) for i in pdf.Observables() )
from ROOT import RooDataSet
data = RooDataSet('data','data', tree, [ i._var for i in pdf.Observables() ], noNAN ) # ' && '.join([ noNAN, 'sel>0.5', '( (triggeredByUnbiasedHlt1AndHlt2>0.5) || (triggeredByBiasedHlt1AndHlt2>0.5) )' ]))
print 'got dataset with %s candidates' % data.numEntries()

# Fit
print 'fitting data'
from P2VVGeneralUtils import numCPU
pdf.fitTo(data, NumCPU = numCPU() , Timer=1, Extended = True, Verbose = False,  Optimize=0)

from ROOT import kDashed, kRed, kGreen, kBlue
from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 1000, 500)
obs = pdf.Observables()
for (p,o) in zip( canvas.pads(len(obs)), obs ) :
    f = o.frame()
    from P2VVGeneralUtils import plot
    plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                        , 'psi*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                                        , 'comb*'    : dict( LineColor = kBlue,  LineStyle = kDashed )
                                        }
                         , dataOpts = dict( MarkerSize = 0.8, MarkerColor = RooFit.kBlack )
                         , pdfOpts  = dict( LineWidth = 2 )
                         )
