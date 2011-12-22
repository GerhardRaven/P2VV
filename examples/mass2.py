from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooMsgService

# RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

ws = RooObject( workspace = 'swimming')

t = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))
m = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5450))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3170))


# now build the actual signal PDF...
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay
from ROOT import RooCBShape as CrystalBall

# Time resolution model
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
tres = LP2011_TimeResolution(time = t)['model']

# Signal time pdf
from P2VVParameterizations.TimePDFs import Single_Exponent_Time
sig_t = Single_Exponent_Time( Name = 'sig_t',time = t, resolutionModel = tres, t_sig_tau = dict( Value = 1.5, Name = 'signal_tau', MinMax=(1,2) ) )

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass
sig_m = LP2011_Signal_Mass( Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5368, MinMax = (5363,5372) ) )
#sig_m.setConstant('.*')

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 14, MinMax = (8, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.90, MinMax = (1, 3))
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 2, MinMax = (0.1, 3))
sig_mpsi = Pdf('sig_mpsi', Type = CrystalBall, Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])


# Create combinatorical background component
from P2VVParameterizations.MassPDFs import LP2011_Background_Mass
bkg_m = LP2011_Background_Mass( Name = 'bkg_m', mass = m )
#bkg_m.setConstant('.*')


psi_c = RealVar( 'mpsi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf('bkg_mpsi', Type = Exponential, Parameters = [mpsi, psi_c])
bkg_tau = RealVar('bkg_tau', Title = 'comb background lifetime', Unit = 'ps', Value = 1, MinMax = (0.0001, 5))
comb_t = Pdf('comb_t', Type = Decay, Parameters = [t, bkg_tau, tres, 'SingleSided'])


# Create psi background component
psi_tau = RealVar('psi_tau',  Unit = 'ps', Value = 0.5, MinMax = (0.001, 1))
psi_t = Pdf('psi_t', Type = Decay, Parameters = [t,psi_tau,tres, 'SingleSided'])

# Create components
signal          = Component('signal',          ( sig_m.pdf(),  sig_mpsi,  sig_t.pdf() ), Yield = ( 32520, 30000, 35000) )
psi_background  = Component('psi_background',  ( bkg_m.pdf(),  sig_mpsi,  comb_t ),      Yield=  ( 28496, 25000, 30000) )
comb_background = Component('comb_background', ( bkg_m.pdf(),  bkg_mpsi,  comb_t ),      Yield = ( 46115, 45000, 50000) )

# Build PDF
pdf = buildPdf((signal, comb_background, psi_background ), Observables = (m,mpsi), Name='pdf')


from ROOT import TFile
#f = TFile.Open('/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root')
f = TFile.Open('/tmp/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root')
tree = f.Get('DecayTree')
noNAN =  ' && '.join( '( %s==%s )' % (i.GetName(),i.GetName()) for i in pdf.Observables() )
from ROOT import RooDataSet
 #= RooArgSet( i._var for i in pdf.Observables() )

data = RooDataSet('data','data', tree, [ i._var for i in pdf.Observables() ])#  , noNAN ) # ' && '.join([ noNAN, 'sel>0.5', '( (triggeredByUnbiasedHlt1AndHlt2>0.5) || (triggeredByBiasedHlt1AndHlt2>0.5) )' ]))
print 'got dataset with %s candidates' % data.numEntries()



# Fit
print 'fitting data'
pdf.fitTo(data, NumCPU = 4 , Timer=1, Extended = True, Verbose = False,  Optimize=1)

from ROOT import kDashed, kRed, kGreen
from ROOT import TCanvas

print 'plotting'

canvas = TCanvas('canvas', 'canvas', 1000, 500)
obs = pdf.Observables()
for (p,o) in zip( canvas.pads(len(obs)), obs ) :
    f = o.frame()
    from P2VVGeneralUtils import plot
    plot( p, o, data, pdf, components = { 'sig*'  : { 'LineColor' : RooFit.kGreen, 'LineStyle' : RooFit.kDashed }
                                        , 'psi*'  : { 'LineColor' : RooFit.kRed,   'LineStyle' : RooFit.kDashed }
                                        , 'comb*' : { 'LineColor' : RooFit.kBlue,  'LineStyle' : RooFit.kDashed }
                                        }
                         , dataOpts = dict( MarkerSize = 0.8, MarkerColor = RooFit.kBlack )
                         , pdfOpts  = dict( LineWidth = 2 )
                         )
