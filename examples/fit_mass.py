from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

dataSample = '2011'
assert(dataSample in ['2011', '2012'])

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
mphi = RealVar('mdau2', Title = 'phi mass', Unit = 'MeV', Observable = True, MinMax = (990, 1050))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

from math import pi
cpsi = RealVar('helcosthetaK', Title = 'cpsi', Observable = True, MinMax = (-1, 1))
ctheta = RealVar('helcosthetaL', Title = 'ctheta', Observable = True, MinMax = (-1, 1))
phi = RealVar('helphi', Title = 'helphi', Observable = True, MinMax = (-pi, pi))
angles = [cpsi, ctheta, phi]

muPlusTrChi2 = RealVar('muplus_track_chi2ndof',  Title = 'mu+ track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
muMinTrChi2  = RealVar('muminus_track_chi2ndof', Title = 'mu- track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
KPlusTrChi2  = RealVar('Kplus_track_chi2ndof',   Title = 'K+ track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )
KMinTrChi2   = RealVar('Kminus_track_chi2ndof',  Title = 'K- track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories
excl_biased = Category('triggerDecisionBiasedExcl', States = {'ExclBiased' : 1, 'Unbiased' : 0})
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mphi, st, excl_biased, selected,
               hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased,
               muPlusTrChi2, muMinTrChi2, KPlusTrChi2, KMinTrChi2] + angles
               
project_vars = [st, excl_biased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 10, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.36, MinMax = (0.5, 3))
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 2, MinMax = (0.1, 4), Constant = True)
psi_m  = Pdf(Name = 'psi_m', Type = CrystalBall, Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# J/psi background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf(Name = 'bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )
background = Component('background', (bkg_m.pdf(), bkg_mpsi), Yield = (20000,2000,100000) )

psi_background = Component('psi_background', (bkg_m.pdf(), psi_m), Yield= (10000,500,50000) )

from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
if dataSample == '2011':
    input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
elif dataSample == '2012':
    input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_2012_ntupleB_20121212.root'
else:
    raise RuntimeError

trackChi2Cuts = ' && '.join('%s < %f' % (trackChi2, trackChi2.getMax()) for trackChi2 in \
                             [muPlusTrChi2, muMinTrChi2, KPlusTrChi2, KMinTrChi2])
## cut = '(runNumber > 0 && sel == 1 && sel_cleantail == 1 && (hlt1_unbiased_dec == 1 || hlt1_biased == 1) && hlt2_biased == 1 && %s)' % trackChi2Cuts
cut = '(runNumber > 0 && sel == 1 && sel_cleantail == 1 && hlt1_unbiased_dec == 1 && hlt1_unbiased_tis == 0 && hlt1_unbiased_tos == 0 && hlt2_biased == 1 && %s)' % trackChi2Cuts
data = readData(input_file, tree_name, ntupleCuts = cut,
                NTuple = True, observables = observables)

# Create signal component
signal = Component('signal', (sig_m.pdf(), psi_m), Yield = (21000,10000,50000))

## Build PDF
mass_pdf = buildPdf(Components = (signal, background), Observables = (m, ), Name='pdf')
mass_pdf.Print("t")

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True,
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
mass_canvas.Print('mass_%s.eps' % dataSample)

for p in mass_pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
splot = SData(Pdf = mass_pdf, Data = data, Name = 'MassSplot')
signal_sdata = splot.data('signal')
bkg_sdata = splot.data('background')

## Fit
from P2VVGeneralUtils import plot
print 'plotting'
canvas = TCanvas('canvas', 'canvas', 1000, 1000)

obs = [t] + angles
for p, o in zip(canvas.pads(2, 2), obs):
    plot(p, o, pdf = None, data = signal_sdata
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , frameOpts = dict(Title = ' ')
         , logy = (True if o == t else False)
         )

canvas.Print('time_angles_%s.eps' % dataSample)
