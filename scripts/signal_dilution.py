from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService
import os
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
hlt1_excl_biased = Category('hlt1_excl_biased_dec', States = {'excl_biased' : 1, 'unbiased' : 0}, Observable = True)
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, mphi, st, hlt1_excl_biased, selected, zerr,
               hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased] + angles
               
project_vars = [st, hlt1_excl_biased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# B mass pdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
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
background = Component('background', (bkg_m.pdf(), bkg_mpsi), Yield = (100000,100,500000) )

psi_background = Component('psi_background', (bkg_m.pdf(), psi_m), Yield= (200000,500,500000) )

from P2VV.Utilities.DataHandling import readData
tree_name = 'DecayTree'
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
if dataSample == '2011':
    input_file = os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root')
elif dataSample == '2012':
    input_file = os.path.join(prefix, '/stuff/PhD/p2vv/data/Bs2JpsiPhi_2012_ntupleB_20121212.root')
else:
    raise RuntimeError

cut = 'runNumber > 0 && sel == 1 && sel_cleantail == 1 && (hlt1_biased == 1 || hlt1_unbiased_dec == 1) && hlt2_biased == 1 && '
## cut = 'runNumber > 0 && sel == 1 && sel_cleantail == 1 && hlt1_unbiased_dec == 1 && hlt2_unbiased == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
data = readData(input_file, tree_name, ntupleCuts = cut,
                NTuple = True, observables = observables)

# Create signal component
signal = Component('signal', (sig_m.pdf(), psi_m), Yield = (27437,10000,50000))

## Build PDF
mass_pdf = buildPdf(Components = (signal, background), Observables = (m, ), Name = 'mass_pdf')
mass_pdf.Print("t")

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True,
               Verbose = False, Optimize = 2, Minimizer = 'Minuit2')

# make sweighted dataset. TODO: use mumu mass as well...
from P2VV.Utilities.SWeights import SData, splot

single_bin_result = mass_pdf.fitTo(data, **fitOpts)

from P2VV.Utilities.General import getSplitPar
# categories for splitting the PDF
# get mass parameters that are split
from ROOT import RooBinning

split_cats = [[hlt1_excl_biased]]
split_pars = [[par for par in mass_pdf.Parameters() if par.getAttribute('Yield')]]

# build simultaneous mass PDF
sWeight_mass_pdf = SimultaneousPdf(mass_pdf.GetName() + '_simul',
                                   MasterPdf       = mass_pdf,
                                   SplitCategories = split_cats,
                                   SplitParameters = split_pars)

from math import sqrt
# set yields for categories
split_cat_pars = sWeight_mass_pdf.getVariables()
for ct in hlt1_excl_biased:
    sig_yield = getSplitPar('N_signal', ct.GetName(), split_cat_pars)
    bkg_yield = getSplitPar('N_background', ct.GetName(), split_cat_pars)
    
    sel_str = '!(%s-%d)' % (hlt1_excl_biased.GetName(), ct.getVal())
    
    nEv    = data.sumEntries()
    nEvBin = data.sumEntries(sel_str)
    
    sig_yield.setVal( sig_yield.getVal() * nEvBin / nEv )
    sig_yield.setError( sqrt( sig_yield.getVal() ) )
    sig_yield.setMin(0.)
    sig_yield.setMax(nEvBin)
    bkg_yield.setVal( bkg_yield.getVal() * nEvBin / nEv )
    bkg_yield.setError( sqrt( bkg_yield.getVal() ) )
    bkg_yield.setMin(0.)
    bkg_yield.setMax(nEvBin)

for i in range(5):
    sWeight_mass_result = sWeight_mass_pdf.fitTo(data, **fitOpts)
    if sWeight_mass_result.status() == 0:
        break
assert(sWeight_mass_result.status() == 0)
sWeight_mass_result.SetName('sWeight_mass_result')

## Canvas for correlation histograms
from ROOT import TCanvas
corr_canvas = TCanvas('corr_canvas', 'corr_canvas', 1000, 500)
corr_canvas.Divide(2, 1)
from ROOT import gStyle
gStyle.SetPalette(1)

## Plot correlation histogram
corr_canvas.cd(1)
corr_hist_sWmass = sWeight_mass_result.correlationHist()
corr_hist_sWmass.GetXaxis().SetLabelSize(0.03)
corr_hist_sWmass.GetYaxis().SetLabelSize(0.03)
corr_hist_sWmass.SetContour(20)
corr_hist_sWmass.Draw('colz')

from P2VV.Utilities.SWeights import SData
sData = SData(Pdf = sWeight_mass_pdf, Data = data, Name = 'SimulMassSPlot')
sig_sdata = sData.data('signal')
bkg_sdata = sData.data('background')

# Plot mass pdf

paper_calibration = FormulaVar('st_calibration', '1.45 * @0', [st])
dg_calibration = FormulaVar('dg_calibration', "-6.81838 * @0 * @0 + 1.72 * @0", [st])
sg_calibration = FormulaVar('sg_calibration', "-0.00455 + 1.92 * @0-9.9352 * @0 * @0", [st])
linear_calibration = FormulaVar('linear_calibration', "0.00555 + 1.267 * @0", [st])
from P2VV.Dilution import dilution
r = signal_dilution(psi_sdata, st, paper_calibration)
r = signal_dilution(psi_sdata, st, dg_calibration)
# r = signal_dilution(sig_sdata, st, sg_calibration)
# r = signal_dilution(sig_sdata, st, linear_calibration)
