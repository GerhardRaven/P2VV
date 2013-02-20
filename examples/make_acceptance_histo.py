from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj = RooObject( workspace = 'w')

from math import pi
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))

# Categories
# Categories
biased = Category('triggerDecision', States = {'Biased' : 1, 'NotBiased' : 0})
unbiased = Category('triggerDecisionUnbiased', States = {'Unbiased' : 1, 'NotUnbiased' : 0})
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})
observables = [t, m, mpsi, unbiased, biased, selected]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# B mass pdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# J/psi mass pdf
mpsi_mean  = RealVar('mpsi_mean',   Unit = 'MeV', Value = 3100, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma',  Unit = 'MeV', Value = 13.6, MinMax = (5, 20))
mpsi_alpha = RealVar('mpsi_alpha',  Unit = '', Value = 1.8, MinMax = (0.5, 3), Constant = True)
mpsi_n = RealVar('mpsi_n',  Unit = '', Value = 2, MinMax = (0.1, 4), Constant = True)
psi_m  = Pdf(Name = 'psi_m', Type = CrystalBall, Parameters = [mpsi, mpsi_mean, mpsi_sigma, mpsi_alpha, mpsi_n])

# J/psi background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_mpsi = Pdf(Name = 'bkg_mpsi',  Type = Exponential, Parameters = [mpsi, psi_c])

# Create signal component
signal = Component('signal', (sig_m.pdf(), psi_m), Yield = (11810,10000,30000))

# Create combinatorical background component
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

# J/psi background component
psi_background = Component('psi_background', (bkg_m.pdf(), psi_m), Yield= (6831,500,50000) )

background = Component('background', (bkg_m.pdf(), bkg_mpsi), Yield = (3054,2000,50000) )

# Apply acceptance
from P2VV.GeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
## input_file = '/stuff/PhD/p2vv/data/B_s0_Output.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
## input_file = '/bfys/raaij/p2vv/data/swimming.root'
data = readData(input_file, tree_name, cuts = '(sel == 1 && triggerDecisionUnbiased == 1)',
                NTuple = True, observables = observables)
## Build PDF
pdf = buildPdf(Components = (signal, background), Observables = (m,), Name='pdf')
pdf.Print("t")

## Fit options
fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Verbose = True, Optimize = 2, Minimizer = 'Minuit2')

## Fit mass pdf
pdf.fitTo(data, **fitOpts)

## from ROOT import TFile
## output = TFile.Open('output.root', 'recreate')
## w = RooWorkspace('workspace')
## w.put(pdf)
## output.WriteTObject(w, w.GetName())
## output.Close()

from P2VV.GeneralUtils import plot
from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
from ROOT import TCanvas
print 'plotting'
canvas = TCanvas('mass_canvas', 'mass_canvas', 1000, 500)
obs = [m, mpsi]
for (p,o) in zip(canvas.pads(len(obs)), obs):
    plot(p, o, pdf = pdf, data = data
         , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
         , pdfOpts  = dict(LineWidth = 2)
         , plotResidHist = True
         , components = { 'psi_*'   : dict( LineColor = kGreen, LineStyle = kDashed )
                          , 'bkg_*' : dict( LineColor = kRed,   LineStyle = kDashed )
                          , 'sig_*' : dict( LineColor = kBlue,  LineStyle = kDashed )
                          }
         )

## sPlot
from P2VV.GeneralUtils import SData
constant = {}
for p in pdf.Parameters() :
    constant[p] = p.isConstant()
    p.setConstant( not p.getAttribute('Yield') )
splot = SData(Pdf = pdf, Data = data, Name = 'MassSPlot')

## Efficiency
from ROOT import TEfficiency
from array import array

bounds = [0.3, 0.33607542514801025, 0.37308794260025024, 0.4110874831676483, 0.4501281678676605, 0.49026864767074585, 0.531572699546814, 0.5741097331047058, 0.61795574426651, 0.6631938219070435, 0.7099152207374573, 0.7582205533981323, 0.8082209825515747, 0.8600398898124695, 0.9138144254684448, 0.9696980714797974, 1.0278630256652832, 1.0885035991668701, 1.1518398523330688, 1.2181226015090942, 1.2876396179199219, 1.360722541809082, 1.437757134437561, 1.5191954374313354, 1.6055715084075928, 1.6975231170654297, 1.7958197593688965, 1.9014027118682861, 2.01543927192688, 2.13940167427063, 2.275183916091919, 2.4252803325653076, 2.5930683612823486, 2.78328275680542, 3.002856731414795, 3.262538194656372, 3.5803275108337402, 3.989957332611084, 4.567111492156982, 5.5529465675354, 14.0]
bins = array('d', [bounds[i] for i in range(0, len(bounds), 2)])

# Create a TEfficiency for each component
efficiencies = {}
for component in (signal, background):
    sdata = splot.data(component.GetName())
    
    observables = sdata.get()
    t_var = observables.find(t.GetName())
    biased_var = observables.find(biased.GetName())
    name = component.GetName() + '_efficiency'
    
    # Use TEfficiency to calculate the efficiency
    efficiency = TEfficiency(name, name, len(bins) - 1, bins)
    # Set the statistic to be used to Bayesian with a uniform prior indicated by the 6
    efficiency.SetStatisticOption(6)
    efficiency.SetUseWeightedEvents()
    for i in range(sdata.numEntries()):
        evt = sdata.get(i)
        efficiency.FillWeighted(biased_var.getIndex(), sdata.weight(), t_var.getVal())
    efficiencies[component] = efficiency

# Draw the efficiencies
ec = TCanvas('efficiency_canvas', 'efficiency_canvas', 1000, 500)
for (p, (c, e)) in zip(ec.pads(len(efficiencies)), efficiencies.items()):
    p.cd()
    e.SetTitle(c.GetName())
    e.Draw()

# Make a TH1F from the efficiencies
from ROOT import TH1F
histos = []
for e in efficiencies.itervalues():
    n = e.GetName() + '_histo'
    h = TH1F(n, n, len(bins) - 1, bins)
    histos.append(h)
    for i in range(1, len(bins)):
        h.SetBinContent(i, e.GetEfficiency(i))
        h.SetBinError(i, e.GetEfficiencyErrorLow(i))

# Save the histograms
from ROOT import TFile
output = TFile('efficiencies.root', 'recreate')
for h in histos:
    r = output.WriteTObject(h, h.GetName()[:-6])
output.Close()
