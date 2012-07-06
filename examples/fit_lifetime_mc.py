from RooFitWrappers import *
from P2VVLoad import P2VVLibrary
from ROOT import RooCBShape as CrystalBall
from ROOT import RooMsgService

## RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Eval))

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.5, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5250, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.02, 0.12))

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories
excl_biased = Category('triggerDecisionBiasedExcl', States = {'Biased' : 1, 'Unbiased' : 0})
unbiased = Category('triggerDecisionUnbiased', States = {'Unbiased' : 1, 'NotUnbiased' : 0})
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, st, excl_biased, unbiased, selected, hlt1_unbiased, hlt2_unbiased]
project_vars = [st, excl_biased, unbiased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

signal_tau = RealVar('signal_tau', Title = 'mean lifetime', Unit = 'ps', Value =  1.43, MinMax = (1., 2.5))

# Time resolution model
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as SignalTimeResolution
sig_tres = SignalTimeResolution(Name = 'sig_tres', time = t, timeResSFConstraint = True, sigmat = st,
                                timeResSF = dict( Name = 'timeResSF', Value = 1.46, MinMax = (0.1,5.),
                                                  Constant = True)
                                )
# Signal time pdf
sig_t = Pdf(Name = 'sig_t', Type = Decay,  Parameters = [t, signal_tau, sig_tres.model(), 'SingleSided'],
            ConditionalObservables = sig_tres.model().ConditionalObservables(),
            ExternalConstraints = sig_tres.model().ExternalConstraints())

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(Name = 'sig_m', mass = m, m_sig_mean = dict(Value = 5365, MinMax = (5363,5372)))

# Apply acceptance
from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
## input_file = '/stuff/PhD/p2vv/data/B_s0_Output.root'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_MC11a_biased_unbiased.root'
## data = readData(input_file, tree_name, cuts = '(sel == 1 && (triggerDecisionBiasedExcl == 1 || triggerDecisionUnbiased == 1))',
##                 NTuple = False, observables = observables)
data = readData(input_file, tree_name, cuts = '(sel == 1 && hlt1_unbiased == 1 && hlt2_unbiased == 1)',
                NTuple = False, observables = observables)

original = data.numEntries()
f = 5e4 / original
# make a NEW dataset with 
import random
new_data = RooDataSet("new_data", "new_data", data.get())
for i, obs in enumerate(data):
    ## b2 = obs.find('hlt2_biased')
    ## ub2 = obs.find('hlt2_unbiased')
    ## if b2.getIndex() == 0:
    ##     pass
    ## elif random.random() < 0.3:
    ##     ub2.setIndex(0)
    if random.random() < f:
        new_data.add(obs)

del data
data = new_data
data.table(RooArgSet(excl_biased, unbiased)).Print('v')

from P2VVParameterizations.TimeAcceptance import Paper2012_TimeAcceptance
sig_acceptance = Paper2012_TimeAcceptance(time = t, Input = '/stuff/PhD/p2vv/data/BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root',
                                          Histograms = {(excl_biased, 'Biased')   : 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1ExclB_40bins',
                                                        (excl_biased, 'Unbiased') : 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_40bins'},
                                          Data = data)
## sig_t = sig_acceptance * sig_t

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True,
               Verbose = True, Optimize = 1, Minimizer = 'Minuit2')

# make sweighted dataset. TODO: use mumu mass as well...
from P2VVGeneralUtils import SData, splot

## Fit
print 'fitting data'
## from profiler import profiler_start, profiler_stop
## profiler_start("acceptance.log")
result = sig_t.fitTo(data, SumW2Error = False, **fitOpts)
result.Print('v')

from P2VVGeneralUtils import plot
print 'plotting'
from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 500, 500)
obs = [t]

for p in sig_t.getParameters(data):
    p.removeMin()
    p.removeMax()

from ROOT import kBlack
pdfOpts = dict(ProjWData = (RooArgSet(*project_vars), data, True))
p = canvas.cd(1)
plot(p, t, pdf = sig_t, data = data
     , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
     , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
     , plotResidHist = True
     , logy = False
     , logx = True
     )
    
