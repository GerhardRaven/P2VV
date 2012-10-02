from RooFitWrappers import *
from P2VVLoad import P2VVLibrary

obj = RooObject( workspace = 'w')
w = obj.ws()

time  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax=(0.3, 14))
mass  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0, 0.12))
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [time, mass, st, hlt1_biased, selected, hlt1_unbiased, hlt2_biased, hlt2_unbiased]

from P2VVGeneralUtils import readData
tree_name = 'DecayTree'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_for_fitting_20120110.root'
## input_file = '/stuff/PhD/p2vv/data/B_s0_Output.root'
## input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
input_file = '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120821_MagDownMagUp.root'
data = readData(input_file, tree_name, cuts = '(sel == 1 && (hlt1_unbiased == 1 || hlt1_biased == 1)' +\
                '&& (hlt2_unbiased == 1 || hlt2_biased == 1))',
                NTuple = True, observables = observables)
