import os
from P2VV.RooFitWrappers import *

obj = RooObject( workspace = 'w')
w = obj.ws()

prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
input_data = {'2011' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_ntupleB_20130207.root'),
                        'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011.root'),
                        'workspace' : 'Bs2JpsiPhiPrescaled_2011_workspace',
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled.root')},
              '2012' : {'data' :os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20121218.root'),
                        'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_2012.root'),
                        'workspace' : 'Bs2JpsiPhiPrescaled_2012_workspace',
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_Prescaled.root')},
              'MC11a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130222.root'),
                         'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC11a.root'),
                         'workspace' : 'Bs2JpsiPhiPrescaled_MC11a_workspace',
                         'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_Prescaled.root')},
              'MC11a_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_ntuple_B_MC11a_incl_Jpsi_20130103.root'),
                                   'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC11a.root'),
                                   'workspace' : 'Bs2JpsiPhiPrescaled_MC11a_workspace',
                                   'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_incl_Jpsi_Prescaled_cache.root')}
              }


t_minmax = (-1.5, 8)
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories needed for selecting events
unbiased = Category('triggerDecisionUnbiasedPrescaled', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
nPV = RealVar('nPV', Title = 'Number of PVs', Observable = True, MinMax = (0, 10))
zerr = RealVar('B_s0_bpv_zerr', Title = 'Best PV Z error', Unit = 'mm', Observable = True, MinMax = (0, 1))
observables = [t, m, mpsi, st, unbiased, nPV, zerr]

cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])
cut += ' && sel_cleantail == 1'

from P2VV.GeneralUtils import readData
tree_name = 'DecayTree'
data = {}
for s in ['2011', 'MC11a_incl_Jpsi']:
    ds = readData(input_data[s]['data'], tree_name, NTuple = True, observables = observables,
                    ntupleCuts = cut, Rename = tree_name + '_' + s)
    ds.SetName(tree_name + '_' + s)
    data[s] = ds

from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

from P2VV.Reweighing import reweigh
reweighed_data, weights = reweigh(data['MC11a_incl_Jpsi'], 'nPV',
                                  data['2011'], nPV, binning = PV_bounds)
