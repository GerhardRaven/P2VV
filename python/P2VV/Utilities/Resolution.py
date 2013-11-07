from P2VV.RooFitWrappers import BinningCategory
from P2VV.RooFitWrappers import MappedCategory
from ROOT import RooBinning
from array import array

class SplitPPT(object):
    def __init__(self, data_type, p, pt):
        self.__p = p
        self.__pt = pt
        self.__bins = {pt : array('d', [0., 1871.9, 3352.0, 1e6]),
                       p : array('d', [0., 58962., 88835., 1e6])}
        self.__format = 'p_pt_{0}bins_simul'.format(len(self.__pt_bins) * len(self.__p_bins))
        
    def split_cats(self, data = None):
        if hasattr(self, '__cats'):
            return self.__cats
        cat = data.get().find(self.__p.GetName() + '_cat')
        if cat:
            p_cat = data.get().find(momentum.GetName() + '_cat')
            p_cat = self.__p.ws().put(p_cat)
            pt_cat = data.get().find(sept.GetName() + '_cat')
            pt_cat = self.__pt.ws().put(pt_cat)
        else:
            pt_bins = self.__bins[self.__pt]
            pt_binning = RooBinning(len(pt_bins) - 1, pt_bins, 'pt_binning')
            self.__pt.setBinning(pt_binning, 'pt_binning')
            pt_cat = BinningCategory(self.__pt.GetName() + '_cat', Observable = self.__pt,
                                     Binning = pt_binning, Fundamental = True, Data = data,
                                     CatTypeName = 'pt_bin_')
            p_bins = self.__bins[self.__p]
            p_binning = RooBinning(len(p_bins) - 1, p_bins, 'momentum_binning')
            self.__p.setBinning(p_binning, 'momentum_binning')
            p_cat = BinningCategory(self.__p.GetName() + '_cat', Observable = self.__p,
                                    Binning = p_binning, Fundamental = True,
                                    Data = data, CatTypeName = 'p_bin_')
        self.__cats = [p_cat, pt_cat]
        return self.__cats
    
    def binning(self, observable):
        return self.__bins[observable]

    def directory(self, hd):
        return self.__format + '/' + hd

    def observables(self):
        return [self.__p, self.__pt]
        
class SplitSigmat(object):
    def __init__(self, data_type, st):
        self.__st = st
        if data_type == 'MC11a':
            self.__bins = array('d', [0.01000, 0.02066, 0.02375, 0.02616, 0.02833,
                                      0.03047, 0.03269, 0.03520, 0.03837, 0.04343, 0.07000])
        elif data_type == 'MC11a_incl_Jpsi':
            self.__bins = array('d', [0.01000, 0.02414, 0.02720, 0.02948, 0.03145,
                                      0.03338, 0.03537, 0.03761, 0.04049, 0.04504, 0.07000])
        else:
            self.__bins = array('d', [0.01000, 0.02410, 0.02727, 0.02969, 0.03186,
                                      0.03398, 0.03632, 0.03923, 0.04378, 0.07000])        
        self.__format = '%sbins_%4.2ffs_simul' % (len(self.__bins) - 1,
                                                  (1000 * (self.__bins[1] - self.__bins[0])))

    def split_cats(self, data):
        if hasattr(self, '__cat'):
            return [self.__cat]
        cat = data.get().find(self.__st.GetName() + '_cat')
        if cat:
            self.__cat = self.__st.ws().put(cat)
        else:
            st_binning = RooBinning(len(self.__bins) - 1, self.__bins, 'st_binning')
            self.__st.setBinning(st_binning, 'st_binning')
            self.__cat = BinningCategory(self.__st.GetName() + '_cat', Observable = self.__st,
                                         Binning = st_binning, Fundamental = True, Data = data,
                                         CatTypeName = 'st_bin_')
        return [self.__cat]
        
    def binning(self, observable):
        assert(observable == self.__st)
        return self.__bins
        
    def directory(self, hd):
        return self.__format + '/' + hd

    def observables(self):
        return [self.__st]

class SplitMomentum(object):
    def __init__(self, data_type, p):
        self.__p = p
        self.__bins = array('d', [0, 42071.68, 49609.56, 56558.79, 63839.62, 71925.85,
                                  81359.74, 92988.52, 1.2e5, 1e6])
        self.__format = 'momentum_{0}bins_simul'.format(len(self.__bins))
        
    def split_cats(self, data = None):
        if hasattr(self, '__cat'):
            return [self.__cat]
        cat = data.get().find(self.__p.GetName() + '_cat')
        if cat:
            self.__cat = self.__p.ws().put(cat)
        else:
            p_binning = RooBinning(len(self.__bins) - 1, self.__bins, 'momentum_binning')
            self.__p.setBinning(p_binning, 'momentum_binning')
            self.__cat = BinningCategory(self.__p.GetName() + '_cat', Observable = self.__p,
                                         Binning = p_binning, Fundamental = True,
                                         Data = data, CatTypeName = 'p_bin_')
        return [self.__cat]
    
    def binning(self, observable):
        assert(observable == self.__p)
        return self.__bins

    def directory(self, hd):
        return self.__format + '/' + hd

    def observables(self):
        return [self.__p]

class SplitPT(object):
    def __init__(self, data_type, pt):
        self.__pt = pt
        self.__bins = array('d', [0, 976.79, 1385.08, 1750.6, 2123.96, 2529.76,
                                  2991.73, 3560.23, 4330.62, 5603.08, 1e5])
        self.__format = 'pt_{0}bins_simul'.format(len(self.__bins))
        
    def split_cats(self, data = None):
        if hasattr(self, '__cat'):
            return [self.__cat]
        cat = data.get().find(self.__pt.GetName() + '_cat')
        if cat:
            self.__cat = self.__pt.ws().put(cat)
        else:
            pt_binning = RooBinning(len(self.__bins) - 1, self.__bins, 'pt_binning')
            self.__pt.setBinning(pt_binning, 'pt_binning')
            self.__cat = BinningCategory(self.__pt.GetName() + '_cat', Observable = self.__pt,
                                         Binning = pt_binning, Fundamental = True, Data = data,
                                         CatTypeName = 'pt_bin_')
        return [self.__cat]
        
    def binning(self, observable):
        assert(observable == self.__pt)
        return self.__bins

    def directory(self, hd):
        return self.__format + '/' + hd

    def observables(self):
        return [self.__pt]
        
class SplitPVZerr(object):
    def __init__(self, data_type, zerr):
        self.__zerr = zerr
        self.__bins = array('d', [0, 0.0237, 0.029, 0.0376, 1])
        self.__format = 'pv_zerr_{0}bins_simul'.format(len(zerr_bins))
        
    def split_cats(self, data = None):
        if hasattr(self, '__cat'):
            return [self.__cat]
        cat = data.get().find(self.__zerr.GetName() + '_cat')
        if cat:
            self.__cat = self.__zerr.ws().put(cat)
        else:
            zerr_binning = RooBinning(len(self.__bins) - 1, self.__bins, 'zerr_binning')
            self.__zerr.setBinning(zerr_binning, 'zerr_binning')
            self.__cat = BinningCategory(self.__zerr.GetName() + '_cat', Observable = self.__zerr,
                                         Binning = zerr_binning, Fundamental = True, Data = data,
                                         CatTypeName = 'zerr_bin_')
        return [self.__cat]
        
    def binning(self, observable):
        assert(observable == self.__zerr)
        return self.__bins

    def directory(self, hd):
        return self.__format + '/' + hd

    def observables(self):
        return [self.__zerr]

class SplitNPV(object):
    def __init__(self, data_type, nPV):
        self.__nPV = nPV
        self.__bins = array('d', [-0.5 + i for i in range(5)] + [12])
        self.__format = 'nPV_{0}bins_simul'.format(len(self.__bins))
        
    def split_cats(self, data = None):
        if hasattr(self, '__cat'):
            return [self.__cat]
        cat = data.get().find(self.__nPV.GetName() + '_cat')
        if cat:
            self.__cat = self.__nPV.ws().put(cat)
        else:
            nPV_binning = RooBinning(len(self.__bins) - 1, self.__bins, 'nPV_binning')
            self.__nPV.setBinning(nPV_binning, 'nPV_binning')
            self.__cat = BinningCategory(self.__nPV.GetName() + '_cat', Observable = self.__nPV,
                                         Binning = nPV_binning, Fundamental = True, Data = data,
                                         CatTypeName = 'nPV_bin_')
        return [self.__cat]
    
    def binning(self, observable):
        assert(observable == self.__nPV)
        return self.__bins

    def directory(self, hd):
        return self.__format + '/' + hd

    def observables(self):
        return [self.__nPV]

parNames = {'N_prompt'      : ('#prompt', '\\# prompt \jpsi'),
            'N_psi_ll'      : ('#longlived', '\\# long--lived \jpsi'),
            'N_bkg'         : ('#background', '\\# background'),
            'N_signal'      : ('#signal', '\\# signal'),
            'N_sig_wpv'     : ('#wpv' , '\\# wrong PV'),
            'psi_t_fml'     : ('frac short lift', 'fraction short lived'),
            'psi_t_ll_tau'  : ('tau long', 'long--lived $\\tau$'),
            'psi_t_ml_tau'  : ('tau short', 'short--lived $\\tau$'),
            'timeResMu'     : ('mean of Gaussians', 'common mean of Gaussians'),
            'timeResComb'   : ('sf comb', '$\\text{sf}_{\\text{comb}}$'),
            'timeResFrac2'  : ('frac G2', 'fraction 2nd Gauss'),
            'timeResSigmaSF_2' : ('sf G2', '$\\text{sf}_{2}$'),
            'timeResSigmaSF_1' : ('sf G1', '$\\text{sf}_{1}$'),
            'timeResSFMean' : ('sf mean', '$\\overline{\\mathrm{sf}}$'),
            'timeResSFSigma' : ('sf sigma', '$\\mathrm{sf}_{\sigma}$'),
            'sf_mean_offset' : ('sf mean offset', '$\\overline{\\sigma}$ offset'),
            'sf_mean_slope' : ('sf mean slope', '$\\overline{\\sigma}$ slope'),
            'sf_sigma_offset' : ('sf sigma offset', '$\\mathrm{sf}_{\\sigma}$ offset'),
            'sf_sigma_slope' : ('sf sigma offset', '$\\mathrm{sf}_{\\sigma}$ slope')
            }

import os
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/bfys/raaij'
input_data = {'2011' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011_ntupleB_20130722.root'),
                        'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_Mixing_2011_DataSet.root'),
                        'workspace' : 'Bs2JpsiPhiPrescaled_2011_workspace',
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled.root')},
              '2011_Reco14' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011_Reco14_Stripv20r1_ntupleB_20131002.root'),
                               'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_Mixing_2011_Reco14_DataSet.root'),
                               'workspace' : 'Bs2JpsiPhi_WPV_2011_Reco14_workspace',
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Reco14_Prescaled.root')},
              '2012' : {'data' :os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20130905.root'),
                        'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_Mixing_2012_DataSet.root'),
                        'workspace' : 'Bs2JpsiPhiPrescaled_WPV_2012_workspace',
                        'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2012_Prescaled.root')},
              'MC11a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130613.root'),
                         'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC11a.root'),
                         'workspace' : 'Bs2JpsiPhiPrescaled_MC11a_workspace',
                         'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_Prescaled.root')},
              'MC11a_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC11a_incl_Jpsi_ntupleB_20130801.root'),
                                   'wpv' : os.path.join(prefix, 'mixing/Bs2JpsiPhiPrescaled_MC11a.root'),
                                   'workspace' : 'Bs2JpsiPhiPrescaled_MC11a_workspace',
                                   'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC11a_incl_Jpsi_Prescaled.root')},
              'MC2012' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_ntupleB_20130904.root'),
                          'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_Mixing_MC2012_DataSet.root'),
                          'workspace' : 'Bs2JpsiPhiPrescaled_WPV_MC2012_workspace',
                          'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_Prescaled.root')},
              'MC2012_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_incl_Jpsi_ntupleB_20130916.root'),
                                    'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_DataSet.root'),
                                    'workspace' : 'Bs2JpsiPhiPrescaled_WPV_MC2012_incl_Jpsi_workspace',
                                    'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_incl_Jpsi_Prescaled.root')},
##              'MC2011_Sim08a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_ntupleB_20130909.root'),
              'MC2011_Sim08a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_PVRefit_ntupleB_20131105.root'),
                                 'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a.root'),
                                 'workspace' : 'Bs2JpsiPhiPrescaled_MC2011_Sim08a_workspace',
                                 'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_Prescaled.root')},
              'MC2011_Sim08a_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_incl_Jpsi_ntupleB_20130909.root'),
                                           'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a.root'),
                                           'workspace' : 'Bs2JpsiPhiPrescaled_MC2011_Sim08a_workspace',
                                           'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_incl_Jpsi_Prescaled.root')}
              }
