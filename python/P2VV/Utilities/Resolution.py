from P2VV.RooFitWrappers import BinningCategory
from P2VV.RooFitWrappers import MappedCategory
from ROOT import RooBinning
from array import array

class SplitUtil(object):
    def __init__(self, observables, binnings, prefixes, fmt):
        self.__observables = observables
        self.__binnings = binnings
        self.__prefixes = prefixes
        assert(all(o in binnings for o in observables))
        assert(all(o in prefixes for o in observables))
        self.__cats = {}
        self.__format = fmt
        
    def observables(self):
        return self.__observables

    def binning(self, observable):
        assert(observable in self.__observables)
        return self.__binnings[o]

    def split_cat(self, observable, prefix, data = None, make_binning = 0):
        assert(observable in self.__observables)

        def __mb(data, observable, make_binning):
            from P2VV.Utilities.General import make_binning as mb
            bins = array('d', mb(data, observable, make_binning))
            print 'Created binning for %s:' % observable.GetName()
            print bins
            return bins

        if observable in self.__cats:
            return self.__cats[observable]
        if data:
            cat = data.get().find(observable.GetName() + '_cat')
        if cat:
            cat = observable.ws().put(cat)
            obs = data.get().find(observable.GetName())
            if make_binning != 0 and (obs.getBinning(prefix + '_binning').numBins() != make_binning):
                return None
            else:
                self.__cats[observable] = cat
        else:
            if make_binning != 0:
                bins = __mb(data, observable, make_binning)
            else:
                bins = self.__binnings[observable]
            binning = RooBinning(len(bins) - 1, bins, prefix + '_binning')
            observable.setBinning(binning, prefix + '_binning')
            args = dict(Observable = observable, Binning = binning,
                        CatTypeName = prefix + '_bin_')
            if data:
                args.update(dict(Data = data, Fundamental = True))
            self.__cats[observable] = BinningCategory(observable.GetName() + '_cat', **args)
        return self.__cats[observable]

    def split_cats(self, data = None, mb = 0):
        sc = [self.split_cat(o, self.__prefixes[o], data, mb) for o in self.observables()]
        if all(e == None for e in sc):
            return []
        else:
            return sc
    
    def directory(self, hd):
        return self.__format + '/' + hd

class SplitPPT(object):
    def __init__(self, data_type, p, pt):
        self.__p = p
        self.__pt = pt
        bins = {pt : array('d', [0., 1871.9, 3352.0, 1e6]),
                p : array('d', [0., 58962., 88835., 1e6])}
        prefixes = {pt : 'pt', p : 'momentum'}
        fmt = 'p_pt_{0}bins_simul'.format(len(bins[pt]) * len(bins[p]))
        SplitUtil.__init__(self, [p, pt], bins, prefixes, fmt)

class SplitSigmat(SplitUtil):
    def __init__(self, data_type, st):
        binnings = {'MC11a' : [0.01000, 0.02066, 0.02375, 0.02616, 0.02833,
                               0.03047, 0.03269, 0.03520, 0.03837, 0.04343, 0.07000],
                    'MC11a_incl_Jpsi' : [0.01000, 0.02414, 0.02720, 0.02948, 0.03145,
                                         0.03338, 0.03537, 0.03761, 0.04049, 0.04504, 0.07000],
                    'MC2011_Sim08a' : [0.01, 0.02147, 0.02473, 0.02728, 0.0296, 0.03186,
                                       0.03423, 0.0369, 0.04029, 0.04556, 0.07],
                    'MC2012' : [0.01, 0.02083, 0.02385, 0.02618, 0.02822, 0.03016, 0.03209,
                                0.03409, 0.03629, 0.03887, 0.04218, 0.04752, 0.07]}
        default = array('d', [0.01000, 0.02410, 0.02727, 0.02969, 0.03186,
                              0.03398, 0.03632, 0.03923, 0.04378, 0.07000])
        bins = binnings.get(data_type, default)
        fmt = '%sbins_%4.2ffs_simul' % (len(bins) - 1, (1000 * (bins[1] - bins[0])))
        SplitUtil.__init__(self, [st], {st : bins}, {st : 'st'}, fmt)
        
class SplitMomentum(SplitUtil):
    def __init__(self, data_type, p):
        bins = array('d', [0, 42071.68, 49609.56, 56558.79, 63839.62, 71925.85,
                           81359.74, 92988.52, 1.2e5, 1e6])
        fmt = 'momentum_{0}bins_simul'.format(len(bins))
        SplitUtil.__init__(self, [p], {p : bins}, {p : 'momentum'}, fmt)
        
class SplitPT(SplitUtil):
    def __init__(self, data_type, pt):
        bins = array('d', [0, 976.79, 1385.08, 1750.6, 2123.96, 2529.76,
                           2991.73, 3560.23, 4330.62, 5603.08, 1e5])
        fmt = 'pt_{0}bins_simul'.format(len(self.__bins))
        SplitUtil.__init__(self, [pt], {pt : bins}, {pt : 'pt'}, fmt)

class SplitPVZerr(SplitUtil):
    def __init__(self, data_type, zerr):
        bins = array('d', [0, 0.0237, 0.029, 0.0376, 1])
        self.__format = 'pv_zerr_{0}bins_simul'.format(len(bins))
        SplitUtil.__init__(self, [zerr], {zerr : bins}, {zerr : 'zerr'}, fmt)

class SplitNPV(SplitUtil):
    def __init__(self, data_type, nPV):
        bins = array('d', [-0.5 + i for i in range(5)] + [12])
        fmt = 'nPV_{0}bins_simul'.format(len(bins))
        SplitUtil.__init__(self, [nPV], {nPV : bins}, {nPV : 'nPV'}, fmt)

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
              'MC2012' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_PVRefit_ntupleB_20131211.root'),
                          'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_Mixing_MC2012_DataSet.root'),
                          'workspace' : 'Bs2JpsiPhiPrescaled_WPV_MC2012_workspace',
                          'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_Prescaled.root')},
              'MC2012_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_incl_Jpsi_ntupleB_20130916.root'),
                                    'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2012_DataSet.root'),
                                    'workspace' : 'Bs2JpsiPhiPrescaled_WPV_MC2012_incl_Jpsi_workspace',
                                    'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2012_incl_Jpsi_Prescaled.root')},
##              'MC2011_Sim08a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_ntupleB_20130909.root'),
              'MC2011_Sim08a' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_PVRefit_ntupleB_20131105.root'),
                                 'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_WPV_MC2011_Sim08a_DataSet.root'),
                                 'workspace' : 'Bs2JpsiPhiPrescaled_WPV_MC2011_Sim08a_workspace',
                                 'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_Prescaled.root')},
              'MC2011_Sim08a_incl_Jpsi' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a_incl_Jpsi_ntupleB_20130909.root'),
                                           'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_MC2011_Sim08a.root'),
                                           'workspace' : 'Bs2JpsiPhiPrescaled_MC2011_Sim08a_workspace',
                                           'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_incl_Jpsi_Prescaled.root')}
              }
