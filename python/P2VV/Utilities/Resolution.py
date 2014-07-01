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
        return self.__binnings[observable]

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
        cat = data.get().find(observable.GetName() + '_cat') if data else None
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
                self.__binnings[observable] = bins
            else:
                bins = self.__binnings[observable]
            from ROOT import RooBinning
            binning = RooBinning(len(bins) - 1, bins, prefix + '_binning')
            observable.setBinning(binning, prefix + '_binning')
            args = dict(Observable = observable, Binning = binning,
                        CatTypeName = prefix + '_bin_')
            if data:
                args.update(dict(Data = data, Fundamental = True))
                from P2VV.RooFitWrappers import BinningCategory
            self.__cats[observable] = BinningCategory(observable.GetName() + '_cat', **args)
        return self.__cats[observable]

    def split_cats(self, data = None, mb = (0,)):
        sc = [self.split_cat(o, self.__prefixes[o], data, n) for o, n in zip(self.observables(), mb)]
        if all(e == None for e in sc):
            return []
        else:
            return sc
    
    def directory(self, hd, make_binning = None):
        return os.path.join(self.__format, hd)

    def sub_dir(self, make_binning = None):
        if make_binning != None and all(i for i in make_binning):
            n_bins = make_binning
        else:
            n_bins = [len(self.binning(o)) - 1 for o in self.observables()]
        return '_'.join([str(i) for i in n_bins] + ['bins'])
    
class SplitPPT(object):
    def __init__(self, data_type, p, pt):
        self.__p = p
        self.__pt = pt
        bins = {pt : array('d', [0., 1871.9, 3352.0, 1e6]),
                p : array('d', [0., 58962., 88835., 1e6])}
        prefixes = {pt : 'pt', p : 'momentum'}
        fmt = 'p_pt_simul'
        SplitUtil.__init__(self, [p, pt], bins, prefixes, fmt)

class SplitSigmat(SplitUtil):
    def __init__(self, data_type, st):
        binnings = {'MC11a' : [0.0, 0.02066, 0.02375, 0.02616, 0.02833,
                               0.03047, 0.03269, 0.03520, 0.03837, 0.04343, 0.12],
                    'MC11a_incl_Jpsi' : [0.0, 0.02414, 0.02720, 0.02948, 0.03145,
                                         0.03338, 0.03537, 0.03761, 0.04049, 0.04504, 0.12],
                    'MC2011_Sim08a' : [0.0, 0.01997, 0.02272, 0.02477, 0.02653, 0.02812,
                                       0.02967, 0.03118, 0.03273, 0.03436, 0.03612, 0.03812,
                                       0.04056, 0.04384, 0.04945, 0.12],
                    'MC2012' : [0.0001, 0.01978, 0.02242, 0.02432, 0.02593, 0.02736, 0.0287,
                                0.02996, 0.0312, 0.03244, 0.03369, 0.03495, 0.03629, 0.0377,
                                0.03924, 0.04095, 0.04297, 0.04548, 0.04894, 0.05502, 0.12],
                    'MC2012_incl_Jpsi' : [0.0, 0.0263, 0.02969, 0.03221, 0.03447, 0.03681,
                                          0.0393, 0.0424, 0.04742, 0.12],
                    '2011' : [0.0, 0.02231, 0.02537, 0.02769, 0.02972, 0.03169, 0.03371,
                              0.03599, 0.03888, 0.04343, 0.12],
                    '2011_Reco14' : [0.0001, 0.02209, 0.0247, 0.0266, 0.02821, 0.02964, 0.03099,
                                     0.03229, 0.03361, 0.03496, 0.0364, 0.03803, 0.03997, 0.04252,
                                     0.04664, 0.12],
                    '2012' : [0.0001, 0.02116, 0.02344, 0.02507, 0.02644, 0.02764, 0.02874,
                              0.02978, 0.03076, 0.03174, 0.03271, 0.0337, 0.03472, 0.0358,
                              0.03697, 0.03827, 0.03976, 0.04161, 0.04413, 0.04836, 0.12]}
        default = [0., 0.02410, 0.02727, 0.02969, 0.03186,
                   0.03398, 0.03632, 0.03923, 0.04378, 0.12]
        bins = array('d', binnings.get(data_type, default))
        fmt = 'sigmat_simul'
        SplitUtil.__init__(self, [st], {st : bins}, {st : 'st'}, fmt)
        
class SplitMomentum(SplitUtil):
    def __init__(self, data_type, p):
        bins = array('d', [0, 42071.68, 49609.56, 56558.79, 63839.62, 71925.85,
                           81359.74, 92988.52, 1.2e5, 1e6])
        fmt = 'momentum_simul'
        SplitUtil.__init__(self, [p], {p : bins}, {p : 'momentum'}, fmt)
        
class SplitPT(SplitUtil):
    def __init__(self, data_type, pt):
        bins = array('d', [0, 976.79, 1385.08, 1750.6, 2123.96, 2529.76,
                           2991.73, 3560.23, 4330.62, 5603.08, 1e5])
        fmt = 'pt_simul'
        SplitUtil.__init__(self, [pt], {pt : bins}, {pt : 'pt'}, fmt)

class SplitPVZerr(SplitUtil):
    def __init__(self, data_type, zerr):
        bins = array('d', [0, 0.0237, 0.029, 0.0376, 1])
        self.__format = 'pv_zerr_simul'
        SplitUtil.__init__(self, [zerr], {zerr : bins}, {zerr : 'zerr'}, fmt)

class SplitNPV(SplitUtil):
    def __init__(self, data_type, nPV):
        bins = array('d', [-0.5 + i for i in range(5)] + [12])
        fmt = 'nPV_simul'
        SplitUtil.__init__(self, [nPV], {nPV : bins}, {nPV : 'nPV'}, fmt)
        
parNames = {'N_prompt'         : ('#prompt',     '$N_{\\text{prompt}}$'),
            'N_psi_ll'         : ('#longlived',  '$N_{\\text{ll}}$'),
            'N_signal'         : ('#signal',     '$N_{\\text{sig}$}'),
            'N_bkg'            : ('#background', '$N_{\\text{bkg}}'),
            'N_sig_wpv'        : ('#wpv' ,       '$N_{\\text{WPV}}$'),
            'psi_t_fml'        : ('frac short lift', '$f_{\\text{sl}}$'),
            'psi_t_ll_tau'     : ('tau long', '$\\tau_{\\text{ll}}$'),
            'psi_t_ml_tau'     : ('tau short', '$\\tau_{\\text{sl}}$'),
            'timeResMu'        : ('mean of Gaussians', '$\\mu$'),
            'timeResComb'      : ('sf comb', '$\\text{sf}_{\\text{comb}}$'),
            'timeResFrac2'     : ('frac G2', '$f_2$'),
            'timeResSigmaSF_2' : ('sf G2', '$\\text{sf}_{2}$'),
            'timeResSigmaSF_1' : ('sf G1', '$\\text{sf}_{1}$'),
            'timeResSFMean'    : ('sf mean', '$\\overline{\\mathrm{sf}}$'),
            'timeResSFSigma'   : ('sf sigma', '$\\mathrm{sf}_{\sigma}$'),
            'sf_mean_offset'   : ('mean width offset', '$\\overline{\\sigma}$ offset'),
            'sf_mean_slope'    : ('mean width slope', '$\\overline{\\sigma}$ $c_0$'),
            'sf_mean_quad'     : ('mean width quad', '$\\overline{\\sigma}$ $c_1$'),
            'sf_sigma_offset'  : ('sigma width offset', '$\\sigma$ offset'),
            'sf_sigma_slope'   : ('sigma width slope', '$\\sigma_\\sigma$ $c_0$'),
            'sf_sigma_quad'    : ('sigma width quad', '$\\sigma_\\sigma$ $c_1$'),
            'sig_wpv_sigma'    : ('sigma WPV', '$\\sigma_{\\text{WPV}}$')
            }

import os
prefix = '/stuff/PhD' if os.path.exists('/stuff') else '/project/bfys/raaij'
input_data = {'2011_Reco14' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011_Reco14_Stripv20r1_ntupleB_20131002.root'),
                               'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_Mixing_2011_Reco14_DataSet.root'),
                               'workspace' : 'Bs2JpsiPhiPrescaled_WPV_2011_Reco14_workspace',
                               'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Reco14_Prescaled.root')},
              '2012' : {'data' :os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2012_ntupleB_20140429.root'),
                        'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_WPV_2012_DataSet.root'),
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
                          'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_WPV_MC2012_DataSet.root'),
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
                                           'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_WPV_MC2011_Sim08a_DataSet.root'),
                                           'workspace' : 'Bs2JpsiPhiPrescaled_WPV_MC2011_Sim08a_workspace',
                                           'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_MC2011_Sim08a_incl_Jpsi_Prescaled.root')}
              }

plot_dir = os.path.join(prefix, 'p2vv/plots')
## {'2011' : {'data' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhiPrescaled_2011_ntupleB_20130722.root'),
##                         'wpv' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_Mixing_2011_DataSet.root'),
##                         'workspace' : 'Bs2JpsiPhiPrescaled_2011_workspace',
##                         'cache' : os.path.join(prefix, 'p2vv/data/Bs2JpsiPhi_2011_Prescaled.root')}}
