from P2VV.RooFitWrappers import BinningCategory
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

class SplitPT(object):
    def __init__(self, data_type, pt):
        self.__pt = pt
        self.__bins = array('d', [0, 976.79, 1385.08, 1750.6, 2123.96, 2529.76,
                                  2991.73, 3560.23, 4330.62, 5603.08, 1e5])
        self.__format = 'pt_{0}bins_simul/{1}'.format(len(self.__bins))
        
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

parNames = {'N_prompt'      : ('#prompt', '\\# prompt \jpsi'),
            'N_psi_ll'      : ('#longlived', '\\# long--lived \jpsi'),
            'N_bkg'         : ('#background', '\\# background'),
            'N_sig_wpv'     : ('#wpv' , '\\# wrong PV'),
            'psi_t_fml'     : ('frac short lift', 'fraction short lived'),
            'psi_t_ll_tau'  : ('tau long', 'long--lived $\\tau$'),
            'psi_t_ml_tau'  : ('tau short', 'short--lived $\\tau$'),
            'timeResMu'     : ('mean of Gaussians', 'common mean of Gaussians'),
            'timeResComb'   : ('sf comb', '$\\text{sf}_{\\text{comb}}$'),
            'timeResFrac2'  : ('frac G2', 'fraction 2nd Gauss'),
            'timeResSigmaSF_2' : ('sf G2', '$\\text{sf}_{\\text{2}}$')}
