from ROOT import RooCategory, RooRealVar
from ROOT import RooBinning
from array import array
from RooFitWrappers import BinningCategory

def reweigh(target, target_cat, source, source_cat, binning = None):
    cats = {}
    for cat, data, name in [(target_cat, target, 'target_cat'), (source_cat, source, 'source_cat')]:
        if hasattr(cat, '_target_'):
            cat = cat._target_()
        elif type(cat) == str:
            cat_name = cat
            cat = data.get().find(cat_name)
            if not cat:
                data.Print()
                raise RuntimeError('observable or category %s is not present in data' % cat_name)
        if not any([isinstance(cat, t) for t in [RooCategory, RooRealVar]]):
            raise RuntimeError('category must be either a RooRealVar or a RooCategory')
        if isinstance(cat, RooRealVar):
            assert(binning)
            if type(binning) == array:
                binning = RooBinning(len(binning) - 1, binning)
                binning.SetName('reweigh')
            cat.setBinning(binning, 'reweigh')
            cat_name = cat.GetName() + '_cat_' + data.GetName()
            test_cat = data.get().find(cat_name)
            if not test_cat:
                cat = BinningCategory(Name = cat_name, Observable = cat,
                                      Binning = binning, Data = data, Fundamental = True)
            else:
                cat = test_cat
        cats[name] = cat

    target_cat = cats['target_cat']
    source_cat = cats['source_cat']

    print target_cat
    print source_cat
    target_bins = dict([(ct.getVal(), ct.GetName()) for ct in target_cat])
    source_bins = dict([(ct.getVal(), ct.GetName()) for ct in source_cat])
    
    print target_bins
    print source_bins
    
    target_table = target.table(target_cat)
    source_table = source.table(source_cat)
    
    target_table.Print('v')
    source_table.Print('v')
    
    from collections import defaultdict
    reweigh_weights = {}
    for i, l in sorted(target_bins.iteritems()):
        sf = source_table.getFrac(source_bins[i])
        if sf == 0:
            w = 0
        else:
            try:
                w = sf / target_table.getFrac(l)
            except ZeroDivisionError:
                print 'Warning bin %s in wpv_data is 0, setting weight to 0' % l
                w = 0.
        reweigh_weights[i] = w
    
    # RooFit infinity
    from ROOT import RooNumber
    RooInf = RooNumber.infinity()
    from RooFitWrappers import RealVar
    weight_var = RealVar('reweigh_var', MinMax = (RooInf, RooInf))
    
    from ROOT import RooDataSet
    data_name = target.GetName() + 'weight_data'
    from ROOT import RooArgSet
    weight_data = RooDataSet(data_name, data_name, RooArgSet(weight_var))
    weight_var = weight_data.get().find(weight_var.GetName())
    
    for i in range(target.numEntries()):
        r = target.get(i)
        n = target_cat.getIndex()
        w = reweigh_weights[n]
        weight_var.setVal(target.weight() * w)
        weight_data.fill()
    
    target.merge(weight_data)
    target = RooDataSet(target.GetName(), target.GetTitle(), target,
                             target.get(), '', weight_var.GetName())
    
    return target, reweigh_weights
