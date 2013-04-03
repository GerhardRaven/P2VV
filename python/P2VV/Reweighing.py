def reweigh(self, target, source, observable, binning, suffix = ''):
    from array import array
    from ROOT import RooBinning

    if type(binning) == array:
        binning = RooBinning(len(binning) - 1, binning)
        binning.SetName('reweigh')
    observable.setBinning(binning, 'reweigh')

    source_cat = source.get().find(observable.GetName())
    if not source_cat:
        cat_name = observable.GetName() + suffix if suffix else observable.GetName()
        source_cat = BinningCategory(Name = cat_name, Observable = observable,
                                     Binning = binning, Data = source, Fundamental = True)
    
    target_cat = target.get().find(observable.GetName())
    print target_obs
    target_bins = dict([(ct.getVal(), ct.GetName()) for ct in target_obs])
    source_bins = dict([(ct.getVal(), ct.GetName()) for ct in source_cat])
    
    print target_bins
    print source_bins
    
    target_table = target.table(target_cat)
    source_table = source.table(source_cat)
    
    target_table.Print('v')
    source_table.Print('v')
    
    from collections import defaultdict
    reweigh_weights = {}
    for i, l in sorted(source_bins.iteritems()):
        try:
            w = source_table.getFrac(l) / target_table.getFrac(target_bins[i])
        except ZeroDivisionError:
            print 'Warning bin %s in wpv_data is 0, setting weight to 0' % l
            w = 0.
        reweigh_weights[i] = w
    
    # RooFit infinity
    from ROOT import RooNumber
    RooInf = RooNumber.infinity()
    weight_var = RealVar('reweigh_var', MinMax = (RooInf, RooInf))
    
    from ROOT import RooDataSet
    data_name = target_data.GetName() + 'weight_data'
    weight_data = RooDataSet(data_name, data_name, RooArgSet(weight_var))
    weight_var = weight_data.get().find(weight_var.GetName())
    
    for i in range(target_data.numEntries()):
        r = target_data.get(i)
        n = target_cat.getIndex()
        w = reweigh_weights[n]
        weight_var.setVal(target_data.weight() * w)
        weight_data.fill()
    
    target_data.merge(weight_data)
    target_data = RooDataSet(target_data.GetName(), target_data.GetTitle(), target_data,
                             target_data.get(), '', weight_var.GetName())
    
    return target_data
