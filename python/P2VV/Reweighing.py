def reweigh(target, target_cat, source, source_cat):
    source_cat = source.get().find(source_cat.GetName())    
    target_cat = target.get().find(target_cat.GetName())
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
