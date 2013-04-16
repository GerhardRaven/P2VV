from array import array
from math import exp, log, sqrt
from RooFitWrappers import FormulaVar

__keep = []
__bin_counter = 0
__histos = []

def dilution_bin(t_var, bin_data, result, bin_name, dilution_binning,
                 signal = [], subtract_pdf = None, subtract = []):
    from ROOT import sigmaFromFT
    dms = -17.768 ** 2 / 2

    n_bins = 512
    neg_range = 0 - t_var.getMin()
    diff = neg_range / float(n_bins)

    from ROOT import TH1D
    data_histo = TH1D('data_histo', 'data_histo', len(dilution_bounds) - 1, dilution_bounds)
    __histos.append(data_histo)
    data_histo.Sumw2()
    time_var = data.get().find(t_var.GetName())

    weighted = data.isWeighted()
    for i in range(data.numEntries()):
        r = data.get(i)
        value = time_var.getVal()
        ## External weight
        if weighted:
            w = data.weight()
        else:
            w = 1.
        r = data_histo.Fill(value, w)

    ft_bounds = array('d', (0. + i * diff for i in range(n_bins + 1)))
    ft_histo = TH1D('ft_histo', 'ft_histo', len(ft_bounds) - 1, ft_bounds)
    if not subtract_pdf:
        # Fill the FT histogram
        total = 0
        for i in range(1, n_bins + 1):
            d = data_histo.GetBinContent(i)
            e = data_histo.GetBinError(i)
            total += d
            ft_histo.SetBinContent(n_bins + 1 - i, d)
            ft_histo.SetBinError(n_bins + 1 - i, e)        
        ft_histo.SetEntries(total)
        return sigmaFromFT(ft_histo, 17.768, 0.024)

    subtract_yields = {}
    for b in subtract:
        y = [p for p in result.floatParsFinal() if p.GetName() == '_'.join(b.getYield().GetName(), bin_name)]
        assert(len(y) == 1)
        subtract_yield[b] = y[0].getVal()

    for c in subtract:
        y = c.getYield()
        result_yield = subtract_yields[c]
        pdf_yield = subtract_pdf.getVariables().find(y.GetName())
        if not pdf_yield:
            break
        print 'setting yield %s of subtract pdf to %f' % (pdf_yield.GetName(), result_yield)
        pdf_yield.setVal(result_yield)
        
    subtract_histo = subtract_pdf.createHistogram(subtract_pdf.GetName().replace('_pdf', '_histo') + '_' + bin_name,
                                                  t_var._target_(), RooFit.Binning(dilution_binning),
                                                  RooFit.Scaling(False))
    n_subtract = sum(subtract_yields.itervalues())
    total = sum([s.getVal() for s in signal_yields] + [n_subtract])


    sig_yield_names = [s.getYield().GetName() for s in signal]
    signal_yields = []
    for p in result.floatParsFinal():
        for s in sig_yield_names:
            if p.GetName() == '_'.join(s, bin_name):
                signal_yields.append(p.getVal())

    total = sum(signal_yields + [n_subtract])
    sub_int = subtract_histo.Integral()
    data_int = data_histo.Integral()
    scale =  (data_int * n_subtract) / (sub_int * total)
    subtract_histo.Scale(scale)

    # Add like this to preserve errors
    tmp_histo = data_histo.Clone('tmp_histo')
    tmp_histo.Add(subtract_histo, -1)

    # Fill the FT histogram
    total = 0
    for i in range(1, n_bins + 1):
        d = tmp_histo.GetBinContent(i)
        total += d
        e = tmp_histo.GetBinError(i)
        ft_histo.SetBinContent(n_bins + 1 - i, d)
        ft_histo.SetBinError(n_bins + 1 - i, e)        
    ft_histo.SetEntries(total)
    return sigmaFromFT(ft_histo, 17.768, 0.024)

# Calculate dilution
def dilution_bins(t_var, data, sigmat, sigmat_cat, result, signal = [], subtract = [], raw = False, calibration = None):
    __bin_counter = 0
    if calibration:
        assert(sigmat.GetName() in [p.GetName() for p in calibration.getVariables()])
        calibration.redirectServers(data.get())

    assert(t_var.getMin() < 0)

    from ROOT import RooBinning
    n_bins = 512
    neg_range = 0 - t_var.getMin()
    diff = neg_range / float(n_bins)
    dilution_bounds = array('d', (t_var.getMin() + i * diff for i in range(n_bins + 1)))
    b = diff
    t_max = t_var.getMax()
    while b < t_max:
        dilution_bounds.append(b)
        b += diff
    if dilution_bounds[-1] != t_max:
        dilution_bounds.append(t_max)
    dilution_binning = RooBinning(len(dilution_bounds) - 1, dilution_bounds)
    dilution_binning.SetName('dilution_binning')

    if not subtract:
        for ct in sigmat_cat:
            bin_data = data.reduce('{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))
            d = dilution_bin(t_var, bin_data, result, signal, ct.GetName(), dilution_binning)


        
    subtract_pdf = buildPdf(Components = subtract, Observables = (t_var,),
                            Name='subtract_%s_pdf' % '_'.join(c.GetName() for c in subtract))    
    sub_pars = subtract_pdf.getVariables()
    for res_par in result.floatParsFinal():
        sub_par = sub_pars.find(res_par.GetName())
        if sub_par:
            print 'setting value of subtract pdf parameter %s' % sub_par.GetName()
            sub_par.setVal(res_par.getVal())
            sub_par.setError(res_par.getError())

    from itertools import product
    sig_yield_names = [s.getYield().GetName() for s in signal]
    from collections import defaultdict
    signal_yields = defaultdict(list)
    bin_names = [ct.GetName() for ct in sigmat_cat]
    for p in result.floatParsFinal():
        for s, b in product(sig_yield_names, bin_names):
            if p.GetName() == '_'.join((s, b)):
                signal_yields[b].append(p.getVal())
    signal_yield = sum(e[0] for e in signal_yields.itervalues())

    total = 0
    for ct in sigmat_cat:
        bin_data = data.reduce('{0} == {0}::{1}'.format(st_cat.GetName(), ct.GetName()))
        bin_name = ct.GetName()
        d = dilution_bin(t_var, bin_data, result, bin_name, dilution_binning,
                         signal, subtract_pdf, subtract)
        total += sum(signal_yields[ct.GetName()]) * d ** 2
    total /= signal_yield
    return total

# Calculate dilution
def dilution(t_var, data, sigmat = None, result = None, signal = [], subtract = [], raw = False, simultaneous = False, calibration = None):
    __bin_counter = 0
    if calibration:
        assert(sigmat.GetName() in [p.GetName() for p in calibration.getVariables()])
        calibration.redirectServers(data.get())

    assert(t_var.getMin() < 0)

    from ROOT import RooBinning
    n_bins = 512
    neg_range = 0 - t_var.getMin()
    diff = neg_range / float(n_bins)
    dilution_bounds = array('d', (t_var.getMin() + i * diff for i in range(n_bins + 1)))
    b = diff
    t_max = t_var.getMax()
    while b < t_max:
        dilution_bounds.append(b)
        b += diff
    if dilution_bounds[-1] != t_max:
        dilution_bounds.append(t_max)
    dilution_binning = RooBinning(len(dilution_bounds) - 1, dilution_bounds)
    dilution_binning.SetName('dilution_binning')
    ## t_var.setBinning(dilution_binning)

    from ROOT import TH1D
    data_histo = TH1D('data_histo', 'data_histo', len(dilution_bounds) - 1, dilution_bounds)
    data_histo.Sumw2()
    time_var = data.get().find(t_var.GetName())
    if sigmat:
        res_var = data.get().find(sigmat.GetName())
    else:
        res_var = None

    dms = -17.768 ** 2 / 2

    weights = []
    weighted = data.isWeighted()
    for i in range(data.numEntries()):
        r = data.get(i)
        value = time_var.getVal()
        ## External weight
        if weighted:
            w = data.weight()
        else:
            w = 1.

        weights.append(w)
        r = data_histo.Fill(value, w)

    if not subtract or raw:
        ft_bounds = array('d', (0. + i * diff for i in range(n_bins + 1)))
        ft_histo = TH1D('ft_histo', 'ft_histo', len(ft_bounds) - 1, ft_bounds)

        # Fill the FT histogram
        total = 0
        for i in range(1, n_bins + 1):
            d = data_histo.GetBinContent(i)
            e = data_histo.GetBinError(i)
            total += d
            ft_histo.SetBinContent(n_bins + 1 - i, d)
            ft_histo.SetBinError(n_bins + 1 - i, e)        
        ft_histo.SetEntries(total)

        from ROOT import TCanvas
        ft_canvas = TCanvas('ft_canvas', 'ft_canvas', 800, 400)
        ft_canvas.Divide(2, 1)
        pad = ft_canvas.cd(1)
        pad.SetLogy()
        data_histo.Draw()
        pad = ft_canvas.cd(2)
        ft_histo.Draw()

        __keep.append(ft_canvas)
        __keep.append(ft_histo)
        __keep.append(data_histo)

        # Calculate the dilution using Wouter's macro
        from ROOT import sigmaFromFT
        D = sigmaFromFT(ft_histo, 17.768, 0.024)
        return D, weights
    else:
        from ROOT import TCanvas
        ft_canvas = TCanvas('ft_canvas', 'ft_canvas', 1200, 400)
        ft_canvas.Divide(3, 1)
        pad = ft_canvas.cd(1)
        pad.SetLogy()
        data_histo.Draw()

        __keep.append(ft_canvas)
        __keep.append(data_histo)

    # Create a histogram of our WPV component
    from ROOT import RooFit
    from P2VV.RooFitWrappers import buildPdf

    sig_yield_names = [s.getYield().GetName() if type(s) != str else s for s in signal]
    signal_yields = [p for p in result.floatParsFinal() if any([p.GetName().startswith(s) for s in sig_yield_names])]
    data_int = data_histo.Integral()

    # Set the correct yields in the pdf
    # Build the PDF used for subtraction.
    subtract_pdf = buildPdf(Components = subtract, Observables = (t_var,),
                            Name='subtract_%s_pdf' % '_'.join(c.GetName() for c in subtract))    
    sub_pars = subtract_pdf.getVariables()
    for res_par in result.floatParsFinal():
        sub_par = sub_pars.find(res_par.GetName())
        if sub_par:
            print 'setting value of subtract pdf parameter %s' % sub_par.GetName()
            sub_par.setVal(res_par.getVal())
            sub_par.setError(res_par.getError())

    subtract_yields = {}
    for b in subtract:
        subtract_yields[b] = [p for p in result.floatParsFinal() if p.GetName().startswith(b.getYield().GetName())]

    for c in subtract:
        y = c.getYield()
        result_yield = sum([cy.getVal() for cy in subtract_yields[c]])
        pdf_yield = subtract_pdf.getVariables().find(y.GetName())
        if not pdf_yield:
            break
        print 'setting yield %s of subtract pdf' % pdf_yield.GetName()
        pdf_yield.setVal(result_yield)
        
    subtract_histo = subtract_pdf.createHistogram(subtract_pdf.GetName().replace('_pdf', '_histo'), t_var._target_(),
                                                  RooFit.Binning(dilution_binning),
                                                  RooFit.Scaling(False))

    # Calculate appropriate scale factor. Scale histograms such that the ration
    # of their integrals match the ratio of wpv / total yields
    n_subtract = sum([b.getVal() for yields in subtract_yields.itervalues() for b in yields])
    total = sum([s.getVal() for s in signal_yields] + [n_subtract])
    sub_int = subtract_histo.Integral()
    scale =  (data_int * n_subtract) / (sub_int * total)
    subtract_histo.Scale(scale)

    # Create the histogram to be transformed
    ft_bounds = array('d', (0. + i * diff for i in range(n_bins + 1)))
    ft_histo = TH1D('ft_histo', 'ft_histo', len(ft_bounds) - 1, ft_bounds)

    # Add like this to preserve errors
    tmp_histo = data_histo.Clone('tmp_histo')
    tmp_histo.Add(subtract_histo, -1)

    # Fill the FT histogram
    total = 0
    for i in range(1, n_bins + 1):
        d = tmp_histo.GetBinContent(i)
        total += d
        e = tmp_histo.GetBinError(i)
        ft_histo.SetBinContent(n_bins + 1 - i, d)
        ft_histo.SetBinError(n_bins + 1 - i, e)        
    ft_histo.SetEntries(total)

    pad = ft_canvas.cd(2)
    pad.SetLogy()
    subtract_histo.Draw()
    pad = ft_canvas.cd(3)
    pad.SetLogy()
    ft_histo.Draw()

    __keep.append(ft_histo)
    __keep.append(data_histo)
    __keep.append(subtract_histo)

    # Calculate the dilution using Wouter's macro
    from ROOT import sigmaFromFT
    D = sigmaFromFT(ft_histo, 17.68, 0.024)
    return D

def signal_dilution(data, sigmat, calibration = None):
    if calibration:
        assert(sigmat.GetName() in [p.GetName() for p in calibration.getVariables()])
        calibration.redirectServers(data.get())

    dms = 17.7 ** 2 / 2
    total = 0
    res_var = data.get().find(sigmat.GetName())
    weighted = data.isWeighted()

    values = []
    total = 0
    for i in range(data.numEntries()):
        r = data.get(i)
        ## External weight
        if weighted:
            w = data.weight()
        else:
            w = 1.

        if calibration:
            res = calibration.getVal()
        else:
            res = res_var.getVal()
        total += w * (exp(- dms * (res ** 2)) ** 2)
    total = sqrt(total / data.sumEntries())
    eff_res = sqrt(-log(total) / dms)
    print 'Dilution = %f' % total
    print 'Effective resolution = %f' % eff_res

    return total, eff_res
