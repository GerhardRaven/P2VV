from array import array
from math import exp, log, sqrt
from RooFitWrappers import FormulaVar

# Calculate dilution
__keep = []
def dilution(t_diff, data, sigmat = None, result = None, signal = [], subtract = [], raw = False, simultaneous = False, calibration = None):
    if calibration:
        calibration.redirectServers(data.get())

    assert(t_diff.getMin() < 0)

    from ROOT import RooBinning
    n_bins = 512
    neg_range = 0 - t_diff.getMin()
    diff = neg_range / float(n_bins)
    dilution_bounds = array('d', (t_diff.getMin() + i * diff for i in range(n_bins + 1)))
    b = diff
    t_max = t_diff.getMax()
    while b < t_max:
        dilution_bounds.append(b)
        b += diff
    if dilution_bounds[-1] != t_max:
        dilution_bounds.append(t_max)
    dilution_binning = RooBinning(len(dilution_bounds) - 1, dilution_bounds)
    dilution_binning.SetName('dilution_binning')
    ## t_diff.setBinning(dilution_binning)
    
    from ROOT import TH1D
    data_histo = TH1D('data_histo', 'data_histo', len(dilution_bounds) - 1, dilution_bounds)
    time_var = data.get().find(t_diff.GetName())
    if sigmat:
        res_var = data.get().find(sigmat.GetName())
    else:
        res_var = None

    dms = -17.7 ** 2 / 2
    
    weighted = data.isWeighted()
    for i in range(data.numEntries()):
        r = data.get(i)
        value = time_var.getVal()
        ## External weight
        if weighted:
            w = data.weight()
        else:
            w = 1.

        ## Weight using sigmat
        if res_var:
            if calibration:
                res = calibration.getVal()
            else:
                res = res_var.getVal()
            w *= (exp(dms * res ** 2) ** 2)

        r = data_histo.Fill(value, w)

    if not subtract or raw:
        # Calculate the dilution using Wouter's macro
        from ROOT import sigmaFromFT
        D = sigmaFromFT(data_histo, 17.7)
        return D

    # Create a histogram of our WPV component
    from ROOT import RooFit
    from P2VV.RooFitWrappers import buildPdf

    signal_yields = [p for p in result.floatParsFinal() if any([p.GetName().startswith(s.getYield().GetName()) for s in signal])]
    data_int = data_histo.Integral()

    # Set the correct yields in the pdf
    # Build the PDF used for subtraction.
    subtract_pdf = buildPdf(Components = subtract, Observables = (t_diff,),
                            Name='subtract_%s_pdf' % '_'.join(c.GetName() for c in subtract))    

    subtract_yields = {}
    for b in subtract:
        subtract_yields[b] = [p for p in result.floatParsFinal() if p.GetName().startswith(b.getYield().GetName())]

    for c in subtract:
        y = c.getYield()
        result_yield = sum([cy.getVal() for cy in subtract_yields[c]])
        pdf_yield = subtract_pdf.getVariables().find(y.GetName())
        if not pdf_yield:
            break
        pdf_yield.setVal(result_yield)
        
    subtract_histo = subtract_pdf.createHistogram(subtract_pdf.GetName().replace('_pdf', '_histo'), t_diff._target_(),
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

    # Fill the FT histogram
    s = 0
    for i in range(1, n_bins + 1):
        d = data_histo.GetBinContent(i)
        sub = subtract_histo.GetBinContent(i)
        diff = d - sub
        s += diff
        ft_histo.SetBinContent(n_bins + 1 - i, diff)
    ft_histo.SetEntries(s)

    from ROOT import TCanvas
    ft_canvas = TCanvas('ft_canvas', 'ft_canvas', 1000, 1000)
    ft_canvas.Divide(2, 2)
    pad = ft_canvas.cd(1)
    pad.SetLogy()
    data_histo.Draw()
    pad = ft_canvas.cd(2)
    pad.SetLogy()
    subtract_histo.Draw()
    pad = ft_canvas.cd(3)
    pad.SetLogy()
    ft_histo.Draw()

    __keep.append(ft_canvas)
    __keep.append(ft_histo)
    __keep.append(data_histo)
    __keep.append(subtract_histo)

    # Calculate the dilution using Wouter's macro
    from ROOT import sigmaFromFT
    D = sigmaFromFT(ft_histo, 17.7)
    return D

def signal_dilution(data, sigmat, calibration = None):
    if calibration:
        calibration.redirectServers(data.get())

    dms = 17.7 ** 2 / 2
    total = 0
    res_var = data.get().find(sigmat.GetName())
    weighted = data.isWeighted()

    values = []
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
        values.append((res, w))
        total += w * exp(- dms * res ** 2)
    total /= data.sumEntries()
    eff_res = sqrt(-log(total) / dms)
    print 'Dilution = %f' % total
    print 'Effective resolution = %f' % eff_res

    return total, eff_res, values
