from array import array
# Calculate dilution

__keep = []
def dilution(t_diff, diff_pdf, result, data, signal, wpv):
    from ROOT import RooBinning
    n_bins = 512
    dilution_bounds = array('d', (-5 + i * 5 / float(n_bins) for i in range(n_bins + 1)))
    dilution_bounds.append(10)
    dilution_binning = RooBinning(len(dilution_bounds) - 1, dilution_bounds)
    dilution_binning.SetName('dilution_binning')
    t_diff.setBinning(dilution_binning)

    # Create a histogram of our WPV component
    from ROOT import RooFit
    wpv_histo = diff_pdf.createHistogram('wpv_histo', t_diff._target_(),
                                         RooFit.Binning(dilution_binning))
    
    from ROOT import TH1D
    data_histo = TH1D('data_histo', 'data_histo', len(dilution_bounds) - 1, dilution_bounds)
    time_var = data.get().find('time_diff')
    for i in range(data.numEntries()):
        r = data.get(i)
        value = time_var.getVal()
        if value > 0.: continue
        r = data_histo.Fill(value)

    # Calculate appropriate scale factor. Scale histograms such that the ration
    # of their integrals match the ratio of wpv / signal yields
    sig_yield = signal.getYield().GetName()
    n_signal = result.floatParsFinal().find(sig_yield).getVal()
    wpv_yield = wpv.getYield().GetName()
    n_wpv = result.floatParsFinal().find(wpv_yield).getVal()
    scale = data_histo.Integral(1, 512) * n_wpv / (wpv_histo.Integral(1, 512) * n_signal)
    wpv_histo.Scale(scale)

    # Create the histogram to be transformed
    ft_bounds = array('d', (0. + i * 5. / n_bins for i in range(n_bins + 1)))
    ft_histo = TH1D('ft_histo', 'ft_histo', len(ft_bounds) - 1, ft_bounds)

    # Fill the FT histogram
    s = 0
    for i in range(1, n_bins + 1):
        d = data_histo.GetBinContent(i)
        wpv = wpv_histo.GetBinContent(i)
        diff = d - wpv
        if diff < 0: continue
        s += diff
        ft_histo.SetBinContent(n_bins + 1 - i, diff if diff > 0 else 0)
    ft_histo.SetEntries(s)

    from ROOT import TCanvas
    ft_canvas = TCanvas('ft_canvas', 'ft_canvas', 500, 500)
    ft_canvas.SetLogy()
    ft_histo.Draw()

    __keep.append(ft_canvas)
    __keep.append(ft_histo)

    # Calculate the dilution using Wouter's macro
    from ROOT import sigmaFromFT
    D = sigmaFromFT(ft_histo, 17.7)
    return D

