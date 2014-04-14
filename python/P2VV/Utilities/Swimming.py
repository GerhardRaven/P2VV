def swimming_histos(bins, t, time_pdf, data, components, debug = False):
    from ROOT import RooArgSet
    iset = RooArgSet(t)
    lifetime_weights = []
    for i in range(1, len(bins)):
        t_min = bins[i - 1]
        t_max = bins[i]
        range_name = 'range_%d' % i
        t.setRange(range_name, (t_min, t_max))
        integral = time_pdf.createIntegral(iset, NormSet = iset, Range = range_name)
        lifetime_weights.append(integral.getVal())

    from ROOT import TH1D
    efficiencies = {}
    denom = {}

    for c in components:
        # Swimming acceptance
        name = 'swimming_acceptance_' + c.GetName()
        h = TH1D(name, name, len(bins) - 1, bins)
        h.SetDirectory(0)
        efficiencies[c] = h

        name = name + '_denom'
        h = TH1D(name, name, len(bins) - 1, bins)
        h.SetDirectory(0)
        denom[c] = h

    for comp in components:
        ranges = []
        for o in data.get():
            if o.GetName().startswith('tp'):
                ranges.append(o)
        ranges = sorted(ranges, key = lambda r: int(r.GetName()[2:]))

        # get the 
        t_data = data.get().find(t.GetName())
        weight_var = data.get().find(comp['Yield'].GetName() + '_sw')
        histo = efficiencies[comp]
        denom_histo = denom[comp]
        n_bins = histo.GetNbinsX()

        for i in range(data.numEntries()):
            data.get(i)
            b = 1
            if weight_var:
                weight = weight_var.getVal()
            else:
                weight = 1.

            # Correct for probability that a particle with time t ends up in a bin.
            t_bin = histo.FindBin(t_data.getVal())
            if t_bin > n_bins:
                continue
            lweight = lifetime_weights[t_bin - 1]

            if debug:
                print ''
                for r in ranges:
                    print r.getVal()
                print ''
            
            # Fill Swimming histograms
            interval_start = 0
            while b <= n_bins:
                lb = histo.GetBinLowEdge(b)
                rb = lb + histo.GetBinWidth(b)
                interval = interval_start
                f = 0.
                while interval < len(ranges) / 2:
                    li = ranges[2 * interval].getVal()
                    ri = ranges[2 * interval + 1].getVal()
                    if li > rb:
                        # bin is completely left of interval
                        if debug:
                            print '{0:8.4f}, {1:8.4f}, {2:8.4f}, {3:8.4f}, {4:4d}, {5:4d}, {6:8.4f}, {7:8.4f}'.format(lb, rb, li, ri, b, interval, f, weight)
                        break
                    elif ri < lb:
                        # bin is completely right of interval
                        pass
                    elif lb > li and rb <= ri:
                        # bin is inside interval
                        f = 1.
                        if debug:
                            print '{0:8.4f}, {1:8.4f}, {2:8.4f}, {3:8.4f}, {4:4d}, {5:4d}, {6:8.4f}, {7:8.4f}'.format(lb, rb, li, ri, b, interval, f, weight)
                        break
                    elif lb > li and rb > ri:
                        f += (ri - lb) / (rb - lb)
                    elif lb <= li and rb > ri:
                        f += (rb - li) / (rb - lb)
                    elif lb <= li and rb <= ri:
                        f += (rb - li) / (rb - lb)
                    if lb > ri:
                        interval_start = interval
                    if debug:
                        print '{0:8.4f}, {1:8.4f}, {2:8.4f}, {3:8.4f}, {4:4d}, {5:4d}, {6:8.4f}, {7:8.4f}'.format(lb, rb, li, ri, b, interval, f, weight)
                    interval += 1

                histo.Fill((lb + rb) / 2., weight * f * lweight)
                b += 1

            for i in range(1, n_bins + 1):
                c = denom_histo.GetBinCenter(i)
                denom_histo.Fill(c, weight * lweight)

    histos = {}
    for i, (c, h) in enumerate(efficiencies.iteritems()):
        # Draw swimming acceptance
        h = h.Clone(h.GetName() + '_clone')
        h.SetDirectory(0)
        d = denom[c]
        h.Divide(d)
        h.SetTitle(c.GetName())
        h.GetYaxis().SetRangeUser(0., 1.0)
        h.GetYaxis().SetTitle("#varepsilon")
        h.GetXaxis().SetTitle("t [ps]")
        h.SetLineWidth(3)
        histos[c] = h
    return histos
