from array import array
from math import exp, log, sqrt

__keep = []
__histos = []

def sigmaFromFT(hist, dMs, dMsErr, quiet = False):
    from math import cos, sqrt, log, sin

    s = 0
    sdt = 0
    sdterr = 0
    scos = 0
    sx2 = 0
    error = 0
    for i in range(1, hist.GetNbinsX() + 1):
        c = hist.GetBinContent(i)
        e = hist.GetBinError(i)
        x = hist.GetBinCenter(i)
        dt = hist.GetBinWidth(i)
        s += c
        sdt += c * dt
        sdterr += e * e * dt * dt

        cosv = cos( - dMs * x )
        scos += c * cosv * dt
        sx2 += c * x * x

        coserr = sin(dMsErr) * dMsErr
        error += (c * c * coserr * coserr +  cosv * cosv * e * e) * dt * dt

    error = sqrt(error / (sdt * sdt) + sdterr * scos * scos / sdt ** 4)

    try:
        rms = sqrt(sx2/s)
    except ValueError:
        rms = None

    D = scos / sdt
    sigma = sqrt( -2*log(D) ) / dMs

    if not quiet:
        print "%f %f %f" % (sdt, scos, sx2)
        r = dict(rms = rms, gd = exp(-0.5*rms*rms*dMs*dMs), ftd = D, ftde = error, cgd = sigma)
        print """
RMS of input histogram: %(rms)f
If distribution were Gaussian, dilution is: %(gd)f
Dilution from FT: %(ftd)f +- %(ftde)f
Corresponding Gaussian resolution: %(cgd)f""" % r

    return D, error

def sigmaFT(data, dMs, dMsErr, quiet = False):
    from math import cos, sqrt, log, sin
    s = 0
    scos = 0
    sx2 = 0
    for t, w in data_list:
        s += w
        scos += w * cos( - dMs * t )
        sx2 += w * t * t

    try:
        rms = sqrt(sx2/s)
    except ValueError:
        rms = None

    D = scos / s
    sigma = sqrt( -2*log(D) ) / dMs

    if not quiet:
        print "%f %f %f" % (s, scos, sx2)
        r = dict(rms = rms, gd = exp(-0.5*rms*rms*dMs*dMs), ftd = D, cgd = sigma)
        print """
RMS of input histogram: %(rms)f
If distribution were Gaussian, dilution is: %(gd)f
Dilution from FT: %(ftd)f
Corresponding Gaussian resolution: %(cgd)f""" % r

    return D, 0

def __make_ft_histo(data_histo, ft_bounds):
    from ROOT import TH1D
    ft_histo = TH1D('ft_histo', 'ft_histo', len(ft_bounds) - 1, ft_bounds)

    # Fill the FT histogram
    total = 0
    n_bins = len(ft_bounds) - 1
    for i in range(1, n_bins + 1):
        d = data_histo.GetBinContent(i)
        e = data_histo.GetBinError(i)
        total += d
        ft_histo.SetBinContent(n_bins + 1 - i, d)
        ft_histo.SetBinError(n_bins + 1 - i, e)
    ft_histo.SetEntries(total)
    return ft_histo

def dilution_bin(data, t_var, dilution_binning, result = None, bin_name = "",
                 signal = [], subtract_pdf = None, subtract = [], n_bins = 512,
                 st = None, mean_param = None):
    from P2VV.Load import P2VVLibrary
    diff = dilution_binning.binWidth(0)
    ft_left = dilution_binning.lowBound()

    dilution_bounds = array('d', [dilution_binning.binLow(i) for i in range(dilution_binning.numBins())])
    dilution_bounds.append(dilution_binning.highBound())
    from ROOT import TH1D
    data_histo = TH1D('data_histo', 'data_histo', len(dilution_bounds) - 1, dilution_bounds)
    data_histo.Sumw2()
    time_var = data.get().find(t_var.GetName())
    st_var = data.get().find(st.GetName())

    weighted = data.isWeighted()
    for i in range(data.numEntries()):
        r = data.get(i)
        value = time_var.getVal()
        if mean_param:
            value -= mean_param(st_var.getVal())
        ## External weight
        if weighted:
            w = data.weight()
        else:
            w = 1.
        r = data_histo.Fill(value, w)

    if not subtract_pdf:
        ft_bounds = array('d', (ft_left + i * diff for i in range(n_bins + 1)))
        ft_histo = __make_ft_histo(data_histo, ft_bounds)
        return sigmaFromFT(ft_histo, 17.768, 0.024, quiet = True)

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
    ft_bounds = array('d', (ft_left + i * diff for i in range(n_bins + 1)))
    ft_histo = __make_ft_histo(tmp_histo, ft_bounds)
    return sigmaFromFT(ft_histo, 17.768, 0.024, quiet = True)

# Calculate dilution
def dilution_bins(data, t_var, sigmat, sigmat_cat, result = None, signal = [], subtract = [], raw = False, calibration = None, n_bins = 512, mean_param = None, t_range = None):
    if calibration:
        assert(sigmat.GetName() in [p.GetName() for p in calibration.getVariables()])
        calibration.redirectServers(data.get())

    assert(t_var.getMin() < 0)

    from ROOT import RooBinning
    if not t_range:
        ft_left = 0
        t_range = ft_left - t_var.getMin()
        diff = t_range / float(n_bins)
        dilution_bounds = array('d', (t_var.getMin() + i * diff for i in range(n_bins + 1)))
        b = diff
        t_max = t_var.getMax()
        while b < t_max:
            dilution_bounds.append(b)
            b += diff
        if dilution_bounds[-1] != t_max:
            dilution_bounds.append(t_max)
    else:
        diff = t_range / float(n_bins)
        ft_left = t_var.getMin()
        dilution_bounds = array('d', (t_var.getMin() + i * diff for i in range(n_bins + 1)))
    dilution_binning = RooBinning(len(dilution_bounds) - 1, dilution_bounds)
    dilution_binning.SetName('dilution_binning')

    if not subtract:
        dilutions = []
        total = 0
        for ct in sigmat_cat:
            bin_data = data.reduce('{0} == {0}::{1}'.format(sigmat_cat.GetName(), ct.GetName()))
            n = bin_data.sumEntries()
            d = dilution_bin(bin_data, t_var, dilution_binning, n_bins = n_bins,
                             st = sigmat, mean_param = mean_param)
            dilutions.append((n, d))
        dt = sum(d[0] * d[1][0] for d in dilutions) / sum(d[0] for d in dilutions)
        de = sum(d[0] * d[1][1] ** 2 for d in dilutions) / sum(d[0] for d in dilutions)
        return dt, sqrt(de)

    subtract_pdf = buildPdf(Components = subtract, Observables = (t_var,),
                            Name='subtract_%s_pdf' % '_'.join(c.GetName() for c in subtract))
    sub_pars = subtract_pdf.getVariables()
    for res_par in result.floatParsFinal():
        sub_par = sub_pars.find(res_par.GetName())
        if sub_par:
            print 'setting value of subtract pdf parameter %s' % sub_par.GetName()
            sub_par.setVal(res_par.getVal())
            sub_par.setError(res_par.getError())

    # Get the signal yields from the fit result
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
                         signal, subtract_pdf, subtract, n_bins = n_bins)
        total += sum(signal_yields[ct.GetName()]) * d ** 2
    total /= signal_yield
    return total

# Calculate dilution
def dilution_ft(data, t_var, t_range = None, parameters = None, signal = [], subtract = [],
                full_pdf = None, raw = False, simultaneous = False, n_bins = 512, st = None,
                mean_param = None, quiet = False):
    if mean_param: assert(st)

    from ROOT import RooBinning
    if not t_range:
        ft_left = 0
        t_range = ft_left - t_var.getMin()
    else:
        ft_left = t_var.getMin()

    diff = t_range / float(n_bins)
    dilution_bounds = array('d', (t_var.getMin() + i * diff for i in range(n_bins + 1)))

    for i in range(len(dilution_bounds) -1):
        if dilution_bounds[i + 1] <= dilution_bounds[i]:
            print i
    dilution_binning = RooBinning(len(dilution_bounds) - 1, dilution_bounds)
    dilution_binning.SetName('dilution_binning')
    ## t_var.setBinning(dilution_binning)

    from ROOT import TH1D
    data_histo = TH1D('data_histo', 'data_histo', len(dilution_bounds) - 1, dilution_bounds)
    data_histo.Sumw2()
    time_var = data.get().find(t_var.GetName())
    if st:
        st_var = data.get().find(st.GetName())

    weights = []
    weighted = data.isWeighted()
    for i in range(data.numEntries()):
        r = data.get(i)
        value = time_var.getVal()
        if mean_param:
            value -= mean_param(st_var.getVal())
        ## External weight
        if weighted:
            w = data.weight()
        else:
            w = 1.

        weights.append(w)
        r = data_histo.Fill(value, w)

    if not subtract or raw:
        ft_bounds = array('d', (ft_left + i * diff for i in range(n_bins + 1)))
        ft_histo = __make_ft_histo(data_histo, ft_bounds)

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
        from P2VV.Load import P2VVLibrary
        return sigmaFromFT(ft_histo, 17.768, 0.024, quiet)
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

    signal_yields = {}
    for c in signal:
        y = filter(lambda p: p.GetName() == c.getYield().GetName(), [p for p in parameters])[0]
        signal_yields[c] = y
    data_int = data_histo.Integral()

    signal_pdfs = {}
    for c in signal:
        pdf = filter(lambda pdf: pdf.GetName() == c[(t_var,)].GetName(), full_pdf.PDFs())[0]
        y = [p for p in parameters if p.GetName() == c.getYield().GetName()][0]
        signal_pdfs[c] = (y, pdf)
    
    subtract_pdfs = {}
    for c in subtract:
        pdf = filter(lambda pdf: pdf.GetName() == c[(t_var,)].GetName(), full_pdf.PDFs())[0]
        y = [p for p in parameters if p.GetName() == c.getYield().GetName()][0]
        subtract_pdfs[c] = (y, pdf)

    from ROOT import RooArgSet
    t_var.setRange('dilution', (t_var.getMin(), t_var.getMin() + t_range))

    # Calculate appropriate scale factor. Scale histograms such that the ratio
    # of their integrals match the ratio of subtract / total yields
    n_subtract = {}
    for c, (y, pdf) in subtract_pdfs.iteritems():
        I = pdf.createIntegral(RooArgSet(t_var), Range = 'dilution')
        n_subtract[c] = I.getVal() * y.getVal()

    total = sum(n_subtract.itervalues())
    for c, (y, pdf) in signal_pdfs.iteritems():
        I = pdf.createIntegral(RooArgSet(t_var), Range = 'dilution')
        total += I.getVal() * y.getVal()

    # Histogram from which stuff will be subtracted
    tmp_histo = data_histo.Clone('tmp_histo')
    subtract_histo = None

    for c, (y, pdf) in subtract_pdfs.iteritems():
        sub_histo = pdf.createHistogram(pdf.GetName() + '_histo', t_var._target_(),
                                        RooFit.Binning(dilution_binning), RooFit.Scaling(False))
        sub_int = sub_histo.Integral()
        scale = (data_int * n_subtract[c]) / (sub_int * total)
        sub_histo.Scale(scale)
        # Add like this to preserve errors
        if not subtract_histo:
            subtract_histo = sub_histo
        else:
            subtract_histo.Add(sub_histo)
        tmp_histo.Add(sub_histo, -1)
        
    # Fill the FT histogram
    ft_bounds = array('d', (ft_left + i * diff for i in range(n_bins + 1)))
    ft_histo = __make_ft_histo(tmp_histo, ft_bounds)

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
    from P2VV.Load import P2VVLibrary
    return sigmaFromFT(ft_histo, 17.768, 0.024, quiet)

def dilution(data, sfs, calib = None, error_fun = None):
    """
    Calculate the dilution of the data.
    data should be [(st, weight), ...]
    sfs should be [(calib_pars, frac), ..., (calib_pars, None)]
    calib_pars = [par0, par1, ...]
    if calib == None, a polynomial with parameters calib_pars is assumed
    else calib = [fun(pars, st), ...]
    error_fun should calculate the error on D^2 and be callable as (st, w, [list of calibration pars])
    """
    from math import pow
    if not calib:
        def __calib(pars, st):
            return sum(pars[i] * pow(st, i) for i in range(len(pars)))
        calib = [__calib for j in range(len(sfs))]

    if sfs[-1][1] == None:
        sfs[-1][1] = 1 - sum(sf[1] for sf in sfs[:-1])

    ## reset the error functor
    if error_fun:
        error_fun.reset()
    dms = 17.768 ** 2 / 2
    total = 0
    err_2 = 0
    values = []
    from math import exp, sqrt
    sw = 0
    for st, w in data:
        d = 0
        for (pars, frac), cal in zip(sfs, calib):
            d += frac * exp(- dms * cal(pars, st) ** 2)
        total += w * d ** 2
        sw += w
        ## Call error calculator to update sums
        if error_fun:
            error_fun(w, st)

    if error_fun:
        ## Error calculator return error on D^2
        err = 2 * total / sw * error_fun.error()
    else:
        err = None

    total = sqrt(total / sw)
    ## eff_res = sqrt(-log(total) / dms)
    return total, err

def signal_dilution(data, sigmat, calibration = None):
    if calibration:
        assert(sigmat.GetName() in [p.GetName() for p in calibration.getVariables()])
        calibration.redirectServers(data.get())

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
        values.append((res, w))
    return dilution(values, [((0, 1), 1)])

def signal_dilution_sfs(data, sigmat, sfm, frac, sfs):
    from math import sqrt
    sf1 = - sqrt(frac / (1 - frac)) * sfs + sfm
    sf2 = sqrt((1 - frac) / frac) * sfs + sfm
    return signal_dilution_dg(data, sigmat, sf1, frac, sf2)

def signal_dilution_dg(data, sigmat, sf1, frac, sf2):
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
        res = res_var.getVal()
        values.append((res, w))

    return dilution(values, [((0., sf1), 1 - frac), ((0., sf2), frac)])

def signal_dilution_sg(data, sigmat, sf1):
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
        res = res_var.getVal()
        values.append((res, w))

    return dilution(values, [((0., sf1), 1)])

class SolveSF(object):
    def __init__(self, data, sigmat):
        from ROOT import RooDataSet
        self.setData(data, sigmat)
        from ROOT import RooRealVar, RooArgList
        self.__sf = 1.

        self.__min = 0.01
        self.__max = 5

        import ROOT
        from ROOT import TF1
        self.__tf1 = TF1('_tf1', self, self.__min, self.__max, 0)
        self.__wtf1 = ROOT.Math.WrappedTF1(self.__tf1)
        self.__D = 1

    def setData(self, data, sigmat):
        from ROOT import RooDataSet
        if isinstance(data, RooDataSet):
            self.__data = []
            st = data.get().find(sigmat.GetName())
            for i in range(data.numEntries()):
                r = data.get(i)
                self.__data.append((st.getVal(), data.weight()))
        else:
            raise TypeError("Data must be RooDataSet")

    def setSF(self, sf):
        self.__sf = sf

    def __call__(self, x, par = []):
        self.__sf = x[0]
        d = dilution(self.__data, [((0, self.__sf), 1)])
        return d[0] - self.__D

    def solve(self, D ):
        self.__D = D
        import ROOT
        rf = ROOT.Math.BrentRootFinder()
        rf.SetFunction(self.__wtf1, self.__min, self.__max)
        rf.SetNpx(20)
        rf.Solve()
        return rf.Root()


from P2VV.PropagateErrors import Parameter

class DilutionCSFS(object):
    def __init__(self, st_mean, result):
        self.__result = result
        def __mkp(result, name):
            p = result.floatParsFinal().find(name)
            return Parameter(name, p.getVal(), p.getError())

        self.__dms = Parameter('dms', 17.768,  0.024)
        self.__st_mean = st_mean
        self.__sfm_offset = __mkp(self.__result, "sf_mean_slope")
        self.__sfm_slope = __mkp(self.__result, "sf_mean_quad")
        self.__sfm_slope.setValue(self.__sfm_slope.value())
        self.__sfm_slope.setError(self.__sfm_slope.error())
        self.__sfs_offset = __mkp(self.__result, "sf_sigma_slope")
        self.__sfs_slope = __mkp(self.__result, "sf_sigma_quad")
        self.__sfs_slope.setValue(self.__sfs_slope.value())
        self.__sfs_slope.setError(self.__sfs_slope.error())
        self.__frac = __mkp(self.__result, "timeResFrac2")
        from ROOT import RooArgList
        self.__cv = self.__result.reducedCovarianceMatrix(RooArgList(*[self.__result.floatParsFinal().find(p.name()) for p in [self.__sfm_offset, self.__sfm_slope, self.__frac, self.__sfs_offset, self.__sfs_slope]]))

        from P2VV.PropagateErrors import ErrorCSFS
        self.__error = ErrorCSFS(st_mean, self.__dms, self.__sfm_offset, self.__sfm_slope,
                                 self.__frac, self.__sfs_offset, self.__sfs_slope, self.__cv)

    def dilution(self, data):
        d = dilution(data, [([(self.__st_mean, self.__sfm_offset.value(), self.__sfm_slope.value()),
                              (self.__st_mean, self.__sfs_offset.value(), self.__sfs_slope.value()),
                              self.__frac.value()], (1 - self.__frac.value())),
                            ([(self.__st_mean, self.__sfm_offset.value(), self.__sfm_slope.value()),
                              (self.__st_mean, self.__sfs_offset.value(), self.__sfs_slope.value()),
                              self.__frac.value()], self.__frac.value())],
                            (self.calib_sf1, self.calib_sf2), self.__error)
        return d

    def calib_sf1(self, pars, st):
        sfm_pars = pars[0]
        sfs_pars = pars[1]
        f = pars[2]
        sfm = self.calib(sfm_pars, st)
        sfs = self.calib(sfs_pars, st)
        return (- sqrt(f / (1 - f)) * sfs + sfm) * st

    def calib_sf2(self, pars, st):
        sfm_pars = pars[0]
        sfs_pars = pars[1]
        f = pars[2]
        sfm = self.calib(sfm_pars, st)
        sfs = self.calib(sfs_pars, st)
        return (sqrt((1 - f) / f) * sfs + sfm) * st

    def calib(self, pars, st):
        return sum(pars[i + 1] * pow(st - pars[0], i) for i in range(len(pars) - 1))

class DilutionCSFC(object):
    def __init__(self, st_mean, result):
        self.__result = result

        def __mkp(result, name):
            p = result.floatParsFinal().find(name)
            return Parameter(name, p.getVal(), p.getError())

        self.__dms = Parameter('dms', 17.768,  0.024)
        self.__st_mean = st_mean
        self.__sfc_offset = __mkp(self.__result, "sfc_slope")
        self.__sfc_slope = __mkp(self.__result, "sfc_quad")
        self.__frac = __mkp(self.__result, "timeResFrac2")
        self.__sf2_offset = __mkp(self.__result, "sf2_slope")
        self.__sf2_slope = __mkp(self.__result, "sf2_quad")
        from ROOT import RooArgList
        self.__cv = self.__result.reducedCovarianceMatrix(RooArgList(*[self.__result.floatParsFinal().find(p.name()) for p in [self.__sfc_offset, self.__sfc_slope, self.__frac, self.__sf2_offset, self.__sf2_slope]]))

        from P2VV.PropagateErrors import ErrorCDG
        self.__error = ErrorCDG(st_mean, self.__dms, self.__sfc_offset, self.__sfc_slope,
                                self.__frac, self.__sf2_offset, self.__sf2_slope, self.__cv)

    def dilution(self, data):
        d = dilution(data, [([(self.__st_mean, self.__sfc_offset.value(), self.__sfc_slope.value()),
                              (self.__st_mean, self.__sf2_offset.value(), self.__sf2_slope.value()), self.__frac.value()],
                             (1 - self.__frac.value())), ((self.__st_mean, self.__sf2_offset.value(),
                                                           self.__sf2_slope.value()), self.__frac.value())],
                     (self.calib_sf1, self.calib_sf2), self.__error)
        return d

    def calib_sf1(self, pars, st):
        sfc_pars = pars[0]
        sf2_pars = pars[1]
        frac = pars[2]
        sfc = self.calib(sfc_pars, st)
        sf2 = self.calib(sf2_pars, st)
        return (sfc - frac * sf2) / (1 - frac) * st

    def calib_sf2(self, pars, st):
        return self.calib(pars, st) * st

    def calib(self, pars, st):
        return sum(pars[i + 1] * pow(st - pars[0], i) for i in range(len(pars) - 1))

class DilutionSG(object):
    def __init__(self):
        self.__sf = Parameter('sf_sg', 1.45, 0.06)

        self.__dms = Parameter('dms', 17.768,  0.024)
        from P2VV.PropagateErrors import ErrorSG
        self.__error = ErrorSG(self.__dms, self.__sf)

    def dilution(self, data):
        return dilution(data, [((self.__sf,), 1)], (self.calib,), self.__error)

    def calib(self, p, st):
        return p[0].value() * st

class DilutionSFC(object):
    def __init__(self, result):
        def __mkp(result, name):
            p = result.floatParsFinal().find(name)
            return Parameter(name, p.getVal(), p.getError())

        self.__result = result
        self.__dms = Parameter('dms', 17.768,  0.024)
        self.__sfc = __mkp(self.__result, "timeResComb")
        self.__frac = __mkp(self.__result, "timeResFrac2")
        self.__sf2 = __mkp(self.__result, "timeResSigmaSF_2")
        from ROOT import RooArgList
        self.__cv = self.__result.reducedCovarianceMatrix(RooArgList(*[self.__result.floatParsFinal().find(p.name()) for p in [self.__sfc, self.__frac, self.__sf2]]))

        from P2VV.PropagateErrors import ErrorSFC
        self.__error = ErrorSFC(self.__dms, self.__sfc, self.__frac, self.__sf2, self.__cv)

    def dilution(self, data):
        return dilution(data, [((self.__sfc, self.__frac, self.__sf2), (1 - self.__frac.value())),
                               ((0, self.__sf2), self.__frac.value())],
                        (self.sf1, self.calib), self.__error)

    def calib(self, p, st):
        return p[1].value() * st

    def sf1(self, p, st):
        return (p[0].value() - p[1].value() * p[2].value()) / (1 - p[1].value()) * st

from math import sqrt

class DilutionSFS(object):
    def __init__(self, result):
        def __mkp(result, name):
            p = result.floatParsFinal().find(name)
            return Parameter(name, p.getVal(), p.getError())

        self.__result = result
        self.__dms = Parameter('dms', 17.768,  0.024)
        self.__sfc = __mkp(self.__result, "sf_mean_slope")
        self.__frac = __mkp(self.__result, "timeResFrac2")
        self.__sfs = __mkp(self.__result, "sf_sigma_slope")
        from ROOT import RooArgList
        self.__cv = self.__result.reducedCovarianceMatrix(RooArgList(*[self.__result.floatParsFinal().find(p.name()) for p in [self.__sfc, self.__frac, self.__sfs]]))

        from P2VV.PropagateErrors import ErrorSFS
        self.__error = ErrorSFS(self.__dms, self.__sfc, self.__frac, self.__sfs, self.__cv)

    def dilution(self, data):
        return dilution(data, [((self.__sfc, self.__frac, self.__sfs), (1 - self.__frac.value())),
                               ((self.__sfc, self.__frac, self.__sfs), self.__frac.value())],
                        (self.sf1, self.sf2), self.__error)

    def sf1(self, p, st):
        f = p[1].value()
        sfc = p[0].value()
        sfs = p[2].value()
        return (- sqrt(f / (1 - f)) * sfs + sfc) * st

    def sf2(self, p, st):
        f = p[1].value()
        sfc = p[0].value()
        sfs = p[2].value()
        return (sqrt((1 - f) / f) * sfs + sfc) * st
