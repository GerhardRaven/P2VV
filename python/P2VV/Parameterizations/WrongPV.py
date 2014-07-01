import os
from P2VV.RooFitWrappers import *

from P2VV.Parameterizations.MassPDFs import DoubleCB_Psi_Mass as PsiMassPdf
from P2VV.Parameterizations.MassPDFs import DoubleGauss_Signal_Mass as BMassPdf
from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BBkgPdf
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf

class ShapeBuilder(object):
    __weights = set(('jpsi', 'B', 'both'))
    __rho = dict(B = 2., jpsi = 2.)

    def __init__(self, time, masses, t_diff = None, MassResult = None,
                 InputFile = "/bfys/raaij/p2vv/data/Bs2JpsiPhiPrescaled_2011.root",
                 Workspace = 'Bs2JpsiPhiPrescaled_2011_workspace', Data = 'data',
                 UseKeysPdf = False, Weights = 'B', Draw = False, Reweigh = {}):
        assert(Weights in ShapeBuilder.__weights)
        self.__weights = Weights
        self.__time = time
        self.__t_diff = t_diff
        pos = Workspace.find('201')
        self.__year = Workspace[pos : pos + 4]
        self.__input_ws = None
        self.__ws = RooObject().ws()
        self.__shapes = {}
        self.__diff__shapes = {}

        self.__masses = masses

        self._sig = Component('wpv_signal', [], Yield = (1000, 10, 5e5))
        self._psi = Component('wpv_jpsi',   [], Yield = (5000, 10, 5e5))
        self._bkg = Component('wpv_bkg',    [], Yield = (5000, 10, 5e5))

        if 'B' in masses:
            ## m_sig_mean  = RealVar('wpv_m_sig_mean',   Unit = 'MeV', Value = 5365, MinMax = (5363, 5372))
            ## m_sig_sigma = RealVar('wpv_m_sig_sigma',  Unit = 'MeV', Value = 10, MinMax = (1, 20))
            ## from ROOT import RooGaussian as Gaussian
            ## self._sig_mass = Pdf(Name = 'wpv_sig_m', Type = Gaussian, Parameters = (masses['B'], m_sig_mean, m_sig_sigma ))
            self._sig_mass = BMassPdf(masses['B'], Name = 'wpv_sig_mass', ParNamePrefix = "wpv",
                                      AvSigParameterisation = True)
            self._bkg_mass = BBkgPdf(masses['B'],  Name = 'wpv_bkg_mass', ParNamePrefix = "wpv",
                                     m_bkg_exp = dict(Name = 'm_bkg_exp', Value = -0.0017,
                                                      MinMax = (-0.01, -0.00001)))
            self._sig[masses['B']] = self._sig_mass.pdf()
            self._psi[masses['B']] = self._bkg_mass.pdf()
            self._bkg[masses['B']] = self._bkg_mass.pdf()
        if 'jpsi' in masses:
            self._sig_mpsi = PsiMassPdf(masses['jpsi'], Name = 'wpv_sig_mpsi', ParNamePrefix = "wpv")
            self._bkg_mpsi = PsiBkgPdf(masses['jpsi'], Name = 'wpv_bkg_mpsi', ParNamePrefix = "wpv")
            self._sig[masses['jpsi']] = self._sig_mpsi.pdf()
            self._psi[masses['jpsi']] = self._sig_mpsi.pdf()
            self._bkg[masses['jpsi']] = self._bkg_mpsi.pdf()

        self.__components = {'jpsi' : dict(jpsi = self._psi, bkg = self._bkg),
                             'B'    : dict(B = self._sig, bkg = self._bkg),
                             'both' : dict(B = self._sig, jpsi = self._psi, bkg = self._bkg)}
        self.__pdf = buildPdf(self.__components[Weights].values(), Observables = masses.values(),
                             Name = 'wpv_mass_pdf')
        if MassResult:
            ## Use the provided mass result to set all the parameter values, only float the yields
            pdf_params = self.__pdf.getParameters(RooArgSet(*masses.values()))
            for p in MassResult.floatParsFinal():
                ## ignore yields
                if any((p.GetName().startswith(n) for n in ['N_', 'mpsi_c'])):
                    continue
                ## Find pdf parameter, add "wpv_" prefix
                pdf_p = pdf_params.find('wpv_' + p.GetName())
                if pdf_p:
                    pdf_p.setVal(p.getVal())
                    pdf_p.setError(p.getError())
                    pdf_p.setConstant(True)

        self.__pdf.Print("t")

        from ROOT import TFile
        input_file = TFile.Open(InputFile)
        if not input_file or not input_file.IsOpen():
            raise OSError

        if Workspace:
            self.__input_ws = input_file.Get(Workspace)
            if not self.__input_ws:
                print 'Cannot find workspace %s in mixing file.' % Workspace
                raise RuntimeError
            self._data = self.__input_ws.data(Data)
            if not self._data:
                print 'Cannot find data in workspace %s.' % Workspace
                raise RuntimeError
        else:
            self._data = input_file.Get(Data)

        if not self._data.get().find(time.GetName()) and "refit" in time.GetName():
            from ROOT import RooFormulaVar, RooArgList
            def __add_alias(name, obs):
                obs_name = obs.GetName()[:-6]
                do = self._data.get().find(obs_name)
                rf = RooFormulaVar(name, name, "@0", RooArgList(do))
                a = self._data.addColumn(rf)
                a.setMin(obs.getMin())
                a.setMax(obs.getMax())
                return a

            ## Add refit observables
            time = __add_alias("time_refit", time)
            self.__time = time
            self.__time.Print()
            if t_diff:
                t_diff = __add_alias("time_diff_refit", t_diff)
                self.__t_diff = t_diff

        if t_diff:
            self._data = self._data.reduce("{0} > {1} && {0} < {2} && {3} > {4} && {3} < {5}".format(time.GetName(), time.getMin(), time.getMax(), t_diff.GetName(), t_diff.getMin(), t_diff.getMax()))
        else:
            self._data = self._data.reduce("{0} > {1} && {0} < {2}".format(time.GetName(), time.getMin(), time.getMax()))

        # self._data = self._data.reduce("mass > 5348 && mass < 5388")
        fitOpts = dict(NumCPU = 4, Save = True, Minimizer = 'Minuit2', Optimize = 2)
        self.__result = self.__pdf.fitTo(self._data, **fitOpts)

        from P2VV.Utilities.SWeights import SData
        for p in self.__pdf.Parameters(): p.setConstant(not p.getAttribute('Yield'))
        splot = SData(Pdf = self.__pdf, Data = self._data, Name = 'MixingMassSplot')
        self.__sdatas = {}
        self.__reweigh_weights = {}
        for key, c in self.__components[Weights].iteritems():
            sdata = splot.data(c.GetName())

            if 'Data' in Reweigh and key in Reweigh['Data']:
                from array import array
                from ROOT import RooBinning

                source = Reweigh['Data'][key]
                binning = Reweigh['Binning']
                if type(binning) == array:
                    binning = RooBinning(len(binning) - 1, binning)
                    binning.SetName('reweigh')
                    Reweigh['DataVar'].setBinning(binning, 'reweigh')

                source_obs = source.get().find('nPV')
                cat_name = 'nPVs_' + key
                source_cat = source.get().find(cat_name)
                if source_cat:
                    # Remove previous weights to make sure we get it right
                    new_vars = source.get()
                    new_vars.remove(source_cat)
                    source = source.reduce(new_vars)
                    source_obs = source.get().find(source_obs.GetName())

                source_obs.setBinning(binning, 'reweigh')                
                source_cat = BinningCategory(Name = cat_name, Observable = source_obs,
                                             Binning = binning, Data = source, Fundamental = True)
                from P2VV.Reweighing import reweigh
                sdata, weights = reweigh(sdata, sdata.get().find('nPV'),
                                         source, source_cat)
                self.__reweigh_weights[key] = weights

            sdata = self.__ws.put(sdata)
            self.__sdatas[c] = sdata

        rho_keys = dict((v, k) for k, v in self.__components[Weights].iteritems())
        self.__shapes = {}
        self.__diff_shapes = {}
        for c, sdata in self.__sdatas.iteritems():
            if UseKeysPdf:
                rk = rho_keys[c]
                rho = ShapeBuilder.__rho[rk] if rk in ShapeBuilder.__rho else 1.
                time_shape = KeysPdf(Name = 'wpv_%s_pdf' % c.GetName(), Observable = time, Data = sdata, Rho = rho)
                if t_diff:
                    diff_shape = KeysPdf(Name = 'wpv_%s_diff_pdf' % c.GetName(), Observable = t_diff, Data = sdata)

            else:
                time_shape = HistPdf(Name = 'wpv_%s_pdf' % c.GetName(), Observables = [time],
                                Data = sdata, Binning = {time : 35})
                if t_diff:
                    diff_shape = HistPdf(Name = 'wpv_%s_diff_pdf' % c.GetName(), Observables = [t_diff],
                                         Data = sdata, Binning = {time : 35})
            self.__shapes[c] = time_shape
            if t_diff:
                self.__diff_shapes[c] = diff_shape
        if Draw:
            self.__draw()

    def __draw(self):
        from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
        from ROOT import TCanvas
        from ROOT import TPaveText
        self.__year_label = TPaveText(0.71, 0.72, 0.89, 0.85, "NDC")
        self.__year_label.SetFillColor(0)
        self.__year_label.AddText(self.__year)
        self.__year_label.SetBorderSize(0)

        obs = self.__masses.values()
        self.__canvas = TCanvas('wpv_canvas', 'wpv_canvas', len(obs) * 600, 533)
        for (p,o) in zip(self.__canvas.pads(len(obs)), obs):
            from P2VV.Utilities.Plotting import plot
            plot(p, o, pdf = self.__pdf, data = self._data
                 , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = 50)
                 , pdfOpts  = dict(LineWidth = 2)
                 , plotResidHist = 'BX'
                 , xTitle = 'M_{#mu^{+}#mu^{-}} [MeV/c^{2}]'
                 , yTitle = 'Candidates / (%2.1f MeV/c^{2})' % ((o.getMax() - o.getMin()) / float(50))
                 , components = { 'wpv_bkg_*'   : dict( LineColor = kRed,   LineStyle = kDashed )
                                  , 'wpv_psi_*' : dict( LineColor = kGreen, LineStyle = kDashed )
                                  , 'wpv_sig_*' : dict( LineColor = kBlue,  LineStyle = kDashed )
                                  }
                 )
            self.__canvas.cd(1)
            self.__year_label.Draw()
            self.__canvas.Update()
        self.__draw_time()

        if self.__t_diff:
            self.__draw_t_diff()

    def __draw_time(self):
        from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
        from ROOT import TCanvas
        from itertools import product
        self.__time_canvases = []
        t = self.__time
        from P2VV.Utilities.Plotting import plot
        from ROOT import RooArgSet
        for logy, (c, shape) in product((True, False), self.__shapes.items()):
            cn = 'wpv_time_%s' % c.GetName()
            if logy: cn += '_log'
            p = TCanvas(cn, cn, 600, 400)
            self.__time_canvases.append(p)
            nBins = 80
            frame = plot(p, t, pdf = shape, data = self.__sdatas[c]
                         , frameOpts = dict(Title = c.GetName())
                         , dataOpts = dict(MarkerSize = 0.8, Binning = nBins, MarkerColor = kBlack)
                         , pdfOpts  = dict(LineWidth = 2)
                         , logy = logy, yScale = (1 if logy else 0, None)
                         , plotResidHist = False)[0]
            frame.GetXaxis().SetTitle('t [ps]')
            if logy:
                frame.GetYaxis().SetRangeUser(1, 7000 if self.__year == '2012' else 5000)
            frame.GetYaxis().SetTitle('Candidates / (%3.2f ps)' % ((t.getMax() - t.getMin()) / float(nBins)))
            self.__year_label.Draw()
            from P2VV.Utilities.Resolution import plot_dir
            p.Print(os.path.join(plot_dir, 'wpv_time_%s_%s_%s.pdf' % (c.GetName()[4:], self.__year, 'log' if logy else 'linear')), EmbedFonts = True)
            
    def __draw_t_diff(self):
        from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
        from ROOT import TCanvas
        self.__diff_canvases = []
        t_diff = self.__t_diff
        from ROOT import RooArgSet
        for c, shape in self.__diff_shapes.iteritems():
            name = 'wpv_diff_%s' % c.GetName()
            canvas = TCanvas(name, name, 600, 400)
            p = canvas.cd()
            self.__diff_canvases.append(canvas)
            from P2VV.Utilities.Plotting import plot
            plot(p, t_diff, pdf = shape, data = self.__sdatas[c]
                 , dataOpts = dict(MarkerSize = 0.8, Binning = 80, MarkerColor = kBlack)
                 , pdfOpts  = dict(LineWidth = 2)
                 , xTitle = 't_{rec} - t_{true} [ps]'
                 , yTitle = 'Candidates / (12.5 fs)'
                 , logy = False
                 , plotResidHist = False)
            self.__year_label.Draw()
            from P2VV.Utilities.Resolution import plot_dir
            canvas.Print(os.path.join(plot_dir, 'wpv_tdiff_%s_%s_linear.pdf' % (c.GetName()[4:], self.__year)), EmbedFonts = True)

    def sdata(self, key):
        c = self.__components[self.__weights][key]
        return self.__sdatas[c]

    def shape(self, key):
        c = self.__components[self.__weights][key]
        return self.__shapes[c]

    def diff_shape(self, key):
        c = self.__components[self.__weights][key]
        return self.__diff_shapes[c]

    def input_ws(self):
        return self.__input_ws

    def reweigh_weights(self, key):
        try:
            return self.__reweigh_weights[key]
        except AttributeError:
            return {}


def add_alias(name, obs):
    obs_name = obs.GetName()[:-6]
    do = self._data.get().find(name)
    rf = RooFormulaVar(name, name, "@0", RooArgList(do))
    a = self._data.addColumn(rf)
    a.setMin(obs.getMin())
    a.setMax(obs.getMax())
    return a
