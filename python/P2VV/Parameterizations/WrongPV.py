from P2VV.RooFitWrappers import *

from P2VV.Parameterizations.MassPDFs import Signal_PsiMass as PsiMassPdf
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as BMassPdf
from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BBkgPdf
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf

class ShapeBuilder(object):
    __weights = set(('jpsi', 'B', 'both'))
    __rho = dict(B = 2., jpsi = 2.)
    
    def __init__(self, time, masses, sigmat = None, t_diff = None,
                 InputFile = "/bfys/raaij/p2vv/data/Bs2JpsiPhiPrescaled_2011.root",
                 Workspace = 'Bs2JpsiPhiPrescaled_2011_workspace', Data = 'data',
                 UseKeysPdf = False, Weights = 'B', Draw = False, Reweigh = {}):
        assert(Weights in ShapeBuilder.__weights)
        self.__weights = Weights
        self.__time = time
        self.__st = sigmat
        self.__t_diff = t_diff

        self.__input_ws = None
        self.__ws = RooObject().ws()
        self.__shapes = {}
        self.__diff__shapes = {}

        self.__masses = masses

        self._sig = Component('wpv_signal', [], Yield = (1000, 10, 50000))
        self._psi = Component('wpv_jpsi',   [], Yield = (5000, 10, 50000))
        self._bkg = Component('wpv_bkg',    [], Yield = (5000, 10, 50000))

        if 'B' in masses:
            m_sig_mean  = RealVar('wpv_m_sig_mean',   Unit = 'MeV', Value = 5365, MinMax = (5363, 5372))
            m_sig_sigma = RealVar('wpv_m_sig_sigma',  Unit = 'MeV', Value = 10, MinMax = (1, 20))
            from ROOT import RooGaussian as Gaussian
            self._sig_mass = Pdf(Name = 'wpv_sig_m', Type = Gaussian, Parameters = (masses['B'], m_sig_mean, m_sig_sigma ))
            ## self._sig_mass = BMassPdf(masses['B'], Name = 'wpv_sig_mass', Prefix = "wpv_")
            self._bkg_mass = BBkgPdf(masses['B'],  Name = 'wpv_bkg_mass', Prefix = "wpv_",
                                     wpv_m_bkg_exp = dict(Name = 'wpv_m_bkg_exp', Value = -0.0017, MinMax = (-0.01, -0.00001)))
            self._sig[masses['B']] = self._sig_mass
            self._psi[masses['B']] = self._bkg_mass.pdf()
            self._bkg[masses['B']] = self._bkg_mass.pdf()
        if 'jpsi' in masses:
            self._sig_mpsi = PsiMassPdf(masses['jpsi'], Name = 'wpv_sig_mpsi', Prefix = "wpv_",
                                        wpv_mpsi_n = dict(Name = 'wpv_mpsi_n', Value = 1,
                                                          Constant = True ),
                                        wpv_mpsi_alpha = dict(Name = 'wpv_mpsi_alpha', Value = 3,
                                                          MinMax = (1, 5)))
            self._bkg_mpsi = PsiBkgPdf(masses['jpsi'], Name = 'wpv_bkg_mpsi', Prefix = "wpv_")
            self._sig[masses['jpsi']] = self._sig_mpsi.pdf()
            self._psi[masses['jpsi']] = self._sig_mpsi.pdf()
            self._bkg[masses['jpsi']] = self._bkg_mpsi.pdf()

        self.__components = {'jpsi' : dict(jpsi = self._psi, bkg = self._bkg),
                             'B'    : dict(B = self._sig, bkg = self._bkg),
                             'both' : dict(B = self._sig, jpsi = self._psi, bkg = self._bkg)}
        self.__pdf = buildPdf(self.__components[Weights].values(), Observables = masses.values(),
                             Name = 'wpv_mass_pdf')
        self.__pdf.Print("t")

        from ROOT import TFile
        input_file = TFile.Open(InputFile)
        if not input_file or not input_file.IsOpen():
            raise OSError
        
        if Workspace:
            self.__input_ws = input_file.Get(Workspace)
            if not self.__input_ws:
                raise RuntimeError
            self._data = self.__input_ws.data(Data)
            if not self._data:
                raise RuntimeError
        else:
            self._data = input_file.Get(Data)

        self._data = self._data.reduce("{0} > {1} && {0} < {2}".format(time.GetName(), time.getMin(),
                                                                       time.getMax()))

        # self._data = self._data.reduce("mass > 5348 && mass < 5388")
        fitOpts = dict(NumCPU = 4, Save = True, Minimizer = 'Minuit2', Optimize = 2)
        self.__result = self.__pdf.fitTo(self._data, **fitOpts)
            
        from P2VV.GeneralUtils import SData
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
        obs = self.__masses.values()
        self.__canvas = TCanvas('wpv_canvas', 'wpv_canvas', len(obs) * 500, 650)
        for (p,o) in zip(self.__canvas.pads(len(obs)), obs):
            from P2VV.GeneralUtils import plot
            plot(p, o, pdf = self.__pdf, data = self._data
                 , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
                 , pdfOpts  = dict(LineWidth = 2)
                 , plotResidHist = True
                 , components = { 'wpv_bkg_*'   : dict( LineColor = kRed,   LineStyle = kDashed )
                                  , 'wpv_psi_*' : dict( LineColor = kGreen, LineStyle = kDashed )
                                  , 'wpv_sig_*' : dict( LineColor = kBlue,  LineStyle = kDashed )
                                  }
                 )
        if self.__st:
            self.__draw_time()

        if self.__t_diff:
            self.__draw_t_diff()

    def __draw_time(self):
        from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
        from ROOT import TCanvas
        self.__time_canvas = TCanvas('wpv_time_canvas', 'wpv_time_canvas', len(self.__shapes) * 500, 500)
        t = self.__time
        st = self.__st
        from ROOT import RooArgSet
        for (p, (c, shape)) in zip(self.__time_canvas.pads(len(self.__shapes)), self.__shapes.items()):
            from P2VV.GeneralUtils import plot
            pdfOpts  = dict(ProjWData = (RooArgSet(st), self.__sdatas[c], True))
            plot(p, t, pdf = shape, data = self.__sdatas[c]
                 , frameOpts = dict(Title = c.GetName())
                 , dataOpts = dict(MarkerSize = 0.8, Binning = 80, MarkerColor = kBlack)
                 , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
                 , logy = False
                 , plotResidHist = False)

    def __draw_t_diff(self):
        from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
        from ROOT import TCanvas
        self.__diff_canvas = TCanvas('wpv_diff_canvas', 'wpv_diff_canvas', len(self.__diff_shapes) * 500, 500)
        t_diff = self.__t_diff
        st = self.__st
        from ROOT import RooArgSet
        for (p, (c, shape)) in zip(self.__diff_canvas.pads(len(self.__diff_shapes)), self.__diff_shapes.items()):
            from P2VV.GeneralUtils import plot
            pdfOpts  = dict(ProjWData = (RooArgSet(st), self.__sdatas[c], True))
            plot(p, t_diff, pdf = shape, data = self.__sdatas[c]
                 , dataOpts = dict(MarkerSize = 0.8, Binning = 80, MarkerColor = kBlack)
                 , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
                 , logy = False
                 , plotResidHist = False)

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
