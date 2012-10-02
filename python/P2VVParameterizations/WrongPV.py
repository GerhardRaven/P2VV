from RooFitWrappers import *

from P2VVParameterizations.MassPDFs import Signal_PsiMass as PsiMassPdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as BMassPdf
from P2VVParameterizations.MassPDFs import LP2011_Background_Mass as BBkgPdf
from P2VVParameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf

class ShapeBuilder(object):
    __weights = set(('jpsi', 'B', 'both'))

    def __init__(self, time, masses, InputFile = "/stuff/PhD/mixing/Bs2JpsiPhiPrescaled.root",
                 WorkSpace = 'Bs2JpsiPhiPrescaled_workspace', Data = 'data', KeysPdf = False,
                 Weights = 'B', Draw = False):
        assert(Weights in ShapeBuilder.__weights)
        self.__weights = Weights

        self.__input_ws = None
        self.__ws = RooObject().ws()
        self.__shapes = {}
        self.__masses = masses
        
        self._sig = Component('wpv_signal', [], Yield = (1000, 100, 50000))
        self._psi = Component('wpv_jpsi',   [], Yield = (5000, 100, 50000))
        self._bkg = Component('wpv_bkg',    [], Yield = (5000, 100, 50000))

            
        if 'B' in masses:
            self._sig_mass = BMassPdf(masses['B'], Name = 'wpv_sig_mass', Prefix = "wpv_")
            self._bkg_mass = BBkgPdf(masses['B'],  Name = 'wpv_bkg_mass', Prefix = "wpv_")
            self._sig[masses['B']] = self._sig_mass.pdf()
            self._psi[masses['B']] = self._bkg_mass.pdf()
            self._bkg[masses['B']] = self._bkg_mass.pdf()
        if 'jpsi' in masses:
            self._sig_mpsi = PsiMassPdf(masses['jpsi'], Name = 'wpv_sig_mpsi', Prefix = "wpv_")
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
        
        if WorkSpace:
            self.__input_ws = input_file.Get(WorkSpace)
            if not self.__input_ws:
                raise RuntimeError
            self._data = self.__input_ws.data(Data)
            if not self._data:
                raise RuntimeError
        else:
            self._data = input_file.Get(Data)
        # self._data = self._data.reduce("mass > 5348 && mass < 5388")
        fitOpts = dict(NumCPU = 4, Save = True, Minimizer = 'Minuit2', Optimize = 2)
        self.__result = self.__pdf.fitTo(self._data, **fitOpts)

        if Draw:
            self.__draw()
            
        from P2VVGeneralUtils import SData
        for p in self.__pdf.Parameters(): p.setConstant(not p.getAttribute('Yield'))
        splot = SData(Pdf = self.__pdf, Data = self._data, Name = 'MixingMassSplot')
        self.__sdatas = {}
        for key, c in self.__components[Weights].iteritems():
            sdata = splot.data(c.GetName())
            self.__sdatas[c] = sdata
            self.__ws.put(sdata)

        self.__shapes = {}
        from ROOT import RooKeysPdf
        for c, sdata in self.__sdatas.iteritems():
            if KeysPdf:
                shape = Pdf(Name = 'wpv_%s_pdf' % c.GetName(), Type = RooKeysPdf,
                            Parameters = (time, sdata))
            else:
                shape = HistPdf(Name = 'wpv_%s_pdf' % c.GetName(), Observables = [time],
                                Data = sdata, Binning = {time : 50})
            self.__shapes[c] = shape

    def __draw(self):
        from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
        from ROOT import TCanvas
        obs = self.__masses.values()
        self.__canvas = TCanvas('wpv_canvas', 'wpv_canvas', 500, 500)
        for (p,o) in zip(wpv_canvas.pads(len(obs)), obs):
            from P2VVGeneralUtils import plot
            pdfOpts  = dict()
            plot(p, o, pdf = self.__pdf, data = self._data
                 , dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack)
                 , pdfOpts  = dict(LineWidth = 2, **pdfOpts)
                 , plotResidHist = True
                 , components = { 'bkg_*'   : dict( LineColor = kRed,   LineStyle = kDashed )
                                  , 'psi_*' : dict( LineColor = kGreen, LineStyle = kDashed )
                                  , 'sig_*' : dict( LineColor = kBlue,  LineStyle = kDashed )
                                  }
                 )

    def sdata(self, key):
        c = self.__components[self.__weights][key]
        return self.__sdatas[c]

    def shape(self, key):
        c = self.__components[self.__weights][key]
        return self.__shapes[c]

    def input_ws(self):
        return self.__input_ws
