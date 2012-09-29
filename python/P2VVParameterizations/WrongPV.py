from RooFitWrappers import *

from P2VVParameterizations.MassPDFs import Signal_PsiMass as PsiMassPdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as BMassPdf
from P2VVParameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf

class ShapeBuilder(object):
    __weights = set(('jpsi', 'B', 'both'))

    def __init__(self, time, mass, InputFile = "/stuff/PhD/mixing/Bs2JpsiPhiPrescaled.root",
                 WorkSpace = 'Bs2JpsiPhiPrescaled_workspace', Data = 'data', KeysPdf = False):

        self.__ws = RooObject().ws()
        self.__shapes = {}

        self._sig_mpsi = PsiMassPdf(mass, Name = 'wpv_sig_mpsi', Prefix = "wpv_")
        self._bkg_mpsi = PsiBkgPdf(mass, Name = 'wpv_bkg_mpsi', Prefix = "wpv_")

        self._psi = Component('jpsi', (self._sig_mpsi.pdf(),), Yield= (20000, 500, 50000))
        self._bkg = Component('bkg', (self._bkg_mpsi.pdf(),), Yield = (20000, 100, 50000) )

        self.__components = dict(jpsi = self._psi, bkg = self._bkg)
        self._pdf = buildPdf(self.__components.values(), Observables = [mass], Name='wpv_mass_pdf')
        self._pdf.Print("t")

        from ROOT import TFile
        input_file = TFile.Open(InputFile)
        if not input_file or not input_file.IsOpen():
            raise OSError
        
        if WorkSpace:
            w = input_file.Get(WorkSpace)
            if not w:
                raise RuntimeError
            self._data = w.data(Data)
            if not self._data:
                raise RuntimeError
        else:
            self._data = input_file.Get(Data)
        self_data = self._data.reduce("mass > 5348 && mass < 5388")

        fitOpts = dict(NumCPU = 4, Save = True, Minimizer = 'Minuit2', Optimize = 2)
        self._pdf.fitTo(self._data, **fitOpts)

        from P2VVGeneralUtils import SData
        for p in self._pdf.Parameters(): p.setConstant(not p.getAttribute('Yield'))
        splot = SData(Pdf = self._pdf, Data = self._data, Name = 'MixingMassSplot')
        self.__sdatas = {}
        for key, c in self.__components.iteritems():
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

    def sdata(self, key):
        c = self.__components[key]
        return self.__sdatas[c]

    def shape(self, key):
        c = self.__components[key]
        return self.__shapes[c]
