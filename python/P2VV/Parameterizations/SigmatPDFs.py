from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin
class SigmatPdf(_util_parse_mixin):
    def __init__(self, **kwargs):
        self._pdf = kwargs.pop('pdf')
        for (k,v) in kwargs.iteritems():
            setattr(self, '_' + k, v)

    def pdf(self):
        return self._pdf

class DoubleLogNormal(SigmatPdf):
    # Double LogNormal PDF
    def __init__(self, sigmat, **kwargs):
        prefix = kwargs.pop( 'Prefix', '' )
        self._median = self._parseArg( '%smedian' % prefix, kwargs, Title = 'median', Unit = 'ps'
                                      , Value = 0.03, Error = 0.003, MinMax = (0.00011, 0.12))

        self._k1 = self._parseArg('%sk1' % prefix, kwargs, Unit = '', Value = 1.5, MinMax = (0.00001, 10))
        self._k2 = self._parseArg('%sk2' % prefix, kwargs, Unit = '', Value = 1.5, MinMax = (0.00001, 10))
        self._frac = self._parseArg('%sfrac_ln2' % prefix, kwargs, Value = 0.5, MinMax = (0.01, 0.99))

        from P2VV.RooFitWrappers import LognormalPdf, SumPdf
        self._ln1 = LognormalPdf('%sln1' % prefix, Observable = sigmat, Median = self._median, Shape = self._k1)
        self._ln2 = LognormalPdf('%sln2' % prefix, Observable = sigmat, Median = self._median, Shape = self._k2)

        # Do our own sum pdf to have a fraction
        SigmatPdf.__init__(self, pdf = SumPdf(Name = '%sln' % prefix, PDFs = [self._ln1, self._ln2],
                                              Yields = {'%sln2' % prefix : self._frac}))
