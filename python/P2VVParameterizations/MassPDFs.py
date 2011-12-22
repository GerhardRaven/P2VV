from P2VVParameterizations.GeneralUtils import _util_parse_mixin

class MassPdf( _util_parse_mixin ) :
    def __init__(self, **kwargs ) :
        for (k,v) in kwargs.iteritems() :
            setattr(self,'_'+k,v)
    def pdf(self) :
        return self._pdf

class LP2011_Signal_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._parseArg('m_sig_mean', kwargs, Title = 'B Mass', Unit = 'MeV/c^2', Value = 5365, MinMax = (5360,5370) )
        self._parseArg('m_sig_sigma_1', kwargs, Title = 'B Mass resolution 1', Unit = 'MeV/c^2', Value = 6, MinMax = (3,20) )
        self._parseArg('m_sig_sigma_2_scale', kwargs, Title = 'B Mass resolution 2 scale factor', Unit = 'MeV/c^2', Value = 2.15, MinMax = (1,5) )
        self._parseArg('m_bkg_f',kwargs, Title = 'B mass fraction 2nd Gaussian', Value = 0.83, MinMax = (0.,1.) )

        from ROOT import RooGaussian as Gaussian
        from RooFitWrappers import Pdf, FormulaVar, SumPdf
        g1 = Pdf( 'm_sig_1', Type = Gaussian
                           , Observables = (mass,)
                           , Parameters = (self._m_sig_mean, self._m_sig_sigma_1 ) 
                           )
        g2 = Pdf( 'm_sig_2', Type = Gaussian
                           , Observables = (mass,)
                           , Parameters = ( self._m_sig_mean
                                          , FormulaVar('m_sig_sigma_2','@0*@1',(self._m_sig_sigma_2_scale, self._m_sig_sigma_1))
                                          )  
                           )
        MassPdf.__init__(self, pdf = SumPdf(Name = kwargs.pop('Name','LP2011_Signal_Mass'),   PDFs = (  g1,  g2)  , Yields = { g1.GetName() : self._m_bkg_f } ) )


class LP2011_Background_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._parseArg('m_bkg_exp', kwargs, Title = 'Mass background slope', Unit = 'MeV/c^2', Value = -0.001, MinMax = (-0.01,-0.0001) )
        from ROOT import RooExponential as Exponential
        from RooFitWrappers import Pdf
        bkg = Pdf( kwargs.pop('Name','LP2011_Background_Mass'), Type = Exponential
                                           , Observables = (mass,)
                                           , Parameters = (self._m_bkg_exp, )
                                           )
        MassPdf.__init__(self, pdf = bkg )



