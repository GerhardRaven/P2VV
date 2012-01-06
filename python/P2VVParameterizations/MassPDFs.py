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
        self._parseArg('m_sig_sigma_1', kwargs, Title = 'B Mass resolution 1', Unit = 'MeV/c^2', Value = 6.447, MinMax = (3,20) )
        self._parseArg('m_sig_sigma_2_scale', kwargs, Title = 'B Mass resolution 2 scale factor', Unit = 'MeV/c^2', Value = 2.14, MinMax = (1,5), Constant = True )
        self._parseArg('m_sig_f',kwargs, Title = 'B mass fraction 2nd Gaussian', Value = 0.83, MinMax = (0.,1.), Constant = True )

        from ROOT import RooGaussian as Gaussian
        from RooFitWrappers import Pdf, FormulaVar, SumPdf
        g1 = Pdf( Name ='m_sig_1'
                , Type = Gaussian
                , Parameters = (mass, self._m_sig_mean, self._m_sig_sigma_1 )
                )
        g2 = Pdf( Name = 'm_sig_2'
                , Type = Gaussian
                , Parameters = ( mass, self._m_sig_mean
                               , FormulaVar('m_sig_sigma_2','@0*@1',(self._m_sig_sigma_2_scale, self._m_sig_sigma_1))
                               )
                )
        MassPdf.__init__(self, pdf = SumPdf(Name = kwargs.pop('Name','LP2011_Signal_Mass'),   PDFs = (  g1,  g2)  , Yields = { g1.GetName() : self._m_sig_f } ) )


class LP2011_Background_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._parseArg('m_bkg_exp', kwargs, Title = 'Mass background slope', Unit = 'MeV/c^2', Value = -0.001, MinMax = (-0.01,-0.0001) )
        from ROOT import RooExponential as Exponential
        from RooFitWrappers import Pdf
        MassPdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name','LP2011_Background_Mass')
                                        , Type = Exponential 
                                        , Parameters = (mass, self._m_bkg_exp, )) )


class Signal_PsiMass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        from ROOT import RooCBShape as CrystalBall
        from RooFitWrappers import Pdf
        self._parseArg( 'mpsi_mean',  kwargs, Title = 'J/psi mass',  Unit = 'MeV', Value = 3097, MinMax = (3090, 3105))
        self._parseArg( 'mpsi_sigma', kwargs, Title = 'J/psi mass resolution',  Unit = 'MeV', Value = 14, MinMax = (8, 20))
        self._parseArg( 'mpsi_alpha', kwargs, Title = 'J/psi mass CB alpha', Unit = '', Value = 1.90, MinMax = (1, 3))
        self._parseArg( 'mpsi_n',     kwargs, Title = 'J/psi mass CB n',  Unit = '', Value = 2, MinMax = (0.1, 5), Constant = True)
        MassPdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name','Signal_PsiMass')
                                        , Type = CrystalBall
                                        , Parameters = [ mass, self._mpsi_mean, self._mpsi_sigma, self._mpsi_alpha, self._mpsi_n]))

class Background_PsiMass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        from ROOT import RooExponential as Exponential
        from RooFitWrappers import Pdf
        self._parseArg( 'mpsi_c', kwargs, Title = 'J/psi mass background slope', Unit = '1/MeV', Value = -0.01, MinMax = (-0.1, -0.0001))
        MassPdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name','Background_PsiMass')
                                        , Type = Exponential
                                        , Parameters = [ mass, self._mpsi_c ]) )
