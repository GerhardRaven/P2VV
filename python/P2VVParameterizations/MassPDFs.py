from P2VVParameterizations.GeneralUtils import _util_parse_mixin

from ROOT import RooNumber
RooInf = RooNumber.infinity()

class MassPdf( _util_parse_mixin ) :
    def __init__(self, **kwargs ) :
        for ( k, v ) in kwargs.iteritems() : setattr( self, '_' + k, v )
    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def pdf(self) : return self._pdf

class Binned_MassPdf( MassPdf ) :
    def __init__( self, Name, Mass, **kwargs ) :
        self._name = Name
        self._mass = Mass

        # get binning
        self._bins = kwargs.pop( 'Binning', None )
        if not self._bins :
            # create binning
            from array import array
            binBounds = kwargs.pop( 'BinBoundaries', [ self._mass.getMin(), self._mass.getMax() ] )
            self._binBounds = array( 'd', binBounds )
            self._numBins = len(binBounds) - 1

            from ROOT import RooBinning
            self._bins = RooBinning( self._numBins, self._binBounds, self._name + '_binning' )
            self._mass.setBinning( self._bins, self._name + '_binning' )

        self._numBins = self._bins.numBins()

        # determine number of events in each bin
        self._data = kwargs.pop( 'Data', None )
        if self._data :
            assert self._mass._var in self._data.get(0),\
                    'Binned_MassPdf.__init__(): %s is not and observable in the provided data set' % self._mass.GetName()
            self._numEvents = self._data.sumEntries()
            self._numEventsBins = self._numBins * [ 0. ]
            for obsSet in self._data :
                bin = self._bins.binNumber( obsSet.getRealValue( self._mass.GetName() ) )
                self._numEventsBins[bin] += self._data.weight()
                

        # create bin coefficients
        from RooFitWrappers import RealVar
        self._coefs = [ RealVar(  '%s_coef%d' % ( self._name, bin )
                                , Title    = '%s bin coefficient %d' % ( self._name, bin )
                                , Value    = self._numEventsBins[bin] / self._numEvents if self._data else 1. / self._numBins
                                , MinMax   = ( 0., 1. )
                                , Constant = False
                               ) if bin != 0 else None for bin in range( self._numBins )
                      ]
        del self._coefs[0]

        # create a BinnedPdf
        from RooFitWrappers import BinnedPdf
        pdf = BinnedPdf(  Name = self._name
                        , Observable = self._mass
                        , Binning = self._bins
                        , Coefficients = self._coefs
                        , BinIntegralCoefs = True
                       )

        # initialize
        MassPdf.__init__( self, pdf = pdf )


class LP2011_Signal_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._prefix = kwargs.pop("Prefix", "")
        self._parseArg( '%sm_sig_mean' % self._prefix,     kwargs, Title = 'B Mass', Unit = 'MeV/c^2'
                       , Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
        self._parseArg( '%sm_sig_sigma_1' % self._prefix,  kwargs, Title = 'B Mass resolution 1', Unit = 'MeV/c^2'
                       , Value = 6.3,   Error = 0.1,  MinMax = ( 0.1, 20. ) )
        self._parseArg( '%sm_sig_sigma_sf' % self._prefix, kwargs, Title = 'B Mass resolution 2:1 scale factor'
                       , Value = 2.3,   Error = 0.1,  MinMax = ( 0.1, 5. ) )
        self._parseArg( '%sm_sig_frac' % self._prefix,     kwargs, Title = 'B mass fraction first Gaussian'
                       , Value = 0.8,   Error = 0.03, MinMax = ( 0., 1. ) )

        from ROOT import RooGaussian as Gaussian
        from RooFitWrappers import Pdf, FormulaVar, SumPdf
        g1 = Pdf(  Name ='%sm_sig_1' % self._prefix
                 , Type = Gaussian
                 , Parameters = (mass, getattr(self, '_%sm_sig_mean' % self._prefix),
                                 getattr(self, '_%sm_sig_sigma_1' % self._prefix))
                )
        g2 = Pdf( Name = '%sm_sig_2' % self._prefix
                 , Type = Gaussian
                 , Parameters = ( mass, getattr(self, '_%sm_sig_mean' % self._prefix)
                                 , FormulaVar('_%sm_sig_sigma_2' % self._prefix, '@0*@1', (getattr(self, '_%sm_sig_sigma_sf' % self._prefix),
                                                                                           getattr(self, '_%sm_sig_sigma_1' % self._prefix))))
                )
        MassPdf.__init__( self, pdf = SumPdf( Name = kwargs.pop( 'Name', 'LP2011_Signal_Mass' ), PDFs = ( g1, g2 )
                         , Yields = { g1.GetName() : getattr(self, '_%sm_sig_frac' % self._prefix) } ) )


class Box_Signal_Mass ( MassPdf ) :
    def __init__( self, mass, **kwargs ) :
        self._prefix = kwargs.pop( 'Prefix', '' )
        self._parseArg( '%sm_sig_mean' % self._prefix,     kwargs, Title = 'B Mass',     Unit = 'MeV/c^2'
                       , Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
        self._parseArg( '%sm_sig_width' % self._prefix,  kwargs, Title = 'B Mass width', Unit = 'MeV/c^2'
                       , Value = 11.,   Error = 0.1,  MinMax = ( 0.1, 35. ) )

        from ROOT import RooBoxPdf as BoxPdf
        from RooFitWrappers import Pdf
        MassPdf.__init__( self, pdf = Pdf(  Name = kwargs.pop( 'Name', 'Box_Signal_Mass' )
                                          , Type = BoxPdf
                                          , Parameters = (  mass, getattr( self, '_%sm_sig_mean' % self._prefix )
                                                          , getattr( self, '_%sm_sig_width' % self._prefix ) )
                                         )
                        )


class LP2011_Background_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._prefix = kwargs.pop("Prefix", "")
        self._parseArg( '%sm_bkg_exp' % self._prefix, kwargs, Title = 'Mass background slope',
                        Unit = 'c^2/MeV', Value = -0.002, Error = 0.0001, MinMax = ( -0.05, 0. ) )

        from ROOT import RooExponential as Exponential
        from RooFitWrappers import Pdf
        MassPdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name','LP2011_Background_Mass')
                                        , Type = Exponential 
                                        , Parameters = (mass, getattr(self, '_%sm_bkg_exp' % self._prefix))))


class Linear_Background_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._prefix = kwargs.pop( 'Prefix', '' )
        constant = kwargs.pop( 'Constant', False )
        self._parseArg( '%sm_bkg_arg' % self._prefix, kwargs, Title = 'Mass background slope',
                        Unit = 'c^2/MeV', Value = -1.7e-4, Error = 2.e-6, MinMax = ( -1.8e-4, 0. ) )
        if constant :
            getattr( self, '_%sm_bkg_arg' % self._prefix ).setVal(0.)
            getattr( self, '_%sm_bkg_arg' % self._prefix ).setConstant(True)

        from RooFitWrappers import GenericPdf
        MassPdf.__init__(self, pdf = GenericPdf(  Name      = kwargs.pop('Name','Linear_Background_Mass')
                                                , Arguments = [ mass, getattr( self, '_%sm_bkg_arg' % self._prefix ) ]
                                                , Formula   = '1+@1*@0'
                                               )
                        )


class Signal_PsiMass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        from ROOT import RooCBShape as CrystalBall
        from RooFitWrappers import Pdf
        self._prefix = kwargs.pop("Prefix", "")

        self._parseArg( '%smpsi_mean' % self._prefix,  kwargs, Title = 'J/psi mass',  Unit = 'MeV', Value = 3097, MinMax = (3090, 3105))
        self._parseArg( '%smpsi_sigma' % self._prefix, kwargs, Title = 'J/psi mass resolution',  Unit = 'MeV', Value = 14, MinMax = (8, 20))
        self._parseArg( '%smpsi_alpha' % self._prefix, kwargs, Title = 'J/psi mass CB alpha', Unit = '', Value = 1.90, MinMax = (1, 3))
        self._parseArg( '%smpsi_n' % self._prefix,     kwargs, Title = 'J/psi mass CB n',  Unit = '', Value = 2, MinMax = (0.1, 5), Constant = True)
        MassPdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name','%sSignal_PsiMass' % self._prefix)
                                        , Type = CrystalBall
                                        , Parameters = [mass, getattr(self, '_%smpsi_mean' % self._prefix),
                                                        getattr(self, '_%smpsi_sigma' % self._prefix),
                                                        getattr(self, '_%smpsi_alpha' % self._prefix),
                                                        getattr(self, '_%smpsi_n' % self._prefix)]))

class Background_PsiMass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._prefix = kwargs.pop("Prefix", "")
        from ROOT import RooExponential as Exponential
        from RooFitWrappers import Pdf
        self._parseArg( '%smpsi_c' % self._prefix, kwargs, Title = 'J/psi mass background slope', Unit = '1/MeV', Value = -0.01, MinMax = (-0.1, -0.0001))
        MassPdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name', '%sBackground_PsiMass' % self._prefix)
                                          , Type = Exponential
                                          , Parameters = [mass, getattr(self, '_%smpsi_c' % self._prefix)]))
