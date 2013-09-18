from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin

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
        self._coefs = [ ]
        for bin in range( 1, self._numBins ) :
            self._parseArg(  '%s_coef%d' % ( self._name, bin ), kwargs
                           , Title    = '%s bin coefficient %d' % ( self._name, bin )
                           , Value    = self._numEventsBins[bin] / self._numEvents if self._data else 1. / self._numBins
                           , MinMax   = ( 0., 1. )
                           , ContainerList = self._coefs
                          )

        # create a BinnedPdf
        from P2VV.RooFitWrappers import BinnedPdf
        pdf = BinnedPdf(  Name = self._name
                        , Observable = self._mass
                        , Binning = self._bins
                        , Coefficients = self._coefs
                        , BinIntegralCoefs = True
                       )

        # initialize
        MassPdf.__init__( self, pdf = pdf )
        self._check_extraneous_kw( kwargs )


class Gauss_Signal_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._parseArg( 'm_sig_mean',  kwargs, Title = 'B Mass', Unit = 'MeV/c^2'
                       , Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
        self._parseArg( 'm_sig_sigma', kwargs, Title = 'B Mass resolution', Unit = 'MeV/c^2'
                       , Value = 7.2,   Error = 0.04, MinMax = ( 0.1, 20. ) )

        from ROOT import RooGaussian as Gauss
        from P2VV.RooFitWrappers import Pdf
        MassPdf.__init__(  self
                         , pdf = Pdf(  Name = kwargs.pop( 'Name', 'Gauss_Signal_Mass' )
                                     , Type = Gauss
                                     , Parameters = ( mass, getattr( self, '_m_sig_mean' ), getattr( self, '_m_sig_sigma' ) )
                                    )
                        )
        self._check_extraneous_kw( kwargs )


class DoubleGauss_Signal_Mass ( MassPdf ) :
    def __init__( self, mass, **kwargs ) :
        namePF = self.getNamePrefix(kwargs)
        self._transWidthPars = kwargs.pop( 'TransformWidthPars', ( ) )
        self._parseArg( 'm_sig_mean', kwargs, Title = 'B Mass', Unit = 'MeV/c^2', Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )

        if self._transWidthPars :
            # transform width parameters
            self._parseArg( 'm_sig_widthPar0', kwargs, Title = 'B mass width parameter 0', Value = 20.,  Error = 1., MinMax = (0., 40.)   )
            self._parseArg( 'm_sig_widthPar1', kwargs, Title = 'B mass width parameter 1', Value = -59., Error = 1., MinMax = (-100., 0.) )
            self._parseArg( 'm_sig_widthPar2', kwargs, Title = 'B mass width parameter 2', Value = 47.,  Error = 1., MinMax = (0., 100.)  )

            for name in ( 'm_sig_frac', 'm_sig_sigma_1', 'm_sig_sigma_2' ) :
                self._parseArg( name, kwargs
                               , Formula = '%.10e*@0+%.10e*@1+%.10e*@2' % ( self._transWidthPars[name][0], self._transWidthPars[name][1]
                                                                           , self._transWidthPars[name][2] )
                               , Arguments = ( self._m_sig_widthPar0, self._m_sig_widthPar1, self._m_sig_widthPar2 )
                               , ObjectType = 'FormulaVar'
                              )

        else :
            # use fraction of first Gaussian and widths directly
            self._parseArg( 'm_sig_frac', kwargs, Title = 'B mass fraction first Gaussian', Value = 0.8, Error = 0.03
                           , MinMax = ( 0., 1. ) )
            self._parseArg( 'm_sig_sigma_1', kwargs, Title = 'B Mass resolution 1', Unit = 'MeV/c^2', Value = 6.3, Error = 0.1
                           , MinMax = ( 0.1, 20. ) )
            self._parseArg( 'm_sig_sigma_2', kwargs, Title = 'B Mass resolution 2', Unit = 'MeV/c^2', Value = 14.5, Error = 0.8
                           , MinMax = ( 0.2, 40. ) )

        from ROOT import RooGaussian as Gaussian
        from P2VV.RooFitWrappers import Pdf, SumPdf
        g1 = Pdf(  Name = '%s_m_sig_1' % namePF
                 , Type = Gaussian
                 , Parameters = ( mass, getattr( self, '_m_sig_mean' ), getattr( self, '_m_sig_sigma_1' ) )
                )
        g2 = Pdf(  Name = '%s_m_sig_2' % namePF
                 , Type = Gaussian
                 , Parameters = ( mass, getattr( self, '_m_sig_mean' ), getattr( self, '_m_sig_sigma_2' ) )
                )
        MassPdf.__init__( self, pdf = SumPdf( Name = kwargs.pop( 'Name', 'DoubleGauss_Signal_Mass' ), PDFs = ( g1, g2 )
                                             , Yields = { g1.GetName() : getattr( self, '_m_sig_frac' ) }
                                            )
                        )
        self._check_extraneous_kw( kwargs )


class LP2011_Signal_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        namePF = self.getNamePrefix(kwargs)
        self._parseArg( 'm_sig_mean', kwargs, Title = 'B Mass', Unit = 'MeV/c^2', Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
        self._parseArg( 'm_sig_sigma_1', kwargs, Title = 'B Mass resolution 1', Unit = 'MeV/c^2', Value = 6.3, Error = 0.1
                       , MinMax = ( 0.1, 20. ) )
        self._parseArg( 'm_sig_sigma_sf', kwargs, Title = 'B Mass resolution 2:1 scale factor', Value = 2.3, Error = 0.1
                       , MinMax = ( 0.1, 5. ) )
        self._parseArg( 'm_sig_frac', kwargs, Title = 'B mass fraction first Gaussian', Value = 0.8, Error = 0.03, MinMax = ( 0., 1. ) )

        from ROOT import RooGaussian as Gaussian
        from P2VV.RooFitWrappers import Pdf, SumPdf
        g1 = Pdf(  Name = '%sm_sig_1' % namePF
                 , Type = Gaussian
                 , Parameters = ( mass, self._m_sig_mean, self._m_sig_sigma_1 )
                )
        g2 = Pdf(  Name = '%sm_sig_2' % namePF
                 , Type = Gaussian
                 , Parameters = ( mass, self._m_sig_mean
                                 , self._parseArg( 'm_sig_sigma_2', kwargs, Formula = '@0*@1'
                                                  , Arguments = ( self._m_sig_sigma_sf, self._m_sig_sigma_1 )
                                                  , ObjectType = 'FormulaVar'
                                                 )
                                )
                )
        MassPdf.__init__( self, pdf = SumPdf( Name = kwargs.pop( 'Name', 'LP2011_Signal_Mass' ), PDFs = ( g1, g2 )
                                             , Yields = { g1.GetName() : self._m_sig_frac }
                                            )
                        )
        self._check_extraneous_kw( kwargs )


class CB_Signal_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._parseArg( 'm_sig_mean', kwargs, Title = 'B Mass', Unit = 'MeV/c^2', Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
        self._parseArg( 'm_sig_sigma', kwargs, Title = 'B Mass resolution', Unit = 'MeV/c^2', Value = 6.3, Error = 0.1
                       , MinMax = ( 0.1, 20. ) )
        self._parseArg( 'm_sig_alpha', kwargs, Title = 'B Mass tail parameter', Value = 9., Error = 15., MinMax = ( 0.1, 30. ) )
        self._parseArg( 'm_sig_n', kwargs, Title = 'B Mass tail order', Value = 1., ObjectType = 'ConstVar' )

        from ROOT import RooCBShape as CrystalBall
        from P2VV.RooFitWrappers import Pdf
        MassPdf.__init__(  self
                         , pdf = Pdf(  Name = kwargs.pop( 'Name', 'CB_Signal_Mass' )
                                     , Type = CrystalBall
                                     , Parameters = ( mass, getattr( self, '_m_sig_mean' ), getattr( self, '_m_sig_sigma' )
                                                     , getattr( self, '_m_sig_alpha' ), getattr( self, '_m_sig_n' )
                                                    )
                                    )
                        )
        self._check_extraneous_kw( kwargs )


class DoubleCB_Signal_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        namePF = self.getNamePrefix(kwargs)
        self._parseArg( 'm_sig_mean', kwargs, Title = 'B Mass', Unit = 'MeV/c^2', Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
        self._parseArg( 'm_sig_sigma_1', kwargs, Title = 'B Mass resolution 1', Unit = 'MeV/c^2', Value = 6.3, Error = 0.1
                       , MinMax = ( 0.1, 20. ) )
        self._parseArg( 'm_sig_sigma_sf', kwargs, Title = 'B Mass resolution 2:1 scale factor', Value = 2.3, Error = 0.1
                       , MinMax = ( 0.1, 5. ) )
        self._parseArg( 'm_sig_alpha_1', kwargs, Title = 'B Mass tail parameter 1', Value = 9., Error = 15., MinMax = ( 0.1, 30. ) )
        self._parseArg( 'm_sig_alpha_sf', kwargs, Title = 'B Mass tail parameter 2:1 scale factor', Value = 1., ObjectType = 'ConstVar' )
        self._parseArg( 'm_sig_n_1', kwargs, Title = 'B Mass tail order 1', Value = 1., ObjectType = 'ConstVar' )
        self._parseArg( 'm_sig_n_2', kwargs, Title = 'B Mass tail order 2', Value = 1., ObjectType = 'ConstVar' )
        self._parseArg( 'm_sig_frac', kwargs, Title = 'B mass fraction first CB', Value = 0.8, Error = 0.03, MinMax = ( 0., 1. ) )

        from ROOT import RooCBShape as CrystalBall
        from P2VV.RooFitWrappers import Pdf, SumPdf
        CB1 = Pdf(  Name ='%sm_sig_1' % namePF
                  , Type = CrystalBall
                  , Parameters = ( mass, self._m_sig_mean, self._m_sig_sigma_1, self._m_sig_alpha_1, self._m_sig_n_1 )
                 )
        CB2 = Pdf( Name = '%sm_sig_2' % namePF
                  , Type = CrystalBall
                  , Parameters = (  mass, self._m_sig_mean
                                  , self._parseArg( '_m_sig_sigma_2', kwargs, Formula = '@0*@1'
                                                   , Arguments = ( self._m_sig_sigma_sf, self._m_sig_sigma_1 ), ObjectType = 'FormulaVar' )
                                  , self._parseArg( '_m_sig_alpha_2', kwargs, Formula = '@0*@1'
                                                   , Arguments = ( self._m_sig_alpha_sf, self._m_sig_alpha_1 ), ObjectType = 'FormulaVar' )
                                  , self._m_sig_n_2
                                 )
                 )
        MassPdf.__init__( self, pdf = SumPdf( Name = kwargs.pop( 'Name', 'DoubleCB_Signal_Mass' ), PDFs = ( CB1, CB2 )
                         , Yields = { CB1.GetName() : self._m_sig_frac } ) )
        self._check_extraneous_kw( kwargs )


class Box_Signal_Mass ( MassPdf ) :
    def __init__( self, mass, **kwargs ) :
        self._parseArg( 'm_sig_mean', kwargs, Title = 'B Mass', Unit = 'MeV/c^2', Value = 5368., Error = 0.05, MinMax = ( 5000., 5700. ) )
        self._parseArg( 'm_sig_width', kwargs, Title = 'B Mass width', Unit = 'MeV/c^2', Value = 11., Error = 0.1,  MinMax = ( 0.1, 35. ) )

        from P2VV.Load import P2VVLibrary
        from ROOT import RooBoxPdf as BoxPdf
        from P2VV.RooFitWrappers import Pdf
        MassPdf.__init__( self, pdf = Pdf( Name = kwargs.pop( 'Name', 'Box_Signal_Mass' ), Type = BoxPdf
                                          , Parameters = ( mass, getattr( self, '_m_sig_mean' ), getattr( self, '_m_sig_width' ) )
                                         )
                        )
        self._check_extraneous_kw( kwargs )


class LP2011_Background_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        self._parseArg( 'm_bkg_exp', kwargs, Title = 'Mass background slope', Unit = 'c^2/MeV', Value = -0.002, Error = 0.0001
                       , MinMax = ( -0.05, 0. ) )

        from ROOT import RooExponential as Exponential
        from P2VV.RooFitWrappers import Pdf
        MassPdf.__init__( self, pdf = Pdf( Name = kwargs.pop('Name','LP2011_Background_Mass'), Type = Exponential 
                                          , Parameters = ( mass, getattr( self, '_m_bkg_exp' ) )
                                         )
                        )
        self._check_extraneous_kw( kwargs )


class Linear_Background_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        constant = kwargs.pop( 'Constant', False )
        self._parseArg( 'm_bkg_arg', kwargs, Title = 'Mass background slope', Unit = 'c^2/MeV', Value = -1.7e-4, Error = 2.e-6
                       , MinMax = ( -1.8e-4, 0. ) )
        if constant :
            getattr( self, '_m_bkg_arg' ).setVal(0.)
            getattr( self, '_m_bkg_arg' ).setConstant(True)

        from P2VV.RooFitWrappers import GenericPdf
        MassPdf.__init__(self, pdf = GenericPdf(  Name      = kwargs.pop('Name','Linear_Background_Mass')
                                                , Arguments = [ mass, getattr( self, '_m_bkg_arg' ) ]
                                                , Formula   = '1+@1*@0'
                                               )
                        )
        self._check_extraneous_kw( kwargs )


class Signal_PsiMass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        from ROOT import RooCBShape as CrystalBall
        from P2VV.RooFitWrappers import Pdf
        namePF = self.getNamePrefix(kwargs)
        self._parseArg( 'mpsi_mean', kwargs, Title = 'J/psi mass', Unit = 'MeV', Value = 3097., MinMax = ( 3090., 3105. ) )
        self._parseArg( 'mpsi_sigma', kwargs, Title = 'J/psi mass resolution',  Unit = 'MeV', Value = 14., MinMax = ( 8., 20. ) )
        self._parseArg( 'mpsi_alpha', kwargs, Title = 'J/psi mass CB alpha', Unit = '', Value = 1.90, MinMax = ( 1., 3. ) )
        self._parseArg( 'mpsi_n', kwargs, Title = 'J/psi mass CB n', Unit = '', Value = 2., MinMax = ( 0.1, 5. ), Constant = True )
        MassPdf.__init__( self, pdf = Pdf(  Name = kwargs.pop( 'Name','%sSignal_PsiMass' % namePF )
                                          , Type = CrystalBall
                                          , Parameters = [ mass, getattr( self, '_mpsi_mean' ), getattr( self, '_mpsi_sigma' )
                                                          , getattr( self, '_mpsi_alpha' ), getattr( self, '_mpsi_n' )
                                                         ]
                                         )
                        )
        self._check_extraneous_kw( kwargs )

class DoubleCB_Psi_Mass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        namePF = self.getNamePrefix(kwargs)
        self._parseArg( 'mpsi_mean', kwargs, Title = 'psi Mass core', Unit = 'MeV/c^2', Value = 3100., Error = 0.05
                       , MinMax = ( 3050., 3150. ) )
        self._parseArg( 'mpsi_sigma_1', kwargs, Title = 'psi Mass resolution 1', Unit = 'MeV/c^2', Value = 6.3, Error = 0.1
                       , MinMax = ( 0.1, 20. ) )
        self._parseArg( 'mpsi_sigma_sf', kwargs, Title = 'psi Mass resolution 2:1 scale factor', Value = 2.3, Error = 0.1
                       , MinMax = ( 0.1, 20. ) )
        self._parseArg( 'mpsi_alpha_1', kwargs, Title = 'psi Mass tail parameter 1', Value = 2., Error = 15., MinMax = ( 0.1, 30. ) )
        self._parseArg( 'mpsi_alpha_sf', kwargs, Title = 'psi Mass tail parameter 2:1 scale factor', Value = 1., ObjectType = 'ConstVar' )
        self._parseArg( 'mpsi_n_1', kwargs, Title = 'psi Mass tail order 1', Value = 2., ObjectType = 'ConstVar' )
        self._parseArg( 'mpsi_n_2', kwargs, Title = 'psi Mass tail order 2', Value = 2.,    ObjectType = 'ConstVar' )
        self._parseArg( 'mpsi_frac', kwargs, Title = 'psi mass fraction first CB', Value = 0.8, Error = 0.03, MinMax = ( 0., 1. ) )

        from ROOT import RooCBShape as CrystalBall
        from P2VV.RooFitWrappers import Pdf, SumPdf
        CB1 = Pdf( Name = '%smpsi_1' % namePF, Type = CrystalBall
                  , Parameters = ( mass, self._mpsi_mean, self._mpsi_sigma_1, self._mpsi_alpha_1, self._mpsi_n_1 )
                 )
        CB2 = Pdf( Name = '%smpsi_2' % namePF, Type = CrystalBall
                  , Parameters = (  mass, self._mpsi_mean
                                  , self._parseArg( 'mpsi_sigma_2', kwargs, Formula = '@0*@1'
                                                   , Arguments = ( self._mpsi_sigma_sf, self._mpsi_sigma_1 ), ObjectType = 'FormulaVar' )
                                  , self._parseArg( 'mpsi_alpha_2', kwargs, Formula = '@0*@1'
                                                   , Arguments = ( self._mpsi_alpha_sf, self._mpsi_alpha_1 ), ObjectType = 'FormulaVar' )
                                  , self._mpsi_n_2
                                 )
                 )
        MassPdf.__init__( self, pdf = SumPdf( Name = kwargs.pop( 'Name', 'DoubleCB_Psi_Mass' ), PDFs = (CB1, CB2)
                                             , Yields = { CB1.GetName() : self._mpsi_frac } )
                        )
        self._check_extraneous_kw( kwargs )

class Background_PsiMass ( MassPdf ) :
    def __init__(self, mass, **kwargs ) :
        namePF = self.getNamePrefix(kwargs)
        from ROOT import RooExponential as Exponential
        from P2VV.RooFitWrappers import Pdf
        self._parseArg( 'mpsi_c', kwargs, Title = 'J/psi mass background slope', Unit = '1/MeV', Value = -0.01
                       , MinMax = ( -0.1, -0.000001 ) )
        MassPdf.__init__( self, pdf = Pdf( Name = kwargs.pop( 'Name', '%sBackground_PsiMass' % namePF ), Type = Exponential
                                          , Parameters = [mass, self._mpsi_c]))
        self._check_extraneous_kw( kwargs )
