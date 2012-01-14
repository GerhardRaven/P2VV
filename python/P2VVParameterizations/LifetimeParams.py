###########################################################################################################################################
## P2VVParameterizations.LifetimeParams: Lifetime parameters                                                                             ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin


class LifetimeParams ( _util_parse_mixin ):
    def __init__( self, **kwargs ) :
        for coef in [ 'MeanLifetime', 'deltaGamma', 'deltaM' ] : setattr( self, '_' + coef, kwargs.pop(coef) )
        self._constraints = kwargs.pop( 'Constraints', None )
        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw )     : return getattr( self, '_' + kw )
    def ExternalConstraints( self ) : return self._constraints

class Gamma_LifetimeParams( LifetimeParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar, ConstVar, Pdf

        self._parseArg( 'Gamma',      kwargs, Title = 'Gamma',       Unit = 'ps^{-1}', Value = 0.68, MinMax = (  0.4,  0.9 ) )
        self._parseArg( 'deltaGamma', kwargs, Title = 'delta Gamma', Unit = 'ps^{-1}', Value = 0.05, MinMax = (- 0.3,  0.3 ) )
        self._parseArg( 'deltaM',     kwargs, Title = 'delta m',     Unit = 'ps^{-1}', Value = 17.8, MinMax = ( 16.5, 18.5 ) )
        
        self._check_extraneous_kw( kwargs )

        from ROOT import RooGaussian as Gaussian
        constraint = Pdf(  Name = self._deltaM.GetName() + '_constraint', Type = Gaussian
                         , Parameters = [  self._deltaM
                                         , ConstVar( Name = 'deltaM_mean',  Value = 17.63 )
                                         , ConstVar( Name = 'deltaM_sigma', Value = 0.11 )
                                        ]
                        )

        LifetimeParams.__init__( self
                                 , MeanLifetime = FormulaVar( 'MeanLifetime', '1. / @0', [self._Gamma], Title = 'B Mean lifetime' )
                                 , deltaGamma  = self._deltaGamma
                                 , deltaM      = self._deltaM
                                 , Constraints = [constraint]
                               )

class Tau_LifetimeParams( LifetimeParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar

        self._parseArg( 'MeanLifetime', kwargs, Title = 'MeanLifetime', Unit = 'ps',      Value = 1.48, MinMax = (  1.0,  2.5 ) )
        self._parseArg( 'deltaGamma',   kwargs, Title = 'delta Gamma',  Unit = 'ps^{-1}', Value = 0.05, MinMax = (- 0.3,  0.3)  )
        self._parseArg( 'deltaM',       kwargs, Title = 'delta m',      Unit = 'ps^{-1}', Value = 17.8, MinMax = ( 13.,  23.)   )

        self._check_extraneous_kw( kwargs )
        LifetimeParams.__init__( self, MeanLifetime = self._MeanLifetime, deltaGamma = self._deltaGamma, deltaM = self._deltaM )

