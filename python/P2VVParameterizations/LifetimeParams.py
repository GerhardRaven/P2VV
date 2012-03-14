###########################################################################################################################################
## P2VVParameterizations.LifetimeParams: Lifetime parameters                                                                             ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin


class LifetimeParams ( _util_parse_mixin, _util_extConstraints_mixin ):
    def __init__( self, **kwargs ) :
        for coef in [ 'MeanLifetime', 'dGamma', 'dM' ] : setattr( self, '_' + coef, kwargs.pop(coef) )
        _util_extConstraints_mixin.__init__( self, kwargs )
        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw )     : return getattr( self, '_' + kw )

class Gamma_LifetimeParams( LifetimeParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar, ConstVar, Pdf

        self._parseArg( 'Gamma',  kwargs, Title = 'Gamma',       Unit = 'ps^{-1}', Value = 0.667,  MinMax = (  0.4,  0.9 ) )
        self._parseArg( 'dGamma', kwargs, Title = 'delta Gamma', Unit = 'ps^{-1}', Value = 0.12,   MinMax = (- 1.,   1.  ) )
        self._parseArg( 'dM',     kwargs, Title = 'delta m',     Unit = 'ps^{-1}', Value = 17.63,  MinMax = ( 16.5, 18.5 ) )
        
        constraints = [  ]
        if kwargs.pop( 'dMConstraint', None ) :
            from ROOT import RooGaussian as Gaussian
            constraints.append( Pdf(  Name = self._dM.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._dM
                                                    , ConstVar( Name = 'dM_mean',  Value = 17.63 )
                                                    , ConstVar( Name = 'dM_sigma', Value = 0.11 )
                                                   ]
                                   )
                              )
            self._dM['Error'] = 0.11

        self._check_extraneous_kw( kwargs )
        LifetimeParams.__init__( self
                                 , MeanLifetime = FormulaVar( 'MeanLifetime', '1. / @0', [self._Gamma], Title = 'B Mean lifetime' )
                                 , dGamma  = self._dGamma
                                 , dM      = self._dM
                                 , Constraints = constraints
                               )

class Tau_LifetimeParams( LifetimeParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar

        self._parseArg( 'MeanLifetime', kwargs, Title = 'MeanLifetime', Unit = 'ps',      Value = 1.48, MinMax = (  1.0,  2.5 ) )
        self._parseArg( 'dGamma',       kwargs, Title = 'delta Gamma',  Unit = 'ps^{-1}', Value = 0.05, MinMax = (- 0.3,  0.3)  )
        self._parseArg( 'dM',           kwargs, Title = 'delta m',      Unit = 'ps^{-1}', Value = 17.8, MinMax = ( 13.,  23.)   )

        self._check_extraneous_kw( kwargs )
        LifetimeParams.__init__( self, MeanLifetime = self._MeanLifetime, dGamma = self._dGamma, dM = self._dM )

