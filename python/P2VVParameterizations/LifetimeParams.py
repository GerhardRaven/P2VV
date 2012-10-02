###########################################################################################################################################
## P2VVParameterizations.LifetimeParams: Lifetime parameters                                                                             ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin

GammaVal  = 0.67
GammaErr  = 0.005
DGammaVal = 0.11
DGammaErr = 0.02
DMVal     = 17.6
DMErr     = 0.1

DMConstrVal = 17.63
DMConstrErr = 0.11

from ROOT import RooNumber
RooInf = RooNumber.infinity()

class LifetimeParams ( _util_parse_mixin, _util_extConstraints_mixin ):
    def __init__( self, **kwargs ) :
        for coef in [ 'MeanLifetime', 'dGamma', 'dM' ] : setattr( self, '_' + coef, kwargs.pop(coef) )
        _util_extConstraints_mixin.__init__( self, kwargs )
        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw )     : return getattr( self, '_' + kw )

class Gamma_LifetimeParams( LifetimeParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar, ConstVar, Pdf

        self._parseArg( 'Gamma', kwargs, Title = 'Gamma', Unit = 'ps^{-1}', Value = GammaVal, Error = GammaErr, MinMax = ( 0.1, 10. ) )
        self._parseArg( 'dGamma', kwargs, Title = 'delta Gamma', Unit = 'ps^{-1}', Value = DGammaVal, Error = DGammaErr
                       , MinMax = ( -RooInf, +RooInf ) )
        self._parseArg( 'dM', kwargs, Title = 'delta m', Unit = 'ps^{-1}', Value = DMVal, Error = DMErr, MinMax = ( -RooInf, +RooInf ) )
        
        constraints = [  ]
        dMConstr = kwargs.pop( 'dMConstraint', None )
        if type(dMConstr) == str and dMConstr == 'fixed' :
            self._dM.setConstant(True)
            self._dM.setVal(DMConstrVal)
            self._dM.setError(DMConstrErr)

        elif dMConstr :
            from ROOT import RooGaussian as Gaussian
            constraints.append( Pdf(  Name = self._dM.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._dM
                                                    , ConstVar( Name = 'dM_mean',  Value = DMConstrVal )
                                                    , ConstVar( Name = 'dM_sigma', Value = DMConstrErr )
                                                   ]
                                   )
                              )

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

        self._parseArg( 'MeanLifetime', kwargs, Title = 'MeanLifetime', Unit = 'ps', Value = 1. / GammaVal, Error = GammaErr / GammaVal**2
                       , MinMax = ( 0.1, 10. ) )
        self._parseArg( 'dGamma', kwargs, Title = 'delta Gamma', Unit = 'ps^{-1}', Value = DGammaVal, Error = DGammaErr
                       , MinMax = ( -RooInf, RooInf ) )
        self._parseArg( 'dM', kwargs, Title = 'delta m', Unit = 'ps^{-1}', Value = DMVal, Error = DMErr, MinMax = ( -RooInf, RooInf ) )

        self._check_extraneous_kw( kwargs )
        LifetimeParams.__init__( self, MeanLifetime = self._MeanLifetime, dGamma = self._dGamma, dM = self._dM )

