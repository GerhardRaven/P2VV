###########################################################################################################################################
## P2VVParameterizations.LifetimeParams: Lifetime parameters                                                                             ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin

GammaVal       =  0.67
GammaErr       =  0.005
DGammaVal      =  0.11
DGammaErr      =  0.02
DMVal          =  17.77 #17.6
DMErr          =  0.1
BetaTimeEffVal = -8.3e-3
BetaTimeEffErr =  4.0e-3

from P2VV.Imports import extConstraintValues
( DMConstrVal,   DMConstrErr   ) = extConstraintValues.getSetVal( 'DM',   (  17.768, 0.024  ) )
( betaConstrVal, betaConstrErr ) = extConstraintValues.getSetVal( 'beta', ( -8.3e-3, 4.0e-3 ) )

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
        # parse arguments
        self._parseArg( 'Gamma', kwargs, Title = 'Gamma', Unit = 'ps^{-1}', Value = GammaVal, Error = GammaErr, MinMax = ( 0.1, 10. ) )
        self._parseArg( 'dGamma', kwargs, Title = 'delta Gamma', Unit = 'ps^{-1}', Value = DGammaVal, Error = DGammaErr
                       , MinMax = ( -RooInf, +RooInf ) )
        self._parseArg( 'dM', kwargs, Title = 'delta m', Unit = 'ps^{-1}', Value = DMVal, Error = DMErr, MinMax = ( -RooInf, +RooInf ) )
        
        dMConstr   = kwargs.pop( 'dMConstraint',   None )
        betaConstr = kwargs.pop( 'betaConstraint', 'noBeta' )
        if betaConstr == 'noBeta' :
            self._timeEffGamma = self._Gamma
        else :
            # include "upper decay-time efficiency" effect in PDF: exp(-Gamma*t) --> exp(+beta*t) * exp(-Gamma*t) = exp(-(Gamma-beta)*t)
            self._parseArg( 'betaTimeEff', kwargs, Title = 'Upper time efficiency factor', Unit = 'ps^{-1}', Value = BetaTimeEffVal
                           , Error = BetaTimeEffErr, MinMax = ( -RooInf, +RooInf ) )
            self._parseArg( 'timeEffGamma', kwargs, Title = 'Efficiency corrected Gamma', Formula = '@0-@1'
                           , Arguments = [ self._Gamma, self._betaTimeEff ], ObjectType = 'FormulaVar' )

        self._check_extraneous_kw( kwargs )

        # build lifetime parmeters object
        LifetimeParams.__init__(  self
                                , MeanLifetime = self._parseArg( 'MeanLifetime', kwargs, Formula = '1. / @0', Title = 'B Mean lifetime'
                                                                , Arguments = [ self._timeEffGamma ], ObjectType = 'FormulaVar' )
                                , dGamma  = self._dGamma
                                , dM      = self._dM
                               )

        # build constraint PDF for Deltam
        if dMConstr == 'fixed' :
            self._dM.setConstant(True)
            self._dM.setVal(DMConstrVal)
            self._dM.setError(DMConstrErr)

        elif dMConstr :
            from P2VV.RooFitWrappers import Pdf
            from ROOT import RooGaussian as Gaussian
            self.addConstraint( Pdf(  Name = self._dM.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._dM
                                                    , self._parseArg( 'dM_mean',  kwargs, Value = DMConstrVal, ObjectType = 'ConstVar' )
                                                    , self._parseArg( 'dM_sigma', kwargs, Value = DMConstrErr, ObjectType = 'ConstVar' )
                                                   ]
                                   )
                              )

        # build constraint for decay-time efficiency factor beta
        if betaConstr == 'fixed' :
            self._betaTimeEff.setConstant(True)
            self._betaTimeEff.setVal(betaConstrVal)
            self._betaTimeEff.setError(betaConstrErr)

        elif betaConstr and betaConstr != 'noBeta' :
            from P2VV.RooFitWrappers import Pdf
            from ROOT import RooGaussian as Gaussian
            self.addConstraint( Pdf(  Name = self._betaTimeEff.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._betaTimeEff
                                                    , self._parseArg('beta_mean',  kwargs, Value = betaConstrVal, ObjectType = 'ConstVar')
                                                    , self._parseArg('beta_sigma', kwargs, Value = betaConstrErr, ObjectType = 'ConstVar')
                                                   ]
                                   )
                              )


class Tau_LifetimeParams( LifetimeParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'MeanLifetime', kwargs, Title = 'MeanLifetime', Unit = 'ps', Value = 1. / GammaVal, Error = GammaErr / GammaVal**2
                       , MinMax = ( 0.1, 10. ) )
        self._parseArg( 'dGamma', kwargs, Title = 'delta Gamma', Unit = 'ps^{-1}', Value = DGammaVal, Error = DGammaErr
                       , MinMax = ( -RooInf, RooInf ) )
        self._parseArg( 'dM', kwargs, Title = 'delta m', Unit = 'ps^{-1}', Value = DMVal, Error = DMErr, MinMax = ( -RooInf, RooInf ) )

        self._check_extraneous_kw( kwargs )
        LifetimeParams.__init__( self, MeanLifetime = self._MeanLifetime, dGamma = self._dGamma, dM = self._dM )

