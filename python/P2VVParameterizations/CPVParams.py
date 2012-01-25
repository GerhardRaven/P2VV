###########################################################################################################################################
## P2VVParameterizations.CPVParams: CP violation parameterizations                                                                       ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin

phiCP      = 0.17
lambdaCPSq = 1.

class CPParam ( _util_parse_mixin ):
    def __init__( self, **kwargs ) :
        for coef in 'CDS' : setattr( self, '_' + coef, kwargs.pop(coef) )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )


class LambdaCarth_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar
        from math import cos, sin

        self._parseArg('ReLambdaCP', kwargs,  Title = 'CPV param. Re(lambda)', Value = cos(phiCP), MinMax = ( -2., 2. ) )
        self._parseArg('ImLambdaCP', kwargs,  Title = 'CPV param. Im(lambda)', Value = sin(phiCP), MinMax = ( -2., 2. ) )

        self._check_extraneous_kw( kwargs )
        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0*@0 - @1*@1) / (1. + @0*@0 + @1*@1)', [ self._ReLambdaCP, self._ImLambdaCP ] )
                             , D = FormulaVar('D', '2. * @0 / (1. + @0*@0 + @1*@1)',              [ self._ReLambdaCP, self._ImLambdaCP ] )
                             , S = FormulaVar('S', '2. * @1 / (1. + @0*@0 + @1*@1)',              [ self._ReLambdaCP, self._ImLambdaCP ] )
                        )


class LambdaSqArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar
        from math import pi

        self._parseArg( 'lambdaCPSq', kwargs,  Title = 'CPV param. lambda^2', Value = lambdaCPSq, MinMax = (  0.,      5.      ) )
        self._parseArg( 'phiCP',      kwargs,  Title = 'CPV param. phi',      Value = phiCP,      MinMax = ( -2. * pi, 2. * pi ) )
        self._check_extraneous_kw( kwargs )
        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0) / (1. + @0)',                [ self._lambdaCPSq              ] )
                             , D = FormulaVar('D', '2. * sqrt(@0) * cos(-@1) / (1. + @0)', [ self._lambdaCPSq, self._phiCP ] )
                             , S = FormulaVar('S', '2. * sqrt(@0) * sin(-@1) / (1. + @0)', [ self._lambdaCPSq, self._phiCP ] )
                        )

class LambdaArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar,ConstVar
        from math import pi

        self._parseArg( 'phiCP',      kwargs,  Title = 'CPV param. phi',      Value = phiCP,      MinMax = ( -2. * pi, 2. * pi ) )
        self._check_extraneous_kw( kwargs )
        CPParam.__init__(self, C = ConstVar( Name = 'C', Value = 0. )
                             , D = FormulaVar('D', ' cos(-@0) ', [  self._phiCP ] )
                             , S = FormulaVar('S', ' sin(-@0) ', [  self._phiCP ] )
                        )
