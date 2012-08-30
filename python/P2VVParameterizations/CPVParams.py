###########################################################################################################################################
## P2VVParameterizations.CPVParams: CP violation parameterizations                                                                       ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin

from math import sqrt, cos, sin
phiVal    = 0.
lambSqVal = 1.
lambVal   = sqrt(lambSqVal)
ReLambVal = lambVal * cos(-phiVal)
ImLambVal = lambVal * sin(-phiVal)
CVal      = ( 1. - lambSqVal ) / ( 1. + lambSqVal )
DVal      = 2. * ReLambVal / ( 1. + lambSqVal )
SVal      = 2. * ImLambVal / ( 1. + lambSqVal )

phiErr    = 0.1
lambSqErr = 0.08
lambErr   = lambSqErr / 2.
ReLambErr = lambSqErr / 2.
ImLambErr = phiErr
CErr      = lambSqErr / 2.
DErr      = lambSqErr / sqrt(2.)
SErr      = phiErr

from ROOT import RooNumber
RooInf = RooNumber.infinity()

class CPParam ( _util_parse_mixin ):
    def __init__( self, **kwargs ) :
        for coef in 'CDS' : setattr( self, '_' + coef, kwargs.pop(coef) )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )

    def CPVInDecay(self) : return False

    def C(self) : return self._C
    def D(self) : return self._D
    def S(self) : return self._S
    def R( self, CPVDecInd, amp0Ind, amp1Ind ) : return None

class CDS_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from math import cos, sin

        self._parseArg('C', kwargs,  Title = 'CPV param. C', Value = CVal, Error = CErr, MinMax = ( -1., 1. ) )
        self._parseArg('D', kwargs,  Title = 'CPV param. D', Value = DVal, Error = DErr, MinMax = ( -2., 2. ) )
        self._parseArg('S', kwargs,  Title = 'CPV param. S', Value = SVal, Error = SErr, MinMax = ( -2., 2. ) )
        self._check_extraneous_kw( kwargs )

class LambdaCarth_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar
        from math import cos, sin

        self._parseArg( 'ReLambdaCP', kwargs,  Title = 'CPV param. Re(lambda)', Value = ReLambVal, Error = ReLambErr
                       , MinMax = ( -RooInf, +RooInf ) )
        self._parseArg( 'ImLambdaCP', kwargs,  Title = 'CPV param. Im(lambda)', Value = ImLambVal, Error = ImLambErr
                       , MinMax = ( -RooInf, +RooInf ) )

        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0*@0 - @1*@1) / (1. + @0*@0 + @1*@1)', [ self._ReLambdaCP, self._ImLambdaCP ] )
                             , D = FormulaVar('D', '2. * @0 / (1. + @0*@0 + @1*@1)',              [ self._ReLambdaCP, self._ImLambdaCP ] )
                             , S = FormulaVar('S', '2. * @1 / (1. + @0*@0 + @1*@1)',              [ self._ReLambdaCP, self._ImLambdaCP ] )
                        )
        self._check_extraneous_kw( kwargs )


class LambdaAbsArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar
        from math import pi

        self._parseArg( 'lambdaCP', kwargs,  Title = 'CPV param. |lambda|', Value = lambVal, Error = lambErr, MinMax = ( 0.,      5.   ) )
        self._parseArg( 'phiCP',    kwargs,  Title = 'CPV param. phi',      Value = phiVal,  Error = phiErr,  MinMax = (-RooInf, +RooInf) )

        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0*@0) / (1. + @0*@0)',       [ self._lambdaCP              ] )
                             , D = FormulaVar('D', '2. * @0 * cos(-@1) / (1. + @0*@0)', [ self._lambdaCP, self._phiCP ] )
                             , S = FormulaVar('S', '2. * @0 * sin(-@1) / (1. + @0*@0)', [ self._lambdaCP, self._phiCP ] )
                        )
        self._check_extraneous_kw( kwargs )


class LambdaSqArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar
        from math import pi

        self._parseArg( 'lambdaCPSq', kwargs,  Title = 'CPV param. |lambda|^2', Value = lambSqVal, Error = lambSqErr
                       , MinMax = ( 0., 25. ) )
        self._parseArg( 'phiCP',      kwargs,  Title = 'CPV param. phi',        Value = phiVal,    Error = phiErr
                       , MinMax = ( -RooInf, +RooInf ) )

        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0) / (1. + @0)',                [ self._lambdaCPSq              ] )
                             , D = FormulaVar('D', '2. * sqrt(@0) * cos(-@1) / (1. + @0)', [ self._lambdaCPSq, self._phiCP ] )
                             , S = FormulaVar('S', '2. * sqrt(@0) * sin(-@1) / (1. + @0)', [ self._lambdaCPSq, self._phiCP ] )
                        )
        self._check_extraneous_kw( kwargs )

class LambdaArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar,ConstVar
        from math import pi

        self._parseArg( 'phiCP', kwargs, Title = 'CPV param. phi', Value = phiVal, Error = phiErr, MinMax = ( -RooInf, +RooInf ) )
        CPParam.__init__(self, C = ConstVar( Name = 'C', Value = 0. )
                             , D = FormulaVar('D', 'cos(-@0) ', [ self._phiCP ] )
                             , S = FormulaVar('S', 'sin(-@0) ', [ self._phiCP ] )
                        )
        self._check_extraneous_kw( kwargs )
