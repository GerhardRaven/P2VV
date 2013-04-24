###########################################################################################################################################
## P2VVParameterizations.CPVParams: CP violation parameterizations                                                                       ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin

###########################################################################################################################################
## follow HFAG conventions:                                                                                                              ##
## lambda_f = q/p * A_f^bar / A_f                                                                                                        ##
##          = eta_f * |lambda_f| * e^(-i * phi_s^f)                                                                                      ##
##          = eta_f * |lambda_f| * [cos(-phi_s^f) + i * sin(-phi_s^f)]                                                                   ##
##          = eta_f * |lambda_f| * [cos(phi_s^f)  - i * sin(phi_s^f) ]                                                                   ##
## C_f = (1 - |lambda_f|^2) / (1 + |lambda_f|^2)                                                                                         ##
## D_f = - Re(lambda_f) / (1 + |lambda_f|^2)                                                                                             ##
##     = -eta_f * |lambda_f| / (1 + |lambda_f|^2) * cos(-phi_s^f)                                                                        ##
##     = -eta_f * |lambda_f| / (1 + |lambda_f|^2) * cos(phi_s^f)                                                                         ##
## S_f = + Im(lambda_f) / (1 + |lambda_f|^2)                                                                                             ##
##     = +eta_f * |lambda_f| / (1 + |lambda_f|^2) * sin(-phi_s^f)                                                                        ##
##     = -eta_f * |lambda_f| / (1 + |lambda_f|^2) * sin(phi_s^f)                                                                         ##
###########################################################################################################################################
from math import sqrt, cos, sin
phiVal    = 0.
lambSqVal = 1.
lambVal   = sqrt(lambSqVal)
ReLambVal = lambVal * cos(-phiVal)
ImLambVal = lambVal * sin(-phiVal)
CVal      = ( 1. - lambSqVal ) / ( 1. + lambSqVal )
DVal      = -2. * ReLambVal / ( 1. + lambSqVal )
SVal      = +2. * ImLambVal / ( 1. + lambSqVal )

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
        self._ampNames = kwargs.pop( 'AmplitudeNames', None )

        if self._ampNames :
            for coefDict in [ 'RPlusDict', 'RMinDict', 'RReDict', 'RImDict' ] : setattr( self, '_' + coefDict, kwargs.pop(coefDict) )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )

    def CPVInDecay(self) : return True if self._ampNames else False
    def AmplitudeNames(self) : return self._ampNames

    def C(self) : return self._C
    def D(self) : return self._D
    def S(self) : return self._S
    def R( self, CPVDecInd, amp0Ind, amp1Ind ) :
        if not self.CPVInDecay() : return None

        assert CPVDecInd in ( 'plus', 'min', 'Re', 'Im' )\
            , 'P2VV - ERROR: CPParam.R(): valid CP violation in decay parameter indices are "plus", "min", "Re" or "Im" (got "%s")'\
              % CPVDecInd
        assert amp0Ind in self._ampNames and amp1Ind in self._ampNames\
            , 'P2VV - ERROR: CPParam.R(): valid amplitude indices are "%s" (got "%s" and "%s")'\
              % ( '", "'.join( amp for amp in self._ampNames ), amp0Ind, amp1Ind )

        return self._RPlusDict[ ( amp0Ind, amp1Ind ) ] if CPVDecInd == 'plus'\
          else self._RMinDict[  ( amp0Ind, amp1Ind ) ] if CPVDecInd == 'min'\
          else self._RReDict[   ( amp0Ind, amp1Ind ) ] if CPVDecInd == 'Re'\
          else self._RImDict[   ( amp0Ind, amp1Ind ) ]

class CDS_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from math import cos, sin

        self._parseArg('C', kwargs,  Title = 'CPV param. C', Value = CVal, Error = CErr, MinMax = ( -1., 1. ) )
        self._parseArg('D', kwargs,  Title = 'CPV param. D', Value = DVal, Error = DErr, MinMax = ( -2., 2. ) )
        self._parseArg('S', kwargs,  Title = 'CPV param. S', Value = SVal, Error = SErr, MinMax = ( -2., 2. ) )
        self._check_extraneous_kw( kwargs )

class LambdaCarth_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import FormulaVar
        from math import cos, sin

        self._parseArg( 'ReLambdaCP', kwargs,  Title = 'CPV param. Re(lambda)', Value = ReLambVal, Error = ReLambErr
                       , MinMax = ( -RooInf, +RooInf ) )
        self._parseArg( 'ImLambdaCP', kwargs,  Title = 'CPV param. Im(lambda)', Value = ImLambVal, Error = ImLambErr
                       , MinMax = ( -RooInf, +RooInf ) )

        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0*@0 - @1*@1) / (1. + @0*@0 + @1*@1)', [ self._ReLambdaCP, self._ImLambdaCP ] )
                             , D = FormulaVar('D', '-2. * @0 / (1. + @0*@0 + @1*@1)',             [ self._ReLambdaCP, self._ImLambdaCP ] )
                             , S = FormulaVar('S', '2. * @1 / (1. + @0*@0 + @1*@1)',              [ self._ReLambdaCP, self._ImLambdaCP ] )
                        )
        self._check_extraneous_kw( kwargs )


class LambdaAbsArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import FormulaVar
        from math import pi

        self._parseArg( 'lambdaCP', kwargs,  Title = 'CPV param. |lambda|', Value = lambVal, Error = lambErr, MinMax = ( 0.,      5.    ) )
        self._parseArg( 'phiCP',    kwargs,  Title = 'CPV param. phi',      Value = phiVal,  Error = phiErr,  MinMax = (-RooInf, +RooInf) )

        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0*@0) / (1. + @0*@0)',       [ self._lambdaCP              ] )
                             , D = FormulaVar('D', '-2. * @0 * cos(@1) / (1. + @0*@0)', [ self._lambdaCP, self._phiCP ] )
                             , S = FormulaVar('S', '-2. * @0 * sin(@1) / (1. + @0*@0)', [ self._lambdaCP, self._phiCP ] )
                        )
        self._check_extraneous_kw( kwargs )


class LambdaSqArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import FormulaVar
        from math import pi

        self._parseArg( 'lambdaCPSq', kwargs,  Title = 'CPV param. |lambda|^2', Value = lambSqVal, Error = lambSqErr
                       , MinMax = ( 0., 25. ) )
        self._parseArg( 'phiCP',      kwargs,  Title = 'CPV param. phi',        Value = phiVal,    Error = phiErr
                       , MinMax = ( -RooInf, +RooInf ) )

        CPParam.__init__(self, C = FormulaVar('C', '(1. - @0) / (1. + @0)',                [ self._lambdaCPSq              ] )
                             , D = FormulaVar('D', '-2. * sqrt(@0) * cos(@1) / (1. + @0)', [ self._lambdaCPSq, self._phiCP ] )
                             , S = FormulaVar('S', '-2. * sqrt(@0) * sin(@1) / (1. + @0)', [ self._lambdaCPSq, self._phiCP ] )
                        )
        self._check_extraneous_kw( kwargs )

class LambdaArg_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import FormulaVar,ConstVar
        from math import pi

        self._parseArg( 'phiCP', kwargs, Title = 'CPV param. phi', Value = phiVal, Error = phiErr, MinMax = ( -RooInf, +RooInf ) )
        CPParam.__init__(self, C = ConstVar( Name = 'C', Value = 0. )
                             , D = FormulaVar('D', '-cos(@0) ', [ self._phiCP ] )
                             , S = FormulaVar('S', '-sin(@0) ', [ self._phiCP ] )
                        )
        self._check_extraneous_kw( kwargs )

class CPVDecay_CPParam( CPParam ) :
    def __init__( self, **kwargs ) :
        ampNames   = kwargs.pop('AmplitudeNames')
        amps       = kwargs.pop('Amplitudes')
        magNames   = kwargs.pop('MagnitudeNames')
        phaseNames = kwargs.pop('PhaseNames')

        RPlusDict = { }
        RMinDict  = { }
        RReDict   = { }
        RImDict   = { }
        try :   from itertools import combinations_with_replacement as cwr  # this requires python 2.7 or later
        except: from P2VV.Compatibility import cwr
        from P2VV.RooFitWrappers import FormulaVar
        for ampComb in cwr( ampNames, 2 ) :
            eta0 = amps[ ampComb[0] ].CP
            eta1 = amps[ ampComb[1] ].CP

            RPlusDict[ampComb] = tuple( [ FormulaVar(  '%sRPlus_%s_%s' % ( comp, ampComb[0], ampComb[1] )
                                                     , '(%s@0*@1*%s(@2-@3)) / 2.'\
                                                       % ( '1. + ' if comp == 'Re' and eta0 == eta1\
                                                           else '1. - ' if comp == 'Re' and eta0 != eta1\
                                                           else '' if eta0 == eta1 else '-'
                                                          , func
                                                         )
                                                     , [  getattr( self, '_%s' % magNames[   ampComb[0] ] )
                                                        , getattr( self, '_%s' % magNames[   ampComb[1] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[0] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[1] ] )
                                                       ]
                                                    ) for ( comp, func ) in [ ( 'Re', 'cos' ), ( 'Im', 'sin' ) ]
                                      ] )
            RMinDict[ampComb]  = tuple( [ FormulaVar(  '%sRMin_%s_%s' % ( comp, ampComb[0], ampComb[1] )
                                                     , '(%s@0*@1*%s(@2-@3)) / 2.'\
                                                       % ( '1. - ' if comp == 'Re' and eta0 == eta1\
                                                           else '1. + ' if comp == 'Re' and eta0 != eta1\
                                                           else '-' if eta0 == eta1 else ''
                                                          , func
                                                         )
                                                     , [  getattr( self, '_%s' % magNames[   ampComb[0] ] )
                                                        , getattr( self, '_%s' % magNames[   ampComb[1] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[0] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[1] ] )
                                                       ]
                                                    ) for ( comp, func ) in [ ( 'Re', 'cos' ), ( 'Im', 'sin' ) ]
                                      ] )
            RReDict[ampComb]   = tuple( [ FormulaVar(  '%sRRe_%s_%s' % ( comp, ampComb[0], ampComb[1] )
                                                     , '(%s@0*%s(@2) %s @1*%s(-@3)) / 2.'\
                                                       % ( '' if eta0 > 0 else '-', func, '+' if eta1 > 0 else '-', func )
                                                     , [  getattr( self, '_%s' % magNames[   ampComb[0] ] )
                                                        , getattr( self, '_%s' % magNames[   ampComb[1] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[0] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[1] ] )
                                                       ]
                                                    ) for ( comp, func ) in [ ( 'Re', 'cos' ), ( 'Im', 'sin' ) ]
                                      ] )
            RImDict[ampComb]   = tuple( [ FormulaVar(  '%sRIm_%s_%s' % ( comp, ampComb[0], ampComb[1] )
                                                     , '(%s@0*%s(@2) %s @1*%s(-@3)) / 2.'\
                                                       % (  '-' if ( eta0 > 0 and comp == 'Re' ) or ( eta0 < 0 and comp == 'Im' ) else ''
                                                          , func
                                                          , '+' if ( eta1 > 0 and comp == 'Re' ) or ( eta1 < 0 and comp == 'Im' ) else '-'
                                                          , func
                                                         )
                                                     , [  getattr( self, '_%s' % magNames[   ampComb[0] ] )
                                                        , getattr( self, '_%s' % magNames[   ampComb[1] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[0] ] )
                                                        , getattr( self, '_%s' % phaseNames[ ampComb[1] ] )
                                                       ]
                                                    ) for ( comp, func ) in [ ( 'Re', 'sin' ), ( 'Im', 'cos' ) ]
                                      ] )

        self._check_extraneous_kw( kwargs )
        CPParam.__init__(  self, AmplitudeNames = ampNames
                         , C = FormulaVar( 'C', '(1. - @0*@0) / (1. + @0*@0)',       [  getattr( self, '_%s' % magNames['mix'] ) ] )
                         , D = FormulaVar( 'D', '-2. * @0 * cos(@1) / (1. + @0*@0)', [  getattr( self, '_%s' % magNames['mix']   )
                                                                                      , getattr( self, '_%s' % phaseNames['mix'] )
                                                                                     ]
                                         )
                         , S = FormulaVar( 'S', '-2. * @0 * sin(@1) / (1. + @0*@0)', [  getattr( self, '_%s' % magNames['mix']   )
                                                                                      , getattr( self, '_%s' % phaseNames['mix'] )
                                                                                     ]
                                         )
                         , RPlusDict = RPlusDict
                         , RMinDict  = RMinDict
                         , RReDict   = RReDict
                         , RImDict   = RImDict
                        )

class LambdaAbsArg_CPVDecay_CPParam( CPVDecay_CPParam ) :
    def __init__( self, **kwargs ) :
        ampNames   = kwargs.pop('AmplitudeNames')
        amps       = kwargs.pop('Amplitudes')
        magNames   = dict( mix = 'rhoCP_m'    )
        phaseNames = dict( mix = 'phiCPRel_m' )
        self._parseArg( 'rhoCP_m', kwargs, Title = 'CPV in mixing param. |rho_m|', Value = 1., Error = lambErr, MinMax = ( 0., 5. )
                       , Constant = True )
        self._parseArg( 'phiCPRel_m', kwargs, Title = 'CPV in mixing param. phi_m - phi_m', Value = 0., ObjectType = 'ConstVar' )

        for amp in ampNames :
            magNames[amp]   = 'rhoCP_%s' % amp
            phaseNames[amp] = 'phiCP_%s' % amp
            self._parseArg( 'rhoCP_%s' % amp, kwargs, Title = 'CPV in decay param. |rho_%s|' % amp, Value = lambVal, Error = lambErr
                           , MinMax = ( 0., 5. ) )
            self._parseArg( 'phiCP_%s' % amp, kwargs, Title = 'CPV in decay param. phi_%s + phi_m' % amp, Value = 0., Error = phiErr
                           , MinMax = ( -RooInf, +RooInf ) )

        self._check_extraneous_kw( kwargs )
        CPVDecay_CPParam.__init__( self, AmplitudeNames = ampNames, Amplitudes = amps, MagnitudeNames = magNames, PhaseNames = phaseNames )

class LambdaAbsArgRel_CPVDecay_CPParam( CPVDecay_CPParam ) :
    def __init__( self, **kwargs ) :
        ampNames   = kwargs.pop('AmplitudeNames')
        amps       = kwargs.pop('Amplitudes')
        magNames   = dict( mix = 'rhoCP_m' )
        phaseNames = dict( mix = 'phiCP_m' )
        self._parseArg( 'rhoCP_m', kwargs, Title = 'CPV in mixing param. |rho_m|', Value = 1., Error = lambErr, MinMax = ( 0., 5. )
                       , Constant = True )
        self._parseArg( 'phiCP_m', kwargs, Title = 'CPV in mixing param. phi_m + phi_0', Value = phiVal,  Error = phiErr
                       , MinMax = ( -RooInf, +RooInf ) )

        for ind, amp in enumerate(ampNames) :
            magNames[amp]   = 'rhoCP_%s'    % amp
            phaseNames[amp] = 'phiCPRel_%s' % amp
            self._parseArg( 'rhoCP_%s' % amp, kwargs, Title = 'CPV in decay param. |rho_%s|' % amp, Value = 1., Error = lambErr
                           , MinMax = ( 0., 5. ) )
            if ind == 0 :
                self._parseArg( 'phiCPRel_%s' % amp, kwargs, Title = 'CPV in decay param. phi_%s - phi_0' % amp, Value = 0.
                               , ObjectType = 'ConstVar' )
            else :
                self._parseArg( 'phiCPRel_%s' % amp, kwargs, Title = 'CPV in decay param. phi_%s - phi_0' % amp, Value = 0., Error = phiErr
                               , MinMax = ( -RooInf, +RooInf ) )

        self._check_extraneous_kw( kwargs )
        CPVDecay_CPParam.__init__( self, AmplitudeNames = ampNames, Amplitudes = amps, MagnitudeNames = magNames, PhaseNames = phaseNames )
