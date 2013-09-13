###########################################################################################################################################
## P2VVParameterizations.DecayAmplitudes: Decay amplitude parameterizations                                                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin, _util_conditionalObs_mixin

# amplitude values
from math import sqrt
A02    = 0.52
Aperp2 = 0.25
Apar2  = 1. - A02 - Aperp2

A0Ph    = 0.
AperpPh = 3.0
AparPh  = 3.3

f_S  = 0.03
AS2  = f_S / ( 1. - f_S )
ASPh = AperpPh

A02Err     = 0.007
Aperp2Err  = 0.01
AperpPhErr = 0.3
AparPhErr  = 0.2

f_SErr     = 0.01
AS2Err     = 0.01
ASPhErr    = 0.2
ASOddPhErr = 0.2

C_SP    = 1.
C_SPErr = 0.5

from ROOT import RooNumber
RooInf = RooNumber.infinity()

class Amplitude :
    def __init__( self, name, Re, Im, Mag2, Phase, CP ) :
        self.name  = name
        self.Re    = Re
        self.Im    = Im
        self.Mag2  = Mag2
        self.Phase = Phase
        self.CP    = CP   # even or odd???

    def __str__(self) : return self.name


# construct amplitudes with carthesian parameters
class Carthesian_Amplitude ( Amplitude ) :
    def __init__( self, name, Re, Im, CP, Mag2 = None, Phase = None, ParNamePrefix = '' ) :
        from P2VV.RooFitWrappers import ConvertPolAmp
        Amplitude.__init__(  self, name
                           , Re
                           , Im
                           , Mag2  if Mag2  else ConvertPolAmp( Name = ParNamePrefix + '%sMag2'  % name, Title = '|%s|^2'  % name
                                                               , Type = 'CarthToMagSq', Arguments = [ Re, Im ] )
                           , Phase if Phase else ConvertPolAmp( Name = ParNamePrefix + '%sPhase' % name, Title = 'arg(%s)' % name
                                                               , Type = 'CarthToPhase', Arguments = [ Re, Im ] )
                           , CP
                          )


# construct amplitudes with polar parameters
class Polar2_Amplitude( Amplitude ) :
    def __init__( self, name, Mag2, Phase, CP, Re = None, Im = None, ParNamePrefix = '' ) :
        from P2VV.RooFitWrappers import ConvertPolAmp
        Amplitude.__init__(  self, name
                           , Re if Re else ConvertPolAmp( Name = ParNamePrefix + 'Re%s' % name, Title = 'Re(%s)' % name, Type = 'PolSqToRe'
                                                         , Arguments = [Mag2, Phase] )
                           , Im if Im else ConvertPolAmp( Name = ParNamePrefix + 'Im%s' % name, Title = 'Im(%s)' % name, Type = 'PolSqToIm'
                                                         , Arguments = [Mag2, Phase] )
                           , Mag2
                           , Phase
                           , CP
                          )


class AmplitudeSet( dict, _util_parse_mixin, _util_conditionalObs_mixin ) :
    def __init__( self, *args, **kwargs ) :
        # maybe make this thing readonly???
        for v in args: 
            self[ v.name ] = v
            assert hasattr( v, 'Re'    )
            assert hasattr( v, 'Im'    )
            assert hasattr( v, 'Mag2'  )
            assert hasattr( v, 'Phase' )
            assert hasattr( v, 'CP'    )
        # require the names in args to be unique...
        assert( len(self)==len(args) )

        _util_conditionalObs_mixin.__init__( self, kwargs )
        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw ) :
        if kw in self : return dict.__getitem__( self, kw )
        return getattr( self, '_' + kw )

    def _createBinnedAmp( self, name, kwargs, title, initVals, minMaxs, category ) :
        # get initial values, minima and maxima for bin coefficients
        if hasattr( initVals, '__iter__' ) and len(initVals) == category.numTypes() : vals = initVals
        else : vals = category.numTypes() * [ initVals ]
        if hasattr( minMaxs, '__iter__' ) and len(minMaxs) == category.numTypes() and type(minMaxs) != tuple : minMax = minMaxs
        else : minMax = category.numTypes() * [ minMaxs ]

        # create bin coefficients
        coefs = [ ]
        binIt = 0
        stateIt = category.typeIterator()
        state = stateIt.Next()
        while state :
            ind = state.getVal()
            self._parseArg( '%s_bin%d' % ( name, ind ), kwargs, Title = '%s %d' % ( title, ind ), Value = vals[binIt]
                           , MinMax = minMax[binIt], ContainerList = coefs )
            binIt += 1
            state = stateIt.Next()
        setattr( self, '_%sCoefs' % name, coefs )

        # create binned amplitude
        self._parseArg( name, kwargs, ObjectType = 'CategoryVar', Category = category, Variables = coefs )


class JpsiVCarthesian_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters', False )
        namePF        = self.getNamePrefix(kwargs)

        from math import sqrt, cos, sin
        ReApar  = sqrt( Apar2  / A02 ) * cos(AparPh  - A0Ph)
        ImApar  = sqrt( Apar2  / A02 ) * sin(AparPh  - A0Ph)
        ReAperp = sqrt( Aperp2 / A02 ) * cos(AperpPh - A0Ph)
        ImAperp = sqrt( Aperp2 / A02 ) * sin(AperpPh - A0Ph)
        ReAS    = sqrt( AS2 / A02 )    * cos(ASPh    - A0Ph)
        ImAS    = sqrt( AS2 / A02 )    * sin(ASPh    - A0Ph)
        if ambiguityPars :
            ImA0    = -ImA0
            ImApar  = -ImApar
            ReAperp = -ReAperp
            ImAS    = -ImAS

        self._parseArg( 'ReA0',    kwargs, Title = 'Re(A_0)',    Value = 1.,      Constant = True      )
        self._parseArg( 'ImA0',    kwargs, Title = 'Im(A_0)',    Value = 0.,      Constant = True      )
        self._parseArg( 'ReApar',  kwargs, Title = 'Re(A_par)',  Value = ReApar,  MinMax = ( -1., 1. ) )
        self._parseArg( 'ImApar',  kwargs, Title = 'Im(A_par)',  Value = ImApar,  MinMax = ( -1., 1. ) )
        self._parseArg( 'ReAperp', kwargs, Title = 'Re(A_perp)', Value = ReAperp, MinMax = ( -1., 1. ) )
        self._parseArg( 'ImAperp', kwargs, Title = 'Im(A_perp)', Value = ImAperp, MinMax = ( -1., 1. ) )

        self._KKMassCat = kwargs.pop( 'KKMassCategory', None )
        if self._KKMassCat :
            self._createBinnedAmp( 'ReAS', kwargs, 'Re(A_S)', ReAS, ( -1., 1. ), self._KKMassCat )
            self._createBinnedAmp( 'ImAS', kwargs, 'Im(A_S)', ImAS, ( -1., 1. ), self._KKMassCat )
        else :
            self._parseArg( 'ReAS', kwargs, Title = 'Re(A_S)', Value = ReAS, MinMax = ( -1., 1. ) )
            self._parseArg( 'ImAS', kwargs, Title = 'Im(A_S)', Value = ImAS, MinMax = ( -1., 1. ) )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, Carthesian_Amplitude( 'A0',    self._ReA0,    self._ImA0,    +1, ParNamePrefix = namePF )
                                   , Carthesian_Amplitude( 'Apar',  self._ReApar,  self._ImApar,  +1, ParNamePrefix = namePF )
                                   , Carthesian_Amplitude( 'Aperp', self._ReAperp, self._ImAperp, -1, ParNamePrefix = namePF )
                                   , Carthesian_Amplitude( 'AS',    self._ReAS,    self._ImAS,    -1, ParNamePrefix = namePF )
                                   , Conditionals = [ self._KKMassCat ] if self._KKMassCat else [ ]
                             )


class JpsiVPolar_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters', False       )
        ASParam       = kwargs.pop( 'ASParameterization',  'deltaPerp' )
        namePF        = self.getNamePrefix(kwargs)

        from math import pi
        deltaPar  = AparPh  - A0Ph
        deltaPerp = AperpPh - A0Ph
        deltaS    = ASPh    - ( AperpPh if ASParam in [ 'deltaPerp', 'Mag2ReImPerp' ] else A0Ph )
        if ambiguityPars :
            deltaPar  = -deltaPar
            deltaPerp = pi - deltaPerp
            deltaS    = ( pi - deltaS ) if ASParam in [ 'deltaPerp', 'Mag2ReImPerp' ] else -deltaS

        # A_0, A_par, A_perp
        self._parseArg( 'A0Mag2',    kwargs, Title = '|A0|^2',     Value = A02,    MinMax = ( 0., 1. ) )
        self._parseArg( 'AperpMag2', kwargs, Title = '|A_perp|^2', Value = Aperp2, MinMax = ( 0., 1. ) )

        PNorm = kwargs.pop( 'PWaveNorm', True )
        if PNorm :
            self._AparMag2 = self._parseArg( 'AparMag2', kwargs, Title = '|A_par|^2', Coefficients = [ self._A0Mag2, self._AperpMag2 ]
                                            , ObjectType = 'ComplementCoef' )
        else :
            self._AparMag2 = self._parseArg( 'AparMag2', kwargs, Title = '|A_par|^2'
                                            , Coefficients = [self._A0Mag2, self._AperpMag2, self._ASMag2], ObjectType = 'ComplementCoef' )

        self._parseArg( 'A0Phase',    kwargs, Title = 'delta_0',    Value = 0.,        Constant = True                )
        self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = deltaPar,  MinMax = ( -2. * pi, 2. * pi ) )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = deltaPerp, MinMax = ( -2. * pi, 2. * pi ) )

        # A_S
        self._KKMassCat = kwargs.pop( 'KKMassCategory', None )
        if self._KKMassCat :
            self._createBinnedAmp( 'ASMag2', kwargs, '|A_S|^2', AS2, ( 0.,  1. ), self._KKMassCat )
            if ASParam == 'deltaPerp' :
                self._createBinnedAmp( 'ASOddPhase', kwargs, 'delta_S - delta_perp', deltaS, ( -2. * pi, 2. * pi ), self._KKMassCat )
            elif ASParam == 'Mag2ReIm' :
                from math import sqrt, cos, sin
                self._createBinnedAmp( 'ReAS', kwargs, 'Re(A_S)', sqrt(AS2) * cos(deltaS), ( -1., 1. ), self._KKMassCat )
                self._createBinnedAmp( 'ImAS', kwargs, 'Im(A_S)', sqrt(AS2) * sin(deltaS), ( -1., 1. ), self._KKMassCat )
            elif ASParam == 'Mag2ReImPerp' :
                from math import sqrt, cos, sin
                self._createBinnedAmp( 'ReASPerp', kwargs, 'Re(A_S*e^(-i*delta_perp))', sqrt(AS2) * cos(deltaS), ( -1., 1. )
                                      , self._KKMassCat )
                self._createBinnedAmp( 'ImASPerp', kwargs, 'Im(A_S*e^(-i*delta_perp))', sqrt(AS2) * sin(deltaS), ( -1., 1. )
                                      , self._KKMassCat )
            else  :
                self._createBinnedAmp( 'ASPhase', kwargs, 'delta_S', deltaS, ( -2. * pi, 2. * pi ), self._KKMassCat )

            if not ASParam.startswith('Mag2ReIm') :
                self._createBinnedAmp( 'C_SP', kwargs, 'S-P wave couping factor', C_SP, ( 0., 2. ), self._KKMassCat )

        else :
            self._parseArg( 'ASMag2', kwargs, Title = '|A_S|^2', Value = AS2, MinMax = ( 0., 1. ) )
            if ASParam == 'deltaPerp' :
                self._parseArg( 'ASOddPhase', kwargs, Title = 'delta_S - delta_perp', Value = deltaS, MinMax = ( -2. * pi, 2. * pi ) )
            elif ASParam == 'Mag2ReIm' :
                from math import sqrt, cos, sin
                self._parseArg( 'ReAS', kwargs, Title = 'Re(A_S)', Value = sqrt(AS2) * cos(deltaS), MinMax = ( -1., 1. ) )
                self._parseArg( 'ImAS', kwargs, Title = 'Im(A_S)', Value = sqrt(AS2) * sin(deltaS), MinMax = ( -1., 1. ) )
            elif ASParam == 'Mag2ReImPerp' :
                from math import sqrt, cos, sin
                self._parseArg( 'ReASPerp', kwargs, Title = 'Re(A_S*e^(-i*delta_perp))', Value = sqrt(AS2) * cos(deltaS)
                               , MinMax = ( -1., 1. ) )
                self._parseArg( 'ImASPerp', kwargs, Title = 'Im(A_S*e^(-i*delta_perp))', Value = sqrt(AS2) * sin(deltaS)
                               , MinMax = ( -1., 1. ) )
            else :
                self._parseArg( 'ASPhase', kwargs, Title = 'delta_S', Value = deltaS, MinMax = ( -2. * pi, 2. * pi ) )

            if not ASParam.startswith('Mag2ReIm') :
                self._parseArg( 'C_SP', kwargs, Title = 'S-P wave couping factor', Value = C_SP, Error = C_SPErr, MinMax = ( 0., 2. ) )

        if ASParam == 'deltaPerp' :
            self._ReAS = self._parseArg( 'ReAS', kwargs, Type = 'PolSqToReRelC', Title = 'Re(A_S)', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._ASMag2, self._ASOddPhase, self._AperpPhase, self._C_SP ] )
            self._ImAS = self._parseArg( 'ImAS', kwargs, Type = 'PolSqToImRelC', Title = 'Im(A_S)', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._ASMag2, self._ASOddPhase, self._AperpPhase, self._C_SP ] )

        elif ASParam == 'Mag2ReImPerp' :
            self._ReAS = self._parseArg( 'ReAS', kwargs, Title = 'Re(A_S)', Type = 'ProdCarthPhaseToRe', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._ReASPerp, self._ImASPerp, self._AperpPhase ] )
            self._ImAS = self._parseArg( 'ImAS', kwargs, Title = 'Im(A_S)', Type = 'ProdCarthPhaseToIm', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._ReASPerp, self._ImASPerp, self._AperpPhase ] )

        elif ASParam != 'Mag2ReIm' :
            self._ReAS = self._parseArg( 'ReAS', kwargs, Title = 'Re(A_S)', Type = 'PolSqToReC', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._ASMag2, self._ASPhase, self._C_SP ] )
            self._ImAS = self._parseArg( 'ImAS', kwargs, Title = 'Im(A_S)', Type = 'PolSqToImC', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._ASMag2, self._ASPhase, self._C_SP ] )

        ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1, getattr( self, '_ASMag2', None ), getattr( self, '_ASPhase', None )
                                     , namePF )

        self._check_extraneous_kw( kwargs ) 
        AmplitudeSet.__init__( self, Polar2_Amplitude( 'A0',    self._A0Mag2,    self._A0Phase,    +1, ParNamePrefix = namePF )
                                   , Polar2_Amplitude( 'Apar',  self._AparMag2,  self._AparPhase,  +1, ParNamePrefix = namePF )
                                   , Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1, ParNamePrefix = namePF )
                                   , ASAmp
                                   , Conditionals = [ self._KKMassCat ] if self._KKMassCat else [ ]
                             )


class JpsiVPolarSWaveFrac_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters',  False       )
        ASParam       = kwargs.pop( 'ASParameterization',   'deltaPerp' )
        AparParam     = kwargs.pop( 'AparParameterization', 'phase'     )
        namePF        = self.getNamePrefix(kwargs)

        from math import pi
        deltaPar   = AparPh - A0Ph
        delParErr  = AparPhErr
        deltaPerp  = AperpPh - A0Ph
        delPerpErr = AperpPhErr
        deltaS     = ASPh - ( AperpPh if ASParam == 'deltaPerp' else A0Ph )
        delSErr    = ASOddPhErr if ASParam == 'deltaPerp' else ASPhErr
        if ambiguityPars :
            deltaPar  = -deltaPar
            deltaPerp = pi - deltaPerp
            deltaS    = ( pi - deltaS ) if ASParam == 'deltaPerp' else -deltaS

        from math import sqrt, sin, cos

        # A_0
        self._parseArg( 'A0Mag2',  kwargs, Title = '|A0|^2',  Value = A02, Error = A02Err, MinMax = ( 0., 1. ) )
        self._parseArg( 'A0Phase', kwargs, Title = 'delta_0', Value = 0.,  Constant = True )
        A0Amp = Polar2_Amplitude( 'A0', self._A0Mag2, self._A0Phase, +1, ParNamePrefix = namePF )

        # A_perp
        self._parseArg( 'AperpMag2',  kwargs, Title = '|A_perp|^2', Value = Aperp2,    Error = Aperp2Err,  MinMax = (  0.,     1.     ) )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = deltaPerp, Error = delPerpErr, MinMax = ( -RooInf, RooInf ) )
        AperpAmp = Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1, ParNamePrefix = namePF )

        # A_par
        self._AparMag2 = self._parseArg( 'AparMag2', kwargs, Title = '|A_par|^2', Coefficients = [ self._A0Mag2, self._AperpMag2 ]
                                        , ObjectType = 'ComplementCoef' )
        if AparParam in [ 'real', 'cos' ] :
            # Re(A_par)
            if AparParam == 'real' :
                self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Value = sqrt(Apar2) * cos(deltaPar), MinMax = ( -1., 1. ) )
            else :
                self._parseArg( 'cosAparPhase', kwargs, Title = 'cos(delta_par)', Value = cos(deltaPar), MinMax = ( -1., 1. ) )
                self._ReApar = self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Type = 'MagSqSinToReIm', ObjectType = 'ConvertPolAmp'
                                              , Arguments = [ self._AparMag2, self._cosAparPhase ] )

            # Im(A_par)
            self._ImAparSign = self._parseArg( 'ImAparSign', kwargs, States = { 'plus' : +1, 'minus' : -1 }, ObjectType = 'Category' )
            if ambiguityPars : self._ImAparSign.setIndex(+1)
            else             : self._ImAparSign.setIndex(-1)
            self._ImApar = self._parseArg( 'ImApar', kwargs, ObjectType = 'FormulaVar', Formula = '@2*sqrt(@0 - @1*@1)'
                                          , Arguments = [ self._AparMag2, self._ReApar, self._ImAparSign ], Title = 'Im(A_par)' )
            AparAmp = Carthesian_Amplitude( 'Apar', self._ReApar, self._ImApar, +1, ParNamePrefix = namePF )

        else :
            # delta_par
            self._parseArg( 'AparPhase', kwargs, Title = 'delta_par', Value = deltaPar, Error = delParErr, MinMax = ( -RooInf, RooInf ) )
            AparAmp = Polar2_Amplitude( 'Apar', self._AparMag2, self._AparPhase, +1, ParNamePrefix = namePF )

        # A_S
        self._KKMassCat = kwargs.pop( 'KKMassCategory', None )
        if self._KKMassCat :
            self._createBinnedAmp( 'C_SP', kwargs, 'S-P wave couping factor', C_SP, ( 0., 1. ), self._KKMassCat )
        else :
            self._parseArg( 'C_SP', kwargs, Title = 'S-P wave couping factor', Value = C_SP, Error = C_SPErr, MinMax = ( 0., 1. ) )

        if ASParam != 'ReIm' :
            if self._KKMassCat :
                self._createBinnedAmp( 'f_S', kwargs, 'S wave fraction', f_S, ( 0., 1. ), self._KKMassCat )
                if ASParam == 'deltaPerp' :
                    self._createBinnedAmp( 'ASOddPhase', kwargs, 'delta_S - delta_perp', deltaS, ( -2. * pi, 2. * pi ), self._KKMassCat )
                else :
                    self._createBinnedAmp( 'ASPhase', kwargs, 'delta_S', deltaS, ( -2. * pi, 2. * pi ), self._KKMassCat )

            else :
                self._parseArg( 'f_S', kwargs, Title = 'S wave fraction', Value = f_S, Error = f_SErr, MinMax = ( 0., 1. ) )
                if ASParam == 'deltaPerp' :
                    self._parseArg( 'ASOddPhase', kwargs, Title = 'delta_S - delta_perp', Value = deltaS, Error = delSErr
                                   , MinMax = ( -RooInf, RooInf ) )
                else :
                    self._parseArg( 'ASPhase', kwargs, Title = 'delta_S', Value = deltaS, Error = delSErr, MinMax = ( -RooInf, RooInf ) )

            self._ASMag2 = self._parseArg( 'ASMag2', kwargs, Title = 'Re(A_S)', Type = 'FracToMagSq', ObjectType = 'ConvertPolAmp'
                                          , Arguments = [ self._f_S ] )
            if ASParam == 'deltaPerp' :
                self._ReAS = self._parseArg( 'ReAS', kwargs, Type = 'FracPolToReRelC', Title = 'Re(A_S)', ObjectType = 'ConvertPolAmp'
                                            , Arguments = [ self._f_S, self._ASOddPhase, self._AperpPhase, self._C_SP ] )
                self._ImAS = self._parseArg( 'ImAS', kwargs, Type = 'FracPolToImRelC', Title = 'Im(A_S)', ObjectType = 'ConvertPolAmp'
                                            , Arguments = [ self._f_S, self._ASOddPhase, self._AperpPhase, self._C_SP ] )

            else :
                self._ReAS = self._parseArg( 'ReAS', kwargs, Title = 'Re(A_S)', Type = 'FracPolToReC', ObjectType = 'ConvertPolAmp'
                                            , Arguments = [ self._f_S, self._ASPhase, self._C_SP ] )
                self._ImAS = self._parseArg( 'ImAS', kwargs, Title = 'Im(A_S)', Type = 'FracPolToImC', ObjectType = 'ConvertPolAmp'
                                            , Arguments = [ self._f_S, self._ASPhase, self._C_SP ] )

        else :
            if self._KKMassCat :
                self._createBinnedAmp( 'sqrtfS_Re', kwargs, 'sqrt(S wave fraction) * cos(delta_S)', sqrt(f_S) * cos(deltaS), ( -1., 1. )
                                      , self._KKMassCat )
                self._createBinnedAmp( 'sqrtfS_Im', kwargs, 'sqrt(S wave fraction) * sin(delta_S)', sqrt(f_S) * sin(deltaS), ( -1., 1. )
                                      , self._KKMassCat )
            else :
                self._parseArg( 'sqrtfS_Re', kwargs, Title = 'sqrt(S wave fraction) * cos(delta_S)', Value = sqrt(f_S) * cos(deltaS)
                               , MinMax = ( -1., 1. ) )
                self._parseArg( 'sqrtfS_Im', kwargs, Title = 'sqrt(S wave fraction) * sin(delta_S)', Value = sqrt(f_S) * sin(deltaS)
                               , MinMax = ( -1., 1. ) )

            self._ASMag2 = self._parseArg( 'ASMag2', kwargs, Title = '|A_S|^2', Type = 'CarthFracToMagSq', ObjectType = 'ConvertPolAmp'
                                          , Arguments = [ self._sqrtfS_Re, self._sqrtfS_Im ] )
            self._ReAS = self._parseArg( 'ReAS', kwargs, Title = 'Re(A_S)', Type = 'FracCarthToReRelC', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._sqrtfS_Re, self._sqrtfS_Im, self._C_SP ] )
            self._ImAS = self._parseArg( 'ImAS', kwargs, Title = 'Im(A_S)', Type = 'FracCarthToImRelC', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._sqrtfS_Re, self._sqrtfS_Im, self._C_SP ] )

        ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1, getattr(self, '_ASMag2', None), getattr(self, '_ASPhase', None)
                                     , namePF )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, A0Amp, AparAmp, AperpAmp, ASAmp, Conditionals = [ self._KKMassCat ] if self._KKMassCat else [ ] )


class JpsiVBank_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters',  False       )
        ASParam       = kwargs.pop( 'ASParameterization',   'deltaPerp' )
        AparParam     = kwargs.pop( 'AparParameterization', 'phase'     )
        namePF        = self.getNamePrefix(kwargs)

        from math import pi
        deltaPar  = AparPh  - A0Ph
        deltaPerp = AperpPh - A0Ph
        deltaS    = ASPh    - AperpPh
        if ambiguityPars :
            deltaPar  = -deltaPar
            deltaPerp = pi - deltaPerp
            deltaS    = pi - deltaS

        from math import sqrt, sin, cos

        # A_0
        self._parseArg( 'ReA0', kwargs, Title = 'Re(A_0)', Value = 1., Constant = True )
        self._parseArg( 'ImA0', kwargs, Title = 'Im(A_0)', Value = 0., Constant = True )
        A0Amp = Carthesian_Amplitude( 'A0', self._ReA0, self._ImA0, +1, ParNamePrefix = namePF )

        # A_perp
        self._parseArg( 'AperpMag2',  kwargs, Title = '|A_perp|^2', Value = Aperp2 / A02, MinMax = ( 0., 1. )            )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = deltaPerp,    MinMax = ( -2. * pi, 2. * pi ) )
        AperpAmp = Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1, ParNamePrefix = namePF )

        # A_par
        if AparParam != 'ReIm' :
            self._parseArg( 'AparMag2', kwargs, Title = '|A_par|^2', Value = Apar2 / A02, MinMax = ( 0., 1. ) )

        if AparParam in [ 'Mag2ReIm', 'ReIm', 'real', 'cos' ] :
            # Re(A_par)
            if AparParam in [ 'Mag2ReIm', 'ReIm', 'real' ] :
                self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Value = sqrt(Apar2 / A02) * cos(deltaPar), MinMax = (-1., 1.) )
            else :
                self._parseArg( 'cosAparPhase', kwargs, Title = 'cos(delta_par)', Value = cos(deltaPar), MinMax = ( -1., 1. ) )
                self._ReApar = self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Type = 'MagSqSinToReIm'
                                              , ObjectType = 'ConvertPolAmp', Arguments = [ self._AparMag2, self._cosAparPhase ] )

            # Im(A_par)
            if AparParam not in [ 'Mag2ReIm', 'ReIm' ] :
                self._ImAparSign = self._parseArg( 'ImAparSign', kwargs, States = { 'plus' : +1, 'minus' : -1 }, ObjectType = 'Category' )
                if ambiguityPars : self._ImAparSign.setIndex(+1)
                else             : self._ImAparSign.setIndex(-1)
                self._ImApar = self._parseArg( 'ImApar', kwargs, ObjectType = 'FormulaVar', Formula = '@2*sqrt(@0 - @1*@1)'
                                              , Arguments = [ self._AparMag2, self._ReApar, self._ImAparSign ], Title = 'Im(A_par)' )
            else :
                self._parseArg( 'ImApar', kwargs, Title = 'Im(A_par)', Value = sqrt(Apar2 / A02) * sin(deltaPar), MinMax = (-1., 1.) )

            if AparParam == 'ReIm' :
                self._AparMag2 = self._parseArg( 'AparMag2', kwargs, Title = '|A_par|^2', Type = 'CarthToMagSq'
                                                , ObjectType = 'ConvertPolAmp', Arguments = [ self._ReApar, self._ImApar ] )

            AparAmp = Carthesian_Amplitude( 'Apar', self._ReApar, self._ImApar, +1, self._AparMag2, ParNamePrefix = namePF )

        else :
            # delta_par
            self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = deltaPar,  MinMax = ( -2. * pi, 2. * pi ) )
            AparAmp = Polar2_Amplitude( 'Apar', self._AparMag2, self._AparPhase, +1, ParNamePrefix = namePF )

        # A_S (--> A_S / A_perp and not A_S / A_0!!!)
        self._KKMassCat = kwargs.pop( 'KKMassCategory', None )
        if self._KKMassCat :
            self._createBinnedAmp( 'C_SP', kwargs, 'S-P wave couping factor', C_SP, ( 0., 2. ), self._KKMassCat )
        else :
            self._parseArg( 'C_SP', kwargs, Title = 'S-P wave couping factor', Value = C_SP, Error = C_SPErr, MinMax = ( 0., 2. ) )

        if ASParam != 'ReIm' :
            if self._KKMassCat :
                self._createBinnedAmp( 'ASOddMag2',  kwargs, '|A_S|^2 / |A_perp|^2', AS2 / Aperp2, (  0.,      1.      ), self._KKMassCat )
                self._createBinnedAmp( 'ASOddPhase', kwargs, 'delta_S - delta_perp', deltaS,       ( -2. * pi, 2. * pi ), self._KKMassCat )
            else :
                self._parseArg( 'ASOddMag2',  kwargs, Title = '|A_S|^2 / |A_perp|^2', Value = AS2 / Aperp2
                               , MinMax = (  0.,      1.      ) )
                self._parseArg( 'ASOddPhase', kwargs, Title = 'delta_S - delta_perp', Value = deltaS
                               , MinMax = ( -2. * pi, 2. * pi ) )

            self._ReAS = self._parseArg( 'ReAS', kwargs, Type = 'PolSqToReRelMagC', Title = 'Re(A_S)', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [self._AperpMag2, self._AperpPhase, self._ASOddMag2, self._ASOddPhase, self._C_SP] )
            self._ImAS = self._parseArg( 'ImAS', kwargs, Type = 'PolSqToImRelMagC', Title = 'Im(A_S)', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [self._AperpMag2, self._AperpPhase, self._ASOddMag2, self._ASOddPhase, self._C_SP] )

        else :
            if self._KKMassCat :
                self._createBinnedAmp( 'ReASOdd', kwargs, 'Re(A_S / A_perp)', sqrt(AS2/Aperp2) * cos(deltaS), (-1., 1.), self._KKMassCat )
                self._createBinnedAmp( 'ImASOdd', kwargs, 'Im(A_S / A_perp)', sqrt(AS2/Aperp2) * sin(deltaS), (-1., 1.), self._KKMassCat )
            else :
                self._parseArg( 'ReASOdd', kwargs, Title = 'Re(A_S / A_perp)', Value = sqrt(AS2 / Aperp2) * cos(deltaS)
                               , MinMax = ( -1., 1. ) )
                self._parseArg( 'ImASOdd', kwargs, Title = 'Im(A_S / A_perp)', Value = sqrt(AS2 / Aperp2) * sin(deltaS)
                               , MinMax = ( -1., 1. ) )

            self._ReAS = self._parseArg( 'ReAS', kwargs, Type = 'MixSqToReRelMagC', Title = 'Re(A_S)', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._AperpMag2, self._AperpPhase, self._ReASOdd, self._ImASOdd, self._C_SP ] )
            self._ImAS = self._parseArg( 'ImAS', kwargs, Type = 'MixSqToImRelMagC', Title = 'Im(A_S)', ObjectType = 'ConvertPolAmp'
                                        , Arguments = [ self._AperpMag2, self._AperpPhase, self._ReASOdd, self._ImASOdd, self._C_SP ] )

        ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1, ParNamePrefix = namePF )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, A0Amp, AparAmp, AperpAmp, ASAmp, Conditionals = [ self._KKMassCat ] if self._KKMassCat else [ ] )
