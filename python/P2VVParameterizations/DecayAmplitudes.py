###########################################################################################################################################
## P2VVParameterizations.DecayAmplitudes: Decay amplitude parameterizations                                                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin

# amplitude values
A02     = 0.522
Aperp2  = 0.249
Apar2   = 1. - A02 - Aperp2

A0Ph    = 0.
AperpPh = 2.77
AparPh  = 3.32

f_S     = 0.016
AS2     = f_S / ( 1. - f_S )
ASPh    = 2.74


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
    def __init__( self, name, Re, Im, CP, Mag2 = None, Phase = None ) :
        from RooFitWrappers import FormulaVar
        Amplitude.__init__(  self, name
                           , Re
                           , Im
                           , Mag2  if Mag2  else FormulaVar( '%sMag2'  % name, '@0*@0+@1*@1',  [ Re, Im ], Title = '|%s|^2'  % name )
                           , Phase if Phase else FormulaVar( '%sPhase' % name, 'atan2(@1,@0)', [ Re, Im ], Title = 'arg(%s)' % name )
                           , CP
                          )


# construct amplitudes with polar parameters
class Polar2_Amplitude( Amplitude ) :
    def __init__( self, name, Mag2, Phase, CP, Re = None, Im = None ) :
        from RooFitWrappers import FormulaVar
        Amplitude.__init__(  self, name
                           , Re if Re else FormulaVar('Re%s' % name, 'sqrt(@0) * cos(@1)', [Mag2, Phase], Title = 'Re(%s)' % name )
                           , Im if Im else FormulaVar('Im%s' % name, 'sqrt(@0) * sin(@1)', [Mag2, Phase], Title = 'Im(%s)' % name )
                           , Mag2
                           , Phase
                           , CP
                          )


class AmplitudeSet( dict, _util_parse_mixin ) :
    def __init__( self, *args ) :
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

    def __getitem__( self, kw ) :
        if kw in self : return dict.__getitem__( self, kw )
        return getattr( self, '_' + kw )


class JpsiVCarthesian_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        from math import sqrt, cos, sin
        self._parseArg( 'ReA0',    kwargs, Title = 'Re(A_0)',    Value = 1.,                                         Constant = True      )
        self._parseArg( 'ImA0',    kwargs, Title = 'Im(A_0)',    Value = 0.,                                         Constant = True      )
        self._parseArg( 'ReApar',  kwargs, Title = 'Re(A_par)',  Value = sqrt( Apar2  / A02 ) * cos(AparPh  - A0Ph), MinMax = ( -1., 1. ) )
        self._parseArg( 'ImApar',  kwargs, Title = 'Im(A_par)',  Value = sqrt( Apar2  / A02 ) * sin(AparPh  - A0Ph), MinMax = ( -1., 1. ) )
        self._parseArg( 'ReAperp', kwargs, Title = 'Re(A_perp)', Value = sqrt( Aperp2 / A02 ) * cos(AperpPh - A0Ph), MinMax = ( -1., 1. ) )
        self._parseArg( 'ImAperp', kwargs, Title = 'Im(A_perp)', Value = sqrt( Aperp2 / A02 ) * sin(AperpPh - A0Ph), MinMax = ( -1., 1. ) )
        self._parseArg( 'ReAS',    kwargs, Title = 'Re(A_S)',    Value = sqrt( AS2    / A02 ) * cos(ASPh    - A0Ph), MinMax = ( -1., 1. ) )
        self._parseArg( 'ImAS',    kwargs, Title = 'Im(A_S)',    Value = sqrt( AS2    / A02 ) * sin(ASPh    - A0Ph), MinMax = ( -1., 1. ) )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, Carthesian_Amplitude( 'A0',    self._ReA0,    self._ImA0,    +1 )
                                   , Carthesian_Amplitude( 'Apar',  self._ReApar,  self._ImApar,  +1 )
                                   , Carthesian_Amplitude( 'Aperp', self._ReAperp, self._ImAperp, -1 )
                                   , Carthesian_Amplitude( 'AS',    self._ReAS,    self._ImAS,    -1 )
                             )


class JpsiVPolar_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        from math import pi
        from RooFitWrappers import FormulaVar
        self._parseArg( 'A0Mag2',    kwargs, Title = '|A0|^2',     Value = A02,    MinMax = ( 0., 1. ) )
        self._parseArg( 'AperpMag2', kwargs, Title = '|A_perp|^2', Value = Aperp2, MinMax = ( 0., 1. ) )
        self._parseArg( 'ASMag2',    kwargs, Title = '|A_S|^2',    Value = AS2,    MinMax = ( 0., 1. ) )

        PNorm = kwargs.pop( 'PWaveNorm', True )
        if PNorm :
            self._AparMag2 = FormulaVar( 'AparMag2', '1.-@0-@1',    [ self._A0Mag2, self._AperpMag2 ],               Title = '|A_par|^2' )
        else :
            self._AparMag2 = FormulaVar( 'AparMag2', '1.-@0-@1-@2', [ self._A0Mag2, self._AperpMag2, self._ASMag2 ], Title = '|A_par|^2' )

        self._parseArg( 'A0Phase',    kwargs, Title = 'delta_0',    Value = 0.,             Constant = True                )
        self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = AparPh  - A0Ph, MinMax = ( -2. * pi, 2. * pi ) )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = AperpPh - A0Ph, MinMax = ( -2. * pi, 2. * pi ) )
        self._parseArg( 'ASPhase',    kwargs, Title = 'delta_S',    Value = ASPh    - A0Ph, MinMax = ( -2. * pi, 2. * pi ) )

        self._check_extraneous_kw( kwargs ) 
        AmplitudeSet.__init__( self, Polar2_Amplitude( 'A0',    self._A0Mag2,    self._A0Phase,    +1 )
                                   , Polar2_Amplitude( 'Apar',  self._AparMag2,  self._AparPhase,  +1 )
                                   , Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )
                                   , Polar2_Amplitude( 'AS',    self._ASMag2,    self._ASPhase,    -1 )
                             )


class JpsiVPolarSWaveFrac_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        from math import pi, sqrt, sin, cos
        from RooFitWrappers import FormulaVar
        polarSWave = kwargs.pop( 'PolarSWave', True )
        AparParam  = kwargs.pop( 'AparParameterization', 'phase' )

        # A_0
        self._parseArg( 'A0Mag2',  kwargs, Title = '|A0|^2',  Value = A02, MinMax = ( 0., 1. ) )
        self._parseArg( 'A0Phase', kwargs, Title = 'delta_0', Value = 0.,  Constant = True     )
        A0Amp = Polar2_Amplitude( 'A0', self._A0Mag2, self._A0Phase, +1 )

        # A_perp
        self._parseArg( 'AperpMag2',  kwargs, Title = '|A_perp|^2', Value = Aperp2,         MinMax = ( 0., 1. )            )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = AperpPh - A0Ph, MinMax = ( -2. * pi, 2. * pi ) )
        AperpAmp = Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )

        # A_par
        self._AparMag2 = FormulaVar( 'AparMag2', '1. - @0 - @1', [ self._A0Mag2, self._AperpMag2 ], Title = '|A_par|^2' )
        if AparParam in [ 'real', 'cos' ] :
            # Re(A_par)
            if AparParam == 'real' :
                self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Value = sqrt(Apar2) * cos(AparPh - A0Ph), MinMax = ( -1., 1. ) )
            else :
                self._parseArg( 'cosAparPhase', kwargs, Title = 'cos(delta_par)', Value = cos(AparPh - A0Ph), MinMax = ( -1., 1. ) )
                self._ReApar = FormulaVar( 'ReApar', 'sqrt(@0)*@1', [ self._AparMag2, self._cosAparPhase ], Title = 'Re(A_par)' )

            # Im(A_par)
            from RooFitWrappers import Category
            self._ImAparSign = Category( 'ImAparSign', States = { 'plus' : +1, 'minus' : -1 } )
            self._ImAparSign.setIndex(-1)
            self._ImApar = FormulaVar( 'ImApar', '@2*sqrt(@0 - @1*@1)', [ self._AparMag2, self._ReApar, self._ImAparSign ]
                                      , Title = 'Im(A_par)' )
            AparAmp = Carthesian_Amplitude( 'Apar', self._ReApar, self._ImApar, +1 )

        else :
            # delta_par
            self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = AparPh - A0Ph, MinMax = ( -2. * pi, 2. * pi ) )
            AparAmp = Polar2_Amplitude( 'Apar', self._AparMag2, self._AparPhase, +1 )

        # A_S
        if polarSWave :
            self._parseArg( 'f_S',     kwargs, Title = 'S wave fraction', Value = f_S,         MinMax = (  0.,      1.      ) )
            self._parseArg( 'ASPhase', kwargs, Title = 'delta_S',         Value = ASPh - A0Ph, MinMax = ( -2. * pi, 2. * pi ) )
            self._ASMag2 = FormulaVar( 'ASMag2', '@0 / (1. - @0)', [ self._f_S ], Title = '|A_S|^2' )
            ASAmp = Polar2_Amplitude( 'AS', self._ASMag2, self._ASPhase, -1 )

        else :
            self._parseArg( 'sqrtfS_Re', kwargs, Title = 'sqrt(S wave fraction) * cos(delta_S)', Value = sqrt(f_S) * cos(ASPh - A0Ph)
                           , MinMax = ( -1., 1. ) )
            self._parseArg( 'sqrtfS_Im', kwargs, Title = 'sqrt(S wave fraction) * sin(delta_S)', Value = sqrt(f_S) * sin(ASPh - A0Ph)
                           , MinMax = ( -1., 1. ) )
            self._ReAS = FormulaVar( 'ReAS', '@0 / sqrt(1. - @0*@0 - @1*@1)', [ self._sqrtfS_Re, self._sqrtfS_Im ], Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@1 / sqrt(1. - @0*@0 - @1*@1)', [ self._sqrtfS_Re, self._sqrtfS_Im ], Title = 'Im(A_S)' )

            ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1 )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, A0Amp, AparAmp, AperpAmp, ASAmp )


class JpsiVBank_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        from math import pi, sqrt, sin, cos
        from RooFitWrappers import FormulaVar
        polarSWave = kwargs.pop( 'PolarSWave',           True    )
        AparParam  = kwargs.pop( 'AparParameterization', 'phase' )

        # A_0
        self._parseArg( 'ReA0', kwargs, Title = 'Re(A_0)', Value = 1., Constant = True )
        self._parseArg( 'ImA0', kwargs, Title = 'Im(A_0)', Value = 0., Constant = True )
        A0Amp = Carthesian_Amplitude( 'A0', self._ReA0, self._ImA0, +1 )

        # A_perp
        self._parseArg( 'AperpMag2',  kwargs, Title = '|A_perp|^2', Value = Aperp2 / A02,   MinMax = ( 0., 1. )            )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = AperpPh - A0Ph, MinMax = ( -2. * pi, 2. * pi ) )
        AperpAmp = Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )

        # A_par
        if AparParam != 'ReIm' :
            self._parseArg( 'AparMag2', kwargs, Title = '|A_par|^2', Value = Apar2 / A02, MinMax = ( 0., 1. ) )

        if AparParam in [ 'ReIm', 'real', 'cos' ] :
            # Re(A_par)
            if AparParam in [ 'ReIm', 'real' ] :
                self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Value = sqrt(Apar2 / A02) * cos(AparPh - A0Ph), MinMax = (-1., 1.) )
            else :
                self._parseArg( 'cosAparPhase', kwargs, Title = 'cos(delta_par)', Value = cos(AparPh - A0Ph), MinMax = ( -1., 1. ) )
                self._ReApar = FormulaVar( 'ReApar', 'sqrt(@0)*@1', [ self._AparMag2, self._cosAparPhase ], Title = 'Re(A_par)' )

            # Im(A_par)
            if AparParam != 'ReIm' :
                from RooFitWrappers import Category
                self._ImAparSign = Category( 'ImAparSign', States = { 'plus' : +1, 'minus' : -1 } )
                self._ImAparSign.setIndex(-1)
                self._ImApar = FormulaVar( 'ImApar', '@2*sqrt(@0 - @1*@1)', [ self._AparMag2, self._ReApar, self._ImAparSign ]
                                          , Title = 'Im(A_par)' )
            else :
                self._parseArg( 'ImApar', kwargs, Title = 'Im(A_par)', Value = sqrt(Apar2 / A02) * sin(AparPh - A0Ph), MinMax = (-1., 1.) )

            AparAmp = Carthesian_Amplitude( 'Apar', self._ReApar, self._ImApar, +1 )

        else :
            # delta_par
            self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = AparPh - A0Ph,  MinMax = ( -2. * pi, 2. * pi ) )
            AparAmp = Polar2_Amplitude( 'Apar', self._AparMag2, self._AparPhase, +1 )

        # A_S (--> A_S / A_perp and not A_S / A_0!!!)
        if polarSWave :
            self._parseArg( 'ASOddMag2',  kwargs, Title = '|A_S|^2 / |A_perp|^2', Value = AS2 / Aperp2,   MinMax = (  0.,      1.      ) )
            self._parseArg( 'ASOddPhase', kwargs, Title = 'delta_S - delta_perp', Value = ASPh - AperpPh, MinMax = ( -2. * pi, 2. * pi ) )
            self._ReAS = FormulaVar(  'ReAS'
                                    , 'sqrt(@0*@2) * cos(@1+@3)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ASOddMag2, self._ASOddPhase ]
                                    , Title = 'Re(A_S)'
                                   )
            self._ImAS = FormulaVar(  'ImAS'
                                    , 'sqrt(@0*@2) * sin(@1+@3)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ASOddMag2, self._ASOddPhase ]
                                    , Title = 'Im(A_S)'
                                   )

        else :
            self._parseArg( 'ReASOdd', kwargs, Title = 'Re(A_S / A_perp)', Value = sqrt(AS2 / Aperp2) * cos(ASPh - AperpPh)
                           , MinMax = ( -1., 1. ) )
            self._parseArg( 'ImASOdd', kwargs, Title = 'Im(A_S / A_perp)', Value = sqrt(AS2 / Aperp2) * sin(ASPh - AperpPh)
                           , MinMax = ( -1., 1. ) )
            self._ReAS = FormulaVar(  'ReAS'
                                    , 'sqrt(@0) * (cos(@1)*@2 - sin(@1)*@3)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ReASOdd, self._ImASOdd ]
                                    , Title = 'Re(A_S)'
                                   )
            self._ImAS = FormulaVar(  'ImAS'
                                    , 'sqrt(@0) * (cos(@1)*@3 + sin(@1)*@2)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ReASOdd, self._ImASOdd ]
                                    , Title = 'Im(A_S)'
                                   )

        ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1 )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, A0Amp, AparAmp, AperpAmp, ASAmp )
