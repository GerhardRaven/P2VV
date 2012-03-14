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


# construct amplitudes with carthesian parameters
class Carthesian_Amplitude :
    def __init__( self,name, Re, Im, CP ) :
        self.name = name
        self.Re = Re
        self.Im = Im
        self.CP = CP # even or odd???
    def __str__(self) : return self.name


# construct amplitudes with polar parameters
class Polar2_Amplitude( Carthesian_Amplitude ) :
    def __init__( self,name, r2, arg, CP ) :
        from RooFitWrappers import FormulaVar
        Carthesian_Amplitude.__init__( self,  name, FormulaVar('Re_%s'%name, 'sqrt(@0) * cos(@1)', [r2,arg], Title = 'Re(%s)'% name )
                                                  , FormulaVar('Im_%s'%name, 'sqrt(@0) * sin(@1)', [r2,arg], Title = 'Im(%s)'% name )
                                                  , CP )


class AmplitudeSet( dict, _util_parse_mixin ) :
    def __init__( self, *args ) :
        # maybe make this thing readonly???
        for v in args: 
            self[ v.name ] = v
            assert hasattr( v, 'Re' )
            assert hasattr( v, 'Im' )
            assert hasattr( v, 'CP' )
        # require the names in args to be unique...
        assert( len(self)==len(args) )

    def __getitem__( self, kw ) :
        if kw in self : return dict.__getitem__( self, kw )
        return getattr( self, '_' + kw )


class JpsiVCarthesian_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        from math import sqrt, cos, sin
        self._parseArg( 'ReA0',    kwargs, Title = 'Re(A_0)',    Value = 1.,                                  Constant = True      )
        self._parseArg( 'ImA0',    kwargs, Title = 'Im(A_0)',    Value = 0.,                                  Constant = True      )
        self._parseArg( 'ReApar',  kwargs, Title = 'Re(A_par)',  Value = sqrt( Apar2  / A02 ) * cos(AparPh ), MinMax = ( -1., 1. ) )
        self._parseArg( 'ImApar',  kwargs, Title = 'Im(A_par)',  Value = sqrt( Apar2  / A02 ) * sin(AparPh ), MinMax = ( -1., 1. ) )
        self._parseArg( 'ReAperp', kwargs, Title = 'Re(A_perp)', Value = sqrt( Aperp2 / A02 ) * cos(AperpPh), MinMax = ( -1., 1. ) )
        self._parseArg( 'ImAperp', kwargs, Title = 'Im(A_perp)', Value = sqrt( Aperp2 / A02 ) * sin(AperpPh), MinMax = ( -1., 1. ) )
        self._parseArg( 'ReAS',    kwargs, Title = 'Re(A_S)',    Value = sqrt( AS2    / A02 ) * cos(ASPh   ), MinMax = ( -1., 1. ) )
        self._parseArg( 'ImAS',    kwargs, Title = 'Im(A_S)',    Value = sqrt( AS2    / A02 ) * sin(ASPh   ), MinMax = ( -1., 1. ) )

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

        self._parseArg( 'A0Phase',    kwargs, Title = 'delta_0',    Value = 0.,      Constant = True                )
        self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = AparPh,  MinMax = ( -2. * pi, 2. * pi ) )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = AperpPh, MinMax = ( -2. * pi, 2. * pi ) )
        self._parseArg( 'ASPhase',    kwargs, Title = 'delta_S',    Value = ASPh,    MinMax = ( -2. * pi, 2. * pi ) )

        self._check_extraneous_kw( kwargs ) 
        AmplitudeSet.__init__( self, Polar2_Amplitude( 'A0',    self._A0Mag2,    self._A0Phase,    +1 )
                                   , Polar2_Amplitude( 'Apar',  self._AparMag2,  self._AparPhase,  +1 )
                                   , Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )
                                   , Polar2_Amplitude( 'AS',    self._ASMag2,    self._ASPhase,    -1 )
                             )


class JpsiVPolarSWaveFrac_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        from math import pi
        from RooFitWrappers import FormulaVar
        polarSWave = kwargs.pop( 'polarSWave', True )

        self._parseArg( 'A0Mag2',    kwargs, Title = '|A0|^2',     Value = A02,    MinMax = ( 0., 1. ) )
        self._parseArg( 'AperpMag2', kwargs, Title = '|A_perp|^2', Value = Aperp2, MinMax = ( 0., 1. ) )

        self._AparMag2 = FormulaVar( 'AparMag2', '1. - @0 - @1',   [ self._A0Mag2, self._AperpMag2 ], Title = '|A_par|^2' )

        self._parseArg( 'A0Phase',    kwargs, Title = 'delta_0',    Value = 0.,      Constant = True                )
        self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = AparPh,  MinMax = ( -2. * pi, 2. * pi ) )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = AperpPh, MinMax = ( -2. * pi, 2. * pi ) )
 
        if polarSWave :
            self._parseArg( 'f_S',     kwargs, Title = 'S wave fraction', Value = f_S,  MinMax = (  0.,      1.      ) )
            self._parseArg( 'ASPhase', kwargs, Title = 'delta_S',         Value = ASPh, MinMax = ( -2. * pi, 2. * pi ) )
            self._ASMag2 = FormulaVar( 'ASMag2', '@0 / (1. - @0)', [ self._f_S ], Title = '|A_S|^2' )

            self._check_extraneous_kw( kwargs ) 
            AmplitudeSet.__init__( self, Polar2_Amplitude( 'A0',    self._A0Mag2,    self._A0Phase,    +1 )
                                       , Polar2_Amplitude( 'Apar',  self._AparMag2,  self._AparPhase,  +1 )
                                       , Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )
                                       , Polar2_Amplitude( 'AS',    self._ASMag2,    self._ASPhase,    -1 )
                                 )

        else :
            from math import sin, cos, sqrt
            self._parseArg( 'sqrtfS_Re', kwargs, Title = 'sqrt(S wave fraction) * cos(delta_S)', Value = sqrt(f_S) * cos(ASPh)
                           , MinMax = ( -1., 1. ) )
            self._parseArg( 'sqrtfS_Im', kwargs, Title = 'sqrt(S wave fraction) * sin(delta_S)', Value = sqrt(f_S) * sin(ASPh)
                           , MinMax = ( -1., 1. ) )
            self._ReAS = FormulaVar( 'ReAS', '@0 / sqrt(1. - @0*@0 - @1*@1)', [ self._sqrtfS_Re, self._sqrtfS_Im ], Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@1 / sqrt(1. - @0*@0 - @1*@1)', [ self._sqrtfS_Re, self._sqrtfS_Im ], Title = 'Im(A_S)' )

            self._check_extraneous_kw( kwargs ) 
            AmplitudeSet.__init__( self, Polar2_Amplitude(     'A0',    self._A0Mag2,    self._A0Phase,    +1 )
                                       , Polar2_Amplitude(     'Apar',  self._AparMag2,  self._AparPhase,  +1 )
                                       , Polar2_Amplitude(     'Aperp', self._AperpMag2, self._AperpPhase, -1 )
                                       , Carthesian_Amplitude( 'AS',    self._ReAS,      self._ImAS,       -1 )
                                 )
