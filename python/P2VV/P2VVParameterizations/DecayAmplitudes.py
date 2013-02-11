###########################################################################################################################################
## P2VVParameterizations.DecayAmplitudes: Decay amplitude parameterizations                                                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin, _util_conditionalObs_mixin

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

    def _createBinnedAmp( self, name, title, initVals, minMaxs, observable, binning ) :
        # get initial values, minima and maxima for bin coefficients
        if hasattr( initVals, '__iter__' ) and len(initVals) == binning.numBins() : vals = initVals
        else : vals = binning.numBins() * [ initVals ]
        if hasattr( minMaxs, '__iter__' ) and len(minMaxs) == binning.numBins() and type(minMaxs) != tuple : minMax = minMaxs
        else : minMax = binning.numBins() * [ minMaxs ]

        # create bin coefficients
        from RooFitWrappers import RealVar
        coefs = [ RealVar( '%s_%d' % ( name, bin ), Title = '%s %d' % ( title, bin ), Value = vals[bin], MinMax = minMax[bin] )\
                  for bin in range( binning.numBins() ) ]
        setattr( self, '_%sCoefs' % name, coefs )

        # abuse BinnedPdf to make a binned function for the amplitude
        from RooFitWrappers import BinnedPdf
        amp = BinnedPdf(  Name = name
                        , Observable = observable
                        , Binning = binning
                        , Coefficients = coefs
                        , BinIntegralCoefs = False
                       )
        setattr( self, '_' + name, amp )

        # put amplitude in parameters
        if not hasattr( self, '_params' ) : self._params = []
        self._params += [ amp ]


class JpsiVCarthesian_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters', False )

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

        if 'KKMass' in kwargs and 'KKMassBinning' in kwargs :
            self._KKMass = kwargs.pop('KKMass')
            self._KKMassBinning = kwargs.pop('KKMassBinning')
            self._createBinnedAmp( 'ReAS', 'Re(A_S)', ReAS, ( -1., 1. ), self._KKMass, self._KKMassBinning )
            self._createBinnedAmp( 'ImAS', 'Im(A_S)', ImAS, ( -1., 1. ), self._KKMass, self._KKMassBinning )
        else :
            self._KKMass = None
            self._KKMassBinning = None
            self._parseArg( 'ReAS', kwargs, Title = 'Re(A_S)', Value = ReAS, MinMax = ( -1., 1. ) )
            self._parseArg( 'ImAS', kwargs, Title = 'Im(A_S)', Value = ImAS, MinMax = ( -1., 1. ) )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, Carthesian_Amplitude( 'A0',    self._ReA0,    self._ImA0,    +1 )
                                   , Carthesian_Amplitude( 'Apar',  self._ReApar,  self._ImApar,  +1 )
                                   , Carthesian_Amplitude( 'Aperp', self._ReAperp, self._ImAperp, -1 )
                                   , Carthesian_Amplitude( 'AS',    self._ReAS,    self._ImAS,    -1 )
                                   , Conditionals = [ self._KKMass ] if self._KKMass else [ ]
                             )


class JpsiVPolar_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters', False       )
        ASParam       = kwargs.pop( 'ASParameterization',  'deltaPerp' )

        from math import pi
        deltaPar  = AparPh  - A0Ph
        deltaPerp = AperpPh - A0Ph
        deltaS    = ASPh    - ( AperpPh if ASParam in [ 'deltaPerp', 'Mag2ReImPerp' ] else A0Ph )
        if ambiguityPars :
            deltaPar  = -deltaPar
            deltaPerp = pi - deltaPerp
            deltaS    = ( pi - deltaS ) if ASParam in [ 'deltaPerp', 'Mag2ReImPerp' ] else -deltaS

        # A_0, A_par, A_perp
        from RooFitWrappers import FormulaVar
        self._parseArg( 'A0Mag2',    kwargs, Title = '|A0|^2',     Value = A02,    MinMax = ( 0., 1. ) )
        self._parseArg( 'AperpMag2', kwargs, Title = '|A_perp|^2', Value = Aperp2, MinMax = ( 0., 1. ) )

        PNorm = kwargs.pop( 'PWaveNorm', True )
        if PNorm :
            self._AparMag2 = FormulaVar( 'AparMag2', '1.-@0-@1',    [ self._A0Mag2, self._AperpMag2 ],               Title = '|A_par|^2' )
        else :
            self._AparMag2 = FormulaVar( 'AparMag2', '1.-@0-@1-@2', [ self._A0Mag2, self._AperpMag2, self._ASMag2 ], Title = '|A_par|^2' )

        self._parseArg( 'A0Phase',    kwargs, Title = 'delta_0',    Value = 0.,        Constant = True                )
        self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = deltaPar,  MinMax = ( -2. * pi, 2. * pi ) )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = deltaPerp, MinMax = ( -2. * pi, 2. * pi ) )

        # A_S
        if not ASParam.startswith('Mag2ReIm') :
            self._parseArg( 'C_SP', kwargs, Title = 'S-P wave couping factor', Value = C_SP, MinMax = ( 0., 2. ) )
        if 'KKMass' in kwargs and 'KKMassBinning' in kwargs :
            self._KKMass = kwargs.pop('KKMass')
            self._KKMassBinning = kwargs.pop('KKMassBinning')
        else :
            self._KKMass = None
            self._KKMassBinning = None

        if self._KKMass and self._KKMassBinning :
            self._createBinnedAmp( 'ASMag2', '|A_S|^2', AS2, ( 0.,  1. ), self._KKMass, self._KKMassBinning )
            if ASParam == 'deltaPerp' :
                self._createBinnedAmp( 'ASOddPhase', 'delta_S - delta_perp', deltaS, ( -2. * pi, 2. * pi )
                                      , self._KKMass, self._KKMassBinning )
            elif ASParam == 'Mag2ReIm' :
                from math import sqrt, cos, sin
                self._createBinnedAmp( 'ReAS', 'Re(A_S)', sqrt(AS2) * cos(deltaS), ( -1., 1. ), self._KKMass, self._KKMassBinning )
                self._createBinnedAmp( 'ImAS', 'Im(A_S)', sqrt(AS2) * sin(deltaS), ( -1., 1. ), self._KKMass, self._KKMassBinning )
            elif ASParam == 'Mag2ReImPerp' :
                from math import sqrt, cos, sin
                self._createBinnedAmp( 'ReASPerp', 'Re(A_S*e^(-i*delta_perp))', sqrt(AS2) * cos(deltaS), ( -1., 1. )
                                      , self._KKMass, self._KKMassBinning )
                self._createBinnedAmp( 'ImASPerp', 'Im(A_S*e^(-i*delta_perp))', sqrt(AS2) * sin(deltaS), ( -1., 1. )
                                      , self._KKMass, self._KKMassBinning )
            else  :
                self._createBinnedAmp( 'ASPhase', 'delta_S', deltaS, ( -2. * pi, 2. * pi ), self._KKMass, self._KKMassBinning )

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


        if ASParam == 'deltaPerp' :
            self._ReAS = FormulaVar( 'ReAS', '@3 * sqrt(@0) * cos(@1+@2)'
                                    , [ self._ASMag2, self._ASOddPhase, self._AperpPhase, self._C_SP ], Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@3 * sqrt(@0) * sin(@1+@2)'
                                    , [ self._ASMag2, self._ASOddPhase, self._AperpPhase, self._C_SP ], Title = 'Im(A_S)' )

        elif ASParam == 'Mag2ReImPerp' :
            self._ReAS = FormulaVar( 'ReAS', '@0 * cos(@2) - @1 * sin(@2)', [ self._ReASPerp, self._ImASPerp, self._AperpPhase ]
                                    , Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@0 * sin(@2) + @1 * cos(@2)', [ self._ImASPerp, self._ImASPerp, self._AperpPhase ]
                                    , Title = 'Im(A_S)' )

        elif ASParam != 'Mag2ReIm' :
            self._ReAS = FormulaVar( 'ReAS', '@2 * sqrt(@0) * cos(@1)', [ self._ASMag2, self._ASPhase, self._C_SP ]
                                    , Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@2 * sqrt(@0) * sin(@1)', [ self._ASMag2, self._ASPhase, self._C_SP ]
                                        , Title = 'Im(A_S)' )

        ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1
                                     , getattr( self, '_ASMag2', None ) 
                                     , getattr( self, '_ASPhase', None )  )

        self._check_extraneous_kw( kwargs ) 
        AmplitudeSet.__init__( self, Polar2_Amplitude( 'A0',    self._A0Mag2,    self._A0Phase,    +1 )
                                   , Polar2_Amplitude( 'Apar',  self._AparMag2,  self._AparPhase,  +1 )
                                   , Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )
                                   , ASAmp
                                   , Conditionals = [ self._KKMass ] if self._KKMass else [ ]
                             )


class JpsiVPolarSWaveFrac_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters',  False       )
        ASParam       = kwargs.pop( 'ASParameterization',   'deltaPerp' )
        AparParam     = kwargs.pop( 'AparParameterization', 'phase'     )

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
        from RooFitWrappers import FormulaVar

        # A_0
        self._parseArg( 'A0Mag2',  kwargs, Title = '|A0|^2',  Value = A02, Error = A02Err, MinMax = ( 0., 1. ) )
        self._parseArg( 'A0Phase', kwargs, Title = 'delta_0', Value = 0.,  Constant = True )
        A0Amp = Polar2_Amplitude( 'A0', self._A0Mag2, self._A0Phase, +1 )

        # A_perp
        self._parseArg( 'AperpMag2',  kwargs, Title = '|A_perp|^2', Value = Aperp2,    Error = Aperp2Err,  MinMax = (  0.,     1.     ) )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = deltaPerp, Error = delPerpErr, MinMax = ( -RooInf, RooInf ) )
        AperpAmp = Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )

        # A_par
        self._AparMag2 = FormulaVar( 'AparMag2', '1. - @0 - @1', [ self._A0Mag2, self._AperpMag2 ], Title = '|A_par|^2' )
        if AparParam in [ 'real', 'cos' ] :
            # Re(A_par)
            if AparParam == 'real' :
                self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Value = sqrt(Apar2) * cos(deltaPar), MinMax = ( -1., 1. ) )
            else :
                self._parseArg( 'cosAparPhase', kwargs, Title = 'cos(delta_par)', Value = cos(deltaPar), MinMax = ( -1., 1. ) )
                self._ReApar = FormulaVar( 'ReApar', 'sqrt(@0)*@1', [ self._AparMag2, self._cosAparPhase ], Title = 'Re(A_par)' )

            # Im(A_par)
            from RooFitWrappers import Category
            self._ImAparSign = Category( 'ImAparSign', States = { 'plus' : +1, 'minus' : -1 } )
            if ambiguityPars : self._ImAparSign.setIndex(+1)
            else             : self._ImAparSign.setIndex(-1)
            self._ImApar = FormulaVar( 'ImApar', '@2*sqrt(@0 - @1*@1)', [ self._AparMag2, self._ReApar, self._ImAparSign ]
                                      , Title = 'Im(A_par)' )
            AparAmp = Carthesian_Amplitude( 'Apar', self._ReApar, self._ImApar, +1 )

        else :
            # delta_par
            self._parseArg( 'AparPhase', kwargs, Title = 'delta_par', Value = deltaPar, Error = delParErr, MinMax = ( -RooInf, RooInf ) )
            AparAmp = Polar2_Amplitude( 'Apar', self._AparMag2, self._AparPhase, +1 )

        # A_S
        self._parseArg( 'C_SP', kwargs, Title = 'S-P wave couping factor', Value = C_SP, Error = C_SPErr, MinMax = ( 0., 1. ) )
        if 'KKMass' in kwargs and 'KKMassBinning' in kwargs :
            self._KKMass = kwargs.pop('KKMass')
            self._KKMassBinning = kwargs.pop('KKMassBinning')
        else :
            self._KKMass = None
            self._KKMassBinning = None

        if ASParam != 'ReIm' :
            if self._KKMass and self._KKMassBinning :
                self._createBinnedAmp( 'f_S', 'S wave fraction', f_S, ( 0., 1. ), self._KKMass, self._KKMassBinning )
                if ASParam == 'deltaPerp' :
                    self._createBinnedAmp( 'ASOddPhase', 'delta_S - delta_perp', deltaS, ( -2. * pi, 2. * pi )
                                          , self._KKMass, self._KKMassBinning )
                else :
                    self._createBinnedAmp( 'ASPhase', 'delta_S', deltaS, ( -2. * pi, 2. * pi ), self._KKMass, self._KKMassBinning )

            else :
                self._parseArg( 'f_S', kwargs, Title = 'S wave fraction', Value = f_S, Error = f_SErr, MinMax = ( 0., 1. ) )
                if ASParam == 'deltaPerp' :
                    self._parseArg( 'ASOddPhase', kwargs, Title = 'delta_S - delta_perp', Value = deltaS, Error = delSErr
                                   , MinMax = ( -RooInf, RooInf ) )
                else :
                    self._parseArg( 'ASPhase', kwargs, Title = 'delta_S', Value = deltaS, Error = delSErr, MinMax = ( -RooInf, RooInf ) )

            self._ASMag2 = FormulaVar( 'ASMag2', '@0 / (1. - @0)', [ self._f_S ], Title = 'Re(A_S)' )
            if ASParam == 'deltaPerp' :
                self._ReAS = FormulaVar( 'ReAS', '@3 * sqrt(@0 / (1. - @0)) * cos(@1+@2)'
                                        , [ self._f_S, self._ASOddPhase, self._AperpPhase, self._C_SP ], Title = 'Re(A_S)' )
                self._ImAS = FormulaVar( 'ImAS', '@3 * sqrt(@0 / (1. - @0)) * sin(@1+@2)'
                                        , [ self._f_S, self._ASOddPhase, self._AperpPhase, self._C_SP ], Title = 'Im(A_S)' )

            else :
                self._ReAS = FormulaVar( 'ReAS', '@2 * sqrt(@0 / (1. - @0)) * cos(@1)', [ self._f_S, self._ASPhase, self._C_SP ]
                                        , Title = 'Re(A_S)' )
                self._ImAS = FormulaVar( 'ImAS', '@2 * sqrt(@0 / (1. - @0)) * sin(@1)', [ self._f_S, self._ASPhase, self._C_SP ]
                                        , Title = 'Im(A_S)' )

        else :
            if self._KKMass and self._KKMassBinning :
                self._createBinnedAmp( 'sqrtfS_Re', 'sqrt(S wave fraction) * cos(delta_S)', sqrt(f_S) * cos(deltaS), ( -1., 1. )
                                      , self._KKMass, self._KKMassBinning )
                self._createBinnedAmp( 'sqrtfS_Im', 'sqrt(S wave fraction) * sin(delta_S)', sqrt(f_S) * sin(deltaS), ( -1., 1. )
                                      , self._KKMass, self._KKMassBinning )
            else :
                self._parseArg( 'sqrtfS_Re', kwargs, Title = 'sqrt(S wave fraction) * cos(delta_S)', Value = sqrt(f_S) * cos(deltaS)
                               , MinMax = ( -1., 1. ) )
                self._parseArg( 'sqrtfS_Im', kwargs, Title = 'sqrt(S wave fraction) * sin(delta_S)', Value = sqrt(f_S) * sin(deltaS)
                               , MinMax = ( -1., 1. ) )

            self._ASMag2 = FormulaVar( 'ASMag2', '@0 / (1. - @0*@0 - @1*@1)', [ self._f_S ], Title = 'Re(A_S)' )
            self._ReAS = FormulaVar( 'ReAS', '@2 * @0 / sqrt(1. - @0*@0 - @1*@1)', [ self._sqrtfS_Re, self._sqrtfS_Im, self._C_SP ]
                                    , Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@2 * @1 / sqrt(1. - @0*@0 - @1*@1)', [ self._sqrtfS_Re, self._sqrtfS_Im, self._C_SP ]
                                    , Title = 'Im(A_S)' )

        ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1
                                     , getattr( self, '_ASMag2', None ) 
                                     , getattr( self, '_ASPhase', None ) )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, A0Amp, AparAmp, AperpAmp, ASAmp, Conditionals = [ self._KKMass ] if self._KKMass else [ ] )


class JpsiVBank_AmplitudeSet( AmplitudeSet ) :
    def __init__( self, **kwargs ) :
        ambiguityPars = kwargs.pop( 'AmbiguityParameters',  False       )
        ASParam       = kwargs.pop( 'ASParameterization',   'deltaPerp' )
        AparParam     = kwargs.pop( 'AparParameterization', 'phase'     )

        from math import pi
        deltaPar  = AparPh  - A0Ph
        deltaPerp = AperpPh - A0Ph
        deltaS    = ASPh    - AperpPh
        if ambiguityPars :
            deltaPar  = -deltaPar
            deltaPerp = pi - deltaPerp
            deltaS    = pi - deltaS

        from math import sqrt, sin, cos
        from RooFitWrappers import FormulaVar

        # A_0
        self._parseArg( 'ReA0', kwargs, Title = 'Re(A_0)', Value = 1., Constant = True )
        self._parseArg( 'ImA0', kwargs, Title = 'Im(A_0)', Value = 0., Constant = True )
        A0Amp = Carthesian_Amplitude( 'A0', self._ReA0, self._ImA0, +1 )

        # A_perp
        self._parseArg( 'AperpMag2',  kwargs, Title = '|A_perp|^2', Value = Aperp2 / A02, MinMax = ( 0., 1. )            )
        self._parseArg( 'AperpPhase', kwargs, Title = 'delta_perp', Value = deltaPerp,    MinMax = ( -2. * pi, 2. * pi ) )
        AperpAmp = Polar2_Amplitude( 'Aperp', self._AperpMag2, self._AperpPhase, -1 )

        # A_par
        if AparParam != 'ReIm' :
            self._parseArg( 'AparMag2', kwargs, Title = '|A_par|^2', Value = Apar2 / A02, MinMax = ( 0., 1. ) )

        if AparParam in [ 'Mag2ReIm', 'ReIm', 'real', 'cos' ] :
            # Re(A_par)
            if AparParam in [ 'Mag2ReIm', 'ReIm', 'real' ] :
                self._parseArg( 'ReApar', kwargs, Title = 'Re(A_par)', Value = sqrt(Apar2 / A02) * cos(deltaPar), MinMax = (-1., 1.) )
            else :
                self._parseArg( 'cosAparPhase', kwargs, Title = 'cos(delta_par)', Value = cos(deltaPar), MinMax = ( -1., 1. ) )
                self._ReApar = FormulaVar( 'ReApar', 'sqrt(@0)*@1', [ self._AparMag2, self._cosAparPhase ], Title = 'Re(A_par)' )

            # Im(A_par)
            if AparParam not in [ 'Mag2ReIm', 'ReIm' ] :
                from RooFitWrappers import Category
                self._ImAparSign = Category( 'ImAparSign', States = { 'plus' : +1, 'minus' : -1 } )
                if ambiguityPars : self._ImAparSign.setIndex(+1)
                else             : self._ImAparSign.setIndex(-1)
                self._ImApar = FormulaVar( 'ImApar', '@2*sqrt(@0 - @1*@1)', [ self._AparMag2, self._ReApar, self._ImAparSign ]
                                          , Title = 'Im(A_par)' )
            else :
                self._parseArg( 'ImApar', kwargs, Title = 'Im(A_par)', Value = sqrt(Apar2 / A02) * sin(deltaPar), MinMax = (-1., 1.) )

            if AparParam == 'ReIm' :
                self._AparMag2 = FormulaVar( 'AparMag2', '@0*@0+@1*@1', [ self._ReApar, self._ImApar ], Title = '|A_par|^2' )

            AparAmp = Carthesian_Amplitude( 'Apar', self._ReApar, self._ImApar, +1, self._AparMag2 )

        else :
            # delta_par
            self._parseArg( 'AparPhase',  kwargs, Title = 'delta_par',  Value = deltaPar,  MinMax = ( -2. * pi, 2. * pi ) )
            AparAmp = Polar2_Amplitude( 'Apar', self._AparMag2, self._AparPhase, +1 )

        # A_S (--> A_S / A_perp and not A_S / A_0!!!)
        self._parseArg( 'C_SP', kwargs, Title = 'S-P wave couping factor', Value = C_SP, MinMax = ( 0., 2. ) )
        if 'KKMass' in kwargs and 'KKMassBinning' in kwargs :
            self._KKMass = kwargs.pop('KKMass')
            self._KKMassBinning = kwargs.pop('KKMassBinning')
        else :
            self._KKMass = None
            self._KKMassBinning = None

        if ASParam != 'ReIm' :
            if self._KKMass and self._KKMassBinning :
                self._createBinnedAmp( 'ASOddMag2',  '|A_S|^2 / |A_perp|^2', AS2 / Aperp2, (  0.,      1.      )
                                      , self._KKMass, self._KKMassBinning )
                self._createBinnedAmp( 'ASOddPhase', 'delta_S - delta_perp', deltaS,       ( -2. * pi, 2. * pi )
                                      , self._KKMass, self._KKMassBinning )
            else :
                self._parseArg( 'ASOddMag2',  kwargs, Title = '|A_S|^2 / |A_perp|^2', Value = AS2 / Aperp2
                               , MinMax = (  0.,      1.      ) )
                self._parseArg( 'ASOddPhase', kwargs, Title = 'delta_S - delta_perp', Value = deltaS
                               , MinMax = ( -2. * pi, 2. * pi ) )

            self._ReAS = FormulaVar( 'ReAS', '@4 * sqrt(@0*@2) * cos(@1+@3)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ASOddMag2, self._ASOddPhase, self._C_SP ]
                                    , Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@4 * sqrt(@0*@2) * sin(@1+@3)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ASOddMag2, self._ASOddPhase, self._C_SP ]
                                    , Title = 'Im(A_S)' )

        else :
            if self._KKMass and self._KKMassBinning :
                self._createBinnedAmp( 'ReASOdd', 'Re(A_S / A_perp)', sqrt(AS2 / Aperp2) * cos(deltaS), ( -1., 1. )
                                     , self._KKMass, self._KKMassBinning )
                self._createBinnedAmp( 'ImASOdd', 'Im(A_S / A_perp)', sqrt(AS2 / Aperp2) * sin(deltaS), ( -1., 1. )
                                     , self._KKMass, self._KKMassBinning )
            else :
                self._parseArg( 'ReASOdd', kwargs, Title = 'Re(A_S / A_perp)', Value = sqrt(AS2 / Aperp2) * cos(deltaS)
                               , MinMax = ( -1., 1. ) )
                self._parseArg( 'ImASOdd', kwargs, Title = 'Im(A_S / A_perp)', Value = sqrt(AS2 / Aperp2) * sin(deltaS)
                               , MinMax = ( -1., 1. ) )

            self._ReAS = FormulaVar( 'ReAS', '@4 * sqrt(@0) * (cos(@1)*@2 - sin(@1)*@3)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ReASOdd, self._ImASOdd, self._C_SP ], Title = 'Re(A_S)' )
            self._ImAS = FormulaVar( 'ImAS', '@4 * sqrt(@0) * (cos(@1)*@3 + sin(@1)*@2)'
                                    , [ self._AperpMag2, self._AperpPhase, self._ReASOdd, self._ImASOdd, self._C_SP ], Title = 'Im(A_S)' )

        ASAmp = Carthesian_Amplitude( 'AS', self._ReAS, self._ImAS, -1 )

        self._check_extraneous_kw( kwargs )
        AmplitudeSet.__init__( self, A0Amp, AparAmp, AperpAmp, ASAmp, Conditionals = [ self._KKMass ] if self._KKMass else [ ] )
