###########################################################################################################################################
## P2VVParameterizations.AngularFunctions: Decay angles and functions of decay angles                                                    ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

#TODO: inherit from UserDict mixin instead of wrapping & forwarding...
class AngularFunctions :
    def __init__(self) :         self._d = dict()
    def __getitem__(self,k) :    return self._d[k]
    def __setitem__(self,k,v) :  self._d[k] = v
    def __contains__(self,k) :   return k in self._d
    def keys(self) :             return self._d.keys()
    def iterkeys(self) :         return self._d.iterkeys()
    def items(self) :            return self._d.items()
    def iteritems(self) :        return self._d.iteritems()
    def values(self) :           return self._d.values()
    def itervalues(self) :       return self._d.itervalues()


from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin
class AngleDefinitions( _util_parse_mixin ) :
    def __init__( self , **kwargs ) :
        self.angles = dict( (i,kwargs.pop(i)) for i in ['cpsi' ,'ctheta','phi' ] )
        self.functions = kwargs.pop('functions')


class JpsiphiHelicityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = {  'cpsi'   : self._parseArg( 'cpsi',   kwargs, SingleArgKey = 'Name', Name = 'helcosthetaK'
                                         , Title = 'cos(#theta_{K})',   MinMax = ( -1,  1 ),  Observable = True )
             , 'ctheta' : self._parseArg( 'ctheta', kwargs, SingleArgKey = 'Name', Name = 'helcosthetaL'
                                         , Title = 'cos(#theta_{#mu})', MinMax = ( -1,  1 ),  Observable = True )
             , 'phi'    : self._parseArg( 'phi',    kwargs, SingleArgKey = 'Name', Name = 'helphi'
                                         , Title = '#phi_{h}',          MinMax = ( -pi, pi ), Observable = True, Unit = 'rad' )
            }
        coef = kwargs.pop( 'DummyCoef', False )
        self._check_extraneous_kw(kwargs)
        d['functions'] = JpsiphiTransversityAmplitudesHelicityAngles( Angles = d, DummyCoef = coef )
        AngleDefinitions.__init__(self, **d )


class JpsiphiTransversityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = { 'cpsi' :   self._parseArg( 'cpsi',   kwargs, SingleArgKey = 'Name', Name='trcospsi'
                                        , Title = 'cos(#psi_{tr})',   MinMax=( -1,  1 ),  Observable = True )
            , 'ctheta' : self._parseArg( 'ctheta', kwargs, SingleArgKey = 'Name', Name='trcostheta'
                                        , Title = 'cos(#theta_{tr})', MinMax=( -1,  1 ),  Observable = True )
            , 'phi'   :  self._parseArg( 'phi',    kwargs, SingleArgKey = 'Name', Name='trphi'
                                        , Title = '#phi_{tr}',        MinMax=( -pi, pi ), Observable = True )
            }
        coef = kwargs.pop( 'DummyCoef', False )
        d['functions'] = JpsiphiTransversityAmplitudesTransversityAngles( Angles = d, DummyCoef = coef )
        self._check_extraneous_kw(kwargs)
        AngleDefinitions.__init__(self, **d )


class JpsiKstHelicityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = {  'cpsi'   : self._parseArg( 'cpsi',   kwargs, SingleArgKey = 'Name', Name = 'helcosthetaK'
                                         , Title = 'cos(#theta_{K})',   MinMax = ( -1,  1 ),  Observable = True )
             , 'ctheta' : self._parseArg( 'ctheta', kwargs, SingleArgKey = 'Name', Name = 'helcosthetaL'
                                         , Title = 'cos(#theta_{#mu})', MinMax = ( -1,  1 ),  Observable = True )
             , 'phi'    : self._parseArg( 'phi',    kwargs, SingleArgKey = 'Name', Name = 'helphi'
                                         , Title = '#phi_{h}',          MinMax = ( -pi, pi ), Observable = True, Unit = 'rad' )
            }
        coef = kwargs.pop( 'DummyCoef', False )
        self._check_extraneous_kw(kwargs)
        d['functions'] = JpsiKstTransversityAmplitudesHelicityAngles( Angles = d, DummyCoef = coef )
        AngleDefinitions.__init__(self, **d )


class JpsiKstTransversityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = { 'cpsi' :   self._parseArg( 'cpsi',   kwargs, SingleArgKey = 'Name', Name='trcospsi'
                                        , Title = 'cos(#psi_{tr})',   MinMax=( -1,  1 ),  Observable = True )
            , 'ctheta' : self._parseArg( 'ctheta', kwargs, SingleArgKey = 'Name', Name='trcostheta'
                                        , Title = 'cos(#theta_{tr})', MinMax=( -1,  1 ),  Observable = True )
            , 'phi'   :  self._parseArg( 'phi',    kwargs, SingleArgKey = 'Name', Name='trphi'
                                        , Title = '#phi_{tr}',        MinMax=( -pi, pi ), Observable = True )
            }
        coef = kwargs.pop( 'DummyCoef', False )
        d['functions'] = JpsiKstTransversityAmplitudesTransversityAngles( Angles = d, DummyCoef = coef )
        self._check_extraneous_kw(kwargs)
        AngleDefinitions.__init__(self, **d )


class JpsiphiTransversityAmplitudesHelicityAngles( AngularFunctions ) :
    def __init__( self, **kwargs ) :
        AngularFunctions.__init__(self)
        angles = kwargs.pop('Angles')
        coef   = kwargs.pop( 'DummyCoef', False )

        from P2VV.Parameterizations.GeneralUtils import getParNamePrefix
        namePF = getParNamePrefix(True)

        # define function to build angular function; create dummy coefficient if we need to split the angular efficiency functions later
        from P2VV.RooFitWrappers import P2VVAngleBasis, Addition, RealVar
        dummy = RealVar( Name = namePF + 'angEffDummyCoef', Value = 1. ) if coef else None
        _ba = lambda name, args :\
            Addition( namePF + name, [ P2VVAngleBasis( Name = namePF + name, Angles = angles, Indices = arg[0], FixedCoef = arg[1]
                                                      , Coefficient = dummy ) for arg in args ] )

        # TODO: generate the following table straight from the physics using PS->(VV,VS) ->ffss  (V=spin 1, f=spin 1/2, PS,S,s = spin 0)
        from math import sqrt
        angFuncs = { ('A0',   'A0')    :  ( _ba('Re_ang_A0_A0',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),       -sqrt( 16. / 5.))
                                                                          , ( ( 2, 0, 0,  0 ),        8.             )
                                                                          , ( ( 2, 0, 2,  0 ),       -sqrt( 64. / 5.))]), None)
                   , ('Apar', 'Apar')  :  ( _ba('Re_ang_Apar_Apar',       [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),        sqrt(  1. / 5.))
                                                                          , ( ( 2, 2, 2,  2 ),       -sqrt(  3. / 5.))]), None)
                   , ('Aperp','Aperp') :  ( _ba('Re_ang_Aperp_Aperp',     [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),        sqrt(  1. / 5.))
                                                                          , ( ( 2, 2, 2,  2 ),        sqrt(  3. / 5.))]), None)
                   , ('A0',   'Apar')  :  ( _ba('Re_ang_A0_Apar',         [ ( ( 2, 1, 2,  1 ),        sqrt( 24. / 5.))]), None)
                   , ('A0',   'Aperp') :  ( None, _ba('Im_ang_A0_Aperp',  [ ( ( 2, 1, 2, -1 ), -1. * -sqrt( 24. / 5.))]))
                   , ('Apar', 'Aperp') :  ( None, _ba('Im_ang_Apar_Aperp',[ ( ( 2, 2, 2, -2 ), -1. *  sqrt( 12. / 5.))]))
                   , ('AS',   'AS')    :  ( _ba('Re_ang_AS_AS',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),       -sqrt( 16. / 5.))]), None)
                   , ('A0',   'AS')    :  ( _ba('Re_ang_A0_AS',           [ ( ( 1, 0, 0,  0 ),        sqrt(192.     ))
                                                                          , ( ( 1, 0, 2,  0 ),       -sqrt(192. / 5.))]), None)
                   , ('Apar', 'AS')    :  ( _ba('Re_ang_Apar_AS',         [ ( ( 1, 1, 2,  1 ),        sqrt( 72. / 5.))]), None)
                   , ('Aperp','AS')    :  ( None, _ba('Im_ang_Aperp_AS',  [ ( ( 1, 1, 2, -1 ), -1. *  sqrt( 72. / 5.))]))
                   }
        for k,v in angFuncs.iteritems() : self[k] = v


class JpsiphiTransversityAmplitudesTransversityAngles( AngularFunctions ) :
    def __init__( self, **kwargs ) :
        AngularFunctions.__init__(self)
        angles = kwargs.pop('Angles')
        coef   = kwargs.pop( 'DummyCoef', False )

        from P2VV.Parameterizations.GeneralUtils import getParNamePrefix
        namePF = getParNamePrefix(True)

        # define function to build angular function; create dummy coefficient if we need to split the angular efficiency functions later
        from P2VV.RooFitWrappers import P2VVAngleBasis, Addition, RealVar
        dummy = RealVar( Name = namePF + 'angEffDummyCoef', Value = 1. ) if coef else None
        _ba = lambda name, args :\
            Addition( namePF + name, [ P2VVAngleBasis( Name = namePF + name, Angles = angles, Indices = arg[0], FixedCoef = arg[1]
                                                      , Coefficient = dummy ) for arg in args ] )

        # TODO: generate the following table straight from the physics using PS->(VV,VS) ->ffss  (V=spin 1, f=spin 1/2, PS,S,s = spin 0)
        from math import sqrt
        angFuncs = { ('A0',   'A0')    :  ( _ba('Re_ang_A0_A0',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),        sqrt(  4. / 5.))
                                                                          , ( ( 0, 0, 2,  2 ),       -sqrt( 12. / 5.))
                                                                          , ( ( 2, 0, 0,  0 ),        8.             )
                                                                          , ( ( 2, 0, 2,  0 ),        sqrt( 16. / 5.))
                                                                          , ( ( 2, 0, 2,  2 ),       -sqrt( 48. / 5.))]), None)
                   , ('Apar', 'Apar')  :  ( _ba('Re_ang_Apar_Apar',       [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),        sqrt(  1. / 5.))
                                                                          , ( ( 2, 2, 2,  2 ),        sqrt(  3. / 5.))]), None)
                   , ('Aperp','Aperp') :  ( _ba('Re_ang_Aperp_Aperp',     [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),       -sqrt(  4. / 5.))]), None)
                   , ('A0',   'Apar')  :  ( _ba('Re_ang_A0_Apar',         [ ( ( 2, 1, 2, -2 ),       -sqrt( 24. / 5.))]), None)
                   , ('A0',   'Aperp') :  ( None, _ba('Im_ang_A0_Aperp',  [ ( ( 2, 1, 2,  1 ),  -1. * sqrt( 24. / 5.))]))
                   , ('Apar', 'Aperp') :  ( None, _ba('Im_ang_Apar_Aperp',[ ( ( 2, 2, 2, -1 ),  -1. * sqrt( 12. / 5.))]))
                   , ('AS',   'AS')    :  ( _ba('Re_ang_AS_AS',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),        sqrt(  4. / 5.))
                                                                          , ( ( 0, 0, 2,  2 ),       -sqrt( 12. / 5.))]),None)
                   , ('A0',   'AS')    :  ( _ba('Re_ang_A0_AS',           [ ( ( 1, 0, 0,  0 ),        sqrt(192.     ))
                                                                          , ( ( 1, 0, 2,  0 ),        sqrt( 48. / 5.))
                                                                          , ( ( 1, 0, 2,  2 ),       -sqrt(144. / 5.))]), None)
                   , ('Apar', 'AS')    :  ( _ba('Re_ang_Apar_AS',         [ ( ( 1, 1, 2, -2 ),       -sqrt( 72. / 5.))]), None)
                   , ('Aperp','AS')    :  ( None, _ba('Im_ang_Aperp_AS',  [ ( ( 1, 1, 2,  1 ), -1. * -sqrt( 72. / 5.))]))
                   }
        for k,v in angFuncs.iteritems() : self[k] = v


class JpsiKstTransversityAmplitudesHelicityAngles( AngularFunctions ) :
    def __init__( self, **kwargs ) :
        AngularFunctions.__init__(self)
        angles = kwargs.pop('Angles')
        coef   = kwargs.pop( 'DummyCoef', False )

        from P2VV.Parameterizations.GeneralUtils import getParNamePrefix
        namePF = getParNamePrefix(True)

        # define function to build angular function; create dummy coefficient if we need to split the angular efficiency functions later
        from P2VV.RooFitWrappers import P2VVAngleBasis, Addition, RealVar
        dummy = RealVar( Name = namePF + 'angEffDummyCoef', Value = 1. ) if coef else None
        _ba = lambda name, args :\
            Addition( namePF + name, [ P2VVAngleBasis( Name = namePF + name, Angles = angles, Indices = arg[0], FixedCoef = arg[1]
                                                      , Coefficient = dummy ) for arg in args ] )

        # TODO: generate the following table straight from the physics using PS->(VV,VS) ->ffss  (V=spin 1, f=spin 1/2, PS,S,s = spin 0)
        from math import sqrt
        angFuncs = { ('A0',   'A0')    :  ( _ba('Re_ang_A0_A0',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),       -sqrt( 16. / 5.))
                                                                          , ( ( 2, 0, 0,  0 ),        8.             )
                                                                          , ( ( 2, 0, 2,  0 ),       -sqrt( 64. / 5.))]), None)
                   , ('Apar', 'Apar')  :  ( _ba('Re_ang_Apar_Apar',       [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),        sqrt(  1. / 5.))
                                                                          , ( ( 2, 2, 2,  2 ),       -sqrt(  3. / 5.))]), None)
                   , ('Aperp','Aperp') :  ( _ba('Re_ang_Aperp_Aperp',     [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),        sqrt(  1. / 5.))
                                                                          , ( ( 2, 2, 2,  2 ),        sqrt(  3. / 5.))]), None)
                   , ('A0',   'Apar')  :  ( _ba('Re_ang_A0_Apar',         [ ( ( 2, 1, 2,  1 ),        sqrt( 24. / 5.))]), None)
                   , ('A0',   'Aperp') :  ( None, _ba('Im_ang_A0_Aperp',  [ ( ( 2, 1, 2, -1 ), -1. * -sqrt( 24. / 5.))]))
                   , ('Apar', 'Aperp') :  ( None, _ba('Im_ang_Apar_Aperp',[ ( ( 2, 2, 2, -2 ), -1. *  sqrt( 12. / 5.))]))
                   , ('AS',   'AS')    :  ( _ba('Re_ang_AS_AS',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),       -sqrt( 16. / 5.))]), None)
                   , ('A0',   'AS')    :  ( _ba('Re_ang_A0_AS',           [ ( ( 1, 0, 0,  0 ),        sqrt(192.     ))
                                                                          , ( ( 1, 0, 2,  0 ),       -sqrt(192. / 5.))]), None)
                   , ('Apar', 'AS')    :  ( _ba('Re_ang_Apar_AS',         [ ( ( 1, 1, 2,  1 ),        sqrt( 72. / 5.))]), None)
                   , ('Aperp','AS')    :  ( None, _ba('Im_ang_Aperp_AS',  [ ( ( 1, 1, 2, -1 ), -1. *  sqrt( 72. / 5.))]))
                   }
        for k,v in angFuncs.iteritems() : self[k] = v


class JpsiKstTransversityAmplitudesTransversityAngles( AngularFunctions ) :
    def __init__( self, **kwargs ) :
        AngularFunctions.__init__(self)
        angles = kwargs.pop('Angles')
        coef   = kwargs.pop( 'DummyCoef', False )

        from P2VV.Parameterizations.GeneralUtils import getParNamePrefix
        namePF = getParNamePrefix(True)

        # define function to build angular function; create dummy coefficient if we need to split the angular efficiency functions later
        from P2VV.RooFitWrappers import P2VVAngleBasis, Addition, RealVar
        dummy = RealVar( Name = namePF + 'angEffDummyCoef', Value = 1. ) if coef else None
        _ba = lambda name, args :\
            Addition( namePF + name, [ P2VVAngleBasis( Name = namePF + name, Angles = angles, Indices = arg[0], FixedCoef = arg[1]
                                                      , Coefficient = dummy ) for arg in args ] )

        # TODO: generate the following table straight from the physics using PS->(VV,VS) ->ffss  (V=spin 1, f=spin 1/2, PS,S,s = spin 0)
        from math import sqrt
        angFuncs = { ('A0',   'A0')    :  ( _ba('Re_ang_A0_A0',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),        sqrt(  4. / 5.))
                                                                          , ( ( 0, 0, 2,  2 ),       -sqrt( 12. / 5.))
                                                                          , ( ( 2, 0, 0,  0 ),        8.             )
                                                                          , ( ( 2, 0, 2,  0 ),        sqrt( 16. / 5.))
                                                                          , ( ( 2, 0, 2,  2 ),       -sqrt( 48. / 5.))]), None)
                   , ('Apar', 'Apar')  :  ( _ba('Re_ang_Apar_Apar',       [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),        sqrt(  1. / 5.))
                                                                          , ( ( 2, 2, 2,  2 ),        sqrt(  3. / 5.))]), None)
                   , ('Aperp','Aperp') :  ( _ba('Re_ang_Aperp_Aperp',     [ ( ( 2, 2, 0,  0 ),        2.             )
                                                                          , ( ( 2, 2, 2,  0 ),       -sqrt(  4. / 5.))]), None)
                   , ('A0',   'Apar')  :  ( _ba('Re_ang_A0_Apar',         [ ( ( 2, 1, 2, -2 ),       -sqrt( 24. / 5.))]), None)
                   , ('A0',   'Aperp') :  ( None, _ba('Im_ang_A0_Aperp',  [ ( ( 2, 1, 2,  1 ),  -1. * sqrt( 24. / 5.))]))
                   , ('Apar', 'Aperp') :  ( None, _ba('Im_ang_Apar_Aperp',[ ( ( 2, 2, 2, -1 ),  -1. * sqrt( 12. / 5.))]))
                   , ('AS',   'AS')    :  ( _ba('Re_ang_AS_AS',           [ ( ( 0, 0, 0,  0 ),        4.             )
                                                                          , ( ( 0, 0, 2,  0 ),        sqrt(  4. / 5.))
                                                                          , ( ( 0, 0, 2,  2 ),       -sqrt( 12. / 5.))]),None)
                   , ('A0',   'AS')    :  ( _ba('Re_ang_A0_AS',           [ ( ( 1, 0, 0,  0 ),        sqrt(192.     ))
                                                                          , ( ( 1, 0, 2,  0 ),        sqrt( 48. / 5.))
                                                                          , ( ( 1, 0, 2,  2 ),       -sqrt(144. / 5.))]), None)
                   , ('Apar', 'AS')    :  ( _ba('Re_ang_Apar_AS',         [ ( ( 1, 1, 2, -2 ),       -sqrt( 72. / 5.))]), None)
                   , ('Aperp','AS')    :  ( None, _ba('Im_ang_Aperp_AS',  [ ( ( 1, 1, 2,  1 ), -1. * -sqrt( 72. / 5.))]))
                   }
        for k,v in angFuncs.iteritems() : self[k] = v

