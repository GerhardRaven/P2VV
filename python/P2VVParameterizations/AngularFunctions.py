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


class AngleDefinitions( object ) :
    def __init__( self , **kwargs ) :
        self.angles = dict( (i,kwargs.pop(i)) for i in ['cpsi' ,'ctheta','phi' ] )
        self.functions = kwargs.pop('functions')

    @staticmethod
    def __intrprt__(name,kwargs,defname,title,minmax) :
        from RooFitWrappers import RealVar
        if name in kwargs and type( kwargs[name] ) == str :
            return RealVar(  Name = kwargs[name], Title = title, Observable = True,  MinMax=minmax )
        elif name not in kwargs :
            return  RealVar(  Name = defname,     Title = title, Observable = True,  MinMax=minmax)
        else  :
            return kwargs[name]


class JpsiphiHelicityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = { 'cpsi'   : AngleDefinitions.__intrprt__( 'cpsi',   kwargs, 'helcthetaK', 'Cosine of kaon polarization angle',   ( -1., 1. ) )
            , 'ctheta' : AngleDefinitions.__intrprt__( 'ctheta', kwargs, 'helcthetaL', 'Cosine of lepton polarization angle', ( -1., 1. ) )
            , 'phi'    : AngleDefinitions.__intrprt__( 'phi',    kwargs, 'helphi',     'Angle between decay planes',          ( -pi, pi ) )
            }
        d['functions'] =  JpsiphiTransversityAmplitudesHelicityAngles( **d )
        AngleDefinitions.__init__(self, **d )


class JpsiphiTransversityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = { 'cpsi' :   AngleDefinitions.__intrprt__( 'cpsi',   kwargs, 'trcpsi',   'Cosine of kaon polarization angle',  ( -1., 1. ) )
            , 'ctheta' : AngleDefinitions.__intrprt__( 'ctheta', kwargs, 'trctheta', 'Cosine of transversity polar angle', ( -1., 1. ) )
            , 'phi'   :  AngleDefinitions.__intrprt__( 'phi',    kwargs, 'trphi',    'Transversity azimuthal angle',       ( -pi, pi ) )
            }
        d['functions'] = JpsiphiTransversityAmplitudesTransversityAngles( **d )
        AngleDefinitions.__init__(self, **d )


class JpsiphiTransversityAmplitudesHelicityAngles( AngularFunctions ) :
    def __init__( self, **kwargs ) :
        AngularFunctions.__init__(self)
        from RooFitWrappers import P2VVAngleBasis, Addition
        from math import sqrt
        _ba = lambda  name,args : Addition(name, [ P2VVAngleBasis(kwargs , *a) for a in args ] )
        # TODO: generate the following table straight from the physics using PS->(VV,VS) ->ffss  (V=spin 1, f=spin 1/2, PS,S,s = spin 0)
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
        from RooFitWrappers import P2VVAngleBasis, Addition
        from math import sqrt
        _ba = lambda  name,args : Addition(name, [ P2VVAngleBasis(kwargs , *a) for a in args ] )
        # TODO: generate the following table straight from the physics using PS->(VV,VS) ->ffss  (V=spin 1, f=spin 1/2, PS,S,s = spin 0)
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

