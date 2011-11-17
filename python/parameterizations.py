
class CPParam :
    def __init__(self,**kwargs) :
        for i in 'CDS' : setattr(self,i,kwargs.pop(i))

class LambdaCarth_CPParam( CPParam ) :
    def __init__(self, **kwargs) :
        from RooFitWrappers import FormulaVar
        CPParam.__init__(self, C = FormulaVar('C', '(1.-@0*@0-@1*@1) / (1. + @0*@0 + @1*@1)', [ kwargs['ReLambda'], kwargs['ImLambda'] ] )
                             , D = FormulaVar('D', '2. * @0 / (1.+@0*@0+@1*@1)',              [ kwargs['ReLambda'], kwargs['ImLambda'] ] )
                             , S = FormulaVar('S', '2. * @1 / (1.+@0*@0+@1*@1)',              [ kwargs['ReLambda'], kwargs['ImLambda'] ] )
                        )

class LambdaSqArg_CPParam( CPParam ) :
    def __init__(self, **kwargs) :
        from RooFitWrappers import FormulaVar
        CPParam.__init__(self, C = FormulaVar('C', '(1.-@0)/(1.+@0)',                  [ kwargs['lambdaSq']              ] )
                             , D = FormulaVar('D', '2 * sqrt(@0) * cos(@1) / (1.+@0)', [ kwargs['lambdaSq'], kwargs['lambdaArg'] ] )
                             , S = FormulaVar('S', '2 * sqrt(@0) * sin(@1) / (1.+@0)', [ kwargs['lambdaSq'], kwargs['lambdaArg'] ] )
                        )

# construct amplitudes with carthesian parameters
class Carthesian_Amplitude :
    def __init__(self,name, Re, Im, CP ) :
        self.name = name
        self.Re = Re
        self.Im = Im
        self.CP = CP # even or odd???
    def __str__(self) : return self.name


# construct amplitudes with polar parameters
class Polar2_Amplitude(Carthesian_Amplitude) :
    def __init__(self,name, r2, arg, CP ) :
        from RooFitWrappers import FormulaVar
        Carthesian_Amplitude.__init__( self,  name, FormulaVar('Re_%s'%name, 'sqrt(@0) * cos(@1)', [r2,arg], Title = 'Re(%s)'% name )
                                                  , FormulaVar('Im_%s'%name, 'sqrt(@0) * sin(@1)', [r2,arg], Title = 'Im(%s)'% name )
                                                  , CP )

class AmplitudeSet ( dict ) :
    def __init__( self, *args) :
        # maybe make this thing readonly???
        for v in args: 
            self[ v.name ] = v
            assert hasattr(v,'Re')
            assert hasattr(v,'Im')
            assert hasattr(v,'CP')
        # require the names in args to be unique...
        assert(len(self)==len(args))

class JpsiphiAmplitudesLP2011 ( AmplitudeSet ) :
    def __init__( self ) :
        from RooFitWrappers import RealVar, FormulaVar
        from math import pi
        _A0Mag2    = RealVar('A0Mag2',    Title = '|A0|^2',      Value = 0.556,   MinMax = (0., 1.))
        _A0Ph      = RealVar('delta0',    Title = 'delta_0',     Value = 0. )
        _AperpMag2 = RealVar('AperpMag2', Title = '|A_perp|^2',  Value = 0.233,   MinMax = ( 0., 1.))
        _AperpPh   = RealVar('deltaPerp', Title = 'delta_perp',  Value = pi-2.91, MinMax=( -2. * pi, 2. * pi))
        _AparMag2  = FormulaVar('AparMag2', '1. - @0 - @1', [_A0Mag2, _AperpMag2],  Title = '|A_par|^2' )
        _AparPh    = RealVar('deltaPar',  Title = 'delta_par',   Value = 2.93,    MinMax = ( -2. * pi, 2. * pi))
        _ASMag2    = RealVar('ASMag2',    Title = '|A_S|^2',     Value = 0.05,    MinMax=( 0., 1.))
        _ASPh      = RealVar('deltaS',    Title = 'delta_S',     Value = 2.2,     MinMax=( -2. * pi, 2. * pi))

        AmplitudeSet.__init__( self, Polar2_Amplitude( 'A0',    _A0Mag2,    _A0Ph,    +1 )
                                   , Polar2_Amplitude( 'Apar',  _AparMag2,  _AparPh,  +1 )
                                   , Polar2_Amplitude( 'Aperp', _AperpMag2, _AperpPh, -1 )
                                   , Polar2_Amplitude( 'AS',    _ASMag2,    _ASPh,    -1 )
                                   )

class CEvenOdd :
    def __init__(self, **kwargs ) :
        for i in ['avgCEven','avgCOdd' ] : setattr(self,i,kwargs.pop(i))
    def __getitem__(self,kw) : return getattr(self,kw)

class Trivial_CEvenOdd( CEvenOdd ) :
    def __init__(self) :
        from RooFitWrappers import ConstVar
        CEvenOdd.__init__(self, avgCEven =  ConstVar('one', Value = 1 )
                              , avgCOdd  =  ConstVar('zero', Value = 0 )
                              )

class ProdTagNorm_CEvenOdd( CEvenOdd ) :
    def __init__(self,**kwargs) :
        _AProd = kwargs.pop('AProd')
        _ANorm = kwargs.pop('ANorm')
        _ATagEff = kwargs.pop('ATagEff')
        if kwargs : raise KeyError('unknown keyword argument: %s' % kwargs )
        from RooFitWrappers import FormulaVar
        CEvenOdd.__init__(self, avgCEven =  FormulaVar( 'avgCEven', '1. + @0*@1 + @0*@2 + @1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average even coefficients') 
                              , avgCOdd  =  FormulaVar( 'avgCOdd',     '@0 + @1 + @2 + @0*@1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average odd coefficients') 
                              )



#TODO: inherit from UserDict mixin instead of wrapping & forwarding...
class AngularFunctions :
    def __init__(self) :         self._d = dict()
    def __getitem__(self,k) :    return self._d[k]
    def __setitem__(self,k,v) :  self._d[k] = v
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
        if name in kwargs and type( kwargs[name] )== str :
            return RealVar(  kwargs[name], Title = title,   Observable = True,  MinMax=minmax )
        elif name not in kwargs :
            return  RealVar(  defname,   Title = title,   Observable = True,  MinMax=minmax)
        else  :
            return kwargs[name]


class JpsiphiHelicityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = { 'cpsi'   : AngleDefinitions.__intrprt__( 'cpsi', kwargs, 'helcthetaK', 'cosine of kaon polarization angle',      (-1., 1.))
            , 'ctheta' : AngleDefinitions.__intrprt__( 'ctheta', kwargs,  'helcthetaL', 'cosine of lepton polarization angle', (-1., 1.))
            , 'phi'    : AngleDefinitions.__intrprt__( 'phi', kwargs,  'helphi', 'angle between decay planes',                 (-pi, pi))
            }
        d['functions'] =  JpsiphiTransversityAmplitudesHelicityAngles( **d )
        AngleDefinitions.__init__(self, **d )

class JpsiphiTransversityAngles( AngleDefinitions ) :
    def __init__( self, **kwargs ) :
        from math import pi
        d = { 'cpsi' :   AngleDefinitions.__intrprt__( 'cpsi', kwargs,  'trcpsi', 'cosine of kaon polarization angle', (-1., 1.))
            , 'ctheta' : AngleDefinitions.__intrprt__( 'ctheta',kwargs, 'trctheta', 'cosine of transversity polar angle',     (-1., 1.))
            , 'phi'   :  AngleDefinitions.__intrprt__( 'phi', kwargs, 'trphi',     'transversity azimuthal angle',          (-pi, pi))
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
        angFuncs = { ('A0',   'A0')    :  ( _ba('Re_ang_A0_A0',           [(0, 0, 0,  0,  4.             )
                                                                          ,(0, 0, 2,  0, -sqrt( 16. / 5.))
                                                                          ,(2, 0, 0,  0,  8.             )
                                                                          ,(2, 0, 2,  0, -sqrt( 64. / 5.))]), None)
                   , ('Apar', 'Apar')  :  ( _ba('Re_ang_Apar_Apar',       [(2, 2, 0,  0,  2.             )
                                                                          ,(2, 2, 2,  0,  sqrt(  1. / 5.))
                                                                          ,(2, 2, 2,  2, -sqrt(  3. / 5.))]), None)
                   , ('Aperp','Aperp') :  ( _ba('Re_ang_Aperp_Aperp',     [(2, 2, 0,  0,  2.             )
                                                                          ,(2, 2, 2,  0,  sqrt(  1. / 5.))
                                                                          ,(2, 2, 2,  2,  sqrt(  3. / 5.))]), None)
                   , ('A0',   'Apar')  :  ( _ba('Re_ang_A0_Apar',         [(2, 1, 2,  1,  sqrt( 24. / 5.))]), None)
                   , ('A0',   'Aperp') :  ( None, _ba('Im_ang_A0_Aperp',  [(2, 1, 2, -1, -sqrt( 24. / 5.))]))
                   , ('Apar', 'Aperp') :  ( None, _ba('Im_ang_Apar_Aperp',[(2, 2, 2, -2,  sqrt( 12. / 5.))])) 
                   , ('AS',   'AS')    :  ( _ba('Re_ang_AS_AS',           [(0, 0, 0,  0,  4.             )
                                                                          ,(0, 0, 2,  0, -sqrt( 16. / 5.))]), None)
                   , ('A0',   'AS')    :  ( _ba('Re_ang_A0_AS',           [(1, 0, 0,  0,  sqrt(192.     ))
                                                                          ,(1, 0, 2,  0, -sqrt(192. / 5.))]), None)
                   , ('Apar', 'AS')    :  ( _ba('Re_ang_Apar_AS',         [(1, 1, 2,  1,  sqrt( 72. / 5.))]), None)
                   , ('Aperp','AS')    :  ( None, _ba('Im_ang_Aperp_AS',  [(1, 1, 2, -1,  sqrt( 72. / 5.))]))
                   }
        for k,v in angFuncs.iteritems() : self[k] = v

class JpsiphiTransversityAmplitudesTransversityAngles( AngularFunctions ) :
    def __init__( self, **kwargs ) :
        AngularFunctions.__init__(self)
        from RooFitWrappers import P2VVAngleBasis, Addition
        from math import sqrt
        _ba = lambda  name,args : Addition(name, [ P2VVAngleBasis(kwargs , *a) for a in args ] )
        # TODO: generate the following table straight from the physics using PS->(VV,VS) ->ffss  (V=spin 1, f=spin 1/2, PS,S,s = spin 0)

        angFuncs = { ('A0',   'A0')    :  ( _ba('Re_ang_A0_A0',           [(0, 0, 0,  0,  4.             )
                                                                          ,(0, 0, 2,  0,  sqrt(  4. / 5.))
                                                                          ,(0, 0, 2,  2, -sqrt( 12. / 5.))
                                                                          ,(2, 0, 0,  0,  8.             )
                                                                          ,(2, 0, 2,  0,  sqrt( 16. / 5.))
                                                                          ,(2, 0, 2,  2, -sqrt( 48. / 5.))]), None)
                   , ('Apar', 'Apar')  :  ( _ba('Re_ang_Apar_Apar',       [(2, 2, 0,  0,  2.             )
                                                                          ,(2, 2, 2,  0,  sqrt(  1. / 5.))
                                                                          ,(2, 2, 2,  2,  sqrt(  3. / 5.))]), None)
                   , ('Aperp','Aperp') :  ( _ba('Re_ang_Aperp_Aperp',     [(2, 2, 0,  0,  2.             )
                                                                          ,(2, 2, 2,  0, -sqrt(  4. / 5.))]), None)

                   , ('A0',   'Apar')  :  ( _ba('Re_ang_A0_Apar',         [(2, 1, 2, -2, -sqrt( 24. / 5.))]), None)
                   , ('A0',   'Aperp') :  ( None, _ba('Im_ang_A0_Aperp',  [(2, 1, 2,  1,  sqrt( 24. / 5.))]))
                   , ('Apar', 'Aperp') :  ( None, _ba('Im_ang_Apar_Aperp',[(2, 2, 2, -1,  sqrt( 12. / 5.))]))

                   , ('AS',   'AS')    :  ( _ba('Re_ang_AS_AS',           [(0, 0, 0,  0,  4.             )
                                                                          ,(0, 0, 2,  0,  sqrt(  4. / 5.))
                                                                          ,(0, 0, 2,  2, -sqrt( 12. / 5.))]),None)
                   , ('A0',   'AS')    :  ( _ba('Re_ang_A0_AS',           [(1, 0, 0,  0,  sqrt(192.     ))
                                                                          ,(1, 0, 2,  0,  sqrt( 48. / 5.))
                                                                          ,(1, 0, 2,  2, -sqrt(144. / 5.))]), None)
                   , ('Apar', 'AS')    :  ( _ba('Re_ang_Apar_AS',         [(1, 1, 2, -2, -sqrt( 72. / 5.))]), None)
                   , ('Aperp','AS')    :  ( None, _ba('Im_ang_Aperp_AS',  [(1, 1, 2,  1, -sqrt( 72. / 5.))]))
                   }
        for k,v in angFuncs.iteritems() : self[k] = v


class BTagDecayBasisCoefficients :
    def __init__(self, **kwargs ) :
        for i in ['sin','cos','sinh','cosh' ] : setattr(self,i,kwargs.pop(i))
        if kwargs : raise KeyError('unknown keyword arguments: %s' % kwargs )
    def __getitem__(self,kw) :
        return getattr(self,kw)

class JpsiphiBTagDecayBasisCoefficients( BTagDecayBasisCoefficients ) :
    def __init__(self,  angFuncs, Amplitudes,CP, order ) :
        def combine( name, afun, A, CPparams, i, j) :
            from RooFitWrappers import ConstVar, FormulaVar, Product
            plus  = ConstVar('plus', Value = 1)
            minus = ConstVar('minus',  Value = -1  )
            one   = plus
            # define functions which return Re(Conj(Ai) Aj), Im( Conj(Ai) Aj)
            # TODO: replace by Addition & Product...
            Re        = lambda ai, aj  : FormulaVar('Re_c_%s_%s'%(ai,aj),'@0*@2+@1*@3',[ai.Re,ai.Im,aj.Re,aj.Im])
            Im        = lambda ai, aj  : FormulaVar('Im_c_%s_%s'%(ai,aj),'@0*@3-@1*@2',[ai.Re,ai.Im,aj.Re,aj.Im])
            # define functions which return the coefficients that define the time-dependence...
            _minus_if = lambda b, x : [ minus, x ] if b else [ x ]
            coef = { 'cosh' : lambda ai,aj,CP : ( one  if ai.CP == aj.CP else CP.C  # make a trivial product just to get the labels right???
                                                , None )
                   , 'cos'  : lambda ai,aj,CP : ( CP.C if ai.CP == aj.CP else one 
                                                , None )
                   , 'sinh' : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP > 0     ,  CP.D ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP < aj.CP ,  CP.S )) )
                   , 'sin'  : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP > 0     ,  CP.S ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP > aj.CP ,  CP.D )) )
                   }
            (c_re,c_im) = coef[name](A[i],A[j],CPparams)
            (f_re,f_im) = afun[(i,j)]
            (a_re,a_im) = ( Re(A[i],A[j]),Im(A[i],A[j]) )
            # this triplet of complex numbers used to be written recursively as a doublet of a single number and another doublet...
            # hence the current structure: Re(xyz) =  Re(x)Re(yz) + Im(x)Im(yz) 
            # TODO: move some minus sign around (ie into afun and coef) so that
            # NOTE: this becomes just the obvious Re(a b c)  = Re(a)Re(b)Re(c) - Re(a)Im(b)Im(c) - Im(a)Re(b)Im(c) - Im(a)Im(b)Re(c)....
            prod = lambda name, args : [ Product(name, args) ] if all(args) else []
            s  = prod('ReReRe_%s_%s_%s'%(name,A[i],A[j]), [        a_re , c_re, f_re ] ) \
               + prod('ImImRe_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_im , c_im, f_re ] ) \
               + prod('ImReIm_%s_%s_%s'%(name,A[i],A[j]), [        a_im , c_re, f_im ] ) \
               + prod('ReImIm_%s_%s_%s'%(name,A[i],A[j]), [        a_re , c_im, f_im ] )
            assert len(s) == 1 # for now, coefficients are either real, or imaginary, but not both... (not true in general, but I'm lazy today ;-)
            return s[0]

        args = dict()
        from RooFitWrappers import Addition
        try : # this requires python 2.7 or later...
            from itertools import combinations_with_replacement as cwr
        except:
            from compatibility import cwr

        for name in [ 'cosh', 'sinh', 'cos', 'sin' ] :
            # NOTE: 'Amplitudes'  must be traversed 'in order' : A0, Apar, Aperp, AS -- so we cannot use Amplitudes.keys() out of the box...
            args[ name ] = Addition( 'a_%s'% name, [ combine(name,angFuncs,Amplitudes,CP,i,j) for (i,j) in cwr( order, 2 ) ] )

        BTagDecayBasisCoefficients.__init__( self, **args )




class ResolutionModelLP2011 :
    def __init__( self, t ) :
        from RooFitWrappers import RealVar, ConstVar, ResolutionModel, AddModel
        mu = ConstVar('tres_mu', Value = -0.0027 )
        SF = RealVar('tres_SF', Value = 1.0, MinMax = (0.5,1.5) )
        from ROOT import RooGaussModel as GaussModel
        sigmas = [ (3,0.513 ), (2,0.0853), (1,0.0434) ]
        frac   = [ (3,0.0017), (2,0.165) ]
        self.Model = AddModel('tres', [ ResolutionModel('tres_%s'%n, Type = GaussModel, Observables = [ t ], Parameters = [ mu, ConstVar('tres_s%s'%n, Value = v  ), SF ] ) for n,v in sigmas ]
                                    , [ ConstVar('tres_f%s'%n,Value = v) for n,v in frac ]
                             )

