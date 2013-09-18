###########################################################################################################################################
## P2VVParameterizations.BBbarAsymmetries: Parameterizations of asymmetries between B and Bbar                                           ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin


class CEvenOdd ( _util_parse_mixin ) :
    def __init__( self, **kwargs ) :
        for i in [ 'avgCEven', 'avgCOdd' ] : setattr( self, i, kwargs.pop(i) )

    def __getitem__( self, kw ) : return getattr( self, kw )


class Trivial_CEvenOdd( CEvenOdd ) :
    def __init__( self ) :
        from P2VV.RooFitWrappers import ConstVar

        self._zero   = ConstVar( Name = 'zero', Value = 0 )
        self._one    = ConstVar( Name =  'one', Value = 1 )
        self._params = [ self._one, self._zero ]

        CEvenOdd.__init__( self, avgCEven = self._one, avgCOdd = self._zero )


class Coefficients_CEvenOdd( CEvenOdd ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'avgCEven', kwargs, Title = 'CP average even coefficients', Value = 1., MinMax = (  0.8, 1.2 ) )
        self._parseArg( 'avgCOdd',  kwargs, Title = 'CP average odd coefficients',  Value = 0., MinMax = ( -1.,  1.  ) )

        self._check_extraneous_kw( kwargs ) 
        CEvenOdd.__init__(self, avgCEven = self._avgCEven, avgCOdd = self._avgCOdd )


class ProdTagNorm_CEvenOdd( CEvenOdd ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import ConstVar

        self._parseArg( 'AProd',   kwargs, Title = 'production asymmetry',         Value = 0., MinMax = ( -1., 1. ) )
        self._parseArg( 'ATagEff', kwargs, Title = 'tagging efficiency asymmetry', Value = 0., MinMax = ( -1., 1. ) )

        if 'C' in kwargs or 'CPParam' in kwargs:
          self._minus = ConstVar( Name = 'minus', Value = -1. )
          if 'C' in kwargs :
            self._parseArg( 'C', kwargs, Title = 'CPV param. C', Value = 0., MinMax = ( -1., 1. ) )
          else :
            CPParam = kwargs.pop('CPParam')
            self._C = CPParam['C']
          self._ANorm = self._parseArg( 'ANorm', kwargs, Arguments = [ self._minus, self._C ], Title = 'normalization asymmetry'
                                       , ObjectType = 'Product' )
          self._params += [ self._minus, self._C ]
        else   :
          self._parseArg( 'ANorm', kwargs, Title = 'normalization asymmetry', Value = 0., MinMax = ( -1., 1. ) )

        self._check_extraneous_kw( kwargs ) 
        CEvenOdd.__init__( self, avgCEven = self._parseArg( 'avgCEven', kwargs, Formula = '1. + @0*@1 + @0*@2 + @1*@2'
                                                           , Arguments =  [ self._AProd, self._ANorm, self._ATagEff ]
                                                           , ObjectType = 'FormulaVar', Title = 'CP average even coefficients' )
                              , avgCOdd   = self._parseArg( 'avgCOdd', kwargs, Formula = '@0 + @1 + @2 + @0*@1*@2'
                                                           , Arguments = [ self._AProd, self._ANorm, self._ATagEff ]
                                                           , ObjectType = 'FormulaVar', Title = 'CP average odd coefficients' )
                         )
