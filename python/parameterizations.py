
class CPParam :
    def __init__(self,**kwargs) :
        for i in 'CDS' : setattr(self,i,kwargs.pop(i))


class LambdaSqArg_CPParam( CPParam ) :
    def __init__(self, **kwargs) :
        CPParam.__init__(self, C = FormulaVar('C', '(1.-@0)/(1.+@0)',                   [ kwargs['lambdaSq']              ] )
                             , D = FormulaVar('D', ' 2 * sqrt(@0) * cos(@1) / (1.+@0)', [ kwargs['lambdaSq'], kwargs['arg'] ] )
                             , S = FormulaVar('S', '-2 * sqrt(@0) * sin(@1) / (1.+@0)', [ kwargs['lambdaSq'], kwargs['arg'] ] )
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
from RooFitWrappers import FormulaVar
class Polar2_Amplitude(Carthesian_Amplitude) :
    def __init__(self,name, r2, arg, CP ) :
        Carthesian_Amplitude.__init__( self,  name, FormulaVar('Re_%s'%name, 'sqrt(@0) * cos(@1)', [r2,arg], Title = 'Re(%s)'% name )
                                                  , FormulaVar('Im_%s'%name, 'sqrt(@0) * sin(@1)', [r2,arg], Title = 'Im(%s)'% name )
                                                  , CP )


class CEvenOdd :
    def __init__(self, **kwargs ) :
        for i in ['avgCEven','avgCOdd' ] : setattr(self,i,kwargs.pop(i))
    def __getitem__(self,kw) :
        return getattr(self,kw)

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
        CEvenOdd.__init__(self, avgCEven =  FormulaVar( 'avgCEven', '1. + @0*@1 + @0*@2 + @1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average even coefficients') 
                              , avgCOdd  =  FormulaVar( 'avgCOdd',     '@0 + @1 + @2 + @0*@1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average odd coefficients') 
                              )

