
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
        self.CP = CP
    def __str__(self) : return self.name


# construct amplitudes with polar parameters
from RooFitWrappers import FormulaVar
class Polar2_Amplitude(Carthesian_Amplitude) :
    def __init__(self,name, r2, arg, CP ) :
        Carthesian_Amplitude.__init__( self,  name, FormulaVar('Re_%s'%name, 'sqrt(@0) * cos(@1)', [r2,arg], Title = 'Re(%s)'% name )
                                                  , FormulaVar('Im_%s'%name, 'sqrt(@0) * sin(@1)', [r2,arg], Title = 'Im(%s)'% name )
                                                  , CP )

