from TimeAnglesPdf import *

class Bs2JpsiPhi(TimeAnglesPdf):
    def __init__(self, Observables, Parametrisation, Parameters, ResolutionModel):
        self._timeComponents = {'|A0|^2'  : {'cosh' : '1', 'cos' : 'C', 'sinh' : '-D', 'sin' : '-S'},
                                '|Apa|^2' : {'cosh' : '1', 'cos' : 'C', 'sinh' : '-D', 'sin' : '-S'},
                                }
