from RooFitWrappers import RealVar

class Base(object):
    """ A base class for all parametrisations. It defines the parameters which
    should be provided by all derived classes."""

    _parameters = set(['1', 'C', 'S','D'])

    def __init__(self, d):
        self._parameters = {}
        for k in Base._parameters:
            self._parameters[k] = d[k]
    
class Trivial(Base):
    def __init__(self):
        One = RealVar('One', Observable = False, Value = 1.)
        C   = RealVar('C', Observable = False, Value = 0., MinMax = (-1, 1))
        S   = RealVar('S', Observable = False, Value = 0., MinMax = (-1, 1))
        D   = RealVar('D', Observable = False, Value = 0., MinMax = (-1, 1))
        
        d = {'1' : One, 'C' : C, 'S' : S, 'D' : D}

        super(Trivial, self).__init__(name, d))

    def __getitem__(self, k):
        return self._parameters[k]
