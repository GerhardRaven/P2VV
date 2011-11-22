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

class Cart_Amplitude(object) :
   def __init__(self,name,CP,x,y) :
     self.name = name
     self.CP = CP
     self.Re = x
     self.Im = y
     self.r2 = '%s^2+%s^2'%(x,y)

class Polar_Amplitude(object) :
   def __init__(self,name,CP,phase,mag) :
     self.name = name
     self.CP = CP
     self.Re = 'cos(%s)*%s'%(phase,mag)
     self.Im = 'sin(%s)*%s'%(phase,mag)
     self.r2 = '%s^2' % mag

class Polar_Amplitude2(object) :
   def __init__(self,name,CP,phase,mag2) :
     self.name = name
     self.CP = CP
     self.Re = 'cos(%s)*sqrt(%s)'%(phase,mag2)
     self.Im = 'sin(%s)*sqrt(%s)'%(phase,mag2)
     self.r2 = '%s' % mag2

## TODO(?): replace the list below with a dedicated class
Amplitudes = [ Polar_Amplitude('Apar',+1,'delta_par','r_par')
             , Polar_Amplitude('Aperp',-1,'delta_perp','r_perp')
             , Polar_Amplitude2('A0',+1,'zero','one')
             , Polar_Amplitude('As',-1,'delta_S','r_S') ]

l = []
from itertools import combinations_with_replacement
for (i,j) in combinations_with_replacement( Amplitudes, 2 ) :
   if i == j : l+=[ i.r2 ]
   else :
      l+= [ 'Re(conj(%s) %s)'%(i.name,j.name) ] # too lazy to write this in terms of {i,j}.{Re,Im} 
      if i.CP != j.CP : l+= [ 'Im(conj(%s) %s)'%(i.name,j.name) ] # too lazy to write this in terms of {i,j}.{Re,Im} 

from pprint import pprint
pprint(l)
