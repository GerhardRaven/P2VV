###### decorate TPad with pads...
from ROOT import TPad
def __pads(self,n=None,m=None) :
    if n : 
        if m : self.Divide(n,m)
        else : self.Divide(n)
    i=1
    while self.GetPad(i) :
        yield self.cd(i) 
        i += 1

def __frames(self) :
    for prim in self.GetListOfPrimitives() :
        if type(prim) == TPad :
            for prim1 in prim.frames() : yield prim1
        elif prim.GetName().startswith('TFrame') :
            yield prim

def __frameHists(self) :
    for prim in self.GetListOfPrimitives() :
        if type(prim) == TPad :
            for prim1 in prim.frameHists() : yield prim1
        elif prim.GetName().startswith('frame') :
            yield prim

TPad.pads       = __pads
TPad.frames     = __frames
TPad.frameHists = __frameHists
