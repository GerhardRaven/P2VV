###### decorate TPad with pads...
from ROOT import TPad
def __pads( self, n = None, m = None, predicate = lambda x : True ) :
    if n : 
        if m : self.Divide( n, m )
        else : self.Divide( n )
    i = 1
    while self.GetPad(i) :
        if predicate(i) : yield self.cd(i)
        i += 1

def __frames(self) :
    for prim in self.GetListOfPrimitives() :
        if isinstance(prim,TPad) :
            for prim1 in prim.frames() : yield prim1
        elif prim.GetName().startswith('TFrame') :
            yield prim

def __frameHists(self) :
    for prim in self.GetListOfPrimitives() :
        if isinstance(prim,TPad) :
            for prim1 in prim.frameHists() : yield prim1
        elif prim.GetName().startswith('frame') :
            yield prim

TPad.pads       = __pads
TPad.frames     = __frames
TPad.frameHists = __frameHists

def __ROOTversion() :
    from ROOT import gROOT
    versionInt = gROOT.GetVersionInt()
    versionMajor = versionInt/10000
    versionMinor = versionInt/100 - versionMajor*100
    versionPatch = versionInt%100 
    return (versionMajor, versionMinor, versionPatch)

ROOTversion = __ROOTversion()


from ROOT import TFile
TFile.__enter__ = lambda self : self
def __TFile__exit__(self,*args) : 
    if self : self.Close()
TFile.__exit__ = __TFile__exit__
