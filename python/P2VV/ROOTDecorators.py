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

__creators = {'TTree'  : ['CopyTree'],
              'TChain' : []}

def __creates(method):
    if hasattr(method, 'im_func') and method.im_func.func_closure:
        method = method.im_func.func_closure[-1].cell_contents
    method._creates = True

def __get_base_methods(cl, methods):
    from ROOT import gInterpreter
    cinfo = gInterpreter.ClassInfo_Factory(cl)
    binfo = gInterpreter.BaseClassInfo_Factory(cinfo)
    ## Loop over base classes
    while gInterpreter.BaseClassInfo_Next(binfo):
        bname = gInterpreter.BaseClassInfo_Name(binfo)
        if bname in __creators:
            methods |= set(__creators[bname])
        __get_base_methods(bname, methods)

def set_class_creates(creators, cl):
    __creators.update(creators)
    _temp = __import__('ROOT', globals(), locals(), [cl], -1)
    cl_type = getattr(_temp, cl)
    bm = set()
    __get_base_methods(cl, bm)
    for method in set(__creators[cl]) | bm:
        method = getattr(cl_type, method)
        __creates(method)
    
def get_creators():
    return __creators

for cl in __creators.keys():
    set_class_creates({}, cl)
