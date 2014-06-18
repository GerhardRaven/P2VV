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

def set_classes_creates(creators, classes = None):
    __creators.update(creators)
    if classes == None:
        classes = creators.iterkeys()
    for cl in classes:    
        _temp = __import__('ROOT', globals(), locals(), [cl], -1)
        cl_type = getattr(_temp, cl)
        bm = set()
        __get_base_methods(cl, bm)
        for method in set(__creators[cl]) | bm:
            method = getattr(cl_type, method)
            __creates(method)
    
def get_creators():
    return __creators

set_classes_creates({}, __creators.iterkeys())

def get_roofit_classes(rootmap):
    classes = []
    with open(rootmap) as rfile:
        for line in rfile:
            cl = line.strip().split()[0].split('.', 1)[1]
            if not cl.startswith('Roo'):
                continue
            
            cl = cl.strip(':')
            cl = cl.replace('-', ' ').replace('@', ':')
            classes.append(cl)
    return classes

def get_preturn_methods(classes):
    from ROOT import gInterpreter
    from ROOT import TDictionary

    ignore = set(['name', 'argName', 'clone', 'Clone', 'Class', 'Class_Name', 'IsA',
                  'DeclFileName', 'ImplFileName', 'GetName'])

    from collections import defaultdict
    p_returns = defaultdict(set)
    for cl in classes:
        cinfo = gInterpreter.ClassInfo_Factory(cl)
        minfo = gInterpreter.MethodInfo_Factory(cinfo)
        while gInterpreter.MethodInfo_Next(minfo):
            if not bool(gInterpreter.MethodInfo_Property(minfo) & TDictionary.kIsPublic):
                continue
            
            tinfo = gInterpreter.MethodInfo_Type(minfo)
            name = gInterpreter.MethodInfo_Name(minfo)
            if ((bool(gInterpreter.TypeInfo_Property(tinfo) & TDictionary.kIsPointer) or
                bool(gInterpreter.TypeInfo_Property(tinfo) & TDictionary.kIsConstPointer))
                and not name in ignore):
                p_returns[cl].add(name)

def __wrap_font_embed(fun):
    import os
    from functools import wraps
    @wraps(fun)
    def _fun(self, *args, **kwargs):
        embed = kwargs.get('EmbedFonts', False)
        font_path = []
        if 'FontPath' in kwargs:
            font_path.append(kwargs['FontPath'])
        elif 'FontPaths' in kwargs:
            font_path.extend(kwargs['FontPaths'])
        from ROOT import gEnv
        root_font_path = gEnv.GetValue('Root.TTFontPath', 'unset')
        if root_font_path != 'unset':
            font_path.append(root_font_path)
        font_path = ':'.join(font_path)
        pdf = None
        if len(args) >= 1 and type(args[0] == str) and '.pdf' in args[0]:
            pdf = os.path.realpath(args[0])
        r = fun(self, *args)
        if pdf and embed:
            import tempfile
            tmp = tempfile.NamedTemporaryFile()
            tmp_name = os.path.realpath(tmp.name)
            tmp.close()
            cmd = ['gs', '-sFONTPATH=%s' % font_path, '-o', tmp_name, '-sDEVICE=pdfwrite',
                   '-dPDFSETTINGS=/prepress', pdf]
            import subprocess
            p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
            o, e = p.communicate()
            os.remove(pdf)
            import shutil
            shutil.move(tmp_name, pdf)
        return r
    return _fun

from ROOT import TPad
TPad.Print = __wrap_font_embed(TPad.Print)
