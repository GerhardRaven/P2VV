###############################################################################
## RooFitDecorators:                                                         ##
##                                                                           ##
## authors:                                                                  ##
##   GR, Gerhard Raven, Nikhef & VU, Gerhard.Raven@nikhef.nl                 ##
##                                                                           ##
###############################################################################

from ROOT import RooFit, RooArgSet, RooArgList, RooDataSet, RooWorkspace

# RooDataSet functions

def _RooDataSetIter(self) :
    for i in range(self.numEntries()) : yield self.get(i)

RooDataSet.__iter__ = _RooDataSetIter


# RooArgSet/RooArgList functions

def _RooArgSetIter(self) :
    z = self.createIterator()
    while True :
        a = z.Next()
        if not a : return
        yield a

RooArgList.__iter__     = _RooArgSetIter
RooArgList.__len__      = lambda s    : s.getSize()
RooArgList.__contains__ = lambda s, i : s.contains(i)
RooArgList.nameList     = lambda s    : [j.GetName() for j in s] 
RooArgList.names        = lambda s    : ','.join(s.nameList())

RooArgSet.__iter__      = _RooArgSetIter
RooArgSet.__len__       = lambda s    : s.getSize()
RooArgSet.__contains__  = lambda s, i : s.contains(i)
RooArgSet.nameList      = lambda s    : [j.GetName() for j in s]
RooArgSet.names         = lambda s    : ','.join(s.nameList())
#RooArgSet.__iadd__      = lambda s,y  : s.add(y)
#RooArgSet.__radd__      = lambda y,s  : s.add(y)
#RooArgSet.__add__       = lambda x,y  : RooArgSet(x).__iadd__(y)


# RooWorkspace functions

def _RooWorkspacePut(self, x) :
  _import = getattr(RooWorkspace, 'import')

  if _import(self, x, RooFit.Silence()) : return None
  return x.GetName()

def _setConstant(ws, pattern, constant = True, value = None):
  import re

  nrexp = re.compile(pattern)

  rc = 0
  for arg in ws.allVars() :
    if not nrexp.match(arg.GetName()) : continue

    arg.setConstant(constant)

    if constant and value :
      if value < arg.getMin() : arg.setMin(value) 
      if value > arg.getMax() : arg.setMax(value) 
      arg.setVal(value) 

    rc += 1

  return rc

RooWorkspace.__getitem__  = lambda s, i    : s.obj(i)
RooWorkspace.__contains__ = lambda s, i    : True if s.obj(i) else False
#RooWorkspace.__setitem__  = lambda s, k, v : s.put('%s[%s]'%(k,v))
RooWorkspace.put          = _RooWorkspacePut
RooWorkspace.setConstant  = _setConstant

