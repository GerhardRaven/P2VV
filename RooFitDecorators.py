from ROOT import RooArgSet

def _RooArgSetIter(self) :
    z = self.createIterator()
    while True :
        a = z.Next()
        if not a : return
        yield a

RooArgSet.__iter__ = _RooArgSetIter
RooArgSet.__len__  = lambda s   : s.getSize()
#RooArgSet.__iadd__ = lambda s,y : s.add(y)
#RooArgSet.__radd__ = lambda y,s : s.add(y)
#RooArgSet.__add__  = lambda x,y : RooArgSet(x).__iadd__(y)

from ROOT import RooWorkspace
RooWorkspace.put = getattr(RooWorkspace,'import')
RooWorkspace.__getitem__ = lambda s,i : s.obj(i)
RooWorkspace.__contains__ = lambda s,i : s.obj(i) is not None
#RooWorkspace.__setitem__ = lambda s,k,v : s.put('%s[%s]'%(k,v))

if False :
    from ROOT import RooWorkspace
    w = RooWorkspace('w')
    w.factory('{x[-1,1],y[0,1],z[-10,10]}')
    z = w.argSet('x,y')
    #for i in z : i.Print()
    #z += w.argSet('z')
    #for i in z : i.Print()


