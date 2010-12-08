from ROOT import RooArgSet, RooArgList

def _RooArgSetIter(self) :
    z = self.createIterator()
    while True :
        a = z.Next()
        if not a : return
        yield a

RooArgList.__iter__ = _RooArgSetIter
RooArgList.__len__  = lambda s   : s.getSize()
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

def setConstant(ws, pattern, constant = True, value = 9999.):
    from ROOT import TRegexp, TString
    rc = int(0)
    rexp = TRegexp(pattern,True)
    for arg in ws.allVars() :
        if TString(arg.GetName()).Index(rexp)>=0 :
            arg.setConstant( constant )
            if constant and value != 9999 :
                if value < arg.getMin() : 
                    arg.setMin(value) 
                elif value > arg.getMax() :
                    arg.setMax(value) 
                arg.setVal(value) 
            rc += 1
    #print 'number of parameters matching ', pattern, rc
    return rc



if False :
    from ROOT import RooWorkspace
    w = RooWorkspace('w')
    w.factory('{x[-1,1],y[0,1],z[-10,10]}')
    z = w.argSet('x,y')
    #for i in z : i.Print()
    #z += w.argSet('z')
    #for i in z : i.Print()


