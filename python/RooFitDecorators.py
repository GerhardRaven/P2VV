from ROOT import (RooArgSet, RooArgList, RooDataSet,
                  RooWorkspace, RooFitResult, RooFit,
                  RooDataHist, RooLinkedList, RooCmdArg)

import ROOTDecorators

def __wrap_kw_subs( fun ) :
    from ROOT import RooFit, RooAbsCollection, TObject
    __doNotConvert = [ RooAbsCollection, TObject ]
    __tbl  = lambda k : getattr(RooFit, k)
    # TODO: anything relying on _target_ or _var should move to RooFitWrappers...
    __dref = lambda o : o._target_() if hasattr(o,'_target_') else o
    def __disp(k, v) :
        if type(v) == type(None) :
            return __tbl(k)()
        elif any( isinstance( v, t ) for t in __doNotConvert ) or not hasattr( v,'__iter__' ) or str(type(v)).find('.Category') != -1 :
            return __tbl(k)( __dref(v) )
        else :
            return __tbl(k)(*v)

    from functools import wraps
    @wraps(fun)
    def _fun(self,*args,**kwargs) :
        args += tuple( RooCmdArg( __disp('Slice',  s) ) for s in kwargs.pop('Slices',  []) )
        args += tuple( RooCmdArg( __disp('Cut',    s) ) for s in kwargs.pop('Cuts',    []) )
        if 'Imports' in kwargs :
            assert 'Index' in kwargs, 'RooFit keyword wrapper: "Imports" argument found without "Index"'
            args += ( RooCmdArg( __disp( 'Index', kwargs.pop('Index') ) ), )  # note order: "Index" before "Import"
            args += tuple( RooCmdArg( __disp('Import', i) ) for i in kwargs.pop('Imports') )

        if 'Minos' in kwargs and hasattr( kwargs['Minos'], '__iter__' ) :
            if kwargs['Minos'] :
                from ROOT import RooArgSet
                RooPars = RooArgSet()
                for par in kwargs['Minos'] : RooPars.add( __dref(par) )
                kwargs['Minos'] = RooPars
            else :
                kwargs.pop('Minos')

        # TODO: this 'if' block needs to move to RooFitWrappers...
        if 'ArgSet' in kwargs or 'ArgList' in kwargs:
            if 'ArgSet' in kwargs :
                from ROOT import RooArgSet
                P2VVVars = kwargs.pop('ArgSet')
                RooVars = RooArgSet()
            else :
                from ROOT import RooArgList
                P2VVVars = kwargs.pop('ArgList')
                RooVars = RooArgList()
            for var in P2VVVars : RooVars.add( __dref(var) ) 
            args += (RooVars,)

        dispatch = ( (k,kwargs.pop(k)) for k in kwargs.keys() if hasattr(RooFit,k) )
        args += tuple(RooCmdArg(__disp(k,v)) for k,v in dispatch )
        try:
            return fun(self, *args, **kwargs)
        except TypeError:
            fun_args = [ a for a in args if not isinstance(a, RooCmdArg) ]
            otr_args = [ a for a in args if     isinstance(a, RooCmdArg) ]
            l = RooLinkedList()
            for a in otr_args: l += a
            return fun(self, *tuple(fun_args + [l]), **kwargs)
    return _fun

def __convert_init_kw_to_setter( cl ) :
    init = cl.__init__
    from functools import wraps
    @wraps(init)
    def _init(self,*args,**kwargs) :
        calls = tuple()
        for i in kwargs.keys() :
            for m in ['Set','set','Add','add'] :
                if hasattr(self, m+i) :
                    calls += (methodcaller( m+i, kwargs.pop(i) ),)
                    break
        init(self,*args,**kwargs)
        for c in calls : c(self)
    cl.__init__ =  _init

from ROOT import TLatex
__convert_init_kw_to_setter( TLatex )
##########################################

from ROOT import RooPlot
def __wrap_RooPlotDraw( __draw ) :
    from functools import wraps
    @wraps(__draw)
    def _RooPlotDraw(self,*args,**kw) :
       pad = kw.pop('pad',None) 
       if kw : raise RuntimeError('unknown keyword arguments: %s' % kw )
       if pad : pad.cd()
       __draw(self,*args)
       if pad : pad.Update()
    return _RooPlotDraw
RooPlot.Draw = __wrap_RooPlotDraw( RooPlot.Draw )

# RooDataSet functions
def __RooDataSetIter(self) :
    for i in range( self.numEntries() ) : yield self.get(i)

from ROOT import RooDataSet
RooDataSet.__iter__ = __RooDataSetIter

def __wrapRooDataSetInit( init ) :
    from ROOT import TObject,TTree,RooDataSet
    __doNotConvert = [ TTree, RooDataSet ]
    def cnvrt(i) :
        if not hasattr(i,'__iter__') or any( isinstance(i, t) for t in __doNotConvert ): return i
        _i = RooArgSet()
        for j in i : 
            from ROOT import RooAbsArg
            if not isinstance(j,RooAbsArg) : return i
            _i.add( j )
        return _i
    from functools import wraps
    @wraps(init)
    def _init( self, *args ) :
        return init(self,*tuple( cnvrt(i) for i in args ))
    return _init
RooDataSet.__init__ = __wrapRooDataSetInit( RooDataSet.__init__ )

def __createRooIterator( create_iterator ) :
    def __iter(self) :
        i = create_iterator(self)
        while True :
            obj = i.Next()
            if not obj : break
            yield obj
    return __iter

def __RooDataSetToTree( self, Name = '', Title = '', BranchList = '', RooFitFormat = True ) :
    from ROOT import RooDataSetToTree
    return RooDataSetToTree( self, Name, Title, BranchList, RooFitFormat )
RooDataSet.buildTree = __RooDataSetToTree

# RooAbsCategory functions
from operator import methodcaller
from ROOT import RooAbsCategory
RooAbsCategory.__iter__ = __createRooIterator( methodcaller('typeIterator') )

from ROOT import RooLinkedList
RooLinkedList.__iadd__ = lambda s,x : s if s.Add(x)    else s  # else None??
RooLinkedList.__isub__ = lambda s,x : s if s.Remove(x) else s  # else None??
# RooAbsCollection/RooArgSet/RooArgList functions
from ROOT import RooAbsCollection
RooAbsCollection.__iter__ = __createRooIterator( methodcaller('createIterator') )
RooAbsCollection.__len__  = lambda s   : s.getSize()
RooAbsCollection.__contains__  = lambda s,i : s.contains(i)
RooAbsCollection.__iadd__ = lambda s,x : s if s.add(x)    else s  # else None??
RooAbsCollection.__isub__ = lambda s,x : s if s.remove(x) else s  # else None??
RooAbsCollection.nameList = lambda s : [ j.GetName() for j in s ] 
RooAbsCollection.names    = lambda s : ','.join( s.nameList() )
RooAbsCollection.__eq__   = lambda s,x : s.equals(x)
RooAbsCollection.__ne__   = lambda s,x : not s.equals(x)
RooAbsCollection.printLatex = __wrap_kw_subs( RooAbsCollection.printLatex )

def __create_RooAbsCollectionInit(t) :
    def cnvrt(i) :
        from ROOT import TObject
        if str(type(i)).find('.Category') != -1:
            return i._target_()
        elif not hasattr(i, '__iter__') and hasattr(i, '_target_'):
            return i._target_()
        elif not hasattr(i, '__iter__') or isinstance(i, TObject):
            return i
        _i = t()
        for j in i : 
            from ROOT import RooAbsArg
            if not isinstance(j,RooAbsArg):
                print "Not a RooAbsArg"
                return i
            _i.add( j._target_() if hasattr(j, '_target_') else j )
        return _i
    __init = t.__init__
    return lambda self,*args : __init(self, *tuple(cnvrt(i) for i in args))

def _RooTypedUnary2Binary( t,op ) :
    return lambda x,y : getattr(t,op)(t(x),y)

from ROOT import RooArgSet, RooArgList
for t in [ RooArgSet,RooArgList ] :
    t.__init__  = __create_RooAbsCollectionInit(t)
    t.__sub__   = _RooTypedUnary2Binary( t, '__isub__' )
    t.__add__   = _RooTypedUnary2Binary( t, '__iadd__' )

# RooWorkspace functions

from ROOT import RooWorkspace, RooFit
RooWorkspace.__getitem__ = lambda s,i : s.obj(i)
RooWorkspace.__contains__ = lambda s,i : bool( s.obj(i) )
#RooWorkspace.__setitem__ = lambda s,k,v : s.put('%s[%s]'%(k,v))

def __RooWorkspacePut( self ,x, **kwargs ) :
    __dref__ = lambda i : i._var if hasattr(i,'_var') else i
    _import = getattr(RooWorkspace,'import')
    if _import(self,__dref__(x),**kwargs) : return None
    return self[kwargs.get('Rename',x.GetName())]

setattr( RooWorkspace, 'import',  __wrap_kw_subs( getattr(RooWorkspace, 'import' ) ) )
RooWorkspace.put = __RooWorkspacePut

def __RooWorkspaceSetConstant(self, pattern, constant = True, value = None):
    import re
    nrexp = re.compile(pattern)
    rc = 0
    for arg in self.allVars() :
        if not nrexp.match(arg.GetName()) : continue
        arg.setConstant( constant )
        if constant and value :
            if value < arg.getMin() : arg.setMin(value) 
            if value > arg.getMax() : arg.setMax(value) 
            arg.setVal(value) 
        rc += 1
    return rc

RooWorkspace.setConstant = __RooWorkspaceSetConstant


# RooFitResult functions

def __RooFitResultResult(self, parList) :
  # get parameter indices in fit result
  indices = {}
  floatPars = self.floatParsFinal()
  for par in parList :
    index = floatPars.index(par)
    if index >= 0 :
      indices[par] = index
    else :
      print 'ERROR: RooFitResult::result(): fit result does not contain parameter', par
      return None
  covMatrix = self.covarianceMatrix()
  values = tuple([floatPars[indices[par]].getVal() for par in parList])
  covariances = tuple([tuple([covMatrix[indices[row]][indices[col]]\
      for col in parList]) for row in parList])

  return (tuple(parList), values, covariances)

RooFitResult.result = __RooFitResultResult


from ROOT import RooFit
from ROOT import RooSimultaneous
RooSimultaneous.plotOn     = __wrap_kw_subs( RooSimultaneous.plotOn )
from ROOT import RooAbsPdf
RooAbsPdf.generate         = __wrap_kw_subs( RooAbsPdf.generate )
RooAbsPdf.fitTo            = __wrap_kw_subs( RooAbsPdf.fitTo )
RooAbsPdf.plotOn           = __wrap_kw_subs( RooAbsPdf.plotOn )
RooAbsPdf.paramOn          = __wrap_kw_subs( RooAbsPdf.paramOn )
RooAbsPdf.createCdf        = __wrap_kw_subs( RooAbsPdf.createCdf )
RooAbsPdf.createNLL        = __wrap_kw_subs( RooAbsPdf.createNLL )
RooAbsPdf.createProjection = __wrap_kw_subs( RooAbsPdf.createProjection )
from ROOT import RooAbsData
RooAbsData.createHistogram = __wrap_kw_subs( RooAbsData.createHistogram )
RooAbsData.reduce          = __wrap_kw_subs( RooAbsData.reduce )
RooAbsData.plotOn          = __wrap_kw_subs( RooAbsData.plotOn )
from ROOT import RooAbsReal
RooAbsReal.plotOn          = __wrap_kw_subs( RooAbsReal.plotOn )
RooAbsReal.fillHistogram   = __wrap_kw_subs( RooAbsReal.fillHistogram )
RooAbsReal.createIntegral  = __wrap_kw_subs( RooAbsReal.createIntegral )
from ROOT import RooRealVar
RooRealVar.format          = __wrap_kw_subs( RooRealVar.format )
from ROOT import RooAbsCollection
RooAbsCollection.printLatex = __wrap_kw_subs( RooAbsCollection.printLatex )
from ROOT import RooMCStudy
RooMCStudy.plotPull = __wrap_kw_subs( RooMCStudy.plotPull)
RooMCStudy.plotError = __wrap_kw_subs( RooMCStudy.plotError)
RooMCStudy.plotNLL = __wrap_kw_subs( RooMCStudy.plotNLL)
RooMCStudy.plotParam = __wrap_kw_subs( RooMCStudy.plotParam)
RooMCStudy.plotParamOn = __wrap_kw_subs( RooMCStudy.plotParamOn)
from ROOT import RooDataSet
RooDataSet.plotOnXY = __wrap_kw_subs( RooDataSet.plotOnXY )

#from ROOT import RooSimCloneTool
#RooSimCloneTool.build = __wrap_kw_subs(RooSimCloneTool.build )
#from ROOT import RooDataHist
from ROOT import RooDataSet, RooChi2Var, RooProdPdf, RooMCStudy
for i in  [ RooDataSet, RooChi2Var, RooProdPdf, RooMCStudy ] :
    i.__init__ = __wrap_kw_subs( i.__init__ )


from ROOT import RooAbsRealLValue
def __convert_kw_to_setter( fun, ret ) :
    from functools import wraps
    @wraps(fun)
    def _fun(self,*args,**kwargs) :
        calls = tuple()
        for i in kwargs.keys() :
            if hasattr(RooFit,i) : continue # leave kw instance as is...
            for m in ['Set','set','Add','add'] :
                if hasattr(ret, m+i) : 
                    calls += (methodcaller( m+i, *kwargs.pop(i) ),)
                    break
        o = fun(self,*args,**kwargs)
        for c in calls : c(o)
        return o
    return _fun

RooAbsRealLValue.frame     = __convert_kw_to_setter( RooAbsRealLValue.frame, RooPlot )
RooAbsRealLValue.frame     = __wrap_kw_subs( RooAbsRealLValue.frame )

from ROOT import RooMsgService
RooMsgService.addStream = __wrap_kw_subs( RooMsgService.addStream )

def __RooMsgService__iter__(self) :
    for i in range(self.numStreams()) : yield self.getStream(i)
RooMsgService.__iter__ = __RooMsgService__iter__

def __wrap_streamconfig_topic_add_sub( fun ) :
    # why are addTopic and removeTopic void -- would be nicer if they were StreamConfig& addTopic( .. ) { ...; return *this; }
    from functools import wraps
    @wraps(fun)
    def _fun( self, arg ) :
        if hasattr(arg,'__iter__') :
            for i in arg : fun(self,i)
        else :
            fun(self,arg)
        return self
    return _fun
RooMsgService.StreamConfig.__iadd__ = __wrap_streamconfig_topic_add_sub( RooMsgService.StreamConfig.addTopic )
RooMsgService.StreamConfig.__isub__ = __wrap_streamconfig_topic_add_sub( RooMsgService.StreamConfig.removeTopic )


def _RooFitResultCorrMatrixLatex(self, name):
    parlist = self.floatParsFinal()
    npar = parlist.getSize()
    corr = []
    for i in range(npar):
        corr.append(self.correlation(parlist[i]))

    layoutstring = '|c'*(npar+1)
    layoutstring += '|}'

    string = '\\documentclass{article}\n'
    string += '\\begin{document}\n'
    string += '\\begin{table}[h]\n'
    string += '\\begin{center}\n'
    string += '\\begin{tabular}{'
    string += layoutstring+'\n' 
    string += '\\hline\n'
    string += 'Parameter '
    for i in range(npar):
        string += ' & %s'%(parlist[i].GetName())
    string += ' \\\\\n'
    string += '\\hline\n\\hline\n'


    for j in range(npar):
        string += parlist[j].GetName() 
        for i in range(npar):
            if i>=j:
                if abs(corr[j][i].getVal()) < 0.005:
                    string += ' & - '
                else:
                    string += ' & ' + str(round(corr[j][i].getVal(),3)) 
            else : 
                string += ' & '
        string += ' \\\\\n'
    string += '\\hline\n'
    string += '\\end{tabular}\n'
    string += '\\end{center}\n'
    string += '\\end{table}\n'
    string += '\\end{document}\n'

    f = open('corrtable_%s.tex'%name, 'w')
    f.write(string)
    f.close

    return

RooFitResult.writecorr = _RooFitResultCorrMatrixLatex

def _RooFitResultParamsLatex(self,name,toys):
    f = open('params_%s.tex'%name,'w')

    string = '\\documentclass{article}\n'
    string += '\\begin{document}\n'

    def __s(n, fmt = '%.3e'):
        if type(n) == float:
            n = fmt % n
        if n.endswith('e+00'):
            return n[:-4]
        else:
            return n

    offset = 0
    for i in self.floatParsFinal():
        name = i.GetName().replace('_', '\_')
        if len(name) > offset:
            offset = len(name)
    fmt = '{0:<' + str(offset) + '} & '
        
    if toys:
        string += '\\begin{tabular}{|c|c|c|c|}\n'
        string += '\\hline\n'
        string += 'parameter & result & original value & $\sigma$ from original \\\\ \n'
        string += '\\hline\n'
        string += '\\hline\n'
        for i,j in zip(self.floatParsFinal(),self.floatParsInit()):
            print i,j
            string += fmt.format(i.GetName().replace('_', '\_'))
            string += '%s $\pm$ %s & ' % (__s(i.getVal()),__s(i.getError()))
            string += '%s & ' % __s(j.getVal())
            string +=  '%s \\\\ \n' % __s((j.getVal() - i.getVal()) / i.getError())
    else:
        string += '\\begin{tabular}{|c|c|}\n'
        string += '\\hline\n'
        string += 'parameter & result \\\\ \n'
        string += '\\hline\n'
        string += '\\hline\n'
        for i,j in zip(self.floatParsFinal(),self.floatParsInit()):
            print i,j
            string += fmt.format(i.GetName().replace('_', '\_'))
            string += '%s $\pm$ %s \\\\ \n' % (__s(i.getVal()), __s(i.getError()))

    string += '\\hline\n'
    string += '\\end{tabular}\n'
    string += '\\end{document}\n'

    f.write(string)
    f.close()

    return

RooFitResult.writepars = _RooFitResultParamsLatex

def _RooFitResultPrint( self, **kwargs ) :
    from math import ceil, log10
    fitPars = self.floatParsFinal()

    text   = kwargs.pop( 'text',   False )
    LaTeX  = kwargs.pop( 'LaTeX',  False )
    normal = kwargs.pop( 'normal', False )
    if not text and not LaTeX : normal = True

    printDecim  = kwargs.pop( 'Precision', 2   )
    parNameDict = kwargs.pop( 'ParNames',  { } )

    print '  fit result:'
    print '  status: %d' % self.status()
    print '  NLL = %f  :  EDM = %f' % ( self.minNll(), self.edm() )
    print

    if text :
        for par in fitPars :
            name = parNameDict[ par.GetName() ][0] if par.GetName() in parNameDict else par.GetName()
            if par.hasAsymError() :
                prec = printDecim - int( ceil( log10( max( -par.getErrorLo(), par.getErrorHi() ) ) ) )
                if prec < 0 : prec = 0
                print ( '  {0:<30s} {1:<+10.%df} {2:<+10.%df} {3:<+10.%df}' % ( prec, prec, prec ) )\
                      .format( name, par.getVal(), par.getErrorHi(), par.getErrorLo() )

            else :
                prec = printDecim - int( ceil( log10( par.getError() ) ) )
                if prec < 0 : prec = 0
                print ( '  {0:<30s} {1:<+10.%df} +/- {2:<10.%df}' % ( prec, prec ) ).format( name, par.getVal(), par.getError() )
        print

    if LaTeX :
        print '  \\begin{tabular}{|c|c|}'
        print '    \\hline'
        print '    parameter  &  value  \\\\'
        print '    \\hline'
        for par in fitPars :
            name = parNameDict[ par.GetName() ][1] if par.GetName() in parNameDict else par.GetName()
            if par.hasAsymError() :
                prec = printDecim - int( ceil( log10( max( -par.getErrorLo(), par.getErrorHi() ) ) ) )
                if prec < 0 : prec = 0
                print ( '    {2:<40s}  &  ${3:<+10.%df}^{0}{4:<+10.%df}{1}_{0}{5:<+10.%df}{1}$  \\\\' % ( prec, prec, prec ) )\
                      .format( '{', '}', name, par.getVal(), par.getErrorHi(), par.getErrorLo() )

            else :
                prec = printDecim - int( ceil( log10( par.getError() ) ) )
                if prec < 0 : prec = 0
                print ( '    {0:<40s}  &  ${1:<+10.%df} \\pm {2:<10.%df}$             \\\\' % ( prec, prec ) )\
                      .format( name, par.getVal(), par.getError() )
        print '    \\hline'
        print '  \\end{tabular}'
        print

    if normal :
        self.Print()

RooFitResult.PrintSpecial = _RooFitResultPrint
