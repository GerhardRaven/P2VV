from ROOT import (RooArgSet, RooArgList, RooDataSet,
                  RooWorkspace, RooFitResult, RooFit,
                  RooDataHist, RooLinkedList, RooCmdArg)
from ROOT import gStyle,gROOT
import ROOTDecorators
gStyle.SetPalette(1)
gROOT.SetStyle("Plain")

# needed to get RooFit.Name, RooFit.Components.... kRed
# how to get just the RooFit namespace ?
#from ROOT import * 

def __wrap_kw_subs( fun ) :
    from ROOT import RooFit, RooAbsCollection, TObject
    __doNotConvert = [ RooAbsCollection, TObject ]
    __tbl  = lambda k : getattr(RooFit, k)
    def __disp(k, v):
        if any( isinstance( v, t ) for t in __doNotConvert ) or not hasattr( v,'__iter__' ) \
           or str(type(v)).find('.Category') != -1:
            return __tbl(k)(v._target_() if hasattr(v, '_target_') else v)
        elif type(v) != type(None):
            return __tbl(k)(*v)
        else:
            __tbl(k)()

    from functools import wraps
    @wraps(fun)
    def _fun(self,*args,**kwargs) :
        if 'Slices' in kwargs:
            args += tuple( RooCmdArg( __disp('Slice', slice) ) for slice in kwargs.pop('Slices') )
        if 'Imports' in kwargs:
            args += tuple( RooCmdArg( __disp('Import', slice) ) for slice in kwargs.pop('Imports') )
        for k,v in kwargs.iteritems():
            args += (RooCmdArg( __disp(k,v)),)
        try:
            return fun(self, *args)
        except TypeError:
            l = RooLinkedList()
            fun_args = [a for a in args if not isinstance(a, RooCmdArg)]
            for a in args:
                if isinstance(a, RooCmdArg):
                    l.Add(a)
            return fun(self, *tuple(fun_args + [l]))
    return _fun

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
            if not obj : raise StopIteration
            yield obj
    return __iter

# RooAbsCategory functions
from operator import methodcaller
from ROOT import RooAbsCategory
RooAbsCategory.__iter__ = __createRooIterator( methodcaller('typeIterator') )

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
    _import = getattr(RooWorkspace,'import')
    if _import(self,x,**kwargs) : return None
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
from ROOT import RooAbsData
RooAbsData.createHistogram = __wrap_kw_subs( RooAbsData.createHistogram )
RooAbsData.reduce          = __wrap_kw_subs( RooAbsData.reduce )
RooAbsData.plotOn          = __wrap_kw_subs( RooAbsData.plotOn )
from ROOT import RooAbsReal
RooAbsReal.plotOn          = __wrap_kw_subs( RooAbsReal.plotOn )
RooAbsReal.fillHistogram   = __wrap_kw_subs( RooAbsReal.fillHistogram )
RooAbsReal.createIntegral  = __wrap_kw_subs( RooAbsReal.createIntegral )
from ROOT import RooAbsRealLValue
RooAbsRealLValue.frame     = __wrap_kw_subs( RooAbsRealLValue.frame )
from ROOT import RooRealVar
RooRealVar.format          = __wrap_kw_subs( RooRealVar.format )
from ROOT import RooAbsCollection
RooAbsCollection.printLatex = __wrap_kw_subs( RooAbsCollection.printLatex )
#from ROOT import RooSimCloneTool
#RooSimCloneTool.build = __wrap_kw_subs(RooSimCloneTool.build )
#from ROOT import RooDataHist
from ROOT import RooDataSet, RooChi2Var, RooProdPdf, RooMCStudy
for i in  [ RooDataSet, RooChi2Var, RooProdPdf, RooMCStudy ] :
    i.__init__ = __wrap_kw_subs( i.__init__ )

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

    if toys:
        string += '\\begin{tabular}{|c|c|c|c|}\n'
        string += '\\hline\n'
        string += 'parameter & result & original value & $\sigma$ from original \\\\ \n'
        string += '\\hline\n'
        string += '\\hline\n'
        for i,j in zip(self.floatParsFinal(),self.floatParsInit()):
            print i,j
            string += '%s & '%i.GetName()
            string += '%s $\pm$ %s & '%(round(i.getVal(),4),round(i.getError(),4))
            string += '%s & '%round(j.getVal(),4)
            string +=  '%s \\\\ \n'%(round((j.getVal()-i.getVal())/i.getError(),4))
    else:
        string += '\\begin{tabular}{|c|c|}\n'
        string += '\\hline\n'
        string += 'parameter & result \\\\ \n'
        string += '\\hline\n'
        string += '\\hline\n'
        for i,j in zip(self.floatParsFinal(),self.floatParsInit()):
            print i,j
            string += '%s & '%i.GetName()
            string += '%s $\pm$ %s \\\\ \n'%(round(i.getVal(),4),round(i.getError(),4))

    string += '\\hline\n'
    string += '\\end{tabular}\n'
    string += '\\end{document}\n'

    f.write(string)
    f.close()

    return

RooFitResult.writepars = _RooFitResultParamsLatex
