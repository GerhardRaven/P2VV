from RooFitDecorators import *
from functools import wraps

def __check_req_kw__( name, kwargs ) :
    if not name in kwargs : raise KeyError( 'Must provide kw argument %s' % name )
def __check_exists_already__( self ) :
    if self._name in self.ws() :
        raise StandardError( 'Recreating %s is not supported atm' % type(self) )

__dref__ = lambda i : i._var if hasattr(i,'_var') else i

def __wrap__dref_var__( fun ) :
    @wraps(fun)
    def _fun(self,*args) :
        return fun(self, *tuple( __dref__(i) for i in args ) )
    return _fun

RooAbsCollection.__contains__ = __wrap__dref_var__( RooAbsCollection.__contains__ )
RooAbsData.table = __wrap__dref_var__( RooAbsData.table )


class RooObject(object) :
    _ws = None
    _setters = {'Title'      : lambda s,v : s.SetTitle(v)
               ,'Observable' : lambda s,v : s.setObservable(v) 
               }
    _getters = {'Name'       : lambda s : s.GetName()
               ,'Title'      : lambda s : s.GetTitle()
               ,'Value'      : lambda s : s.getVal()
               ,'Type'       : lambda s : s.Type()
               ,'Observable' : lambda s : s.observable() 
               ,'Constant'   : lambda s : s.isConstant() 
               }
    def _factory(self,spec) :
        return self.ws().factory(spec)
        
    def __setitem__(self,k,v):
        from itertools import ifilter, imap
        for setters in imap( lambda x: x._setters, ifilter( lambda x : hasattr(x,'_setters'), type(self).__mro__) ): 
            if k in setters :  return setters[k](self,v )
        raise KeyError('\'%s\' is not known for class %s' % (k, type(self) ) )
    def __getitem__(self,k):
        from itertools import ifilter, imap
        for getters in imap( lambda x: x._getters, ifilter( lambda x : hasattr(x,'_getters'), type(self).__mro__) ): 
            if k in getters :  return getters[k](self)
        raise KeyError('\'%s\' is not known for class %s' % (k, type(self) ) )


    def __init__(self,workspace = None) :
        if workspace :
            from ROOT import RooWorkspace
            if type(workspace) != RooWorkspace : workspace = RooWorkspace(workspace)
            self.setWorkspace(workspace)

    def ws(self) : 
        if not RooObject._ws : raise RuntimeError('No workspace defined!')
        return RooObject._ws

    def setWorkspace(self,ws):
        RooObject._ws = ws
        if not hasattr(ws, '_objects') :   ws._objects     = {} # name -> object
        if not hasattr(ws, '_rooobject') : ws._rooobjects  = {} # name -> rooobject
        if not hasattr(ws, '_mappings') :  ws._mappings    = {}
        if not hasattr(ws, '_spec') :      ws._spec        = {} # factory string -> object
        
    def _rooobject(self,Name) :
        if type(Name)!=str : Name = Name.GetName()
        return self.ws()._rooobject[Name]

    def _declare(self,spec):
        """
        create the underlying C++ object in the workspace
        """
        # canonicalize 'spec' a bit by getting rid of spaces
        spec = spec.strip()
        # TODO: Wouter promised to add a method that, given the factory 'spec' above returns 
        #       the value of 'factory_tag' which is used internally in the conflict resolution
        #       and which is the 'canonical' recipe to build an object
        #       That way, we can improve the check for re-use, as currently it _is_ possible
        #       to get collisions, as two different 'spec' may end up trying to build the same object
        #       where 'same' is defined as having the same name.
        #       For now, we deal with this by raising an exception when the factory call encounters 
        #       a conflict.
        if spec not in self.ws()._spec : 
            x = self._factory(spec)
            if not x: raise NameError("workspace factory failed to return an object for factory string '%s' "%spec)
            if hasattr(x,'setStringAttribute') : x.setStringAttribute('RooFitWrappers.RooObject::spec',spec) 
            #  
            # Keep the PyROOT objects in a container so they don't get garbage
            # collected.
            self.ws()._objects[x.GetName()] = x # Note: use explicit GetName, not str, as x is a 'bare' PyROOT object!!!
            # and keep track what we made 
            self.ws()._spec[ spec ] = x
        else :
            x = self.ws()._spec[ spec ] 
            # print 'INFO: spec not unique, returning pre-existing object: %s -> %s' %( spec, x.GetName() )
        return x

    def _init(self,name,type_) :
        """
        match ourselves to the underlying C++ object in the workspace
        This is done by assigning _var 
        """
        # If the object was mapped to something from the dataset, put the right
        # variable behind it.
        if name in self.ws()._mappings:
            name = self.ws()._mappings[name]

        # Get the right object from our own cache, KeyError is raised correctly.
        x = self.ws()._objects[name]
        if not x.InheritsFrom(type_): 
            raise KeyError('%s is %s, not %s' % (name, x.ClassName(), type_))
        self._var = x
        self.ws()._rooobjects[x.GetName()] = self

    def __getattr__(self, name):
        return getattr(self._target_(), name)      

    def _target_(self) :
        return self._var

    def __cmp__(self, other):
        o = other if type(other) == str else other.GetName()
        return self.GetName() < o

    def __eq__(self, other):
        o = other if type(other) == str else other.GetName()
        return self.GetName() == o

    def __ne__(self, other):
        o = other if type(other) == str else other.GetName()
        return self.GetName() != o

    def __hash__(self):
        return self.GetName().__hash__()

    def __str__(self):
        return self.GetName()

    def Observables(self) :
        return set( self.ws()._rooobjects[i.GetName()] for i in self._var.getVariables() if i.getAttribute('Observable') )
    def Parameters(self) :
        return set( self.ws()._rooobjects[i.GetName()] for i in self._var.getVariables() if not i.getAttribute('Observable') )

    ## FIXME: Should these be in RooObject? Do all RooObjects always have a non-empty _dict???
    def Type(self) :
        _t = self._dict['Type']
        return _t if type(_t)==str else _t.__name__

        
    ## FIXME: Should these be in RooObject?? Do we need an LValue wrapper and move these there?
    def observable(self) : 
        return self._var.getAttribute('Observable')
    def setObservable(self, observable) :
        from ROOT import RooAbsLValue
        assert isinstance(self._var,RooAbsLValue) # if we're not an LValue, we cannot be Observable!!!
        assert self._var.isLValue() # not sure the best way to check for LValue...
        self._var.setAttribute('Observable',observable)

    def mappings(self):
        return self.ws()._mappings
# TODO: make this more of a 'borg' by overloading __new__ instead of __init__
#       otherwise properties of the proxy not in the 'target' are not shared
#       across multiple instances which all defer to the same 'target'
# NOTE: we push property values into the underlying target to share them, so this is OK ;-)

# TODO?: Instead of using GetName everywhere, perhaps an appropriate __hash__ function is more
##       elegant?


class ArgSet(RooObject) :
    def _factory(self,spec):
        import re
        match = re.match('set::([^(]+)\(.*\)',spec)
        assert match
        name =  match.groups(1)[0]
        super(ArgSet,self)._factory(spec) # this will return None for us...
        x = self.ws().set( name )
        x.setName( name ) # argsets come out of the ws as brand new clones without name...
        return x
    def __init__(self,name,args) :
        spec = 'set::%s(%s)' % (name, ','.join( i['Name'] for i in args) )
        self._declare(spec)
        self._init(name,'RooArgSet')



class Category (RooObject): 
    _getters = {'Index'      : lambda s : s.getIndex() 
               ,'Label'      : lambda s : s.getLabel()
               ,'States'     : lambda s : s.states()
               } 
    _setters = {'Index'      : lambda s,v : s.setIndex(v)
               ,'Label'      : lambda s,v : s.setLabel(v)
               ,'Constant'   : lambda s,v : s.setConstant(v)
               }

    def __init__(self,Name,**kwargs):
        # construct factory string on the fly...
        __check_req_kw__( 'States', kwargs )
        states = kwargs.pop('States')
        if   type(states) == list:
            states = ','.join(states)
        elif type(states) == dict:
            states = ','.join(['%s=%d' % i for i in states.iteritems()])
        else:
            raise KeyError('%s does not exist yet, bad States %s defined.' % (Name,states))

        # Create the category and extract states into storage
        self._declare("%s[%s]"%(Name,states))
        self._init(Name,'RooCategory')
        self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
            
    def states(self):
        return self._states

class SuperCategory( Category ) :
    def __init__(self,Name,cats,**kwargs):
        self._declare("SuperCategory::%s({%s})"%(Name,','.join( [ c.GetName() for c in cats ] ) ) )
        self._init(Name,'RooSuperCategory')
        self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class MappedCategory( Category ) :
    def __init__(self,Name,cat,mapper,**kwargs):
        if Name not in self.ws():
            # construct factory string on the fly: no factory string for MappedCategory???
            from ROOT import RooMappedCategory
            obj =  RooMappedCategory(Name,Name,__dref__(cat) )
            for k,vs in mapper.iteritems() : 
                for v in vs : obj.map( v, k )
            self.ws()._objects[ Name ] = self.ws().put( obj )
            self._init(Name,'RooMappedCategory')
            self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(Name,'RooMappedCategory')
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                # Skip these to avoid failure in case we were loaded from a
                # DataSet in the mean time
                if k in ['Index', 'Label']: continue
                assert v == self[k]


class Product(RooObject) :
    def __init__(self,Name,fargs,**kwargs) :
        spec =  "prod::%s(%s)"%(Name,','.join(i['Name'] for i in fargs)) 
        self._declare( spec )
        self._init(Name,'RooProduct')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
            

class Addition(RooObject) :
    def __init__(self,Name,fargs,**kwargs) :
        # construct factory string on the fly...
        self._declare("sum::%s(%s)"%(Name,','.join(i.GetName() for i in fargs)) )
        self._init(Name,'RooAddition')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class FormulaVar (RooObject): 
    _getters = {'Formula'    : lambda s : s.formula()
               ,'Dependents' : lambda s : s.dependents()
               ,'Value'      : lambda s : s.getVal()
               }
    def __init__(self, Name, formula, fargs, **kwargs):
        # construct factory string on the fly...
        spec = "expr::%s('%s',{%s})" % (Name, formula, ','.join(i.GetName() for i in fargs)) 
        self._declare(spec)
        self._init(Name, 'RooFormulaVar')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class ConstVar (RooObject): 
    def __init__(self,**kwargs):
         # construct factory string on the fly...
         __check_req_kw__( 'Value', kwargs )
         __check_req_kw__( 'Name', kwargs )
         self._declare("ConstVar::%(Name)s(%(Value)s)" % kwargs )
         (Name,value) = (kwargs.pop('Name'),kwargs.pop('Value'))
         self._init(Name,'RooConstVar')
         for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        
class P2VVAngleBasis (RooObject) : 
    # TODO: move 'c' out of this class (and into an explicit product), 
    #       which will allow more re-use of existing objects, and hence 
    #       make things faster
    def __init__(self, angles, ind,c=1,ind2 = None) :
        assert c!=0
        name =      'p2vvab_%d%d%d%d' % ind
        if ind2 : name += '_%d%d%d%d' % ind2
        if c!=1 : name += '_%3.2f' % c
        name = name.replace('.','d').replace('-','m')
        spec = "RooP2VVAngleBasis::%s(" % name
        # WARNING: angles may contain barebones PyROOT objects!!!
        spec += "%s, %s, %s, " % tuple( angles[i].GetName() for i in [ 'cpsi', 'ctheta', 'phi'] )
        spec += "%d,%d,%d,%d, " % ind
        if ind2 : spec += "%d,%d,%d,%d, " % ind2
        spec += " %f)"% c
        #NOTE: this requires libP2VV.so to be loaded 
        from P2VVLoad import P2VVLibrary
        self._declare( spec )
        self._init(name,'RooP2VVAngleBasis')
 

class AbsRealMoment( object ):
    def __init__( self, moment )  : self._var = moment
    def __getattr__( self, name ) : return getattr(self._var, name)      
    def GetName( self )           : return self.basisFunc().GetName()
 
class RealMoment( AbsRealMoment ):
    def __init__( self, BasisFunc, Norm ) :
        # get arguments
        self._basisFunc = BasisFunc
        self._norm      = Norm
        # create moment
        from P2VVLoad import P2VVLibrary
        from ROOT import RooRealMoment
        AbsRealMoment.__init__( self, RooRealMoment( __dref__(self._basisFunc), self._norm ) )
          
class RealEffMoment( AbsRealMoment ):
    def __init__( self, BasisFunc, Norm, PDF, NormSet ) :
        # get arguments
        self._basisFunc = BasisFunc
        self._norm      = Norm
        self._pdf       = PDF
        self._normSet   = NormSet

        # build a RooFit normalisation set
        from ROOT import RooArgSet
        self._rooNormSet = RooArgSet( __dref__(var) for var in self._normSet )

        # create efficiency moment
        from P2VVLoad import P2VVLibrary
        from ROOT import RooRealEffMoment
        AbsRealMoment.__init__( self, RooRealEffMoment( __dref__(self._basisFunc), self._norm, __dref__(self._pdf), self._rooNormSet ) )


class RealVar (RooObject): 
    # WARNING: multiple instances don't share proxy state at this time...
    # TODO: move common things like Name and Title in RooObject...
    # TODO: provide scaffolding in RooObject to extend getters & setters on a class-by-class basis
    _getters = {'Unit'       : lambda s : s.getUnit() 
               ,'Value'      : lambda s : s.getVal()
               ,'MinMax'     : lambda s : s.getRange()
               ,'nBins'      : lambda s : s.getBins()
               }
    _setters = {'Unit'       : lambda s,v : s.setUnit(v) 
               ,'Value'      : lambda s,v : s.setVal(v)
               ,'MinMax'     : lambda s,v : s.setRange(v)
               ,'Constant'   : lambda s,v : s.setConstant(v)
               ,'nBins'      : lambda s,v : s.setBins(v)
               ,'Ranges'     : lambda s,v : s.setRanges(v)
               }
    # TODO: provide a copy constructor
    def __init__(self,Name ,**kwargs):
        if 'name' in kwargs : raise RuntimeError('Please replace name argument with Name = xyz' )
        if Name not in self.ws():
            # construct factory string on the fly...
            if 'Value'  not in kwargs  and 'MinMax' not in kwargs and 'Formula' not in kwargs:
                raise KeyError('%s does not exist yet, neither Value nor MinMax nor Formula specified'%Name)
            if 'Formula' not in kwargs :
                if 'Value' not in kwargs:
                    (mi,ma) = kwargs.pop('MinMax')
                    self._declare("%s[%s,%s]"%(Name,mi,ma))
                elif 'MinMax' not in kwargs:
                    self._declare("%s[%s]"%(Name,kwargs.pop('Value')))
                else:
                    (mi,ma) = kwargs.pop('MinMax')
                    val = kwargs.pop('Value')
                    if val < mi or val > ma : raise RuntimeError('Specified Value %s not contained in MinMax (%s,%s)' % ( val,mi,ma))
                    self._declare("%s[%s,%s,%s]"%(Name,val,mi,ma))
            else :
                assert 'Value' not in kwargs
                (formula, args, dset) = kwargs.pop('Formula')
                fname = '_%s__formula'%Name
                dummy  = self._declare("expr::%s('%s',{%s})"%(fname,formula,','.join(i.GetName() for i in args) ) )
                self._declare("dataobs::%s(%s,%s)"%(Name,dset.GetName(),dummy.GetName()))

            if 'Blind' in kwargs: # wrap the blinding class around us...
                b = kwargs.pop('Blind')
                _type = b[0] if type(b[0])==str else b[0].__name__ 
                _bs   = b[1]
                _args = b[2:]
                self._declare("%s::%s_blind('%s',%s,%s)"%(_type,Name,_bs,','.join('%s'%i for i in _args),Name))
            self._init(Name,'RooRealVar')
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(Name,'RooRealVar')
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                # Skip these to avoid failure in case we were loaded from a
                # DataSet in the mean time
                if k == 'Value': continue
                assert v == self[k], '\'%s\' is not the same for %s' % ( k, Name )
            
    # overrule RooRealVar.setRange
    @wraps(RooRealVar.setRange)
    def setRange(self, *args):
        assert args
        if type(args[0]) == str :
            assert len(args)==2
            (mi,ma) = args[1]
            if mi == None : mi = self._var.getMin()
            if ma == None : ma = self._var.getMax()
            self._var.setRange( args[0], mi, ma )
        else  :
            assert len(args)==1
            (mi,ma) = args[0]
            if mi == None : mi = self._var.getMin()
            if ma == None : ma = self._var.getMax()
            self._var.setRange(mi,ma)
            if self.getVal() < mi or self.getVal() > ma : self.setVal( 0.5*(ma+mi) )

    def setRanges(self,arg) :
        for k,v in arg.iteritems() : self.setRange(k,v)

    def getRange(self):
        return self._target_().getMin(), self._target_().getMax()

##TODO, factor out common code in Pdf and ResolutionModel

class Pdf(RooObject):
    # TODO: support conditional observables!!!
    #       multiply automatically with flat PDF for conditional observables in no PDF for them is given??
    #       or intercept fitTo and add ConditionalObservables = ( ... ) ??
    _getters = {'Type'        : lambda s : s._get('Type')
               ,'Parameters'  : lambda s : s._get('Parameters')
               ,'Name'        : lambda s : s._get('Name')
               ,'ConditionalObservables' : lambda s : s.ConditionalObservables()
               }
    _setters = { 'ConditionalObservables' : lambda s,v : s.setConditionalObservables(v)
               }

    ## TODO: define operators
    def __init__(self, **kwargs):
        __check_req_kw__( 'Type', kwargs )
        __check_req_kw__( 'Name', kwargs )

        # Save the keyword args as properties
        self._dict = kwargs
        self._make_pdf()

    def __str__(self):
        d = dict([(a, self[a]) for a in Pdf._getters if hasattr(self, a)])
        return '%s' % d

    def _get(self, name):
        return getattr(self._target_(), '_' + name.lower())
    
    def __getitem__(self, k):
        if hasattr(self, '_dict') and self._dict and k in self._dict:
            return self._dict[k]
        else:
            try:
                return RooObject.__getitem__( self, k )
            except AttributeError as error:
                raise KeyError(str(error))

    def _make_pdf(self):
        if self._dict['Name'] not in self.ws():
            v = list(self._dict['Parameters'])
            deps = ','.join([i.GetName() if type(i) != str else i for i in v])
            x = self._declare( '%s::%s(%s)' % (self.Type(), self._dict['Name'], deps) )
            from ROOT import RooAbsPdf
            assert isinstance(x,RooAbsPdf)
            self._init(self._dict['Name'], x.ClassName())

            # Change self._dict into attributes. Cannot be done before since the
            # underlying object does only exists at this point.
            for k, v in self._dict.iteritems():
                # change into a frozenset AFTER _makeRecipe has been invoked
                # as _makeRecipe relies on the a-priori given order so it can
                # match the c'tor arguments properly!!!!
                if k in [ 'Parameters']: v = frozenset(v)
                attr = '_' + k.lower()
                setattr(self._target_(), attr, v)
            del self._dict # no longer needed
        else:
            self._init(self._dict['Name'], 'RooAbsPdf')
            # Make sure we are the same as last time
            for k, v in self._dict.iteritems():
                if k in [ 'Parameters']: v = frozenset(v)
                if v != self._get(k) : print k,v,self._get(k)
                assert v == self._get(k)

    def _separator(self):
        return '_'

    def ConditionalObservables(self) :
        if not hasattr(self,'_conditionals') : return set()
        return set( i for i in self.Observables() if i.GetName() in self._conditionals )
    def setConditionalObservables(self, obs ) :
        # TODO: check if we have a conditional request for something which isn't one of 
        #       our observables and provide a warning...
        self._conditionals = set( o if type(o)==str else o.GetName() for o in obs )
        

    @wraps(RooAbsPdf.fitTo)
    def fitTo( self, data, **kwargs ) :
        if 'ConditionalObservables' in kwargs :
            kwargs['ConditionalObservables'] = RooArgSet( __dref__(var) for var in kwargs.pop('ConditionalObservables') )
        return self._var.fitTo( data, **kwargs )

    @wraps(RooAbsPdf.generate)
    def generate(self, whatvars, *args, **kwargs):
        #if not whatvars : whatvars = [ i for i in self._var.getVariables() if i.getAttribute('Observable') ]
        return self._var.generate(RooArgSet([ __dref__(i) for i in whatvars] if not isinstance(__dref__(whatvars),RooAbsCategory) else __dref__(whatvars)), *args,**kwargs)

    @wraps(RooAbsPdf.plotOn)
    def plotOn( self, frame, **kwargs ) :
        if 'Slice' in kwargs :
            sl = kwargs.pop('Slice')
            kwargs['Slice'] = ( __dref__(sl[0]), sl[1] )
        if 'Slices' in kwargs :
            sls = kwargs.pop('Slices')
            kwargs['Slices'] = [ ( __dref__(sl[0]), sl[1] ) for sl in sls ]
        return self._var.plotOn( frame, **kwargs )

    def __rmul__(self, eff):
        if not isinstance(eff, FormulaVar) and not isinstance(eff, HistFunc):
            raise RuntimeError, 'trying to multiply a %s with %s; this is not supported!' % (type(self), type(eff))
        name = eff.GetName() + '_X_' + self['Name']
        return EffProd(name, Original = self, Efficiency = eff)

class ProdPdf(Pdf):
    def __init__(self, Name, PDFs, **kwargs):
        self._dict = { 'PDFs' : frozenset(PDFs)
                     , 'Name' : Name + '_' + self._separator().join([i.GetName() for i in PDFs])
                     }
        self._make_pdf()
        del self._dict
        
    def _make_pdf(self):
        if self._dict['Name'] not in self.ws():
            self._declare(self._makeRecipe())
            self._init(self._dict['Name'], 'RooProdPdf')

            # Change self._dict into attributes. Cannot be done before since the
            # underlying object does only exists at this point.
            for k, v in self._dict.iteritems():
                attr = '_' + k.lower()
                setattr(self._target_(), attr, v)
        else:
            self._init(self._dict['Name'], 'RooProdPdf')
            # Make sure we are the same as last time
            for k, v in self._dict.iteritems():
                assert v == self._get(k)

    def _makeRecipe(self):
        def _handleConditional( pdf ) :
            name = pdf.GetName()
            cond = pdf.ConditionalObservables()
            if cond : 
                # first define a named set for our conditionals, as the
                # parsing of this is utterly borken in case we try PROD::name( f | { x,y }, g )
                # as it interprets this as f|x ,y, g due to not recognizing the closing } as
                # the token being parsed is already split on the ',', hence asSET gets "{x" and
                # is very happy with that... and the closing } silenty disappears as well...
                __borken_parser_workaround = ArgSet( name+'_conditional_obs', cond )
                name += '|%s'% __borken_parser_workaround.GetName()
            return name
        pdfs = ','.join( _handleConditional(p) for p in self._dict['PDFs'])
        return 'PROD::%s(%s)' % (self._dict['Name'], pdfs)

    def _separator(self):
        return '_X_'

class SumPdf(Pdf):
    def __init__(self, Name, PDFs, Yields):
        self._yields = {}
        self._dict = {'Name'  : Name,
                      'Yields': Yields,
                      'PDFs'  : PDFs,
                      'Type'  : 'RooAddPdf'}
        pdfs = list(PDFs)
        diff = set([p.GetName() for p in pdfs]).symmetric_difference(set(Yields.keys()))
        if len(diff) not in [0, 1]:
            raise StandardError('The number of yield variables must be equal to or 1'
                                + 'less then the number of PDFs.')
        # self._dict['Name'] = self._separator().join([p.GetName() for p in pdfs])
        self._make_pdf()
        del self._dict
        
    def _make_pdf(self):
        if self._dict['Name'] not in self.ws():
            self._declare(self._makeRecipe())
            self._init(self._dict['Name'], self.Type())

            # Change self._dict into attributes. Cannot be done before since the
            # underlying object does only exists at this point.
            for k, v in self._dict.iteritems():
                attr = '_' + k.lower()
                setattr(self._target_(), attr, v)
        else:
            self._init(self._dict['Name'], 'RooAddPdf') # single component won't be an RooAddPdf
            # Make sure we are the same as last time
            for k, v in self._dict.iteritems():
                print k,v,self._get(k)
                assert v == self._get(k)

    def _makeRecipe(self):
        yields = self._dict['Yields']
        pdfs = ','.join(['%s * %s' % (yields[p.GetName()], p.GetName())
                           if p.GetName() in yields else p.GetName()
                           for p in self._dict['PDFs']])
        return 'SUM::%s(%s)' % (self._dict['Name'], pdfs)

    def _separator(self):
        return '_P_'

class SimultaneousPdf(Pdf):
    def __init__(self, Name, **kwargs):
        d = { 'Name' : Name 
            , 'States' : kwargs.pop('States')
            , 'Cat' : kwargs.pop('SplitCategory')['Name']
            }
        # construct factory string on the fly...
        ## pdfs = sorted([(s, pdf) for s, pdf in d['States'].iteritems()], key = lambda (s, pdf): d['Cat'].lookupType(s).getVal())
        ## pdfs = [e[1] for e in pdfs]
        d['States'] = ','.join(['%s = %s' % (s, pdf['Name']) for s, pdf in d['States'].iteritems()])
        s = "SIMUL::%(Name)s(%(Cat)s,%(States)s)" % d
        self._declare(s)
        self._init(Name,'RooSimultaneous')
        Pdf.__init__(self , Name = Name , Type = 'RooSimultaneous')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
    def _make_pdf(self) : pass

class RealSumPdf( Pdf ):
    def __init__( self, name, functions, **kwargs ) :
        # get the name of the PDF, its functions and its coefficients
        self._dict = { 'Name' : name }
        self._dict['Functions'] = functions
        from itertools import repeat
        self._dict['Coefficients'] = kwargs.pop('coefficients',repeat(ConstVar( Name = 'one', Value = 1 ), len(self._dict['Functions'])))

        # make pdf
        self._make_pdf()
        del self._dict

    def _make_pdf(self):
        if self._dict['Name'] not in self.ws() :
            # create PDF in workspace and initialize
            self._declare( self._makeRecipe() )
            self._init( self._dict['Name'], 'RooRealSumPdf' )

            # change self._dict into attributes
            # (cannot be done before since the underlying object does only exists at this point)
            for k, v in self._dict.iteritems() :
                attr = '_' + k.lower()
                setattr( self._target_(), attr, v )
        else:
            # initialize
            self._init( self._dict['Name'], 'RooRealSumPdf' )

            # make sure we are the same as last time
            for k, v in self._dict.iteritems() : assert v == self._get(k)

    def _makeRecipe(self):
        # create string for declaration of PDF in workspace
        coefficients = ','.join( [ coef.GetName() for coef in self._dict['Coefficients'] ] )
        functions    = ','.join( [ func.GetName() for func in self._dict['Functions'] ] )
        return 'RealSumPdf::%s({%s}, {%s})' % ( self._dict['Name'], functions, coefficients )

class HistPdf( Pdf ) :
    def _make_pdf(self) : pass
    def __init__(self,Name,**kwargs) :
        d = { 'Name' : Name 
            , 'Observables' : list(kwargs.pop('Observables') )
            , 'Data' : kwargs.pop('Data')
            }
        binning = kwargs.pop('Binning', None )
        nbins = dict()
        if binning :
            for o,n in binning.iteritems() :
                if o not in d['Observables'] or o.getBins()==n : continue
                nbins[o] = o.getBins()
                o.setBins(n)
        dhs_name =  Name + '_' + '_'.join( i.GetName() for i in d['Observables'] )
        rdh = self.ws().put(RooDataHist( dhs_name, dhs_name,RooArgSet( i._var for i in d['Observables'] ), d['Data']))

        # construct factory string on the fly...
        self._declare("HistPdf::%s( { %s }, %s )" % (Name, ','.join( i.GetName() for i in d['Observables'] ), dhs_name  )  )
        self._init(Name,'RooHistPdf')
        Pdf.__init__(self , Name = Name , Type = 'RooHistPdf')
        for o,n in nbins.iteritems() : o.setBins(n) 
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class EditPdf( Pdf ) :
    def _make_pdf(self) : pass
    def __init__(self,Name,**kwargs) :
        d = { 'Name' : Name 
            , 'Original' : kwargs.pop('Original')
            , 'Rules' : kwargs.pop('Rules')
            }
        # construct factory string on the fly...
        self._declare("EDIT::%s( %s, %s )" % ( Name, d['Original'].GetName(), ','.join([ '%s=%s'%(k.GetName(),v.GetName()) for k,v in d['Rules'].iteritems()])  ) )
        self._init(Name,type(__dref__(d['Original'])).__name__)
        extraOpts = dict()
        cond =  d['Original'].ConditionalObservables()
        if cond : extraOpts['ConditionalObservables'] = cond
        Pdf.__init__(self , Name = Name , Type = type(__dref__(d['Original'])).__name__,**extraOpts)
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class HistFunc(RooObject):
    _getters = {'Histogram'   : lambda s: s.dataHist() 
               ,'Value'       : lambda s: s.getVal()
               ,'Observables' : lambda s: s.getDependents()
               } 
    def __init__(self, Name, **kwargs):
        __check_req_kw__( 'Histogram', kwargs )
        __check_req_kw__( 'Observables', kwargs )
        _hist = kwargs.pop('Histogram')
        if str(type(_hist)).find('TH1') == -1:
            raise TypeError, "HistFunc can only handle 1D historgrams"
        _dn = Name + '_data_hist'
        # Create Datahist and Import with density set to false
        _data = RooDataHist(_dn, _dn, RooArgList(*[__dref__(o) for o in kwargs['Observables']]),
                            RooFit.Import(_hist, False))
        self.ws().put(_data)
        self._declare('RooHistFunc::%s({%s}, %s)' % (Name, ','.join([o.GetName() for o in kwargs.pop('Observables')]), _dn))
        self._init(Name, 'RooHistFunc')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        
class EffProd(Pdf):
    def _make_pdf(self) : pass
    def __init__(self, Name, **kwargs):
        d = { 'Name' : Name 
            , 'Original' : kwargs.pop('Original')
            , 'Efficiency' : kwargs.pop('Efficiency')
            }
        if 'Type' in kwargs:
            d['Type'] = kwargs.pop('Type').__name__
        elif type(__dref__(d['Efficiency'])).__name__.find('Hist') != -1:
            d['Type'] = 'RooEffHistProd'
        else:
            d['Type'] = 'RooEffProd'
        if d['Type'] not in ['RooEffProd', 'RooEffHistProd']:
            raise TypeError, "An efficiency can only be of type RooEffProd or RooEffHistProd"
        # construct factory string
        self._declare("%s::%s(%s, %s)" % ( d['Type'], Name, d['Original'].GetName(),
                                           d['Efficiency'].GetName()))
        self._init(Name, d['Type'])
        extraOpts = dict()
        cond =  d['Original'].ConditionalObservables()
        if cond : extraOpts['ConditionalObservables'] = cond
        Pdf.__init__(self , Name = Name , Type = type(__dref__(d['Original'])).__name__,**extraOpts)
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class GenericPdf( Pdf ) :
    def _make_pdf(self) : pass
    def __init__(self,Name,**kwargs) :
        d = { 'Name' : Name 
            , 'Args' : ','.join( '%s'%i for i in kwargs.pop('Arguments') )
            , 'Formula' : kwargs.pop('Formula')
            }
        # construct factory string on the fly...
        self._declare("GenericPdf::%(Name)s( '%(Formula)s', { %(Args)s } )" % d )
        self._init(Name,'RooGenericPdf')
        Pdf.__init__(self , Name = Name , Type = 'RooGenericPdf')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class UniformPdf( Pdf ) :
    def _make_pdf(self) : pass
    def __init__(self,Name,**kwargs) :
        d = { 'Name' : Name 
            , 'Args' : ','.join( '%s'%i for i in kwargs.pop('Arguments') )
            }
        # construct factory string on the fly...
        self._declare("Uniform::%(Name)s( { %(Args)s } )" % d )
        self._init(Name,'RooUniform')
        Pdf.__init__(self , Name = Name , Type = 'RooUniform')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)


class BDecay( Pdf ) :
    def __init__(self,Name, **kwargs) :
        d = dict( name = Name
                , time = kwargs.pop('time')
                , tau = kwargs.pop('tau')
                , dGamma = kwargs.pop('dGamma')
                , dm = kwargs.pop('dm')
                , resolutionModel = kwargs.pop('resolutionModel')
                , coshCoef = kwargs.pop('coshCoef') if 'coshCoef' in kwargs else kwargs.pop('cosh')
                , sinhCoef = kwargs.pop('sinhCoef') if 'sinhCoef' in kwargs else kwargs.pop('sinh')
                , cosCoef = kwargs.pop('cosCoef') if 'cosCoef' in kwargs else kwargs.pop('cos')
                , sinCoef = kwargs.pop('sinCoef') if 'sinCoef' in kwargs else kwargs.pop('sin')
                , decayType = kwargs.pop( 'decayType', 'SingleSided' )
                )
        self._declare("BDecay::%(name)s( %(time)s, %(tau)s, %(dGamma)s, "\
                                          " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s, "\
                                          " %(dm)s, %(resolutionModel)s, %(decayType)s  )" % d )
        self._init(Name,'RooBDecay')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
            
 
class BTagDecay( Pdf ) :
    def _make_pdf(self) : pass
    def __init__( self, Name, **kwargs ) :
        from P2VVLoad import P2VVLibrary
        argDict = { 'Name' : Name, 'checkVars' : '1', 'decayType' : 'SingleSided' }

        # construct factory string on the fly...
        convert = lambda arg : str(arg) if type(arg) != list else '{%s}' % ','.join( str(listItem) for listItem in arg )
        if 'tagCat' in kwargs :
            for argName in [  'time', 'iTag', 'tagCat', 'tau', 'dGamma', 'dm'
                            , 'dilutions', 'ADilWTags', 'avgCEvens', 'avgCOdds', 'tagCatCoefs'
                            , 'coshCoef', 'sinhCoef', 'cosCoef', 'sinCoef'
                            , 'resolutionModel', 'decayType', 'checkVars'
                           ] :
                if argName not in argDict or argName in kwargs : argDict[argName] = convert(kwargs.pop(argName))

            self._declare("BTagDecay::%(Name)s( %(time)s, %(iTag)s, %(tagCat)s, %(tau)s, %(dGamma)s, %(dm)s, "\
                                              " %(dilutions)s, %(ADilWTags)s, %(avgCEvens)s, %(avgCOdds)s, %(tagCatCoefs)s,"\
                                              " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s, "\
                                              " %(resolutionModel)s, %(decayType)s, %(checkVars)s )" % argDict
                         )
        else :
            for argName in [  'time', 'iTag', 'tau', 'dGamma', 'dm'
                            , 'dilution', 'ADilWTag', 'avgCEven', 'avgCOdd'
                            , 'coshCoef', 'sinhCoef', 'cosCoef', 'sinCoef'
                            , 'resolutionModel', 'decayType', 'checkVars'
                           ] :
                if argName not in argDict or argName in kwargs : argDict[argName] = convert(kwargs.pop(argName))

            self._declare("BTagDecay::%(Name)s( %(time)s, %(iTag)s, %(tau)s, %(dGamma)s, %(dm)s, "\
                                              " %(dilution)s, %(ADilWTag)s, %(avgCEven)s, %(avgCOdd)s, "\
                                              " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s, "\
                                              " %(resolutionModel)s, %(decayType)s, %(checkVars)s )" % argDict
                         )

        self._init( Name, 'RooBTagDecay' )
        Pdf.__init__(  self
                     , Name = Name
                     , Type = 'RooBTagDecay'
                    )
        for ( k, v ) in kwargs.iteritems() : self.__setitem__( k, v )


class BinnedPdf( Pdf ) :
    def __init__( self, Name, **kwargs ) :
        from P2VVLoad import P2VVLibrary
        argDict = { 'Name' : Name }

        # declare PDF in workspace
        if 'baseCat' in kwargs :
            # single category dependence
            argDict['baseCat']  = str(kwargs.pop('baseCat'))
            argDict['coefList'] = '{%s}' % ','.join( str(listItem) for listItem in kwargs.pop('coefList') )
            self._declare( "BinnedPdf::%(Name)s( %(baseCat)s, %(coefList)s )" % argDict )

        elif 'baseCats' in kwargs :
            # multiple category dependence
            listArrayName = Name + '_coefLists'
            argDict['ignoreFirstBin'] = kwargs.pop( 'ignoreFirstBin', 0 )
            argDict['baseCats']       = '{%s}' % ','.join( str(listItem) for listItem in kwargs.pop('baseCats') )
            argDict['coefLists']      = listArrayName

            # create an array for the coefficient lists
            from ROOT import TObjArray
            wsListArray = TObjArray()
            wsListArray.SetName(listArrayName)

            # create coefficient lists
            from ROOT import RooArgList
            for varNum, coefList in enumerate( kwargs.pop('coefLists') ) :
                listName = Name + '_coefList%d' % varNum
                wsList = RooArgList(listName)
                for coef in coefList : wsList.add( self.ws().arg( str(coef) ) )
                wsListArray.Add(wsList)

            self.ws().put(wsListArray)
            self._declare( "BinnedPdf::%(Name)s( %(baseCats)s, %(coefLists)s, %(ignoreFirstBin)s )" % argDict )

        elif 'baseVar' in kwargs or 'baseVars' in kwargs :
            # continuous variable(s) dependence
            raise KeyError('P2VV - ERROR: BinnedPdf: dependence on continuous variables not (yet) implemented')

        else :
            raise KeyError('P2VV - ERROR: BinnedPdf: please specify variable(s)')

        # initialize PDF
        self._init( Name, 'RooBinnedPdf' )
        Pdf.__init__(  self
                     , Name = Name
                     , Type = 'RooBinnedPdf'
                    )
        for ( k, v ) in kwargs.iteritems() : self.__setitem__( k, v )

    def _make_pdf(self) : pass


class ResolutionModel(RooObject):
    _getters = { 'Parameters'  : lambda s : s._get('Parameters') }

    def __init__(self, **kwargs):
        __check_req_kw__( 'Type', kwargs )
        __check_req_kw__( 'Name', kwargs )

        # Save the keyword args as properties
        _dict = {}
        for k, v in kwargs.iteritems(): _dict[k] = v
        variables = list()
        if 'Parameters' in _dict: variables += list( _dict['Parameters'])
        deps = ','.join([v.GetName() for v in variables])
        _t = kwargs.pop('Type')
        if type(_t) != str : _t = _t.__name__
        name = kwargs.pop('Name')
        self._declare(  '%s::%s(%s)' % ( _t, name, deps) )
        self._init(name, 'RooResolutionModel')
        for k, v in _dict.iteritems():
            attr = '_' + k.lower()
            if hasattr(v, '__iter__'): v = frozenset(v)
            setattr(self._target_(), attr, v)

    def _get(self, name):
        return getattr(self,'_' + name.lower())
            

class AddModel(ResolutionModel) :
    def __init__(self,name,models,fractions,**kwargs) :
        # construct factory string on the fly...
        self._declare("AddModel::%s({%s},{%s})"%(name,','.join(i.GetName() for i in models),','.join(j.GetName() for j in fractions) ) )
        self._init(name,'RooAddModel')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
            


class Component(object):
    _d = {}
    def __init__(self,Name,*args,**kw) :
        # TODO: make things singletons, indexed by 'Name'
        if Name in Component._d : raise KeyError('Name %s is not unique'%name)
        self.name = Name
        Component._d[Name] = dict()
        Component._d[Name]['Name'] = Name
        if len(args) >1: raise IndexError('too many arguments %s' % args )
        if len(args) == 1:
            if type(args[0]) == dict :
                for i,j in args[0].iteritems() : self[i] = j
            else :
                for j in args[0] : self.append(j)
        if 'Yield' in kw : self.setYield( *kw.pop('Yield') )
        if kw : raise IndexError('unknown keyword arguments %s' % kw.keys() )
    def _yieldName(self) : return 'N_%s' % self.name
    def setYield(self, *args):
        y = None
        if len(args) == 1 and type(args[0]) == RealVar:
            y = args[0]
        elif len(args) == 3:
            n, nlo, nhi = args
            assert n>=nlo
            assert n<=nhi
            y = RealVar(self._yieldName(), MinMax=(nlo,nhi), Value=n)
        Component._d[self.name]['Yield'] = y
        Component._d[self.name]['Yield'].setAttribute('Yield',True)
    def __iadd__(self,pdf) :
        self.append(pdf)
        return self
    def append(self,pdf ) :
        if hasattr(pdf,'iteritems'):
            for obs,pdf in pdf.iteritems() : 
                self[obs] = pdf
        else :
            if not hasattr(pdf,'__iter__') : pdf = (pdf,)
            for p in pdf :
                obs = ( o for o in p.Observables() if o not in p.ConditionalObservables() ) 
                self[ obs ] = p
        
    def __setitem__(self, observable, pdf) :
        if not hasattr(observable,'__iter__') : observable = (observable,)

        # create a set of incoming observables
        k = set(o if type(o)==str else o.GetName() for o in observable )

        if type(pdf)!=type(None) : # allow 'None' as a placeholder for 'no PDF', which in RooFit turns into an implicit uniform PDF
            assert k == set( i.GetName() for i in pdf.Observables() if i not in pdf.ConditionalObservables() )
        ####
        # do NOT allow overlaps with already registered observables!!!!!! (maybe allow in future....)
        present = set()
        for i in  Component._d[self.name].iterkeys() :
            if type(i) != frozenset : continue # not an observable, but either name or yield
            for j in i : present.add(j)
                
        if not k.isdisjoint(present) : raise KeyError('sets are not disjoint, overlap: %s' % k.intersection(present))
        # TODO: check whether observable exists in the workspace...
        # TODO: and check it has its observable flag set

        ## Get the right sub-pdf from the Pdf object
        Component._d[self.name][frozenset(k)] = pdf

    def __eq__(self, other):
        o = other if type(other) == str else other.GetName()
        return self.GetName() == o

    def __eq__(self, other):
        return not self == other

    def __hash__(self):
        return self.GetName().__hash__()

    def GetName(self) : return self.name

    def __getitem__(self,k) :
        # TODO: if we return one a-priori build PDF, rename it properly??
        d = Component._d[self.name]
        
        # Catch yields and name here
        if type(k)==str and k in d : return d[k]

        if type(k) != frozenset : k = frozenset(k)
        if k not in d: 
            # try to build product -- note that d.keys() are non-overlapping by requirement
            # first, find the entry with the largest overlap, which is a subset (otherwise we'd have to marginalize)
            terms = []
            nk = k
            while len(nk) :
                overlap = lambda i : len(i.intersection(nk)) if type(i)==frozenset and i.issubset(nk) else 0  # avoid marginalization for now...
                from operator import itemgetter
                (kmax,mo) = max(((i,overlap(i)) for i in d.iterkeys()), key = itemgetter(1))
                if not mo : break # no overlap left...
                terms.append(kmax)
                nk = nk - kmax 
            if nk : raise IndexError('could not construct matching product -- no PDF for observables: %s' % [i for i in nk ])
            nk = frozenset.union(*terms)
            pdfs = [self[i] for i in terms if type(self[i])!=type(None)]
            #TODO: check that the conditional handling is internally consistent...
            d[nk] = ProdPdf(self.name, PDFs = pdfs)
        return d[k]

def buildPdf(Components, Observables, Name) :
    # multiply PDFs for observables (for each component)
    if not Observables : raise RuntimeError('no Observables??')
    obs = [o if type(o)==str else o.GetName() for o in Observables]
    args = {'Yields' : {}, 'PDFs'   : [] }
    for c in Components:
        pdf = c[obs]
        args['Yields'][pdf.GetName()] = c['Yield']
        args['PDFs'].append(pdf)
    if len(Components)==1 : return args['PDFs'][0] # TODO: how to change the name?
    # and sum components (inputs should already be extended)
    return SumPdf(Name,**args)

def buildSimultaneousPdf(Components, Observables, Spec, Name) :
    # multiply PDFs for observables (for each component)
    if not Observables : raise RuntimeError('no Observables??')
    if not Spec : raise RuntimeError('no Spec??')

    if len(Components)==1 : return args['PDFs'][0] # TODO: how to change the name?

    obs = [o if type(o)==str else o.GetName() for o in Observables]
    if len(Spec) != 1:
        raise ValueError, "Cannot handle more than one split category."
    states = {}
    key = Spec.keys()[0]
    observables, split_cat = key

    # - Figure out which components to clone
    # - Create new components where needed (don't forget yields)
    # - Make appropriate PDFs using new compoments when needed
    # - Make simultaneous PDF
    yield_names = {}
    for state, split_def in Spec[key].iteritems():
        args = {'Yields' : {}, 'PDFs' : []}
        suffix = '_'.join((split_cat['Name'], state))
        for c in Components:
            y = c['Yield']
            if not c in yield_names:
                yield_names[c] = y['Name']
            name = c['Name']
            rest = list(set(Observables) - set(observables))
            assert(set(rest).union(set(observables)) == set(Observables))
            rest_pdf = c[rest]
            if c in split_def:
                pdfs = []
                comp = Component('_'.join((name, suffix)), (split_def[c], rest_pdf),
                                 Yield = [y['Value']] + list(y['MinMax']))
                pdf = comp[obs]
                y = comp['Yield']
                yield_names[c]
            else:
                obs_pdf = c[observables]
                comp = Component('_'.join((name, suffix)), (obs_pdf, rest_pdf),
                                 Yield = [y['Value']] + list(y['MinMax']))
                pdf = comp[obs]
                y = comp['Yield']
            args['Yields'][pdf.GetName()] = y
            args['PDFs'].append(pdf)
        states[state] = SumPdf('_'.join((Name, suffix)) , **args)
    ## return states
    return SimultaneousPdf(Name, SplitCategory = split_cat, States = states)

                                 
