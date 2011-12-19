from RooFitDecorators import *
from functools import wraps

def __check_req_kw__( name, kwargs ) :
    if not name in kwargs : raise KeyError( 'Must provide kw argument %s' % name )
def __check_exists_already__( self ) :
    if self._name in self.ws() :
        raise StandardError( 'Recreating %s is not supported atm' % type(self) )

def __wrap__dref_var__( fun ) :
    @wraps(fun)
    def _fun(self,*args) :
        return fun(self, *tuple( i._var if hasattr(i,'_var') else i for i in args ) )
    return _fun

RooAbsCollection.__contains__ = __wrap__dref_var__( RooAbsCollection.__contains__ )
RooAbsData.table = __wrap__dref_var__( RooAbsData.table )


class RooObject(object) :
    # TODO: add functionality to ask for RooObject instances by name from a 'static' dictionary
    #       that way, when given eg. a RooRealVar, we ask for its name, and can look up the
    #       the corresponding RealVar:  x = RooObjects[ rrv.GetName() ]
    # NOTE: this looks like _objects, it has the same key, but the value isn't the underlying PyROOT
    #       object but it's wrapper...
    _ws = None
    _dict = None
    _setters = {'Title'      : lambda s,v : s.SetTitle(v)
               }
    _getters = {'Name'       : lambda s : s.GetName()
               ,'Title'      : lambda s : s.GetTitle()
               ,'Value'      : lambda s : s.getVal()
               }
    def _factory(self,spec) :
        return self.ws().factory(spec)
        
    def __setitem__(self,k,v):
        from itertools import ifilter, imap
        for setters in imap( lambda x: x._setters, ifilter( lambda x : hasattr(x,'_setters'), type(self).__mro__) ): 
            if k in setters :  return setters[k](self,v )
        raise KeyError('%s is not known for class %s' % (k, type(self) ) )
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
        if not hasattr(ws, '_objects') : ws._objects  = {}  # name -> object
        if not hasattr(ws, '_mappings'): ws._mappings = {}
        if not hasattr(ws, '_spec'):     ws._spec = {}      # factory string -> object
        
    def _declare(self,spec):
        """
        create the underlying C++ object in the workspace
        """
        # TODO: add mapping of name -> spec so we can gate identical invocations.
        #       problem is that we don't know 'name' a-priori....
        #       so maybe we either require a name, and verify what we got back has
        #       the matching name, or, in case name is 'None' then it has to be unique
        #       and not yet used???
        #   No! we turn the mapping around, from spec -> ( object -> ) name
        #
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
        return set( i.GetName() for i in self._var.getVariables() if i.getAttribute('Observable') )
    def Parameters(self) :
        return set( i.GetName() for i in self._var.getVariables() if not i.getAttribute('Observable') )

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
    _getters = {'Observable' : lambda s : s.observable() 
               ,'Index'      : lambda s : s.getIndex() 
               ,'Label'      : lambda s : s.getLabel()
               ,'States'     : lambda s : s.states()
               } 
    _setters = {'Observable' : lambda s,v : s.setObservable(v) 
               ,'Index'      : lambda s,v : s.setIndex(v)
               ,'Label'      : lambda s,v : s.setLabel(v)
               ,'Constant'   : lambda s,v : s.setConstant(v)
               }

    def __init__(self,name,**kwargs):
        if name not in self.ws():
            # construct factory string on the fly...
            __check_req_kw__( 'States', kwargs )
            states = None
            if   type(kwargs['States']) == list:
                states = ','.join(kwargs.pop('States'))
            elif type(kwargs['States']) == dict:
                states = ','.join(['%s=%d' % i for i in kwargs.pop('States').items()])
            else:
                raise KeyError('%s does not exist yet, bad States defined.' % name)

            # Create the category and extract states into storage
            self._declare("%s[%s]"%(name,states))
            self._init(name,'RooCategory')
            self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(name,'RooCategory')
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                # Skip these to avoid failure in case we were loaded from a
                # DataSet in the mean time
                if k in ['Index', 'Label']: continue
                assert v == self[k]
            
    def states(self):
        return self._states

class SuperCategory( Category ) :
    def __init__(self,name,cats,**kwargs):
        if name not in self.ws():
            self._declare("SuperCategory::%s({%s})"%(name,','.join( [ c.GetName() for c in cats ] ) ) )
            self._init(name,'RooSuperCategory')
            self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(name,'RooSuperCategory')
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                # Skip these to avoid failure in case we were loaded from a
                # DataSet in the mean time
                if k in ['Index', 'Label']: continue
                assert v == self[k]

class MappedCategory( Category ) :
    def __init__(self,name,cat,mapper,**kwargs):
        if name not in self.ws():
            # construct factory string on the fly: no factory string for SuperCategory???
            from ROOT import RooMappedCategory
            cast = lambda i : i._target_() if hasattr(i,'_target_') else i
            obj =  RooMappedCategory(name,name,cast(cat) )
            for k,vs in mapper.iteritems() : 
                for v in vs : obj.map( v, k )
            self.ws()._objects[ name ] = self.ws().put( obj )
            self._init(name,'RooMappedCategory')
            self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(name,'RooMappedCategory')
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                # Skip these to avoid failure in case we were loaded from a
                # DataSet in the mean time
                if k in ['Index', 'Label']: continue
                assert v == self[k]


class Product(RooObject) :
    def __init__(self,name,fargs,**kwargs) :
        spec =  "prod::%s(%s)"%(name,','.join(i['Name'] for i in fargs)) 
        self._declare( spec )
        self._init(name,'RooProduct')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
            

class Addition(RooObject) :
    def __init__(self,name,fargs,**kwargs) :
        # construct factory string on the fly...
        self._declare("sum::%s(%s)"%(name,','.join(i.GetName() for i in fargs)) )
        self._init(name,'RooAddition')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class FormulaVar (RooObject): 
    _getters = {'Formula'    : lambda s : s.formula()
               ,'Dependents' : lambda s : s.dependents()
               ,'Value'      : lambda s : s.getVal()
               }
    def __init__(self, name, formula, fargs, **kwargs):
        # construct factory string on the fly...
        spec = "expr::%s('%s',{%s})" % (name, formula, ','.join(i.GetName() for i in fargs)) 
        present = name in self.ws() 
        match = present and spec == self.ws()[name].getStringAttribute('RooFitWrappers.RooObject::spec')
        if not present or not match :
            self._declare(spec)
            self._init(name, 'RooFormulaVar')
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(name,'RooFormulaVar')
            for k, v in kwargs.iteritems():
                # Skip these to avoid failure in case we were loaded from a
                # DataSet in the mean time
                if k == 'Value':
                    continue
                assert v == self[k]

class ConstVar (RooObject): 
    def __init__(self,name,**kwargs):
        if name not in self.ws():
            # construct factory string on the fly...
            __check_req_kw__( 'Value', kwargs )
            self._declare("ConstVar::%s(%s)"%(name,kwargs.pop('Value')))
            self._init(name,'RooConstVar')
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(name,'RooConstVar')
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                # Skip these to avoid failure in case we were loaded from a
                # DataSet in the mean time
                if k == 'Value':
                    continue
                assert v == self[k]
        
class P2VVAngleBasis (RooObject) : 
    # TODO: make a 'borg' out of this which avoids re-creating ourselves by construction...
    def __init__(self, angles, i,j,k,l,c=1,i2=None,j2=None,k2=None,l2=None) :
        assert c!=0
        # compute name, given angles,i,j,k,l,c!
        # WARNING: angles may contain barebones PyROOT objects!!!
        name = '_'.join( angles[a].GetName()  for a in ['cpsi','ctheta','phi']) # aargh... too long for RooFit workspace parsing  code....
        name = ''
        second = (i2 or j2 or k2 or l2 )
        if second : assert i2!=None and j2!=None and k2!=None and l2!=None
        name = 'p2vvab_%s%d%d%d%d' % (name, i,j,k,l) 
        if second : name = '%s%d%d%d%d' % (name, i2,j2,k2,l2)
        if c!=1 :   name = '%s%3.2f'%(name,c) # truncate printing of 'c' to 3 decimals?
        name = name.replace('-', 'm')
        name = name.replace('.', 'd')
        if second :
            spec = "RooP2VVAngleBasis::%s(%s, %s, %s, %d,%d,%d,%d, %d,%d,%d,%d, %f)" % (name, angles['cpsi'].GetName(),angles['ctheta'].GetName(),angles['phi'].GetName(), i, j, k, l, i2, j2, k2, l2, c) 
        else :
            spec = "RooP2VVAngleBasis::%s(%s, %s, %s, %d,%d,%d,%d, %f)"              % (name, angles['cpsi'].GetName(),angles['ctheta'].GetName(),angles['phi'].GetName(), i, j, k, l, c) 
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

        cast = lambda var : var._target_() if hasattr( var, '_target_' ) else var

        # create moment
        from P2VVLoad import P2VVLibrary
        from ROOT import RooRealMoment
        AbsRealMoment.__init__( self, RooRealMoment( cast(self._basisFunc), self._norm ) )
          
class RealEffMoment( AbsRealMoment ):
    def __init__( self, BasisFunc, Norm, PDF, NormSet ) :
        # get arguments
        self._basisFunc = BasisFunc
        self._norm      = Norm
        self._pdf       = PDF
        self._normSet   = NormSet

        cast = lambda var : var._target_() if hasattr( var, '_target_' ) else var

        # build a RooFit normalisation set
        from ROOT import RooArgSet
        self._rooNormSet = RooArgSet( cast(var) for var in self._normSet )

        # create efficiency moment
        from P2VVLoad import P2VVLibrary
        from ROOT import RooRealEffMoment
        AbsRealMoment.__init__( self, RooRealEffMoment( cast(self._basisFunc), self._norm, cast(self._pdf), self._rooNormSet ) )


class RealVar (RooObject): 
    # WARNING: multiple instances don't share proxy state at this time...
    # TODO: move common things like Name and Title in RooObject...
    # TODO: provide scaffolding in RooObject to extend getters & setters on a class-by-class basis
    _getters = {'Observable' : lambda s : s.observable() 
               ,'Unit'       : lambda s : s.getUnit() 
               ,'Value'      : lambda s : s.getVal()
               ,'MinMax'     : lambda s : s.getRange()
               }
    _setters = {'Observable' : lambda s,v : s.setObservable(v) 
               ,'Unit'       : lambda s,v : s.setUnit(v) 
               ,'Value'      : lambda s,v : s.setVal(v)
               ,'MinMax'     : lambda s,v : s.setRange(v)
               ,'Constant'   : lambda s,v : s.setConstant(v)
               }

    def __init__(self,Name ,**kwargs):
        if 'name' in kwargs : raise RuntimeError('Please replace name argument with Name = xyz' )
        if Name not in self.ws():
            # construct factory string on the fly...
            if 'Value'  not in kwargs  and 'MinMax' not in kwargs:
                raise KeyError('%s does not exist yet, neither Value nor MinMax specified'%Name)
            if 'Value' not in kwargs:
                (mi,ma) = kwargs.pop('MinMax')
                self._declare("%s[%s,%s]"%(Name,mi,ma))
            elif 'MinMax' not in kwargs:
                self._declare("%s[%s]"%(Name,kwargs.pop('Value')))
            else:
                (mi,ma) = kwargs.pop('MinMax')
                val = kwargs.pop('Value')
                if val < mi or val > ma : raise RuntimeError('Specified Value not contained in MinMax')
                self._declare("%s[%s,%s,%s]"%(Name,val,mi,ma))
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
                try:
                    assert v == self[k]
                except AssertionError:
                    print "%s is not the same" % k
                    raise
            
    # overrule RooRealVar.setRange
    @wraps(RooRealVar.setRange)
    def setRange(self, *args):
        assert args
        if type(args[0]) == str :
            assert len(args)==2
            self._var.setRange( args[0], *args[1] )
        else  :
            assert len(args)==1
            (mi,ma) = args[0]
            self._var.setRange(mi,ma)
            if self.getVal() < mi or self.getVal() > ma : self.setVal(0.5*(ma+mi) )

    def getRange(self):
        return self._target_().getMin(), self._target_().getMax()

##TODO, factor out common code in Pdf and ResolutionModel

class Pdf(RooObject):
    # TODO: support conditional observables!!!
    #       multiply automatically with flat PDF for conditional observables in no PDF for them is given??
    #       or intercept fitTo and add ConditionalObservables = ( ... ) ??
    _getters = {'Observables' : lambda s : s._get('Observables') 
               ,'Type'        : lambda s : s._get('Type')
               ,'Parameters'  : lambda s : s._get('Parameters')
               ,'Name'        : lambda s : s._get('Name')
               }

    ## TODO: define operators
    def __init__(self, Name, **kwargs):
        __check_req_kw__( 'Observables', kwargs )
        __check_req_kw__( 'Type', kwargs )

        # Save the keyword args as properties
        self._dict = kwargs
        self._dict['Name'] = Name
        self._make_pdf()
        print 'Pdf(%s): deduced observables: %s' %  (Name,self.Observables() )
        print 'Pdf(%s): specified observables : %s' % (Name,set( i.GetName() for i in kwargs['Observables'] ))

    def __str__(self):
        d = dict([(a, self[a]) for a in Pdf._getters if hasattr(self, a)])
        return '%s' % d

    def _get(self, name):
        return getattr(self._target_(), '_' + name.lower())
    
    def __getitem__(self, k):
        if self._dict and k in self._dict:
            return self._dict[k]
        else:
            try:
                return RooObject.__getitem__( self, k )
            except AttributeError as error:
                raise KeyError(str(error))

    def _make_pdf(self):
        if self._dict['Name'] not in self.ws():
            v = list(self._dict['Observables'])  
            if 'Parameters' in self._dict: v += list(self._dict['Parameters'])
            if 'ResolutionModel' in self._dict: v.append(self._dict['ResolutionModel'])
            if 'Options' in self._dict: v += list(self._dict['Options'])
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
                if k in ['Observables', 'Parameters']: v = frozenset(v)
                attr = '_' + k.lower()
                setattr(self._target_(), attr, v)
            del self._dict # no longer needed
        else:
            self._init(self._dict['Name'], 'RooAbsPdf')
            # Make sure we are the same as last time
            for k, v in self._dict.iteritems():
                if k in ['Observables', 'Parameters']: v = frozenset(v)
                if v != self._get(k) : print k,v,self._get(k)
                assert v == self._get(k)

    def _separator(self):
        return '_'


    @wraps(RooAbsPdf.fitTo)
    def fitTo( self, data, **kwargs ) :
        if 'ConditionalObservables' in kwargs :
            cvrt = lambda i : i._target_() if hasattr(i,'_target_') else i
            kwargs['ConditionalObservables'] = RooArgSet( cvrt(var) for var in kwargs.pop('ConditionalObservables') )
        return self._var.fitTo( data, **kwargs )

    @wraps(RooAbsPdf.generate)
    def generate(self, whatvars, *args, **kwargs):
        #if not whatvars : whatvars = [ i for i in self._var.getVariables() if i.getAttribute('Observable') ]
        cvrt = lambda i : i._target_() if hasattr(i,'_target_') else i
        return self._var.generate(RooArgSet([ cvrt(i) for i in whatvars] if not isinstance(cvrt(whatvars),RooAbsCategory) else cvrt(whatvars)), *args,**kwargs)

    @wraps(RooAbsPdf.plotOn)
    def plotOn( self, frame, **kwargs ) :
        if 'Slice' in kwargs :
            cvrt = lambda i : i._target_() if hasattr(i,'_target_') else i
            sl = kwargs.pop('Slice')
            kwargs['Slice'] = ( cvrt(sl[0]), sl[1] )
        return self._var.plotOn( frame, **kwargs )


    def edit(self, Name, rule ) :
        #TODO: wrap result in Pdf!!!
        return self.ws().factory('EDIT::%s(%s,%s)' % ( Name, self.GetName(), ','.join( '%s=%s'%(k.GetName(),v.GetName()) for k,v in rule.iteritems() ) ) )


class ProdPdf(Pdf):
    # TODO: support conditional terms, use 'Conditional' key word for that...
    # ProdPdf( 'foo', [a,b], Conditional = [c,d] ) -> PROD::foo(a,b|c,d)
    def __init__(self, Name, PDFs, **kwargs):
        self._dict = { 'PDFs' : frozenset(PDFs)
                     , 'Name' : Name + '_' + self._separator().join([i.GetName() for i in PDFs])
                     , 'Observables' : frozenset( i for p in PDFs for i in p['Observables'] )
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
        pdfs = ','.join([p.GetName() for p in self._dict['PDFs']])
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
        self._dict['Observables'] = frozenset( i for p in pdfs for i in p['Observables' ] )
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

class RealSumPdf( Pdf ):
    def __init__( self, name, functions, **kwargs ) :
        # get the name of the PDF, its functions and its coefficients
        self._dict = { 'Name' : name }
        self._dict['Functions'] = functions
        from itertools import repeat
        self._dict['Coefficients'] = kwargs.pop('coefficients',repeat(ConstVar( 'one', Value = 1. ), len(self._dict['Functions'])))

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
        Pdf.__init__(self
                    , Name = Name
                    , Type = 'RooGenericPdf'
                    , Observables = () # let Pdf figure this out itself...
                    )
        kwargs.pop('Observables',None) # TODO: remove this bad hack!!! Observables are automatically determined
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class BTagDecay( Pdf ) :
    def _make_pdf(self) : pass
    def __init__(self,Name,**kwargs) :
        d = { 'Name' : Name, 'checkVars' : 1, 'decayType' : 'SingleSided', 'tagCat' : None }
        for i in [ 'time','iTag','tau', 'dGamma', 'dm'
                 , 'dilution', 'ADilWTag', 'avgCEven', 'avgCOdd'
                 , 'coshCoef', 'sinhCoef', 'cosCoef', 'sinCoef'
                 , 'resolutionModel', 'decayType', 'checkVars' ] :
           if i not in d or i in kwargs : d[i] = kwargs.pop(i) 

        from P2VVLoad import P2VVLibrary
        # construct factory string on the fly...
        if d['tagCat'] :
            self._declare("BTagDecay::%(Name)s( %(time)s, %(iTag)s, %(tagCat)s, %(tau)s, %(dGamma)s, %(dm)s, "\
                                              " %(dilutions)s, %(ADilWTags)s, %(avgCEvens)s, %(avgCOdds)s, %(tagCatCoefs)s"\
                                              " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s, "\
                                              " %(resolutionModel)s, %(decayType)s, %(checkVars)s )" % d
                         )
        else :
            self._declare("BTagDecay::%(Name)s( %(time)s, %(iTag)s, %(tau)s, %(dGamma)s, %(dm)s, "\
                                              " %(dilution)s, %(ADilWTag)s, %(avgCEven)s, %(avgCOdd)s, "\
                                              " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s, "\
                                              " %(resolutionModel)s, %(decayType)s, %(checkVars)s )" % d
                         )

        self._init(Name,'RooBTagDecay')
        Pdf.__init__(self
                    , Name = Name
                    , Type = 'RooBTagDecay'
                    , Observables = () # let Pdf figure this out itself...
                    )
        kwargs.pop('Observables',None) # TODO: remove this bad hack!!! Observables are automatically determined
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class ResolutionModel(RooObject):
    _getters = {'Observables' : lambda s : s._get('Observables')
               ,'Parameters'  : lambda s : s._get('Parameters')
               ,'Type'        : lambda s : s._get('Type')
               }

    def __init__(self, name, **kwargs):
        __check_req_kw__( 'Type', kwargs )
        __check_req_kw__( 'Observables', kwargs )

        if name not in self.ws():
            # Save the keyword args as properties
            _dict = {}
            _dict['Name'] = name
            for k, v in kwargs.iteritems(): _dict[k] = v
            variables = list(_dict['Observables'])
            if 'Parameters' in _dict: variables += list( _dict['Parameters'])
            deps = ','.join([v.GetName() for v in variables])
            _t = kwargs.pop('Type')
            if type(_t) != str : _t = _t.__name__
            self._declare(  '%s::%s(%s)' % ( _t, name, deps) )
            self._init(name, 'RooResolutionModel')
            for k, v in _dict.iteritems():
                attr = '_' + k.lower()
                if hasattr(v, '__iter__'): v = frozenset(v)
                setattr(self._target_(), attr, v)
        else:
            self._init(name, 'RooResolutionModel')
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                if hasattr(v, '__iter__'): v = frozenset(v)
                assert v == self._get(k)

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
                for j in args[0] :
                    i = tuple( o.GetName() for o in j.getVariables() if o.getAttribute('Observable') )
                    self[i] = j 
        if 'Yield' in kw : self.setYield( *kw.pop('Yield') )
        if kw : raise IndexError('unknown keyword arguments %s' % kw.keys() )
    def _yieldName(self) : return 'N_%s' % self.name
    def setYield(self, n, nlo, nhi) :
        Component._d[self.name]['Yield'] = RealVar(self._yieldName(), MinMax=(nlo,nhi), Value=n).GetName()
    def __iadd__(self,pdf) :
        self.append(pdf)
        return self
    def append(self,pdf ) :
        self[ pdf.Observables() ] = pdf
        
    def __setitem__(self, observable, pdf) :
        if not hasattr(observable,'__iter__') : observable = (observable,)

        # create a set of incoming observables
        k = set(o if type(o)==str else o.GetName() for o in observable )
        #### 
        print 'Component: specified observables: %s' % k
        print 'Component: deduced observables: %s' % pdf.Observables()
        assert k == pdf.Observables()
        ####
        # do NOT allow overlaps with already registered observables!!!!!! (maybe allow in future....)
        present = set()
        for i in  Component._d[self.name].iterkeys() :
            if type(i) != frozenset : continue # not an observable, but either name or yield
            for j in i : present.add(j)
                
        if not k.isdisjoint(present) : raise KeyError('sets are not disjoint, overlap: %s' % k.intersection(present))
        # TODO: check whether observable exists in the workspace...
        # TODO: and check it has its observable flag set
        ### parse recipe: x(a,b[1,2,3],c(4,5,6)) -> ('x' , 'a,b[1,2,3],c(4,5,6)')
        ### i.e. such that %s(%s)'%(x[0], x[1]) recovers the result...
        ## import re
        ## ex = '^([^(]+)\((.*)\)$'
        ## r = re.match(ex,pdfrecipe)
        ## if not r : raise KeyError('could not parse recipe %s'%pdfrecipe)
        ## (typ,args) = r.groups()
        ## if ':' in typ : 
        ##     raise KeyError('please do not explicitly name pdfs -- %s' % typ)
        #TODO: should build the PDF at this point!

        ## Get the right sub-pdf from the Pdf object
        Component._d[self.name][frozenset(k)] = pdf


    def __getitem__(self,k) :
        # TODO: if we return one a-priori build PDF, rename it properly??
        d = Component._d[self.name]
        
        # Catch yields and name here
        if type(k)==str and k in d : return d[k]

        if type(k) != frozenset : k = frozenset(k)
        if k not in d : 
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
            pdfs = [self[i] for i in terms]
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
    # and sum components (inputs should already be extended)
    return SumPdf(Name,**args)
