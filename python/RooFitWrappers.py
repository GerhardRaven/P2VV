from RooFitDecorators import *
from functools import wraps

def __check_req_kw__( name, kwargs ) :
    if not name in kwargs : raise KeyError( 'Must provide kw argument %s' % name )
def __check_exists_already__( self ) :
    if self._name in self.ws() :
        raise StandardError( 'Recreating %s is not supported atm' % type(self) )
def __check_name_syntax__( name ) :
    from string import whitespace
    if any( whiteChar in name for whiteChar in whitespace ) :
        raise KeyError( 'Whitespace in variable names not supported: "%s"' % name )

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
        for setters in imap( lambda x: x._setters, ifilter( lambda x : hasattr(x,'_setters'), type(self).__mro__) ) :
            if k in setters :  return setters[k](self,v )
        raise KeyError('\'%s\' is not known for class %s' % (k, type(self) ) )
    def __getitem__(self,k):
        from itertools import ifilter, imap
        for getters in imap( lambda x: x._getters, ifilter( lambda x : hasattr(x,'_getters'), type(self).__mro__) ) :
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
        # get name string
        if type(Name) != str : Name = Name.GetName()

        if Name not in self.ws()._rooobjects :
            # object is not in work space dictionary of RooObjects
            if Name in self.ws() :
                # try to create a RooObject wrapper if the requested object exists in work space
                import ROOT
                if   isinstance( self.ws()[Name], ROOT.RooRealVar  ) : return RealVar(Name)
                elif isinstance( self.ws()[Name], ROOT.RooCategory ) : return Category(Name)
            else :
                # object does not exist
                raise KeyError, 'P2VV - ERROR: RooObject._rooobject(): object does not exist (%s)' % Name

        # return object
        return self.ws()._rooobjects[Name]

    # WARNING: the object 'o' given to _addObject should NEVER be used again
    # instead, use the item returned by _addObject
    def _addObject(self, o):
        if o.GetName() not in self.ws(): 
            # print 'P2VV: WARNING: importing %s into workspace.... ' % o.GetName() 
            # assert o.GetName() != 'sig_t_angles'
            self.ws().put(o)
        o = self.ws()[o.GetName()]
        if o.GetName() not in self.ws()._objects:
            self.ws()._objects[o.GetName()] = o
        else:
            raise TypeError, "Adding the same object twice should not happen! %s" % o.GetName()
        return o

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
            # print "Spec: %s" % spec
            x = self._factory(spec)
            if not x: raise NameError("workspace factory failed to return an object for factory string '%s' "%spec)
            if hasattr(x,'setStringAttribute') : x.setStringAttribute('RooFitWrappers.RooObject::spec',spec)
            #
            # Keep the PyROOT objects in a container so they don't get garbage
            # collected.
            # Note: use explicit GetName, not str, as x is a 'bare' PyROOT object!!!
            x = self._addObject(x)
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
        if name in self.ws()._objects :
            x = self.ws()._objects[name]
        elif name in self.ws() :
            x = self._addObject( self.ws()[name] )
        else :
            raise KeyError, 'P2VV - ERROR: RooObject._init(): object not found in work space (%s)' % name

        if not x.InheritsFrom(type_) :
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
        return set( self._rooobject(i) for i in self._var.getVariables() if i.getAttribute('Observable') )
    def Parameters(self) :
        return set( self._rooobject(i) for i in self._var.getVariables() if not i.getAttribute('Observable') )

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
    def _addObject(self,o) :
        if o.GetName() not in self.ws()._objects :
            self.ws()._objects[o.GetName()] = o
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
        __check_name_syntax__(name)
        spec = 'set::%s(%s)' % (name, ','.join( i.GetName() for i in args) )
        self._declare(spec)
        self._init(name,'RooArgSet')




class Category (RooObject) :
    _getters = {  'Index'    : lambda s   : s.getIndex()
                , 'Label'    : lambda s   : s.getLabel()
                , 'States'   : lambda s   : s.states()
               }
    _setters = { 'Index'     : lambda s,v : s.setIndex(v)
                , 'Label'    : lambda s,v : s.setLabel(v)
                , 'Constant' : lambda s,v : s.setConstant(v)
               }

    def __init__( self, Name, **kwargs ) :
        __check_name_syntax__(Name)

        if Name not in self.ws() :
            states = kwargs.pop( 'States', None )
            data   = kwargs.pop( 'Data', None )
            if   type(states) == list : states = ','.join(states)
            elif type(states) == dict : states = ','.join( [ '%s=%d' % i for i in states.iteritems() ] )

            if states and not data :
                # Create the category and extract states into storage
                obj = self._declare('%s[%s]' % ( Name, states ) )

            elif states and data :
                index = kwargs.pop('DataIndex')
                from ROOT import RooCategory
                obj = RooCategory( Name, Name )
                for l, i in states.iteritems() : obj.defineType( l, i )
                obj.setIndex(index)
                obj = data.addColumn(obj)
                obj = self._addObject(obj)

            elif not states and not data:
                from ROOT import RooCategory
                obj = RooCategory( Name, Name )
                obj = self._addObject(obj)

            else:
                raise RuntimeError('Passing data and not states is illegal')

            self._init( Name,'RooCategory' )
            self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

        else :
            self._init( Name, 'RooCategory' )
            # Make sure we are the same as last time
            for k, v in kwargs.iteritems():
                assert v == self[k], '\'%s\' is not the same for %s' % ( k, Name )

    def states(self):
        return self._states

    def __iter__(self):
        return self._target_().__iter__()

    def __contains__(self, v):
        for cat in self:
            if cat.getVal() == v or cat.GetName() == v:
                return True
        else:
            return False

class BinningCategory( Category ) :
    def __init__( self, Name, **kwargs ) :
        __check_req_kw__( 'Observable', kwargs )
        __check_req_kw__( 'Binning', kwargs )

        obs = __dref__( kwargs.pop('Observable') )
        binning = kwargs.pop('Binning')
        if type(binning) != str : binning = binning.GetName()

        from ROOT import RooBinningCategory
        if 'CatTypeName' in kwargs : binCat = RooBinningCategory( Name, Name, obs, binning, kwargs.pop('CatTypeName') )
        else                       : binCat = RooBinningCategory( Name, Name, obs, binning                            )

        if kwargs.pop( 'Fundamental', False ) :
            __check_req_kw__( 'Data', kwargs )
            data = kwargs.pop('Data')
            if type(data) not in [ list, tuple ] : data = [ data ]

            cat = data[0].addColumn(binCat)
            for dataSet in data[ 1 : ] : dataSet.addColumn(binCat)
            cat = self._addObject( __dref__(cat) )
            self._init( Name, 'RooCategory' )

        else :
            binCat = self._addObject(binCat)
            self._init( Name, 'RooBinningCategory' )

        for ( k, v ) in kwargs.iteritems() : self.__setitem__( k, v )

class ThresholdCategory( Category ) :
    def __init__(self,Name,**kwargs):
        __check_req_kw__( 'Observable', kwargs )
        __check_req_kw__( 'Data', kwargs )
        obs = __dref__(kwargs.pop('Observable'))
        data = kwargs.pop('Data', None)
        boundlist = kwargs.pop('Boundaries')
        defaultstring = kwargs.pop('Default')
        from ROOT import RooThresholdCategory
        obj = RooThresholdCategory(Name, Name, obs, defaultstring)
        for i,value in [(i,v) for i,v in enumerate(boundlist)][1:]:
            obj.addThreshold(value,"Bin%s"%(i))
        obj = data.addColumn(__dref__(obj))
        obj = self._addObject(obj)
        t = 'RooCategory'
        self._init(Name, t)
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class SuperCategory( Category ) :
    def __init__(self,Name,cats,**kwargs):
        __check_name_syntax__(Name)
        data = kwargs.pop('Data', None)
        if data:
            dummy = self._declare('SuperCategory::%s_sc({%s})'%(Name,','.join( c.GetName() for c in cats )))
            self._declare( "dataobs::%s(%s,%s)"%(Name,data.GetName(),dummy.GetName()))
            t = 'RooCategory'
        else:
            self._declare("SuperCategory::%s({%s})"%(Name,','.join( c.GetName() for c in cats ) ) )
            t = 'RooSuperCategory'
        self._init(Name, t)
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
            data = kwargs.pop('Data', None)
            def init(o, t):
                obj = self._addObject(o)
                self._init(obj.GetName(), t)
            if data:
                obj = data.addColumn(__dref__(obj))
                init(obj, 'RooCategory')
            else:
                init(obj, 'RooMappedCategory')
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

## class MultiCategory( Category ) :
##     def __init__(self,Name,cats,**kwargs):
##         __check_name_syntax__(Name)
##         data = kwargs.pop('Data', None)
##         t = 'RooMultiCategory'
##         if data:
##             from ROOT import RooMultiCategory
##             obj = RooMultiCategory(Name, Name, RooArgSet(*cats))
##             obj = data.addColumn(__dref__(obj))
##             obj = self._addObject(obj)
##             t = 'RooCategory'
##         else:
##             self._declare("MultiCategory::%s({%s})"%(Name,','.join( [ c.GetName() for c in cats ] ) ) )
##         self._init(Name, t)
##         self._target_()._states = dict( ( s.GetName(), s.getVal()) for s in self._target_() )
##         for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class Product(RooObject) :
    def __init__(self,Name,fargs,**kwargs) :
        __check_name_syntax__(Name)
        spec =  "prod::%s(%s)"%(Name,','.join(i.GetName() for i in fargs))
        self._declare( spec )
        self._init(Name,'RooProduct')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)


class Addition(RooObject) :
    def __init__(self,Name,fargs,**kwargs) :
        # construct factory string on the fly...
        __check_name_syntax__(Name)
        def cn( x ) :
            try :
                return x.GetName()
            except :
                (a,b) = x
                return '%s*%s' % ( a.GetName(),b.GetName() )
        self._declare( "sum::%s(%s)" % ( Name,','.join( cn(i) for i in fargs ) ) )
        self._init(Name,'RooAddition')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class FormulaVar (RooObject) :
    _getters = {'Formula'    : lambda s : s.formula()
               ,'Dependents' : lambda s : s.dependents()
               ,'Value'      : lambda s : s.getVal()
               }
    def __init__(self, Name, formula, fargs, **kwargs):
        # construct factory string on the fly...
        __check_name_syntax__(Name)
        data = kwargs.pop('data', None)
        if not data:
            spec = "expr::%s('%s',{%s})" % (Name, formula, ','.join(i.GetName() for i in fargs))
            self._declare(spec)
            self._init(Name, 'RooFormulaVar')
        else:
            l = RooArgList()
            for arg in fargs:
                l.add(__dref__(arg))
            from ROOT import RooFormulaVar
            form = RooFormulaVar(Name, Name, formula, l)
            form = data.addColumn(form)
            form = self._addObject(form)
            self._init(Name, 'RooRealVar')
            self.setObservable(True)
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class RealCategory(RooObject) :
    def __init__(self, Name,Category ) :
        __check_name_syntax__(Name)
        spec = 'RooRealCategory::%s(%s)'%(Name, Category)
        self._declare(spec)
        self._init(Name,'RooRealCategory')

class ConstVar(RooObject) :
    def __init__(self,**kwargs):
        # construct factory string on the fly...
        __check_req_kw__( 'Value', kwargs )
        __check_req_kw__( 'Name', kwargs )
        __check_name_syntax__( kwargs['Name'] )
        self._declare("ConstVar::%(Name)s(%(Value)s)" % kwargs )
        (Name, value) = (kwargs.pop('Name'),kwargs.pop('Value'))
        self._init(Name,'RooConstVar')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class LinearVar(RooObject) :
    def __init__(self,**kwargs):
        # construct factory string on the fly...
        __check_req_kw__('Name', kwargs)
        __check_req_kw__('Observable', kwargs )
        __check_req_kw__('Slope', kwargs)
        __check_req_kw__('Offset', kwargs)
        __check_name_syntax__(kwargs['Name'])
        args = {}
        for k in ['Name', 'Observable', 'Slope', 'Offset']:
            v = kwargs.pop(k)
            args[k] = v if type(v) == str else v.GetName()
        self._declare("LinearVar::%(Name)s(%(Observable)s,%(Slope)s,%(Offset)s )" % args )
        self._init(args['Name'], 'RooLinearVar')
        for (k, v) in kwargs.iteritems() : self.__setitem__(k, v)

class P2VVAngleBasis (RooObject) :
    # TODO: replace use of RooP2VVAngleBasis with an explicit product with
    #       some attribute set so we can recognize it by attribute instead
    #       of by type...
    # TODO: move 'c' out of this class (and into an explicit product),
    #       which will allow more re-use of existing objects, and hence
    #       make things faster
    def __init__(self, angles, ind,c=1,ind2 = None,**kwargs) :
        assert c!=0
        namePostFix = kwargs.pop( 'NamePostFix', '')
        if namePostFix : namePostFix += '_'
        name = 'p2vvab_%s%d%d%d%d' % ( namePostFix, ind[0], ind[1], ind[2], ind[3] )
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
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)


class AbsRealMoment( object ):
    def __init__( self, moment )  : self._var = moment
    def __getattr__( self, name ) : return getattr(self._var, name)
    def GetName( self )           : return self.basisFunc().GetName()
    def basisFunc( self )         : return self._basisFunc
    def pdf( self )               : return self._pdf
    def normSet( self )           : return self._normSet

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
    def __init__( self, BasisFunc, Norm, PDF, IntSet, NormSet ) :
        # get arguments
        self._basisFunc = BasisFunc
        self._norm      = Norm
        self._pdf       = PDF
        self._intSet    = IntSet
        self._normSet   = NormSet

        # build RooFit integration and normalisation sets
        self._rooIntSet  = ArgSet( self._basisFunc.GetName() + '_intSet',  ( var for var in self._intSet  ) )
        self._rooNormSet = ArgSet( self._basisFunc.GetName() + '_normSet', ( var for var in self._normSet ) )

        # create efficiency moment
        from P2VVLoad import P2VVLibrary
        from ROOT import RooRealEffMoment
        AbsRealMoment.__init__( self, RooRealEffMoment(  __dref__(self._basisFunc)
                                                       , self._norm
                                                       , __dref__(self._pdf)
                                                       , __dref__(self._rooIntSet)
                                                       , __dref__(self._rooNormSet)
                                                      )
                              )

class CalibratedDilution( RooObject ) :
    def __init__( self, **kwargs ) :
        __check_req_kw__( 'Name', kwargs )
        Name = kwargs.pop('Name')
        __check_name_syntax__(Name)

        if 'WTag' in kwargs :
            arguments = [ 'WTag' ]
            if 'AWTag' in kwargs : arguments.append('AWTag')

        else :
            arguments = [ 'EstWTag', 'AvgEstWTag', 'P0', 'P1' ]
            for arg in arguments : __check_req_kw__( arg, kwargs )
            if all( kw in kwargs for kw in [ 'AP0', 'AP1' ] ) : arguments += [ 'AP0', 'AP1' ]

        self._declare( 'CalibratedDilution::%s(%s)' % ( Name, ','.join( kwargs.pop(arg).GetName() for arg in arguments ) ) )
        for arg in [ 'WTag', 'AWTag', 'EstWTag', 'AvgEstWTag', 'P0', 'P1', 'AP0', 'AP1' ] :
            assert arg not in kwargs, 'CalibratedDilution: superfluous argument found: %s' % arg

        self._init( Name, 'RooCalibratedDilution' )
        for k, v in kwargs.iteritems() : self.__setitem__( k, v )

class ComplementCoef( RooObject ) :
    def __init__( self, **kwargs ) :
        __check_req_kw__( 'Name', kwargs )
        __check_req_kw__( 'Coefficients', kwargs )
        name = kwargs.pop('Name')
        __check_name_syntax__(name)

        # build a RooComplementCoef (no workspace declaration, since factory string has limited length!!!)
        from ROOT import RooArgList, RooComplementCoef
        coefList = RooArgList( __dref__(coef) for coef in kwargs.pop('Coefficients') )
        complCoef = RooComplementCoef( name, name, coefList )
        self._addObject(complCoef)
        complCoef.IsA().Destructor(complCoef)

        # initialize
        self._init( name, 'RooComplementCoef' )
        for ( k, v ) in kwargs.iteritems() : self.__setitem__( k, v )

class EfficiencyBin(RooObject):
    def __init__( self, **kwargs ) :
        __check_req_kw__('Name', kwargs)
        __check_req_kw__('Bins', kwargs)
        Name = kwargs.pop('Name')
        __check_name_syntax__(Name)
        Bins = kwargs.pop('Bins')        

        from ROOT import RooEfficiencyBin
        b = RooEfficiencyBin(Name, Name)
        for v, f in Bins.iteritems():
            b.addEntry(__dref__(v), f)
        b = self._addObject(b)
        self._init( Name, 'RooEfficiencyBin' )
        for k, v in kwargs.iteritems() : self.__setitem__( k, v )
    
class RealVar (RooObject) :
    # WARNING: multiple instances don't share proxy state at this time...
    # TODO: move common things like Name and Title in RooObject...
    # TODO: provide scaffolding in RooObject to extend getters & setters on a class-by-class basis
    _getters = {'Unit'       : lambda s : s.getUnit()
               ,'Value'      : lambda s : s.getVal()
               ,'MinMax'     : lambda s : s.getRange()
               ,'nBins'      : lambda s : s.getBins()
               ,'Error'      : lambda s : s.getError()
               }
    _setters = {'Unit'       : lambda s,v : s.setUnit(v)
               ,'Value'      : lambda s,v : s.setVal(v)
               ,'MinMax'     : lambda s,v : s.setRange(v)
               ,'Constant'   : lambda s,v : s.setConstant(v)
               ,'nBins'      : lambda s,v : s.setBins(v)
               ,'Ranges'     : lambda s,v : s.setRanges(v)
               ,'Error'      : lambda s,v : s.setError(v)
               }
    # TODO: provide a copy constructor
    def __init__(self,Name ,**kwargs):
        if 'name' in kwargs : raise RuntimeError('Please replace name argument with Name = xyz' )
        __check_name_syntax__(Name)

        self.__blind = 'Blind' in kwargs

        blindName = '__' + Name + '__' if 'Blind' in kwargs else Name
        if Name not in self.ws():
            # construct factory string on the fly...
            if not any( kw in kwargs for kw in [ 'Value', 'MinMax', 'Formula' ] ) :
                raise KeyError('%s does not exist yet, neither Value nor MinMax nor Formula specified'%Name)

            if 'Formula' not in kwargs :
                if 'Value' not in kwargs:
                    (mi,ma) = kwargs.pop('MinMax')
                    self._declare("%s[%s,%s]"%(blindName,mi,ma))
                elif 'MinMax' not in kwargs:
                    self._declare("%s[%s]"%(blindName,kwargs.pop('Value')))
                else:
                    (mi,ma) = kwargs.pop('MinMax')
                    val = kwargs.pop('Value')
                    if val < mi or val > ma : raise RuntimeError('Specified Value %s not contained in MinMax (%s,%s)' % ( val,mi,ma))
                    self.__original = self._declare("%s[%s,%s,%s]"%(blindName,val,mi,ma))
            else :
                assert 'Value' not in kwargs

                ( formula, args, data ) = kwargs.pop('Formula')
                formName = '_%s__formula'%Name
                dummyForm = self._declare( "expr::%s('%s',{%s})" % ( formName,formula,','.join( arg.GetName() for arg in args ) ) )
                self._declare( "dataobs::%s(%s,%s)" % ( blindName, data.GetName(), dummyForm.GetName() ) )

                varMin = float('+infinity')
                varMax = float('-infinity')
                for varSet in data :
                    varVal = varSet.getRealValue(blindName)
                    if varVal < varMin : varMin = varVal
                    if varVal > varMax : varMax = varVal
                varMin -= 0.0001 * abs(varMin)
                varMax += 0.0001 * abs(varMax)

                if 'MinMax' in kwargs :
                    if kwargs['MinMax'][0] > varMin or kwargs['MinMax'][1] < varMax :
                        raise RuntimeError( 'RealVar(%s): Not all values in data set "%s" (%s, %s) contained within MinMax (%s, %s)'\
                                            % (  blindName, data.GetName(), varMin, varMax, kwargs['MinMax'][0], kwargs['MinMax'][1] ) )
                else :
                    kwargs['MinMax'] = ( varMin, varMax )

                if not 'Constant' in kwargs : kwargs['Constant'] = False

            if 'Blind' in kwargs: # wrap the blinding class around us...
                b = kwargs.pop('Blind')
                _type = b[0] if type(b[0])==str else b[0].__name__
                _bs   = b[1]
                _args = b[2:]
                self._declare("%s::%s('%s',%s,%s)"%(_type,Name,_bs,','.join('%s'%i for i in _args),blindName))
                self._init(Name,_type if _type.startswith('Roo') else 'Roo' + _type)
            else :
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

    def __setitem__(self,k,v):
        from itertools import ifilter, imap
        for setters in imap( lambda x: x._setters, ifilter( lambda x : hasattr(x,'_setters'), type(self).__mro__) ) :
            if k in setters :  return setters[k](self if not self.__blind else self.__original, v)
        raise KeyError('\'%s\' is not known for class %s' % (k, type(self) ) )

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
    _getters = {'Type'        : lambda s : s._get('Type')
               ,'Parameters'  : lambda s : s._get('Parameters')
               ,'Name'        : lambda s : s._get('Name')
               ,'ConditionalObservables' : lambda s : s.ConditionalObservables()
               ,'ExternalConstraints'    : lambda s : s.ExternalConstraints()
               }
    _setters = { 'ConditionalObservables' : lambda s,v : s.setConditionalObservables(v)
               , 'ExternalConstraints'    : lambda s,v : s.setExternalConstraints(v)
               }

    ## TODO: define operators
    def __init__(self, **kwargs):
        __check_req_kw__( 'Type', kwargs )
        __check_req_kw__( 'Name', kwargs )
        __check_name_syntax__( kwargs['Name'] )

        # Save the keyword args as properties
        self._dict = kwargs
        self._make_pdf()
        for d in set(('ConditionalObservables','ExternalConstraints')).intersection( kwargs ) :
            self[d] = kwargs.pop(d)

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
                if v != self._get(k) : print k, v, self._get(k)
                assert v == self._get(k)

    def _separator(self):
        return '_'

    def ConditionalObservables(self) :
        if not hasattr(self, '_conditionals'): return set()
        return set( i for i in self.Observables() if i.GetName() in self._conditionals )
    def setConditionalObservables(self, obs ) :
        # TODO: check if we have a conditional request for something which isn't one of
        #       our observables and provide a warning...
        self._conditionals = set( o if type(o)==str else o.GetName() for o in obs )
    def ExternalConstraints(self) :
        if not hasattr(self, '_externalConstraints'): return list()
        return self._externalConstraints
    def setExternalConstraints(self, constraints ) :
        self._externalConstraints = constraints
    def GlobalObservables(self):
        if not hasattr(self,'_globalObservables') : return set()
        return set( i for i in self.Observables() if i.GetName() in self._globalObservables )
    def setGlobalObservables(self, observables ) :
        self._globalObservables = observables

    @wraps(RooAbsPdf.createNLL)
    def createNLL( self, data, **kwargs ) :
        condObs  = self.ConditionalObservables()
        if condObs :
            assert 'ConditionalObservables' not in kwargs or condObs == set(kwargs['ConditionalObservables']) , 'Inconsistent Conditional Observables'
            print 'INFO: adding ConditionalObservables: %s' % [ i.GetName() for i in  condObs ]
            kwargs['ConditionalObservables'] = condObs
        extConst = self.ExternalConstraints()
        if extConst :
            assert 'ExternalConstraints' not in kwargs or extConst== kwargs['ExternalConstraints'] , 'Inconsistent External Constraints'
            print 'INFO: adding ExternalConstraints: %s' % [ i.GetName() for i in extConst ]
            kwargs['ExternalConstraints'] = extConst
        for d in set(('ConditionalObservables','ExternalConstraints')).intersection( kwargs ) :
            kwargs[d] = RooArgSet( __dref__(var) for var in kwargs.pop(d) )
        print kwargs
        return self._var.createNLL( data, **kwargs )

    @wraps(RooAbsPdf.fitTo)
    def fitTo( self, data, **kwargs ) :
        condObs  = self.ConditionalObservables()
        if condObs :
            assert 'ConditionalObservables' not in kwargs or condObs == set(kwargs['ConditionalObservables']) , 'Inconsistent Conditional Observables'
            print 'INFO: adding ConditionalObservables: %s' % [ i.GetName() for i in  condObs ]
            kwargs['ConditionalObservables'] = condObs
        extConst = self.ExternalConstraints()
        if extConst :
            assert 'ExternalConstraints' not in kwargs or extConst== kwargs['ExternalConstraints'] , 'Inconsistent External Constraints'
            print 'INFO: adding ExternalConstraints: %s' % [ i.GetName() for i in extConst ]
            kwargs['ExternalConstraints'] = extConst
        globalObs = self.GlobalObservables()
        if globalObs:
            assert 'GlobalObservables' not in kwargs or extConst== kwargs['GlobalObservables'] , 'Inconsistent Global Observables'
            print 'INFO: adding GlobalObservables: %s' % [ i.GetName() for i in globalObs ]
            kwargs['GlobalObservables'] = globalObs
        for d in set(('ConditionalObservables','ExternalConstraints', 'GlobalObservables')).intersection( kwargs ) :
            kwargs[d] = RooArgSet( __dref__(var) for var in kwargs.pop(d) )
        print kwargs
        return self._var.fitTo( data, **kwargs )

    @wraps(RooAbsPdf.generate)
    def generate(self, whatvars, *args, **kwargs):
        #if not whatvars : whatvars = [ i for i in self._var.getVariables() if i.getAttribute('Observable') ]
        whatvars = RooArgSet([ __dref__(i) for i in whatvars] if not isinstance(__dref__(whatvars),RooAbsCategory) else __dref__(whatvars))
        conditionals = set(o.GetName() for o in self.ConditionalObservables())
        pdfVars = self._var.getVariables()
        for v in pdfVars:
            if v.GetName() in conditionals:
                v.setAttribute("GenerateConditional", True)
        data = self._var.generate(whatvars, *args,**kwargs)
        for v in pdfVars:
            if v.GetName() in conditionals:
                v.setAttribute("GenerateConditional", False)
        return data

    @wraps(RooAbsPdf.plotOn)
    def plotOn( self, frame, **kwargs ) :
        if 'Slice' in kwargs :
            sl = kwargs.pop('Slice')
            sld = __dref__(sl[0])
            if isinstance(sld, RooArgSet):
                kwargs['Slice'] = sld
            else:
                kwargs['Slice'] = ( sld, sl[1] )
        if 'Slices' in kwargs :
            sls = kwargs.pop('Slices')
            kwargs['Slices'] = [ ( __dref__(sl[0]), sl[1] ) for sl in sls ]
        return self._var.plotOn( frame, **kwargs )

class ProdPdf(Pdf):
    def __init__(self, Name, PDFs, **kwargs):
        assert 'ConditionalObservables' not in kwargs, 'ConditionalObservables are done automatically. Do not specify'
        conds = set()
        obs = set()
        ec = set()
        for pdf in PDFs:
            obs   |=  pdf.Observables() - pdf.ConditionalObservables()
            conds |=  pdf.ConditionalObservables()
            ec    |=  set( pdf.ExternalConstraints() )
        #print "Conditional Observables %s: %s" % (Name, [o['Name'] for o in conds])
        #print "Free Observables %s: %s" % (Name, [o['Name'] for o in obs])

        d = {  'PDFs'                   : frozenset(PDFs)
             , 'Name'                   : Name + '_' + self._separator().join( [ i.GetName() for i in PDFs ] )
             , 'ConditionalObservables' : list( conds - obs )
             , 'ExternalConstraints'    : list(ec)
            }
        Pdf.__init__(self, Type = 'RooProdPdf', **d)

        ##### FIXME/BUG/WORKAROUND #####
        # in RooVectorDataStore::cacheArgs, leafs are cached with a (forced) normalization
        # set corresponding to the normalization set of the top level PDF, unless the leaf
        # has the NOCacheAndTrack attribute set.
        # So here we walk through the leafs of our PDF, and add NOCacheAndTrack for all those
        # PDFs which are dependent on a conditional observable which appears in the toplevel...
        #print 'ProdPDF wrapper -- checking whether RooVectorDataStore::cacheArgs normalization bug workaround is required'
        ## create RooArgSet of conditional observables
        #cset = RooArgSet(__dref__(c) for c in conds )
        #for pdf in PDFs :
        #    for c in pdf.getComponents() :
        #        if c.getAttribute('NOCacheAndTrack') : continue
        #        depset = c.getObservables( cset )
        #        if depset.getSize() > 0 :
        #            print 'setting NOCacheAndTrack for %s as it depends on conditional observable %s ' % (c.GetName(),[i.GetName() for i in depset])
        #            c.setAttribute('NOCacheAndTrack')
        #print 'ProdPDF wrapper -- done checking for RooVectorDataStore::cacheArgs normalization bug '


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
                # parsing of this is utterly broken in case we try PROD::name( f | { x,y }, g )
                # as it interprets this as f|x ,y, g due to not recognizing the closing } as
                # the token being parsed is already split on the ',', hence asSET gets "{x" and
                # is very happy with that... and the closing } silenty disappears as well...
                __broken_parser_workaround = ArgSet( name+'_conditional_obs', cond )
                print 'P2VV - INFO: ProdPdf: Adding conditional observables [%s] to %s ' % ( ', '.join( i.GetName() for i in cond ), name )
                name += '|%s'% __broken_parser_workaround.GetName()
            return name
        # NOTE: this construction 'absorbs' all conditional observables, so the output
        #       does not have any explicit conditional observables anymore...
        pdfs = ','.join( _handleConditional(p) for p in self._dict['PDFs'])
        return 'PROD::%s(%s)' % (self._dict['Name'], pdfs)
        # TODO: check wheter the result has explicitly conditional observables,
        #       or whether they have been 'eaten' by the other PDFs in the product

    def _separator(self):
        return '_X_'

class SumPdf(Pdf):
    def __init__(self, **kwargs) :
        self._yields = {}
        pdfs = list(kwargs['PDFs'])
        co = set([ i for pdf in pdfs for i in pdf.ConditionalObservables() ])
        if 'ConditionalObservables' in kwargs :
            if co != set(kwargs['ConditionalObservables']):
                print 'WARNING: inconsistent conditional observables: %s vs %s' % ( co, kwargs['ConditionalObservables'] )
        elif co :
            kwargs['ConditionalObservables'] = list(co)

        ec = set([ i for pdf in pdfs for i in pdf.ExternalConstraints() ])
        if 'ExternalConstraints' in kwargs:
            if ec != set(kwargs['ExternalConstraints']):
                print 'WARNING: inconsistent external constraints: %s vs %s' % ( ec, kwargs['ExternalConstraints'] )
        elif ec:
            kwargs['ExternalConstraints'] = list(ec)

        diff = set([p.GetName() for p in pdfs]).symmetric_difference(set(kwargs['Yields'].keys()))
        if len(diff) not in [0, 1]:
            raise StandardError('The number of yield variables must be equal to or 1'
                                + 'less then the number of PDFs.')
        # self._dict['Name'] = self._separator().join([p.GetName() for p in pdfs])
        Pdf.__init__( self, Type = 'RooAddPdf', **kwargs )
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

class SimultaneousPdf( Pdf ) :
    def __init__( self, Name, **kwargs ) :
        args = { 'Name' : Name }
        pdfOpts = { }
        if 'States' in kwargs :
            ## pdfs = sorted([(s, pdf) for s, pdf in args['States'].iteritems()], key = lambda (s, pdf): args['Cat'].lookupType(s).getVal())
            ## pdfs = [e[1] for e in pdfs]
            simul = RooSimultaneous(Name, Name, __dref__(kwargs.pop('SplitCategory')))
            for s, pdf in kwargs.pop('States').iteritems():
                simul.addPdf(pdf, s)
            self._addObject(simul)
            simul.IsA().Destructor(simul)
        elif 'SplitParameters' in kwargs :
            args['Master']     = kwargs.pop('MasterPdf')
            args['SplitCats']  = [ kwargs.pop('SplitCategory').GetName() ] if 'SplitCategory' in kwargs\
                                 else [ ','.join( cat.GetName() for cat in cats ) for cats in kwargs.pop('SplitCategories') ]
            args['SplitPars']  = [ ','.join( par.GetName() for par in kwargs.pop('SplitParameters') ) ]\
                                 if len( args['SplitCats'] ) == 1 and not hasattr( kwargs['SplitParameters'][0], '__iter__' )\
                                 else [ ','.join( par.GetName() for par in pars ) for pars in kwargs.pop('SplitParameters') ]
            args['SplitSpecs'] = ','.join( '$SplitParam({%s},{%s})' % parsCats for parsCats in zip(args['SplitPars'], args['SplitCats']) )
            spec = 'SIMCLONE::%(Name)s(%(Master)s,%(SplitSpecs)s)' % args

            cond = args['Master'].ConditionalObservables()
            if cond : pdfOpts['ConditionalObservables'] = cond
            extCon = args['Master'].ExternalConstraints()
            if extCon : pdfOpts['ExternalConstraints' ] = extCon
            self._declare(spec)

        else :
            raise KeyError, 'P2VV - ERROR: SimultaneousdPdf: Must specify either SplitParameters or States'

        self._init( Name, 'RooSimultaneous' )
        Pdf.__init__( self , Name = Name , Type = 'RooSimultaneous', **pdfOpts )

        for ( k, v ) in kwargs.iteritems() : self.__setitem__( k, v )

    def _make_pdf(self) : pass

class RealSumPdf( Pdf ):
    def __init__( self, name, functions, **kwargs ) :
        # get the name of the PDF, its functions and its coefficients
        __check_name_syntax__(name)
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
        exCon = d['Original'].ExternalConstraints()
        if exCon : extraOpts['ExternalConstraints' ] = exCon
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
        __check_name_syntax__(Name)
        _hist = kwargs.pop('Histogram')
        if str(type(_hist)).find('TH1') == -1:
            raise TypeError, "HistFunc can only handle 1D historgrams"
        _dn = Name + '_data_hist'
        # Create Datahist and Import with density set to false
        _data = RooDataHist(_dn, _dn, RooArgList(*[__dref__(o) for o in kwargs['Observables']]),
                            RooFit.Import(_hist, False))
        self.ws().put(_data)
        self._declare('RooHistFunc::%s({%s}, %s)' % (Name, ','.join([o.GetName() for o in kwargs.pop('Observables')]), _data.GetName()))
        self._init(Name, 'RooHistFunc')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

class Projection(Pdf):
    def __init__(self,**kwargs) :
        name = kwargs.pop('Name')
        original = kwargs.pop('Original')
        projectVars = kwargs.pop('ProjectVars')

        # construct factory string on the fly...
        projection = original.createProjection(RooArgSet(*list(projectVars)))
        projection.SetName(name)
        self._addObject(projection)
        self._init(name, 'RooProjectedPdf')
        extraOpts = dict()
        cond =  original.ConditionalObservables()
        if cond : extraOpts['ConditionalObservables'] = cond
        exCon = original.ExternalConstraints()
        if exCon : extraOpts['ExternalConstraints' ] = exCon
        Pdf.__init__(self , Name = name , Type = 'RooProjectedPdf', **extraOpts)
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

    def _make_pdf(self):
        pass

class KeysPdf(Pdf):
    def __init__(self,**kwargs) :
        name = kwargs.pop('Name')
        observable = kwargs.pop('Observable')
        data = kwargs.pop('Data')
        rho = kwargs.pop('Rho', 1.)

        from ROOT import RooKeysPdf
        keysPdf = RooKeysPdf(name, name, __dref__(observable), data, RooKeysPdf.NoMirror, rho)

        keysPdf = self._addObject(keysPdf)
        self._init(name, 'RooKeysPdf')
        Pdf.__init__(self , Name = name , Type = 'RooKeysPdf')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)

    def _make_pdf(self):
        pass

class MultiHistEfficiencyModel(Pdf):
    def _make_pdf(self) : pass
    def __init__(self, **kwargs):

        ## jpsiphi_1TrackMuon,1TrackAllL0: 8.285e-01 +- 1.633e-02
        ## jpsiphi_bdi2                  : 8.330e-01 +- 1.402e-02

        self.__pdf_name = kwargs.pop('Name')
        self.__original = kwargs.pop('Original')
        self.__resolution_model = kwargs.pop('ResolutionModel')
        self.__bins = kwargs.pop('Bins')
        self.__fit_bins = kwargs.pop('FitBins', True)
        relative = kwargs.pop('Relative')
        self.__observable = kwargs.pop('Observable')
        self.__cc = kwargs.pop('ConditionalCategories', False)
        self.__conditionals = self.__original.ConditionalObservables() | \
                              set(kwargs.pop('ConditionalObservables', [])) | \
                              self.__resolution_model.ConditionalObservables()
        self.__fit_bins = kwargs.pop('FitAcceptance', True)
        self.__use_bin_constraint = kwargs.pop('UseSingleBinConstraint', True) and self.__fit_bins
        self.__heights = {}
        self.__shapes = []
        self.__coefficients = {}
        self.__base_bounds = None
        self.__constraints = []
        self.__binHeightMinMax = kwargs.pop('BinHeightMinMax', None)
        if not self.__binHeightMinMax : self.__binHeightMinMax = (0.01, 0.999)

        from copy import copy
        from ROOT import RooBinning        
        from ROOT import RooArgList

        obs_set = RooArgSet()
        for o in self.__original.Observables().intersection(self.__original.ConditionalObservables()):
            obs_set.add(__dref__(o))

        for (category, entries) in self.__bins.iteritems():
            states = set([s.GetName() for s in category])
            coef_info = {}

            for state, state_info in [(s, entries[s]) for s in states if s in entries]:
                heights = state_info['heights']
                heights = [RealVar('%s_%s_bin_%03d' % (category.GetName(), state, i + 1),
                                   Observable = False, Value = v,
                                   MinMax = self.__binHeightMinMax) for i, v in enumerate(heights)]
                self.__heights[(category, state)] = heights
                # Add a binning for this category and state
                bounds = state_info['bins']
                binning_name = '_'.join([category.GetName(), state, 'binning'])
                shape_binning = RooBinning(len(bounds) - 1, bounds)
                shape_binning.SetName(binning_name)
                self.__observable.setBinning(shape_binning, binning_name)

                # Calculate a first scale factor according to the constraint
                if not self.__fit_bins:
                    # If we're not fitting set all bins constant
                    for h in heights:
                        h.setConstant(True)
                elif len(heights) == 1 or self.__use_bin_constraint:
                    # Fix the bin if there is only one.
                    heights[0].setConstant(True)
                else:
                    # We're fitting and using the average constraint, first
                    # create a shape and then the constraint.
                    shape_name = '_'.join([category.GetName(), state, 'shape'])
                    shape = BinnedPdf(shape_name, Observable = self.__observable,
                                      Binning = binning_name, Coefficients = heights)
                    shape.setForceUnitIntegral(True)

                    eff_model = EffResModel(Efficiency = shape,
                                            ResolutionModel = self.__resolution_model)

                    # Set all observables constant for the shape to work around a
                    # RooFit limitation
                    observables = eff_model.getObservables(obs_set)
                    for o in observables:
                        o.setConstant(True)

                    self.__shapes.append(shape)

                    heights_list = RooArgList(__dref__(h) for h in heights)
                    av_name = "%s_%s_average" % (category, state)
                    value, error = state_info['average']
                    mean  = RealVar(Name = av_name + '_constraint_mean', Value = value, Constant = True)
                    sigma = RealVar(Name = av_name + '_constraint_sigma', Value = error, Constant = True)
                    from ROOT import RooAvEffConstraint
                    constraint = Pdf(Name = av_name + '_constraint', Type = RooAvEffConstraint,
                                     Parameters = [self.__original, eff_model, mean, sigma])
                    self.__constraints.append(constraint)
                if not self.__base_bounds or len(bounds) > len(self.__base_bounds):
                    self.__base_bounds = bounds
                    self.__base_binning = shape_binning
                coef_info[state] = copy(state_info)
                coef_info[state].update({'heights' : heights})
            self.__coefficients[category] = coef_info

        # Set the binning on the observable
        self.__observable.setBinning(self.__base_binning, self.__base_binning.GetName())

        # Build relative efficiencies
        self.__relative_efficiencies = {}
        remaining = None

        for categories, re in relative.iteritems():
            state_name = '__'.join(['%s_%s' % (c.GetName(), s) for c, s in categories])

            # Make realvars for relative efficiencies
            if re != None:
                efficiency = RealVar('%s_efficiency' % state_name, Observable = False,
                                     **re)
                self.__relative_efficiencies[state_name] = efficiency
            elif remaining == None:
                remaining = state_name
            else:
                raise RuntimeError("More than one relative efficiency is None")
        # FIXME: perhaps this should be a dedicated class too
        form = '-'.join(['1'] + [e.GetName() for e in self.__relative_efficiencies.itervalues()])
        self.__relative_efficiencies[remaining] = FormulaVar("remaining_efficiency", form, self.__relative_efficiencies.values())

        efficiency_entries = self.__build_shapes(relative)

        from ROOT import RooMultiEffResModel
        mhe = RooMultiEffResModel(self.__pdf_name, self.__pdf_name, efficiency_entries)
        self._addObject(mhe)
        mhe.IsA().Destructor(mhe)

        extraOpts = dict()
        if self.__conditionals:
            extraOpts['ConditionalObservables'] = self.__conditionals

        if self.__fit_bins and not self.__use_bin_constraint:
            self.__add_constraints()

        constraints = list(set(kwargs.pop('ExternalConstraints', [])) | \
                           set(self.__original.ExternalConstraints()) | \
                           set(self.__resolution_model.ExternalConstraints()) | \
                           set(self.ExternalConstraints()))
        if constraints : extraOpts['ExternalConstraints' ] = constraints
        self._init(self.__pdf_name, 'RooMultiEffResModel')
        Pdf.__init__(self, Name = self.__pdf_name , Type = 'RooMultiEffResModel', **extraOpts)

    def binning(self):
        return self.__base_binning

    def bounds(self):
        return self.__base_bounds

    def shapes(self):
        return self.__shapes

    def heights(self):
        return self.__heights

    def __build_shapes(self, relative):
        from ROOT import std
        from ROOT import MultiHistEntry

        efficiency_entries = std.vector('MultiHistEntry*')()

        for categories, relative_efficiency in relative.iteritems():
            # Make EfficiencyBins for the bin values
            heights = []
            bin_vars = [{} for i in range(len(self.__base_bounds) - 1)]
            state_name = '__'.join(['%s_%s' % (c.GetName(), s) for c, s in categories])
            prefix = self.__pdf_name + '_' + state_name
            for category, state in categories:
                if self.__cc: self.__conditionals.add(category)
                category_info = self.__coefficients[category]
                states = [s.GetName() for s in category if s.GetName() in category_info]
                if len(states) == 1:
                    category_heights = category_info[states[0]]['heights']
                    category_bounds  = category_info[states[0]]['bins']
                elif len(states) > 1:
                    category_heights = category_info[state]['heights']
                    category_bounds  = category_info[state]['bins']
                else:
                    raise ValueError("Number of states must not be 0")
                for i in range(len(self.__base_bounds) - 1):
                    val = (self.__base_bounds[i] + self.__base_bounds[i + 1]) / 2
                    coefficient = self.__find_coefficient(val, category_bounds, category_heights)
                    bin_vars[i][__dref__(coefficient)] = state in category_info
            for i, d in enumerate(bin_vars):
                name = '%s_%d' % (prefix, i)
                heights.append(EfficiencyBin(Name = name, Bins = d))

            # BinnedPdf for the shape
            binned_pdf = BinnedPdf(Name = '%s_shape' % prefix, Observable = self.__observable,
                                   Binning = self.__base_binning.GetName(), Coefficients = heights)
            # note: constant optimization WILL evaluate RooBinnedPdf as a PDF, and thus normalize it...
            binned_pdf.setForceUnitIntegral(True)

            # EffResModel to combine shape with PDF
            eff_model = EffResModel(Name = '%s_efficiency' % prefix,
                                    ResolutionModel = self.__resolution_model, Efficiency = binned_pdf)

            # MultiHistEntry
            category_map = std.map('RooAbsCategory*', 'string')
            category_pair = std.pair('RooAbsCategory*', 'string')
            cm = category_map()
            for category, state in categories:
                # cp = category_pair(__dref__(category), state)
                cm[__dref__(category)] = state
            efficiency = self.__relative_efficiencies[state_name]
            entry = MultiHistEntry(cm, __dref__(eff_model), __dref__(efficiency))
            efficiency_entries.push_back(entry)
        return efficiency_entries

    def __add_constraints(self):
        self.setExternalConstraints(self.__constraints)
        
    def __find_coefficient(self, val, bounds, coefficients):
        for i in range(len(bounds) - 1):
            if val > bounds[i] and val < bounds[i + 1]:
                break
        else:
            raise RuntimeError;
        return coefficients[i]


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
        __check_name_syntax__(Name)
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
        cstr = lambda arg : arg if type(arg) == str else arg.GetName() if hasattr(arg,'GetName') else str(arg)
        for k in d.keys() : d[k] = cstr(d[k])
        if type(d['resolutionModel'])!=str : d['resolutionModel'] = d['resolutionModel'].GetName()
        self._declare("BDecay::%(name)s( %(time)s, %(tau)s, %(dGamma)s, "\
                                          " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s, "\
                                          " %(dm)s, %(resolutionModel)s, %(decayType)s  )" % d )
        self._init(Name,'RooBDecay')
        for (k,v) in kwargs.iteritems() :
            self.__setitem__(k,v)


class BTagDecay( Pdf ) :
    def _make_pdf(self) : pass
    def __init__( self, Name, **kwargs ) :
        from P2VVLoad import P2VVLibrary
        argDict = { 'Name' : Name, 'checkVars' : '1', 'decayType' : 'SingleSided' }

        cstr = lambda arg : arg if type(arg) == str else arg.GetName() if hasattr(arg,'GetName') else str(arg)
        convert = lambda arg : cstr(arg) if type(arg) != list else '{%s}' % ','.join( str(listItem) for listItem in arg )
        if 'tagCat0' in kwargs and 'tagCat1' in kwargs :
            # two tagging categories
            for argName in [  'time', 'iTag0', 'iTag1', 'tagCat0', 'tagCat1', 'tau', 'dGamma', 'dm'
                            , 'dilutions0', 'dilutions1', 'ADilWTags0', 'ADilWTags1'
                            , 'coshCoef', 'sinhCoef', 'cosCoef', 'sinCoef'
                            , 'resolutionModel', 'decayType', 'checkVars'
                           ] :
                if argName not in argDict or argName in kwargs : argDict[argName] = convert( kwargs.pop(argName) )

            # put tagging category coefficients and average even and odd coefficients in TObjArrays of RooArgLists
            _a  = lambda x : self.ws()[ argDict[ x ] ]
            _l2 = lambda x : self.ws().factory().asLIST(  x  )
            _l  = lambda x : _l2( argDict[x] ) 
            def _oa( x ) : 
                from ROOT import TObjArray
                z = TObjArray()
                for i in kwargs.pop(x) : z.Add( _l2( convert( i ) ) )
                return z

            from ROOT import RooBTagDecay
            decay = RooBTagDecay(  argDict['Name'], argDict['Name'], _a('time')
                                 , _a('iTag0'), _a('iTag1')
                                 , _a('tagCat0'), _a('tagCat1')
                                 , _a('tau'), _a('dGamma'), _a('dm')
                                 , _l('dilutions0'), _l('dilutions1'), _l('ADilWTags0'), _l('ADilWTags1')
                                 , _oa( 'avgCEvens' ),  _oa( 'avgCOdds' ), _oa( 'tagCatCoefs' )
                                 , _a('coshCoef'), _a('sinhCoef')
                                 , _a('cosCoef'), _a('sinCoef')
                                 , _a('resolutionModel')
                                 , 1 if argDict['decayType'] == 'DoubleSided' else 2 if argDict['decayType'] == 'Flipped' else 0
                                 , int( argDict['checkVars'] )
                                )
            self._addObject(decay)
            decay.IsA().Destructor(decay)

        elif 'tagCat' in kwargs :
            # one tagging category
            for argName in [  'time', 'iTag', 'tagCat', 'tau', 'dGamma', 'dm'
                            , 'dilutions', 'ADilWTags', 'avgCEvens', 'avgCOdds', 'tagCatCoefs'
                            , 'coshCoef', 'sinhCoef', 'cosCoef', 'sinCoef'
                            , 'resolutionModel', 'decayType', 'checkVars'
                           ] :
                if argName not in argDict or argName in kwargs : argDict[argName] = convert( kwargs.pop(argName) )

            self._declare("BTagDecay::%(Name)s( %(time)s, %(iTag)s, %(tagCat)s, %(tau)s, %(dGamma)s, %(dm)s,"\
                                              " %(dilutions)s, %(ADilWTags)s, %(avgCEvens)s, %(avgCOdds)s, %(tagCatCoefs)s,"\
                                              " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s,"\
                                              " %(resolutionModel)s, %(decayType)s, %(checkVars)s )" % argDict
                         )

        else :
            # no tagging categories
            for argName in [  'time', 'iTag', 'tau', 'dGamma', 'dm'
                            , 'dilution', 'ADilWTag', 'avgCEven', 'avgCOdd'
                            , 'coshCoef', 'sinhCoef', 'cosCoef', 'sinCoef'
                            , 'resolutionModel', 'decayType', 'checkVars'
                           ] :
                if argName not in argDict or argName in kwargs : argDict[argName] = convert( kwargs.pop(argName) )

            self._declare("BTagDecay::%(Name)s( %(time)s, %(iTag)s, %(tau)s, %(dGamma)s, %(dm)s,"\
                                              " %(dilution)s, %(ADilWTag)s, %(avgCEven)s, %(avgCOdd)s,"\
                                              " %(coshCoef)s, %(sinhCoef)s, %(cosCoef)s, %(sinCoef)s,"\
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
	# !!! Since the workspace factory doesn't know about RooBinnedPdf and its constructors, the default approach doesn't seem to
	# !!! work very well. We create the object directly and then add it to RooObject and the workspace.
        from P2VVLoad import P2VVLibrary
        argDict = { 'Name' : Name }

        # declare PDF in workspace
        if 'Category' in kwargs :
            # single category dependence
            argDict['cat']   = str(kwargs.pop('Category'))
            argDict['coefs'] = '{%s}' % ','.join( str(listItem) for listItem in kwargs.pop('Coefficients') )
            self._declare( "BinnedPdf::%(Name)s(%(cat)s, %(coefs)s)" % argDict )
            bPdf = None

        elif 'Categories' in kwargs :
            # multiple category dependence

            # build list of base categories
            from ROOT import RooArgList
            categories = kwargs.pop('Categories')
            varList = RooArgList(__dref__(var) for var in categories)

            if hasattr( kwargs['Coefficients'][0], '__iter__' ) :
                # coefficients for different variables factorize

                # build coefficients lists
                assert len(kwargs['Coefficients']) == len(categories),\
                        'P2VV - ERROR: BinnedPdf: number of specified coefficient lists is not equal to the number of specified variables'
                from ROOT import TObjArray, RooArgList
                coefLists = TObjArray()
                for coefficients in kwargs.pop('Coefficients') :
                    coefLists.Add(RooArgList( __dref__(coef) for coef in coefficients ))

                from ROOT import RooBinnedPdf
                bPdf = RooBinnedPdf( argDict['Name'], argDict['Name'], varList, coefLists, int( kwargs.pop( 'IgnoreFirstBin', 0 ) ) )

            else :
                # coefficients for different variables don't factorize

                # build coefficients list
                from ROOT import RooArgList
                coefList = RooArgList( __dref__(coef) for coef in kwargs.pop('Coefficients') ) 

                from ROOT import RooBinnedPdf
                bPdf = RooBinnedPdf( argDict['Name'], argDict['Name'], varList, coefList )

        elif 'Observable' in kwargs :
            # single continuous variable dependence
            var = kwargs.pop('Observable')
            binning = kwargs.pop('Binning')
            binning = var.getBinning(binning).GetName() if type(binning) == str else binning.GetName()

            if 'Function' in kwargs :
                # bin coefficients are given by a function
                from ROOT import RooBinnedPdf
                bPdf = RooBinnedPdf( argDict['Name'], argDict['Name'], __dref__(var), binning, __dref__( kwargs.pop('Function') ) )

            else:
                # independent bin coefficients are specified

                # build coefficients list
                from ROOT import RooArgList
                coefList = RooArgList( __dref__(coef) for coef in kwargs.pop('Coefficients') )

                from ROOT import RooBinnedPdf
                bPdf = RooBinnedPdf( argDict['Name'], argDict['Name'], __dref__(var), binning, coefList
                                    , int( kwargs.pop( 'BinIntegralCoefs', 0 ) ) )

        elif 'Observables' in kwargs :
            # multiple continuous variable dependence

            # build list of base variables
            from ROOT import RooArgList
            observables = kwargs.pop('Observables')
            varList = RooArgList(__dref__(var) for var in observables )

            # build list of binning names
            assert len(kwargs['Binnings']) == len(observables),\
                    'P2VV - ERROR: BinnedPdf: number of specified binnings is not equal to the number of specified variables'
            from ROOT import TObjArray, TObjString
            binningList = TObjArray()
            for binning, var in zip( kwargs.pop('Binnings'), observables ) :
                binningList.Add( TObjString( var.getBinning(binning).GetName() if type(binning) == str else binning.GetName() ) )

            if 'Function' in kwargs:
                # bin coefficients are given by a function
                from ROOT import RooBinnedPdf
                bPdf = RooBinnedPdf( argDict['Name'], argDict['Name'], varList, binningList, __dref__( kwargs.pop('Function') ) )

            else:
                # independent bin coefficients are specified
                if kwargs['Coefficients'] and hasattr( kwargs['Coefficients'][0], '__iter__' ) :
                    # coefficients for different variables factorize

                    # build coefficients lists
                    assert len(kwargs['Coefficients']) == len(observables),\
                            'P2VV - ERROR: BinnedPdf: number of specified coefficient lists is not equal to the number of specified variables'
                    from ROOT import TObjArray, RooArgList
                    coefLists = TObjArray()
                    for coefficients in kwargs.pop('Coefficients') :
                        coefLists.Add(RooArgList(__dref__(coef) for coef in coefficients))

                    from ROOT import RooBinnedPdf
                    bPdf = RooBinnedPdf( argDict['Name'], argDict['Name'], varList, binningList, coefLists
                                        , kwargs.pop( 'BinIntegralCoefs', 0 ), int( kwargs.pop( 'IgnoreFirstBin', 0 ) ) )

                else :
                    # coefficients for different variables don't factorize

                    # build coefficients list
                    from ROOT import RooArgList
                    coefList = RooArgList(__dref__(coef) for coef in kwargs.pop('Coefficients') )

                    from ROOT import RooBinnedPdf
                    bPdf = RooBinnedPdf( argDict['Name'], argDict['Name'], varList, binningList, coefList
                                        , int( kwargs.pop( 'BinIntegralCoefs', 0 ) ) )

        else :
            raise KeyError('P2VV - ERROR: BinnedPdf: please specify variable(s)')

        if bPdf :
            # import the BinnedPdf in the workspace
            self._addObject(bPdf)
            bPdf.IsA().Destructor(bPdf)

        # initialize PDF
        self._init( Name, 'RooBinnedPdf' )
        Pdf.__init__(  self
                     , Name = Name
                     , Type = 'RooBinnedPdf'
                    )
        for ( k, v ) in kwargs.iteritems() : self.__setitem__( k, v )

    def _make_pdf(self) : pass


class Customizer(Pdf) :
    def __init__( self, **kwargs ) :
        __check_req_kw__( 'Pdf',            kwargs )
        __check_req_kw__( 'OriginalArgs',   kwargs )
        __check_req_kw__( 'SubstituteArgs', kwargs )
        pdf     = kwargs.pop('Pdf')
        origSet = kwargs.pop('OriginalArgs')
        subsSet = kwargs.pop('SubstituteArgs')

        byName = kwargs.pop( 'ReplaceByName', False )
        if byName : __check_req_kw__( 'ArgumentSuffix', kwargs )
        argSuff  = kwargs.pop( 'ArgumentSuffix', '' )
        nameSuff = kwargs.pop( 'NameSuffix', argSuff )

        from ROOT import RooCustomizer
        customizer = RooCustomizer( __dref__(pdf), nameSuff )

        if byName :
            origDict = dict( ( origItem.GetName(), origItem ) for origItem in origSet )
            subsDict = dict( ( subsItem.GetName(), subsItem ) for subsItem in subsSet )
            for item in origSet:
                customizer.replaceArg( __dref__(item), __dref__( subsDict[ item.GetName() + argSuff ] ) )
        else :
            for origItem, subsItem in zip( origSet, subsSet ) :
                customizer.replaceArg( __dref__(origItem), __dref__(subsItem) )

        custom = customizer.build()

        self._addObject(custom)
        self._init( custom.GetName(), pdf.Type() )
        Pdf.__init__(  self
                     , Name = custom.GetName()
                     , Type = pdf.Type()
                     , ConditionalObservables = pdf.ConditionalObservables()
                     , ExternalConstraints = pdf.ExternalConstraints()
                    )

    def _make_pdf(self) : pass


class ResolutionModel(Pdf):
    def __init__(self, **kwargs):
        if type(kwargs['Type']) != str : kwargs['Type'] = kwargs['Type'].__name__
        Pdf.__init__(self,**kwargs)

class AddModel(ResolutionModel) :
    def __init__(self, name, **kwargs) :
        #TODO: forward conditionalObservables and ExternalConstraints from dependents...
        # construct factory string on the fly...
        __check_name_syntax__(name)

        models = kwargs.pop('Models')
        conditionals = set()
        externals = set()
        for model in models:
            conditionals |= model.ConditionalObservables()
            externals |= set(model.ExternalConstraints())
        ResolutionModel.__init__(self, Name = name, Type = 'RooAddModel',
                                 Models = models, Fractions = kwargs.pop('Fractions'),
                                 ConditionalObservables = conditionals,
                                 ExternalConstraints = list(externals))

    def _make_pdf(self):
        if self._dict['Name'] not in self.ws():
            self._declare(self._makeRecipe())
            self._init(self._dict['Name'], 'RooAddModel')

            # Change self._dict into attributes. Cannot be done before since the
            # underlying object does only exists at this point.
            for k, v in self._dict.iteritems():
                attr = '_' + k.lower()
                setattr(self._target_(), attr, v)
        else:
            self._init(self._dict['Name'], 'RooAddModel')
            # Make sure we are the same as last time
            for k, v in self._dict.iteritems():
                assert v == self._get(k)

    def _makeRecipe(self):
        models = self._dict['Models']
        fractions = self._dict['Fractions']
        return "AddModel::%s({%s},{%s})"%(self._dict['Name'],','.join(i.GetName() for i in models),','.join(j.GetName() for j in fractions) ) 

class EffResModel(ResolutionModel) :
    def __init__(self,**kwargs) :
        #TODO: forward conditionalObservables and ExternalConstraints from dependents...
        # construct factory string on the fly...
        self.__eff = kwargs.pop('Efficiency')
        self.__res = kwargs.pop('ResolutionModel')
        name = kwargs.pop('Name', '%s_x_%s' % ( self.__eff.GetName() , self.__res.GetName() ) )
        __check_name_syntax__(name)

        conditionals = self.__eff.ConditionalObservables() | self.__res.ConditionalObservables()
        externals = list(set(self.__eff.ExternalConstraints()) | set(self.__res.ExternalConstraints()))
        ResolutionModel.__init__(self, Name = name, Type = 'RooEffResModel',
                                 Parameters = [self.__res, self.__eff],
                                 ConditionalObservables = conditionals,
                                 ExternalConstraints = externals)

        def efficiency(self):
            return self.__eff

class Component(object):
    _d = {}
    def __init__(self,Name,*args,**kw) :
        # TODO: make things singletons, indexed by 'Name'
        if Name in Component._d : raise KeyError('Name %s is not unique' % Name)
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
    def getYield(self):
        if 'Yield' in Component._d[self.name] : return Component._d[self.name]['Yield']
        else : return None
    def setYield(self, *args):
        y = None
        if len(args) == 1 and type(args[0]) == RealVar:
            y = args[0]
        else :
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
                obs = [o for o in p.Observables() if o not in p.ConditionalObservables()]
                self[ obs ] = p

    def __setitem__(self, observable, pdf) :
        from ROOT import RooAbsCategory
        #TODO: Need to deal with fact the categories have iterators
        if not hasattr(observable,'__iter__') or isinstance(observable._var if hasattr(observable,'_var') else observable ,RooAbsCategory) : observable = (observable,)

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
            pdfs = [self[i] for i in terms if type(self[i])!=type(None)]
            externalConstraints = set([ec for pdf in pdfs for ec in pdf.ExternalConstraints()])
            d[nk] = ProdPdf(self.name, PDFs = pdfs, ExternalConstraints = externalConstraints )
        return d[k]

def buildPdf( Components, Observables, Name ) :
    # multiply PDFs for observables (for each component)

    # get observables
    if not Observables : raise RuntimeError('no Observables??')
    obsList = [ obs if type(obs) == str else obs.GetName() for obs in Observables ]

    # loop over components
    args = { 'Yields' : {}, 'PDFs' : [], 'ExternalConstraints' : set() }
    for comp in Components:
        # build PDF
        pdf = comp[obsList]
        if len(Components) > 1 : args['Yields'][pdf.GetName()] = comp['Yield']
        args['PDFs'].append(pdf)

        # add external constraints
        for ec in pdf.ExternalConstraints(): args['ExternalConstraints'].add(ec)

    args['ExternalConstraints'] = list(args['ExternalConstraints'])

    # return product directly if PDF consists of only one component
    if len(Components) == 1 : return args['PDFs'][0] # TODO: how to change the name?

    # add sum components (inputs should already be extended)
    return SumPdf( Name = Name, **args )

def buildSimultaneousPdf(Components, Observables, Spec, Name) :
    # multiply PDFs for observables (for each component)
    if not Observables : raise RuntimeError('no Observables??')
    if not Spec : raise RuntimeError('no Spec??')

    if len(Components)==1 :
        raise RuntimeError("Not implemented")

    obs = [o if type(o)==str else o.GetName() for o in Observables]
    if len(Spec) != 1:
        raise ValueError, "Cannot handle more than one split category."
    states = {}
    key = Spec.keys()[0]
    observables, split_cat = key
    s_obs = set(observables)
    s_Obs = set(Observables)
    # Either full overlap or no overlap
    assert(s_obs.issubset(s_Obs) or s_obs.isdisjoint(s_Obs))

    # - Figure out which components to clone
    # - Create new components where needed (don't forget yields)
    # - Make appropriate PDFs using new compoments when needed
    # - Make simultaneous PDF
    for state, split_def in Spec[key].iteritems():
        args = {'Yields' : {}, 'PDFs' : []}
        suffix = '_'.join((split_cat['Name'], state))
        for c in Components:
            y = c['Yield']
            name = c['Name']
            rest = list(set(Observables) - set(observables))
            if c in split_def:
                if 'Yield' in split_def[c]:
                    yld = split_def[c]['Yield']
                else:
                    yld = [y['Value']] + list(y['MinMax'])
                if s_obs.issubset(s_Obs) and 'PDF' in split_def[c]:
                    pdfs = (split_def[c]['PDF'], c[rest])
                else:
                    pdfs = c[obs],
                comp = Component('_'.join((name, suffix)), pdfs, Yield = yld)
            else:
                if s_obs.issubset(s_Obs):
                    pdfs = (c[observables], c[rest])
                else:
                    pdfs = c[obs],
                comp = Component('_'.join((name, suffix)), pdfs,
                                 Yield = [y['Value']] + list(y['MinMax']))
            pdf = comp[obs]
            args['Yields'][pdf.GetName()] = comp['Yield']
            args['PDFs'].append(pdf)
        states[state] = SumPdf(Name = '_'.join((Name, suffix)) , **args)
    ## return states
    return SimultaneousPdf(Name, SplitCategory = split_cat, States = states)
