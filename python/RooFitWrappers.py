from RooFitDecorators import *
from copy import copy

class RooObject( object ) :
    _ws = None
    _dict = None

    def ws(self) : 
        if not RooObject._ws : raise RuntimeError('No workspace defined!')
        return RooObject._ws
    def setWorkspace(self,ws) : RooObject._ws = ws
    def _declare(self,spec) :
        x = self.ws().factory(spec)
        if not x : raise NameError('failed to _declare %s to workspace factory'%spec)
        x._observable = False
        return x
    def create(self,spec) : return self._declare(spec)
    def _init(self,name,type) :
        x = self.ws()[name]
        if not x : raise KeyError('could not locate %s'%name)
        if not x.InheritsFrom(type): 
            raise KeyError('%s is %s, not %s'%(name,x.ClassName(),type))
        self._var = x
    def __getattr__( self, name ):
        return getattr( self._target_(), name )      
    def _target_( self ) :
        return self._var

    def __cmp__( self, other ):
        o = other if type( other ) == str else other.GetName()
        return self.GetName() < other

    def __eq__( self, other ):
        o = other if type( other ) == str else other.GetName()
        return self.GetName() == other

    def __ne__( self, other ):
        o = other if type( other ) == str else other.GetName()
        return self.GetName() != other

    def __hash__( self ):
        return self.GetName().__hash__()

    def __str__( self ):
        return self.GetName()
        
    ## FIXME: Should these be in RooObject??
    def observable(self) : 
        return self._var._observable
    def setObservable(self, observable) :
        self._var._observable = bool( observable )

# TODO: make this more of a 'borg' by overloading __new__ instead of __init__
#       otherwise properties of the proxy not in the 'target' are not shared
#       across multiple instances which all defer to the same 'target'
# NOTE: we push property values into the underlying target to share them, so this is OK ;-)
# NOTE: Actually, not really since not all information we want to store has a c++ counterpart

# TODO?: Instead of using GetName everywhere, perhaps an appropriate __hash__ function is more
##       elegant?

class Category (RooObject): 
    def __init__(self,name) :
        self._init(name,'RooCategory')

class FormulaVar (RooObject): 
    def __init__(self,name) :
        self._init(name,'RooFormulaVar')

class RealVar (RooObject): 
    # WARNING: multiple instances don't share proxy state at this time...
    _setters = { 'Observable' : lambda s,v : s.setObservable(v) 
               , 'Unit'       : lambda s,v : s.setUnit(v) 
               , 'Value'      : lambda s,v : s.setVal(v)
               , 'MinMax'     : lambda s,v : s.setRange(v)
               }
    _getters = { 'Observable' : lambda s : s.observable() 
               , 'Unit'       : lambda s : s.getUnit() 
               , 'Value'      : lambda s : s.getVal()
               , 'MinMax'     : lambda s : s.getRange()
               }

    def __init__(self,name,**kwargs) :
        # TODO: add blinding support to kwargs
        if name not in self.ws() : 
            # construct factory string on the fly...
            if 'Value'  not in kwargs  and 'MinMax' not in kwargs : 
                raise KeyError('%s does not exist yet, neither Value nor MinMax specified'%name)
            if 'Value' not in kwargs : 
                (mi,ma) = kwargs.pop('MinMax')
                self._declare("%s[%s,%s]"%(name,mi,ma))
            elif 'MinMax' not in kwargs : 
                self._declare("%s[%s]"%(name,kwargs.pop('Value')))
            else :
                (mi,ma) = kwargs.pop('MinMax')
                self._declare("%s[%s,%s,%s]"%(name,kwargs.pop('Value'),mi,ma))
            self._init(name,'RooRealVar')
        for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
    def __setitem__(self,k,v) :
        return RealVar._setters[k]( self, v )
    def __getitem__(self,k) :
        return RealVar._getters[k]( self )
    # overrule RooRealVar.setRange
    def setRange(self, v) :
        (mi,ma) = v
        self._target_().setRange( mi,ma )
        if self.getVal() < mi or self.getVal() > ma : self.setVal( 0.5*(ma+mi)  )


##TODO, factor out common code in Pdf and ResolutionModel

class Pdf( RooObject ):
    ## TODO: define operators
    def __init__( self, name, *args, **kwargs ):
        if 'Observables' not in kwargs:
            raise KeyError( 'Must provide observables %s' % kwargs )
        if 'Type' not in kwargs:
            raise KeyError( 'Must provide type %s' % kwargs )

        # Save the keyword args as properties
        self._dict = {}
        for k, v in kwargs:
            self._dict[ k.lower() ] = v
        self._dict[ 'observables' ] = frozenset( o )

        o = self.__getitem__( 'observables' )
        self._dict[ 'name' ]  = self._separator().join( [ name ] + [ i.GetName() for i in o ] )
        self.__make_pdf()

    def __str__( self ):
        s = '%s, %d\n' % ( self._name, self.single() )
        for o, sp in self._pdfs.iteritems():
            n = sp.GetName() if sp == self else sp.__str__()
            o = '{' + ','.join( [ i.GetName() for i in o ] ) + '}'
            s += '%s : %s\n' % ( o, n )
        return s

    def __getitem__( self, k ):
        if k in self._dict:
            return self._dict[ k ]
        else:
            return None

    def __make_pdf( self ):
        if self._dict[ 'name' ] not in self.ws():
            v = list( self._dict[ 'observables' ] )
            if 'parameters' in self._dict:
                v += list( self._dict[ 'parameters' ] )
            if 'resolutionmodel' in self._dict:
                rm = self._dict[ 'resolutionmodel' ]
                if rm.observables() != frozenset( observables ):
                    raise StandardError ( 'Error, resolution model must have the same '
                                          + 'observables as PDF' )
                ## TODO add conditional observables
                v += [ rm ]
            if 'options' in self._dict:
                v += list( self._dict[ 'options' ] )
            self._declare( self._makeRecipe() )
            self._init( self._dict[ 'name' ], 'RooAbsPdf' )
        else:
            raise StandardError( "Recreating the same Pdf is not supported atm" )

    def _separator( self ):
        return '_'

    def _makeRecipe( self ):
        if 'variables' not in self._dict:
            raise StandardError( "This method requires a single 'variables' to be set" )
        k = self._dict[ 'variables' ]
        deps = ','.join( [ v.GetName() if type( v ) != str else v for v in k ] )
        return '%s::%s(%s)' % ( self._type, self._dict[ 'name' ], deps )

    def generate( self, whatvars, *args ):
        s = RooArgSet()
        for i in whatvars :
            s.add( i._target_() if hasattr( i,'_target_' ) else i )
        return self._var.generate( s, *args )

class ProdPdf( Pdf ):
    def __init__( self, name, PDFs, **kwargs ):
        self._dict = { 'pdfs' = frozenset( PDFs ),
                       'name' = name }
        self._dict[ 'name' ] = name + '_' + self._separator().join( [ i.GetName() for i
                                                                      in self._dict[ 'pdfs' ] ] )
        o = set()
        for p in pdfs:
            for i in p[ 'observables' ]:
                o.add( i )
        self._dict[ 'observables' ] = frozenset( o )
        self.__make_pdf()
        
    def __make_pdf( self ):
        if name not in self.ws():
            self._declare( self._makeRecipe( name ) )
            self._init( self._dict[ 'name' ], 'RooProdPdf' )
        else:
            raise StandardError( "Recreating the same Pdf is not supported atm" )
        return self

    def _makeRecipe( self ):
        pdfs = ','.join( [ p.GetName() for p in self._dict[ 'pdfs' ] ] )
        return 'PROD::%s(%s)' % ( self._dict[ 'name' ], pdfs )

    def _separator( self ):
        return '_X_'

class SumPdf( Pdf ):
    def __init__( self, name, PDFs, Yields ):
        self._yields = {}
        self._dict[ 'name' ] = name
        pdfs = list( PDFs )
        if ( len( pdfs ) - len( Yields ) ) not in [ 0, 1 ]:
            raise StandardError( 'The number of yield variables must be equal to or 1'
                                 + 'less then the number of PDFs.' )
            
        ## WORKING HERE

        o = set()
        for p in pdfs:
            if not p.single():
                raise StandardError( 'Can only multiply single PDFs' )
            if p.GetName() not in yields:
                raise KeyError( 'Must provide a yield for each PDF' )
            for i in p.observables().front():
                o.add( i )
        self._name = self._separator().join( [ p.GetName() for p in pdfs ] )
        self._pdfs[ frozenset( o ) ] = self.__make_pdf( self._name, **kwargs )
        
    def __make_pdf( self, name, **kwargs ):
        if name not in self.ws():
            if not 'PDFs' in kwargs:
                raise StandardError( "Must provide PDFs to make a SumPDF" )
            if not 'Yields' in kwargs:
                raise StandardError( "Must provide PDFs to make a SumPDF" )
            self._declare( self._makeRecipe( name, **kwargs ) )
            self._init( name, 'RooAddPdf' )
        else:
            raise StandardError( "Recreating the same Pdf is not supported atm" )
        return self

    def _makeRecipe( self, name, **kwargs ):
        if not 'PDFs' in kwargs:
            raise StandardError( "Must provide PDFs for recipe" )
        if not 'Yields' in kwargs:
            raise StandardError( "Must provide Yields to make a SumPDF" )
        pdfs = ','.join( [ '%s * %s' % ( kwargs[ 'Yields' ][ p.GetName() ], p.GetName() ) \
                           for p in kwargs[ 'PDFs' ] ] )
        return 'SUM::%s(%s)' % ( name, pdfs )

    def _separator( self ):
        return '_P_'

class ResolutionModel( RooObject ):
    _getters = { "Observables" : lambda s : s.observables()
               , "Parameters"  : lambda s : s.parameters()
               ## , "Options"     : lambda s : s.options()
               }

    def __init__( self, name, **kwargs ):
        self._name = name
        if name not in self.ws():
            if 'Type'  not in kwargs: 
                raise KeyError( 'Must specify type' )
            self._type = kwargs.pop( 'Type' )
            if 'Observables'  not in kwargs: 
                raise KeyError( 'Must specify observables' )
            self._observables = frozenset( kwargs.pop( 'Observables' ) )
            variables = list( self._observables )
            if 'Parameters' in kwargs:
                self._parameters = frozenset( kwargs.pop( 'Parameters' ) )
                variables += list( self._parameters )
            ## if 'Options' in kwargs:
            ##     self._options = kwargs.pop( 'Options' )
            ##     variables += self._options
            self._declare( self._makeRecipe( name, Variables = variables ) )
            self._init( name, 'RooResolutionModel' )
        else:
            raise StandardError( 'Recreating the same Resolution Model is not supported atm' )

    def observables( self ):
        return self._observables

    def parameters( self ):
        if hasattr( self, "_parameters" ):
            return self._parameters
        else:
            return None

    ## def options():
    ##     if hasattr( self, "_options" ):
    ##         return self._options 
    ##     else:
    ##         return None

    def __getitem__( self, k ):
        return ResolutionModel._getters[ k ]( self )
    
    def _makeRecipe( self, name, **kwargs ):
        if len( kwargs ) != 1 or 'Variables' not in kwargs:
            raise StandardError( "This method requires a single 'Variables' argument" )
        k = kwargs[ 'Variables' ]
        deps = ','.join( [ v.GetName() for v in k ] )
        return '%s::%s(%s)' % ( self._type, name, deps )
        
class Component( object ):
    _d = {}
    def __init__(self,name) :
        if name in Component._d : 
            # TODO: make things singletons, indexed by 'name'
            raise KeyError('Name %s is not unique'%name )
        self.name = name
        Component._d[name] = dict()
        Component._d[name]['name'] = name
    def _yieldName(self) : return 'N_%s' % self.name
    def setYield( self, n, nlo, nhi ) :
        Component._d[self.name]['yield'] = RealVar( self._yieldName(), MinMax=(nlo,nhi), Value=n).GetName()
    def __setitem__(self, observable, pdf ) :
        if type( observable ) is not tuple : observable = (observable,)

        # create a set of incoming observables
        k = set()
        for o in observable : 
            if type(o) is not str: o = o.GetName()
            k.add(o)
        # do NOT allow overlaps with already registered observables!!!!!! (maybe allow in future....)
        present = set()
        for i in  Component._d[self.name].iterkeys() :
            if type(i) != frozenset : continue # not an observable, but either name or yield
            for j in i : present.add(j)
                
        if not k.isdisjoint( present ) : raise KeyError('sets are not disjoint, overlap: %s' % k.intersection(present) )
        # TODO: check whether observable exists in the workspace...
        # TODO: and check it has its observable flag set
        ### parse recipe: x(a,b[1,2,3],c(4,5,6)) -> ( 'x' , 'a,b[1,2,3],c(4,5,6)' )
        ### i.e. such that %s(%s)'%( x[0], x[1]) recovers the result...
        ## import re
        ## ex = '^([^(]+)\((.*)\)$'
        ## r = re.match(ex,pdfrecipe)
        ## if not r : raise KeyError('could not parse recipe %s'%pdfrecipe)
        ## (typ,args) = r.groups()
        ## if ':' in typ : 
        ##     raise KeyError('please do not explicitly name pdfs -- %s' % typ )
        #TODO: should build the PDF at this point!

        ## Get the right sub-pdf from the Pdf object
        Component._d[ self.name ][ frozenset( k ) ] = pdf[ k ]

    def __iadd__(self,item) :
        self.__setitem__( item.observables(), item )
        return self
        
    def __getitem__(self,k) :
        # TODO: build PDF, and return it -- instead of returning the recipe...
        print 'request from %s: ' % self.name, k
        # TODO: if we return one a-priori build PDF, rename it properly??
        d = Component._d[self.name]
        print d
        
        # TODO: catch yields and name here
        if type(k)==str and k in d : return d[k]

        if type(k) != frozenset : k = frozenset(k)
        if k not in d : 
            # try to build product -- note that d.keys() are non-overlapping by requirement
            # first, find the entry with the largest overlap, which is a subset (otherwise we'd have to marginalize)
            terms = [ ]
            nk = k
            while len(nk) :
                overlap = lambda i : len( i.intersection(nk) ) if type(i)==frozenset and i.issubset(nk) else 0  # avoid marginalization for now...
                from operator import itemgetter
                (kmax,mo) = max( ( (i,overlap(i)) for i in d.iterkeys() ), key = itemgetter(1) )
                if not mo : break # no overlap left...
                terms.append( kmax )
                nk = nk - kmax 
            if len(nk) : raise IndexError( 'could not construct matching product' )
            nk = frozenset.union( *terms )
            print '%s %s' % ( nk, terms )
            pdfs = [ self[ i ][ i.intersection( k ) ] for i in terms ]
            print pdfs
            d[ nk ] = ProdPdf( self.name, PDFs = pdfs )
        print 'returning %s' % d[ k ][ k ]
        return d[ k ][ k ]

def buildPdf( components, observables, name ) :
    # multiply PDFs for observables (for each component)
    if not observables : raise RuntimeError('no observables??')
    obs = [ o if type(o)==str else o.GetName() for o in observables ]
    args = { 'Yields' : {},
             'PDFs'   : [] }
    for c in components:
        pdf = c[ obs ]
        args[ 'Yields' ][ pdf.GetName() ] = c[ 'yield' ]
        args[ 'PDFs'   ].append( pdf )
    # and sum components (inputs should already be extended)
    return SumPdf( name, **args )
