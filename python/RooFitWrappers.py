from RooFitDecorators import *
from copy import copy

class RooObject( object ) :
    _ws = None
    _dict = None

    def ws(self) : 
        if not RooObject._ws : raise RuntimeError('No workspace defined!')
        return RooObject._ws

    def setWorkspace(self,ws):
        RooObject._ws = ws
        if not hasattr(ws, '_objects'):
            ws._objects = {}
        
    def _declare(self,spec):
        x = self.ws().factory(spec)
        if not x:
            raise NameError('failed to _declare %s to workspace factory'%spec)
        self.ws()._objects[x.GetName()] = x
        x._observable = False
        return x

    def _init(self,name,type) :
        x = self.ws()[name]
        if not x:
            raise KeyError('could not locate %s' % name)
        if not x.InheritsFrom(type): 
            raise KeyError('%s is %s, not %s' % (name, x.ClassName(), type))
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
               , 'Name'       : lambda s : s.GetName()
               }

    def __init__(self,name,**kwargs):
        # TODO: add blinding support to kwargs
        if name not in self.ws():
            # construct factory string on the fly...
            if 'Value'  not in kwargs  and 'MinMax' not in kwargs:
                raise KeyError('%s does not exist yet, neither Value nor MinMax specified'%name)
            if 'Value' not in kwargs:
                (mi,ma) = kwargs.pop('MinMax')
                self._declare("%s[%s,%s]"%(name,mi,ma))
            elif 'MinMax' not in kwargs:
                self._declare("%s[%s]"%(name,kwargs.pop('Value')))
            else:
                (mi,ma) = kwargs.pop('MinMax')
                self._declare("%s[%s,%s,%s]"%(name,kwargs.pop('Value'),mi,ma))
            self._init(name,'RooRealVar')
            for (k,v) in kwargs.iteritems() : self.__setitem__(k,v)
        else:
            self._init(name,'RooRealVar')
            
    def __setitem__(self,k,v):
        return RealVar._setters[k]( self, v )
    def __getitem__(self,k):
        return RealVar._getters[k]( self )
    # overrule RooRealVar.setRange
    def setRange(self, v):
        (mi,ma) = v
        self._var.setRange( mi,ma )
        if self.getVal() < mi or self.getVal() > ma : self.setVal( 0.5*(ma+mi)  )

    def getRange(self):
        return self._target_().getMin(), self._target_().getMax()

##TODO, factor out common code in Pdf and ResolutionModel

class Pdf( RooObject ):
    _getters = { 'Observables' : lambda s : s._get( 'Observables' ) 
               , 'Type'        : lambda s : s._get( 'Type' )
               , 'Options'     : lambda s : s._get( 'Options' )
               , 'Parameters'  : lambda s : s._get( 'Parameters' )
               , 'Name'        : lambda s : s._get( 'Name' )
               }

    ## TODO: define operators
    def __init__(self, name, **kwargs):
        if 'Observables' not in kwargs:
            raise KeyError('Must provide observables %s' % kwargs)
        if 'Type' not in kwargs:
            raise KeyError('Must provide type %s' % kwargs)

        # Save the keyword args as properties
        self._dict = {}
        for k, v in kwargs.iteritems():
            self._dict[k] = v
        for k in ['Observables', 'Parameters']:
            self._dict[k] = frozenset(self._dict[k])
        self._dict['Name'] = self._separator().join([name] + [i.GetName() for i in self._dict['Observables']])

        self.__make_pdf()
        del self._dict

    def __str__( self ):
        d = dict( [ ( a, self[ a ] ) for a in Pdf._getters if hasattr( self, a ) ] )
        return '%s' % d

    def _get(self, name):
        attr = '_' + name.lower()
        return getattr(self, attr)
    
    def __getitem__( self, k ):
        if self._dict and k in self._dict:
            return self._dict[k]
        else:
            try:
                return Pdf._getters[k](self)
            except AttributeError as error:
                raise KeyError(str(error))

    def __make_pdf( self ):
        if self._dict['Name'] not in self.ws():
            v = list(self._dict['Observables'])
            if 'Parameters' in self._dict:
                v += list(self._dict['Parameters'])
            if 'ResolutionModel' in self._dict:
                rm = self._dict['ResolutionModel']
                if rm['Observables'] != frozenset(self._dict['Observables']):
                    raise StandardError ('Error, resolution model must have the same '
                                         + 'observables as PDF')
                ## TODO add conditional observables
                v += [rm]
            if 'Options' in self._dict:
                v += list(self._dict['Options'])
            self._declare(self._makeRecipe(v))
            self._init(self._dict['Name'], 'RooAbsPdf')

            # Change self._dict into attributes. Cannot be done before since the
            # underlying object does only exists at this point.
            for k, v in self._dict.iteritems():
                attr = '_' + k.lower()
                setattr(self._target_(), attr, v)
        else:
            self._init(self._dict['Name'], 'RooAbsPdf')

    def _separator( self ):
        return '_'

    def _makeRecipe( self, variables ):
        deps = ','.join( [ v.GetName() if type( v ) != str else v for v in variables ] )
        return '%s::%s(%s)' % ( self._dict[ 'Type' ], self._dict[ 'Name' ], deps )

    def generate( self, whatvars, *args ):
        s = RooArgSet()
        for i in whatvars :
            s.add( i._target_() if hasattr( i,'_target_' ) else i )
        return self._var.generate( s, *args )

class ProdPdf(Pdf):
    def __init__(self, name, PDFs, **kwargs):
        self._dict = {'PDFs' : frozenset(PDFs)}
        self._dict['Name'] = name + '_' + self._separator().join([i.GetName() for i
                                                                  in self._dict['PDFs']])
        o = set()
        for p in self._dict['PDFs']:
            for i in p['Observables']:
                o.add( i )
        self._dict['Observables'] = frozenset(o)
        self.__make_pdf()
        del self._dict
        
    def __make_pdf(self):
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
        return self

    def _makeRecipe(self):
        pdfs = ','.join( [ p.GetName() for p in self._dict[ 'PDFs' ] ] )
        return 'PROD::%s(%s)' % ( self._dict[ 'Name' ], pdfs )

    def _separator( self ):
        return '_X_'

class SumPdf( Pdf ):
    def __init__( self, name, PDFs, Yields ):
        self._yields = {}
        self._dict = {'Name'  : name,
                      'Yields': Yields,
                      'PDFs'  : PDFs}
        pdfs = list(PDFs)
        diff = set([p.GetName() for p in pdfs]).symmetric_difference(set(Yields.keys()))
        if len(diff) not in [0, 1]:
            raise StandardError('The number of yield variables must be equal to or 1'
                                + 'less then the number of PDFs.')

        o = set()
        for p in pdfs:
            for i in p['Observables']:
                o.add(i)
        self._dict['Name'] = self._separator().join([p.GetName() for p in pdfs])
        self.__make_pdf()
        del self._dict
        
    def __make_pdf( self ):
        if self._dict[ 'Name' ] not in self.ws():
            self._declare( self._makeRecipe() )
            self._init( self._dict[ 'Name' ], 'RooAddPdf' )

            # Change self._dict into attributes. Cannot be done before since the
            # underlying object does only exists at this point.
            for k, v in self._dict.iteritems():
                attr = '_' + k.lower()
                setattr(self._target_(), attr, v)
        else:
            self._init( self._dict['Name'], 'RooAddPdf')

    def _makeRecipe( self ):
        yields = self._dict[ 'Yields' ]
        pdfs = ','.join( [ '%s * %s' % ( yields[ p.GetName() ], p.GetName() )
                           if p.GetName() in yields else p.GetName()
                           for p in self._dict[ 'PDFs' ] ] )
        return 'SUM::%s(%s)' % ( self._dict[ 'Name' ], pdfs )

    def _separator( self ):
        return '_P_'

class ResolutionModel( RooObject ):
    _getters = { 'Observables' : lambda s : s._get( 'Observables' )
               , 'Parameters'  : lambda s : s._get( 'Parameters' )
               , 'Type'        : lambda s : s._get( 'Type' )
               }

    def __init__( self, name, **kwargs ):
        if 'Type'  not in kwargs: 
            raise KeyError( 'Must specify type' )
        if 'Observables'  not in kwargs: 
            raise KeyError( 'Must specify observables' )

        if name not in self.ws():
            # Save the keyword args as properties
            self._dict = {}
            self._dict[ 'Name' ] = name
            def _fs( n ):
                self._dict[ n ] = frozenset( self._dict[ n ] )
                
            for k, v in kwargs.iteritems():
                self._dict[ k ] = v

            for a in [ 'Observables', 'Parameters' ]:
                if a in self._dict: _fs( a )

            self._declare( self._makeRecipe() )
            self._init( name, 'RooResolutionModel' )
            for k, v in self._dict.iteritems():
                attr = '_' + k.lower()
                setattr( self._target_(), attr, v )
            del self._dict
        else:
            self._init( name, 'RooResolutionModel' )

    def _get( self, name ):
        attr = '_' + name.lower()
        return getattr( self, attr )
            
    def __getitem__( self, k ):
        return ResolutionModel._getters[ k ]( self )
    
    def _makeRecipe( self ):
        variables = list( self._dict[ 'Observables' ] )
        if 'Parameters' in self._dict:
            variables += list( self._dict[ 'Parameters' ] )
        deps = ','.join( [ v.GetName() for v in variables ] )
        return '%s::%s(%s)' % ( self._dict[ 'Type' ], self._dict[ 'Name' ], deps )
        
class Component( object ):
    _d = {}
    def __init__(self,name) :
        if name in Component._d : 
            # TODO: make things singletons, indexed by 'Name'
            raise KeyError('Name %s is not unique'%name )
        self.name = name
        Component._d[name] = dict()
        Component._d[name]['Name'] = name
    def _yieldName(self) : return 'N_%s' % self.name
    def setYield( self, n, nlo, nhi ) :
        Component._d[self.name]['Yield'] = RealVar( self._yieldName(), MinMax=(nlo,nhi), Value=n).GetName()
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
        Component._d[ self.name ][ frozenset( k ) ] = pdf

    def __iadd__(self,item) :
        self.__setitem__( item.observables(), item )
        return self
        
    def __getitem__(self,k) :
        # TODO: if we return one a-priori build PDF, rename it properly??
        d = Component._d[self.name]
        
        # Catch yields and name here
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
            pdfs = [ self[ i ] for i in terms ]
            d[ nk ] = ProdPdf( self.name, PDFs = pdfs )
        return d[ k ]

def buildPdf( components, observables, name ) :
    # multiply PDFs for observables (for each component)
    if not observables : raise RuntimeError('no observables??')
    obs = [ o if type(o)==str else o.GetName() for o in observables ]
    args = { 'Yields' : {},
             'PDFs'   : [] }
    for c in components:
        pdf = c[obs]
        args['Yields'][pdf.GetName()] = c['Yield']
        args['PDFs'].append(pdf)
    # and sum components (inputs should already be extended)
    return SumPdf(name,**args)
