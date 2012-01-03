###########################################################################################################################################
## P2VVParameterizations.GeneralUtils: General P2VV parameterization utilities                                                           ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

class _util_parse_mixin( object ) :
    def parameters( self ) :
        return self._params

    def _parseArg( self, arg, kwargs, ContainerList = None, SingleArgKey = 'Value', **d ) :
        def _create( arg,kwargs, **d ) :
            from RooFitWrappers import RealVar, RooObject
            from copy import copy
            _d = copy(d) # make sure we do not modify the input!!!
            if arg in kwargs :
                a = kwargs.pop(arg)
                if isinstance( a, RooObject ) : return a
                _d.update( a if type(a) == dict else { SingleArgKey : a } )
            if 'Name' not in _d : _d[ 'Name' ] = arg
            return RealVar(**_d)

        # create object
        obj = _create( arg, kwargs, **d )

        # either put object in container list or set it as attribute
        if ContainerList != None : ContainerList.append(obj)
        else : setattr( self, '_%s' % arg, obj )

        # put object in parameters
        if not hasattr( self, '_params' ) : self._params = []
        self._params += [ obj ]

        return obj

    def _check_extraneous_kw( self, kwargs ) :
        if kwargs: 
            raise KeyError('got unknown keywords %s for %s' % ( kwargs, type(self) ) )

    def setValues( self, **kwargs ) :
        for ( k, v ) in kwargs.iteritems() : 
          arg = getattr( self, '_' + k )
          if v < arg.getMin() : arg.setMin(v) 
          if v > arg.getMax() : arg.setMax(v) 
          arg['Value'] = v

    def setConstant( self, pattern, constant = True ) :
        import re
        rc = 0
        nrexp = re.compile(pattern)
        for i in self.parameters(): 
            if not nrexp.match( i.GetName() ) : continue
            from ROOT import RooAbsLValue
            if not isinstance( i._var, RooAbsLValue) : continue
            i.setConstant (constant)
            rc += 1
        return rc



#def normalize_individual( name, pdf, tag ) :
#    pl = RooArgList()
#    for t in tag._var :
#        tr = ConstVar('const_%s_%s'%(tag,t.GetName(),Value = t.getVal() )
#        from ROOT import RooCustomizer
#        customizer = RooCustomizer( pdf._var, '%s_%s'%(name,t.GetName() )
#        customizer.replaceArg( tag._var, tr._var )
#        pl += customizer.build(True)
#    # TODO: wrap RooSimultaneous in a RooObject...
#    return RooSimultaneous( name, name, pl, tag )


