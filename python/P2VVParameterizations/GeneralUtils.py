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

    def _parseArg( self, arg, kwargs, **d ) : 
        def _create( arg,kwargs, **d ) :
            from RooFitWrappers import RealVar, RooObject
            from copy import copy
            _d = copy(d) # make sure we do not modify the input!!!
            if arg in kwargs :
                a = kwargs.pop(arg)
                if isinstance( a, RooObject ) : return a
                _d.update( a if type(a) == dict else { 'Value' : a } ) 
            if 'Name' not in _d : _d[ 'Name' ] = arg
            return RealVar(**_d)
        obj = _create( arg, kwargs, **d )
        setattr(self,'_%s'%arg,obj)
        if not hasattr( self, '_params' ) : self._params = []
        self._params += [ obj ]
        return obj

    def _check_extraneous_kw( self, kwargs ) :
        if len(kwargs): 
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


#def buildEff_x_PDF(name,pdf,eff) :
#   if not eff : return pdf
#   # now we need to multiply all relevant components (i.e. all RooP2VVAngleBasis ones) 
#   # of "pdf" with their efficiency corrected versions, multiply them with the right basis fcn & coefficent
#   # those are assumed to be in eff....
#   from ROOT import RooCustomizer, RooP2VVAngleBasis, RooArgSet
#   customizer = RooCustomizer(pdf,name)
#   for c in pdf.getComponents() :
#        if type(c) is not RooP2VVAngleBasis : continue  # TODO: do not use type to recognize, but name??
#        n = "%s_%s_eff" % (name,c.GetName())
#        s = RooArgSet()
#        from RooFitWrappers import Addition
#        a = [ c.createProduct( fijk, cijk ) for (fijk,cijk) in eff ]
#        # put a in ws...
#        rep = Addition( n, [ c.createProduct( fijk, cijk ) for (fijk,cijk) in eff ] )
#        customizer.replaceArg( c, rep._var() )
#   return customizer.build(True)

