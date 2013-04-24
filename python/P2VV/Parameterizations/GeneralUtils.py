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

    def parameter( self, name ) :
        for par in self._params :
            if par.GetName() == name : return par
        return None

    def _parseArg( self, arg, kwargs, **parsDict ) :
        def _create( arg, kwargs, **parsDict ) :
            import P2VV.RooFitWrappers
            from copy import copy
            _parsDict = copy(parsDict) # make sure we do not modify the input!!!
            objType      = _parsDict.pop( 'ObjectType',   'RealVar' )
            singleArgKey = _parsDict.pop( 'SingleArgKey', 'Value'   )
            if arg in kwargs or arg in _parsDict :
                argPars = kwargs.pop(arg) if arg in kwargs else _parsDict.pop(arg)
                if isinstance( argPars, P2VV.RooFitWrappers.RooObject ) : return argPars
                if type(argPars) == dict and 'ObjectType' in argPars : objType = argPars.pop('ObjectType')
                _parsDict.update( argPars if type(argPars) == dict else { singleArgKey : argPars } )
            if 'Name' not in _parsDict : _parsDict[ 'Name' ] = arg

            return vars(P2VV.RooFitWrappers)[objType](**_parsDict)

        # get object containter list
        contList = parsDict.pop( 'ContainerList', None )

        # create object
        obj = _create( arg, kwargs, **parsDict )

        # either put object in container list or set it as attribute
        if contList != None : contList.append(obj)
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

class _util_extConstraints_mixin( object ) :
    def __init__( self, kwargs ) :
        if   'Constraints' in kwargs : self._constraints = [ constraint for constraint in kwargs.pop('Constraints') ]
        elif 'Constraint'  in kwargs : self._constraints = [ kwargs.pop('Constraint') ]
        else                         : self._constraints = [  ]

    def externalConstraints( self ) : return self._constraints
    def hasExtConstraints( self )   : return len(self._constraints) > 0

class _util_conditionalObs_mixin( object ) :
    def __init__( self, kwargs ) :
        if   'Conditionals' in kwargs : self._conditionals = [ conditional for conditional in kwargs.pop('Conditionals') ]
        elif 'Conditional'  in kwargs : self._conditionals = [ kwargs.pop('Conditional') ]
        else                          : self._conditionals = [  ]

    def conditionalObservables( self ) : return self._conditionals
    def hasCondObservables( self )     : return len(self._conditionals) > 0

    def addConditional( self, condObs ) :
        if not condObs in self._conditionals : self._conditionals.append(condObs)
        else : print 'P2VV - WARNING: _util_conditionalObs_mixin.addConditional(): observable "%s" is already conditional' % condObs

    def addConditionals( self, condObsList ) :
        for obs in condObsList : self.addConditional(obs)


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

from itertools import product
def valid_combinations(states):
    all_states = []
    for level in states:
        all_states.extend([e[0] for e in level])
    all_states = list(set(all_states))
    labels = [[(state, label.GetName()) for label in state] for state in all_states]
    all_combinations = list(product(*labels))
    valid = []
    def good(combination):
        s = set(combination)
        for level in states:
            level_good = False
            for entry in level:
                if entry in s:
                    level_good = True
                    break
            if not level_good:
                return level_good
        return True
    return filter(good, all_combinations)

def exclusive_combinations(states):
    all_states = [e[0] for e in states]
    labels = [[(state, label.GetName()) for label in state] for state in all_states]
    all_combinations = list(product(*labels))
    valid = []
    def good(combination):
        s = set(combination)
        r = set(states)
        return len(s & r) == 1
    return filter(good, all_combinations)
