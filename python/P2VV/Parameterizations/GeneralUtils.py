###########################################################################################################################################
## P2VVParameterizations.GeneralUtils: General P2VV parameterization utilities                                                           ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

_parNamePrefix = ''

def setParNamePrefix(prefix) :
    global _parNamePrefix
    if prefix : _parNamePrefix = prefix
    else      : _parNamePrefix = ''

def getParNamePrefix( fullPrefix = False ) :
    global _parNamePrefix
    if not _parNamePrefix : return ''
    if fullPrefix : return _parNamePrefix + '_'
    return _parNamePrefix

class _util_parse_mixin( object ) :
    def parameters( self ) :
        return self._params

    def parameter( self, name ) :
        for par in self._params :
            if par.GetName() == name : return par
        return None

    def getNamePrefix( self, kwargs ) :
        if hasattr( self, '_namePF' ) :
            assert 'ParNamePrefix' not in kwargs or kwargs['ParNamePrefix'] == self._namePF\
                , 'P2VV -- ERROR: _util_parse_mixin.getNamePrefix: parameter name prefix from arguments is not equal to existing prefix'
        else :
            self._namePF = kwargs.pop( 'ParNamePrefix', None )

        if self._namePF : return self._namePF + '_'
        else :            return getParNamePrefix(True)

    def _parseArg( self, argName, kwargs, **parsDict ) :
        def _create( argName, kwargs, **parsDict ) :
            # get dictionary of parameters to construct the variable
            if argName in kwargs :
                argPars = kwargs.pop(argName)
                parsDict.pop( argName, None )
            else :
                argPars = parsDict.pop( argName, { } )

            # parse parameter dictionary
            import P2VV.RooFitWrappers
            singleArgKey = parsDict.pop( 'SingleArgKey', 'Value' )
            if argPars :
                if isinstance( argPars, P2VV.RooFitWrappers.RooObject ) : return argPars
                parsDict.update( argPars if type(argPars) == dict else { singleArgKey : argPars } )

            # construct variable name
            if 'Name' not in parsDict : parsDict['Name'] = argName
            namePF = self.getNamePrefix(kwargs)
            if 'NamePrefix' in parsDict :
                namePF = parsDict.pop('NamePrefix')
                if namePF : namePF += '_'
            if not parsDict.get( 'Observable', False ) and namePF : parsDict['Name'] = namePF + parsDict['Name']

            # create variable
            objType = parsDict.pop( 'ObjectType', 'RealVar' )
            return vars(P2VV.RooFitWrappers)[objType](**parsDict)

        # create object
        contList = parsDict.pop( 'ContainerList', None )
        obj = _create( argName, kwargs, **parsDict )

        # either put object in container list or set it as attribute
        if contList != None : contList.append(obj)
        else : setattr( self, '_%s' % argName, obj )

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
        assert not hasattr(self,'_constraints')
        if   'Constraints' in kwargs : self._constraints = kwargs.pop('Constraints') 
        elif 'Constraint'  in kwargs : self._constraints = set(kwargs.pop('Constraint') )
        else                         : self._constraints = set()
        assert hasattr(self,'_constraints')

    def ExternalConstraints( self ) : return self._constraints
    def hasExtConstraints( self )   : return len(self._constraints) > 0

    def addConstraint( self, constr ) :
        self._constraints.add(constr)

    def addConstraints( self, constrList ) :
        for constr in constrList : self.addConstraint(constr)

class _util_conditionalObs_mixin( object ) :
    def __init__( self, kwargs ) :
        assert not hasattr(self,'_conditionals')
        if   'Conditionals' in kwargs : self._conditionals = set( kwargs.pop('Conditionals') )
        elif 'Conditional'  in kwargs : self._conditionals = set( kwargs.pop('Conditional') )
        elif 'ConditionalObservables'  in kwargs : self._conditionals = set( kwargs.pop('ConditionalObservables') )
        else                          : self._conditionals = set()
        assert hasattr(self,'_conditionals')

    def ConditionalObservables( self ) : return self._conditionals

    def addConditional( self, condObs ) :
        self._conditionals.add(condObs)

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
