__all__ = ['TimeAnglesPdf']

from RooFitWrappers import Pdf

class TimeAnglesPdf(Pdf):
    def __init__(self, name, Observables, TimeComponents, AngularIndeces, Parameters, ResolutionModel):
        if len(Observables) != 4:
            raise StandardError('Need four observables: time and three angles')

        diff = set(TimeComponents.keys()).symmetric_difference(set(AngularIndeces.keys()))
        if len(diff) != 0:
            print 'Error, time and angular components do not contain the same labels: %s' % diff
            raise RuntimeError()
        
        self._angles = Observables[1:]

        self._angle_functions = {}
        for l, (n, i) in AngularIndeces.iteritems():
            self._angle_functions[l] = self.__angular_function(n, i)

        d = {}
        d['Observables'] = Observables
        # Put the names of the terms in the Parameters
        d['Parameters'] = Parameters[:-1] + TimeComponents[TimeComponents.keys()[0]].keys() \
                          + Parameters[-1:]
        d['ResolutionModel'] = ResolutionModel
        
        super(TimeAnglesPdf, self).__init__(name, **d))

    def __angular_function(self, n, d):
        if n in self.ws()._objects:
            return self.ws()._objects[n]
        functions = [self.__basis_function( n, *c ) for c in d]
        return self._declare('sum::%s(%s)' % tuple(n, ','.join(functions)))

    def __basis_function(self, name, i, j, k, l, c):
        # construct name
        name = '%s_P%d%dY%d%d' % tuple([name, i, j, k, l])
        # build basis function
        if name not in self.ws()._objects:
            bf = self._declare('P2VVAngleBasis::%s(%s, %s, %s, %d, %d, %d, %d, %f)' \
                               % tuple([name] + [a.GetName() for a in self._angles]
                                       + [i, j, k, l, c]))
            self.ws()._objects[name] = bf
        return name

    def _makeRecipe(self, variables):
        terms = self._time_components[self._time_components.keys()[0]].keys()
        sums = {}
        for term in terms:
            s = Zero
            for label in self._time_components.iterkeys():
                s += self._time_components[label][term] * self._angle_functions[label][term]
            sums[term] = s
        variables = variables[:-1] + [s['Name'] for s in sums] + variables[-1:]
        deps = ','.join([v.GetName() if type(v) != str else v for v in variables])
        return "BDecay::%s(%s, %s, %s, %s, %s, %s, %s, %s,\
        %s, SingleSided)" % deps

## Definition of angular functions in helicity basis.
## angFuncs.append(('J_0020x0020_0', [ ( 0, 0, 0,  0,  4.               ),
##                                     ( 0, 0, 2,  0, -sqrt( 16. / 5. ) ),
##                                     ( 2, 0, 0,  0,  8.               ),
##                                     ( 2, 0, 2,  0, -sqrt( 64. / 5. ) ) ] ) )
## angFuncs.append(('J_22x002022_0', [ ( 2, 2, 0,  0,  2.               ),
##                                     ( 2, 2, 2,  0,  sqrt(  1. / 5. ) ),
##                                     ( 2, 2, 2,  2, -sqrt(  3. / 5. ) ) ] ) )
## angFuncs.append(('J_22x002022_1', [ ( 2, 2, 0,  0,  2.               ),
##                                     ( 2, 2, 2,  0,  sqrt(  1. / 5. ) ),
##                                     ( 2, 2, 2,  2,  sqrt(  3. / 5. ) ) ] ) )
## angFuncs.append(('J_21x21_0',     [ ( 2, 1, 2,  1,  sqrt( 24. / 5. ) ) ] ) )
## angFuncs.append(('J_21x2m1_0',    [ ( 2, 1, 2, -1, -sqrt( 24. / 5. ) ) ] ) )
## angFuncs.append(('J_22x2m2_0',    [ ( 2, 2, 2, -2,  sqrt( 12. / 5. ) ) ] ) )
