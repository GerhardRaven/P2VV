from RooFitWrappers import Pdf

class TimeAnglesPdf( Pdf ):
    def __init__(self, name, Observables, AngularIndeces, ResolutionModel):
        if len(Observables) != 4:
            raise StandardError('Need four observables: time and three angles')
        
        self._angles = Observables[1:]

        for n, d in AngularIndeces:
            self.__angular_function(n, d)

        super(TimeAnglesPdf, self).__init__(name, o, **d))

    def __angular_function(self, n, d):
        if n in self.ws(): return self.ws()[ n ]
        functions = [self.__basis_function( n, *c ) for c in d]
        return self._declare('sum::%s(%s)' % tuple(n, ','.join( functions )))

    def __basis_function(self, name, i, j, k, l, c):
        # construct name
        name = '%s_P%d%dY%d%d' % tuple([name, i, j, k, l])
        # build basis function
        if name not in self.ws():
            self._declare('P2VVAngleBasis::%s(%s, %s, %s, %d, %d, %d, %d, %f)' \
                          % tuple([name] + [a.GetName() for a in self._angles]
                                  + [i, j, k, l, c]))
        return name

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
