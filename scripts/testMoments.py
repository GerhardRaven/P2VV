nEvents = 20000
nIters  = 1

# create workspace and observable
from P2VV.RooFitWrappers import RooObject, RealVar
ws = RooObject(workspace = 'momentsWS').ws()
obs = RealVar( Name = 'x', Value = 0., MinMax = ( -1., +1. ) )
weight = RealVar( Name = 'weight', Value = 0. )
ws.defineSet( 'obsSet', 'x' )
ws.defineSet( 'obsWSet', 'x,weight' )
obsSet = ws.set('obsSet')
obsWSet = ws.set('obsWSet')

# build some functions
from math import pi
funcNames = [ 'one', 'cos', 'sin', 'sin2' ]
funcMeans = { 'one' : 1., 'cos' : 1. / 3., 'sin' : pi / 4., 'sin2' : pi / 8. }

ws.factory('Legendre::one(x,0,0)')                        # 1
ws.factory('Legendre::cos(x,1,0)')                        # x = cos(theta)
ws.factory('prod::sin(minus[-1.],Legendre::P11(x,1,1))')  # sqrt(1 - x^2) = sin(theta)
ws.factory('prod::sin2(min23[0.],Legendre::P21(x,2,1))')  # 2x * sqrt(1 - x^2) = sin(2*theta)
ws['min23'].setVal( -2. / 3. )

# create data sets
from ROOT import RooDataSet
noWeightData = RooDataSet( 'noWeightData', 'noWeightData', obsSet            )
weightData   = RooDataSet( 'weightData',   'weightData',   obsWSet, 'weight' )

# build moments builders
from P2VV.Utilities.DataMoments import RealMomentsBuilder
from P2VV.RooFitWrappers import RealMoment
moms  = RealMomentsBuilder( Moments = [ RealMoment( Name = ws[func].GetName(), BasisFunc = ws[func], Norm = 1. ) for func in funcNames ] )
momsW = RealMomentsBuilder( Moments = [ RealMoment( Name = ws[func].GetName(), BasisFunc = ws[func], Norm = 1. ) for func in funcNames ] )
moms.initCovariances()
momsW.initCovariances()

# initialize variables for calculation of covariances
covs  = dict( [ ( name1, dict( [ ( name2, 0. ) for name2 in funcNames ] ) ) for name1 in funcNames ] )
covsW = dict( [ ( name1, dict( [ ( name2, 0. ) for name2 in funcNames ] ) ) for name1 in funcNames ] )

# generate data and compute moments in loop
import sys
from ROOT import RooRandom
for it in range(nIters) :
    if it % 100 == 0 :
        print 'iteration %d' % it
        sys.stdout.flush()

    # generate data with (1 + x) / 2
    RooRandom.randomGenerator().SetSeed( 100000 + it )
    noWeightData.reset()
    weightData.reset()
    for evIt in range(nEvents) :
        weight = RooRandom.uniform()
        obs.setVal( 2. * weight - 1. )
        if RooRandom.uniform() <= weight :
            noWeightData.add(obsSet)
        weightData.add( obsSet, weight )

    # compute moments of functions
    moms.compute(  noWeightData, ResetFirst = True, Verbose = False )
    momsW.compute( weightData,   ResetFirst = True, Verbose = False )

    # update covariance sums
    momCoefs  = moms.coefficients()
    momCoefsW = momsW.coefficients()
    for name1 in funcNames :
        for name2 in funcNames :
            covs[name1][name2]  += ( momCoefs[name1][0]  - funcMeans[name1] ) * ( momCoefs[name2][0]  - funcMeans[name2] )
            covsW[name1][name2] += ( momCoefsW[name1][0] - funcMeans[name1] ) * ( momCoefsW[name2][0] - funcMeans[name2] )

# divide covariance sums by number of entries to get covariances
for name1 in funcNames :
    for name2 in funcNames :
        covs[name1][name2]  /= float(nIters)
        covsW[name1][name2] /= float(nIters)

# print data sets from last iteration
noWeightData.Print()
weightData.Print()
for evIt in range(10) :
    print 'not weighted: %+.2f (%.2f)   weighted: %+.2f (%.2f)'\
         % ( noWeightData.get(evIt).getRealValue('x'), noWeightData.weight(), weightData.get(evIt).getRealValue('x'), weightData.weight() )
print

# print moments from last iteration and computed covariances
maxLenName = max( [ len(name) for name in funcNames ] )
from math import sqrt
def printCovs(covMat) :
    print '  errors:'
    for name in funcNames :
        print ( '  {0:>%ds}:  {1:>7.5f}' % maxLenName ).format( name, sqrt( covMat[name][name] ) )
    print

    print '  correlations:'
    print ' ' * ( maxLenName + 5 ) + '   '.join( '{0:>6s}'.format(name) for name in funcNames )
    for name1 in funcNames :
        corrs = [ ( covMat[name1][name2] / sqrt( covMat[name1][name1] * covMat[name2][name2] ) )\
                  if covMat[name1][name1] > 0. and covMat[name2][name2] > 0.\
                  else float( cmp( covMat[name1][name2], 0. ) ) * float('inf') if covMat[name1][name2] != 0.\
                  else 1. if name1 == name2\
                  else 0.\
                  for name2 in funcNames
                ]
        print ( '  {0:>%ds}   ' % maxLenName ).format(name1) + '   '.join( '{0:>+6.3f}'.format(corr) for corr in corrs )
    print

moms.Print()
printCovs(covs)
momsW.Print()
printCovs(covsW)
