from RooFitDecorators import *

ws = RooWorkspace('ws')
ws.factory('{m[-10,10],t[0,10]}')
ws.defineSet('obs','m,t')
ws.factory('Gaussian::gauss_m(m,m_mean[0],m_width[5])')
ws.factory('Exponential::exp_t(t,t_tau[-0.5])')
ws.factory('PROD::pdf(gauss_m,exp_t)')

pdf = ws['pdf']
obs = ws.set('obs')

data = pdf.generate( obs,  1000)
#pdf.fitTo(data)
ref  = pdf.generate( obs, 10000)

#ws.factory("expr::i2('pdf*pdf',{pdf})")
#ip = pdf.createIntegral(obs)
#ip2 = ws['i2'].createIntegral(obs)
#print ip2.getVal(),ip.getVal(),ip2.getVal()/ip.getVal()

# scramble dataset, and create a new one
def mix( data, ref) :
    if data.get() != ref.get() : raise LogicError('incompatible datasets')
    nd = data.numEntries()
    nr = ref.numEntries()
    d = RooDataSet('d','d',data.get())
    from random import sample
    for i in sample(xrange(nd+nr),nd) : d.add( data.get(i) if i<nd else ref.get(i-nd) )
    return d

def psi(data,ref,sigma=1.0) :
    from math import sqrt,exp
    dsi = data['pdf']/sigma
    rsi = ref['pdf']/sigma
    sqr = lambda x : x*x
    z =   sum( sqr((data[i]-ref[i])) for i in data.iterkeys() if i != 'pdf' ) 
    return exp( -0.5 * z * dsi * rsi )

def normalize(pdf, data, m, s ) :
    ds = list()
    for event in data :
        dx = dict( (x.GetName(),(x.getVal()-m[x.GetName()])/s[x.GetName()]) for x in event )
        pdf.getObservables( event ).assignValueOnly( event )
        # TODO: add phase space volume!!!
        dx['pdf'] = pdf.getVal() # / normalize by norm set volume 
        ds.append( dx ) 
    return ds

def PPD(data,ref,pdf,obs,sigma=0.01) :
    rm = dict( (i.GetName(), ref.mean(i))  for i in obs ) # don't really need this...
    rs = dict( (i.GetName(), ref.sigma(i)) for i in obs )

    # create a secondary dataset for both data and ref, which contains 
    #   1) PDF value for the current observables 
    #   2) normalized values of observables (i.e. corrected for mean & sigma)
    # for each event, add PDF evaluated at that point

    xdata = normalize( pdf, data, rm, rs )
    xref  = normalize( pdf, ref,  rm, rs )

    apdf = sum( x['pdf'] for x in xref ) / len(xref)
    print apdf # this is equal to int dx f(x)^2 / int dx f(x)
    

    from itertools import combinations,product
    x1 = sum( psi(i,j,sigma) for (i,j) in combinations( xdata, 2) ) / (len(xdata)*len(xdata))
    x2 = sum( psi(i,j,sigma) for (i,j) in product( xdata, xref  ) ) / (len(xdata)*len(xref) )
    return x1 - x2


print '     data : %s' %  PPD( data, ref, pdf, obs, 0.01 )
for i in xrange(10) : print 'admixture: %s' % PPD( data, mix(data,ref), pdf,obs, 0.01 )

