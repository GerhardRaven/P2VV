import RooFitDecorators
#from itertools import product,count
from itertools import count
from ModelBuilders import *
from ROOT import *
gSystem.Load("libp2vv.so")
from math import sqrt,pi

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

def compute_moments( c ) :
    return { 'AparApar_basis'   :   ( c[(0,0,0)]-  c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]-  c[(2,2,0)]/5) + sqrt(3./20)*(c[(0,2,2)]-  c[(2,2,2)]/5)  )
           , 'AzAz_basis'       :   ( c[(0,0,0)]+2*c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]+2*c[(2,2,0)]/5) - sqrt(3./20)*(c[(0,2,2)]+2*c[(2,2,2)]/5)  )
           , 'AperpAperp_basis' :   ( c[(0,0,0)]-  c[(2,0,0)]/5 - sqrt(1./ 5)*( c[(0,2,0)]-  c[(2,2,0)]/5 ) )
           , 'AparAperp_basis'  :   sqrt(3./5.)*( c[(0,2,-1)] - c[(2,2,-1)]/5 )
           , 'AzAperp_basis'    : - sqrt(6./5.)* 3*pi/32 *( c[(1,2, 1)] - c[(3,2, 1)]/4) #   - 5*c[(5,2, 1)]/128  - 7*c[(7,2, 1)]/512 - 105*c[(9,2, 1)]/16384) 
           , 'AzApar_basis'     :   sqrt(6./5.)* 3*pi/32 *( c[(1,2,-2)] - c[(3,2,-2)]/4) #   - 5*c[(5,2,-2)]/128  - 7*c[(7,2,-2)]/512 - 105*c[(9,2,-2)]/16384)
           }

class efficiency :
    def __init__(self,*args) : 
        (self.cpsi,self.ctheta,self.phi) = args[0] if len(args)==1 else args
    def accept(self) :
        x = self.cpsi.getVal()
        y = self.ctheta.getVal()
        from math import cos
        z = -1+2*cos( self.phi.getVal() )
        from random import random
        return random() > ( x*x*y )/2 # /( 1 if z<0 else 1-z ) )

fname="p2vv_10.root"
pdfName = "jpsiphipdf"
dataName = "jpsiphipdfData"
workspaceName = "w"

f = TFile(fname)
w = f.Get(workspaceName) 
pdf = w[pdfName]
data = w[dataName] 

allObs = pdf.getObservables( data.get() )

angles = w.set("transversityangles")
#angles = w.set("helicityangles")
angles.Print()
angles.remove( w['tagdecision'] )
angles.remove( w['t'] )

angles.Print()

#replace input by inefficient data
eps = efficiency( angles )
inEffData = RooDataSet( "inEffData","inEffData", allObs )
for event in  data :
    allObs.assignValueOnly( event )
    if eps.accept() : inEffData.add( allObs )
origdata = data
data = inEffData;

### make some plots...
c = TCanvas()
c.Divide(3,2)
for (i,var) in enumerate(angles) :
     c.cd(1+i) # start argument of enumerate is python >= 2.6...
     plot = var.frame()
     origdata.plotOn(plot)
     plot.Draw()
for (i,var) in enumerate(angles) :
    c.cd(4+i)
    plot2 = var.frame()
    data.plotOn(plot2)
    plot2.Draw()
c.Flush()

# define the moments used to describe the efficiency
# for this, we need the PDF used to generate the data
marginalObs = pdf.getObservables( data.get() )
marginalObs.remove( angles )
# marginalize pdf over 'the rest' so we get the normalization of the moments right...
pdf_marginal = pdf.createProjection(marginalObs)

# compute the 'canonical' six moments
bnames = [ 'AzAz','AparApar','AperpAperp','AparAperp','AzAperp','AzApar' ]
sixmom = [ EffMoment( w['%s_basis'%n], 1., pdf_marginal, allObs ) for n in bnames ]

xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in sixmom ] )
print 'eerste keer: xi_m =', xi_m
computeMoments(data,sixmom)
xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in sixmom ] )
print 'tweede keer: xi_m =', xi_m

# compute the Fourier series for the efficiency
moments = []
ab = abasis(w,angles)
for (i,l) in product(range(4),range(4)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, pdf_marginal, allObs ) for m in range(-l,l+1) ]

print 'eerste keer:'
for m in moments:
    print m.coefficient()

# loop over all data, determine moments
computeMoments(data,moments)

print 'tweede keer:'
for m in moments:
    print m.coefficient()

#################################################
### Plot the efficiency as function of angles ###
#################################################
#
addlist = []
for m in moments :
    string = str(m.basis().i())+str(m.basis().l())+str(m.basis().m())
    basis = m.basis()
    coeffrrv = w.put(RooRealVar('coeffrrv_%s'%(string),'coeffrrv_%s'%(string),m.coefficient()))
    prodset = RooArgSet(basis,coeffrrv)
    prod = w.put(RooProduct('prod_%s'%(string),'prod_%s'%(string),prodset))
    addlist.append(prod)

addset = RooArgSet()
for item in addlist:
    addset.add(item)

eff_function = RooAddition('eff_function','eff_function',addset)

from ROOT import gStyle
gStyle.SetOptStat(0)

ceff = TCanvas()
ceff.Divide(3,2)
for (f,v) in enumerate( angles ) :
    ceff.cd(1+f)
    frame = v.frame()
    eff_function.plotOn(frame,RooFit.LineColor(kBlue))
    frame.Draw()
    
    ceff.cd(4+f)
    others = RooArgList( angles )
    others.remove( v )
    hist = eff_function.createHistogram( others.names() )
    eff_function.fillHistogram( hist,others,1., RooArgSet(v))
    hist.Draw('COLZ')
        # create residuals in 2D
        #datahist = data.createHistogram( others.name() )
        #self._dataw[i].fillHistogram( datahist )

##########################################################################

# compute the 'canonical' moments given the Fourier series
c = dict()
for m in moments : c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()

print c

#this calls the reverse engineering function compute_moments, that transforms the efficiency weights (Fourier coeff) into 6 xi's ('moments')
xi_c = compute_moments( c )

print 'xi_c =', xi_c
# normalize moments and compare
def norm_xi( d ) :
    n =  (d['AparApar_basis'] + d['AzAz_basis'] + d['AperpAperp_basis'])/3
    for i in d.iterkeys() : d[i] = d[i]/n

norm_xi(xi_c)
norm_xi(xi_m)

for name in xi_c.iterkeys() :
    print '%s : direct moment: %s ;  moment computed from Fourier series: %s ; ratio = %s ' % ( name, xi_m[name], xi_c[name], xi_m[name]/xi_c[name])

## build PDF using the Fourier series efficiency...
## TODO: normalize relative to 0000 so that c_000 = 1
print 'using the following terms in Fourier expansion: '
for n,c in [ ( m.basis().GetName() , m.coefficient() ) for m in moments if m.significance()>2] :
    print '%s : %s ' % (n,c)

# now compute the six moments with the terms we'll use only (i.e. put the other coefficients to zero)
c_sup = dict()
for m in moments :
    print m.significance()
    c_sup[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient() if m.significance()>2 else 0
   
#print c_sup
xi_c_sup = compute_moments( c_sup )

print 'xi_c_sup =', xi_c_sup

pdf_eff = buildEff_x_PDF(w,'fourier_eff',pdf,[ ( m.basis() , m.coefficient() ) for m in moments if m.significance()>2] )

## build PDF using Fourier coefficients reverse engineered from the six moments...
#
coef = [ (0,0,0,   ( xi_m['AparApar_basis']+xi_m['AzAz_basis']+xi_m['AperpAperp_basis'])/3 )
       , (2,0,0,   ( xi_m['AzAz_basis']-xi_m['AparApar_basis'] )*float(5)/3      )
       , (0,2,0,   ( xi_m['AparApar_basis']-xi_m['AperpAperp_basis'] ) * sqrt(float(20)/9) )  # 1/[sqrt(1/20) + sqrt(1/5) ] = 1/[ 3/2 sqrt(1/5) ] = 2*sqrt(5)/3
       , (0,2,-1,  ( xi_m['AparAperp_basis'] ) * sqrt(float(5)/3) ) 
       , (1,2,1,   ( xi_m['AzAperp_basis'] ) * sqrt(float(5)/6) * float(32)/(3*pi) )
       , (1,2,-1 , ( xi_m['AzApar_basis'] ) * sqrt(float(5)/6)*float(32)/(3*pi) )
       ]
#
print 'reverse engineerd c_ijk:'
print coef
pdf_mom = buildEff_x_PDF(w,'reve_mom_eff',pdf,[ ( ab.build('reve_mom_eff',c[0],0,c[1],c[2],1.), c[3] ) for c in coef ] )

print '*'*80
print '*'*20 + ' fitting original PDF on original data'
print '*'*80
pdf.fitTo(origdata,RooFit.NumCPU(2))
print '*'*80
print '*'*20 + ' fitting original PDF on inefficient data'
print '*'*80
pdf.fitTo(data,RooFit.NumCPU(2))
print '*'*80
print '*'*20 + ' fitting reverse engineerd moment PDF on inefficient data'
print '*'*80
pdf_mom.fitTo(data,RooFit.NumCPU(2))
print '*'*80
print '*'*20 + ' fitting Fourier expansion efficiency PDF on inefficient data'
print '*'*80
pdf_eff.fitTo(data,RooFit.NumCPU(2))

### make some plots...
c2 = TCanvas()
c2.Divide(3,1)
for (i,var) in enumerate(angles) :
     c2.cd(1+i) # start argument of enumerate is python >= 2.6...
     plot = var.frame()
     data.plotOn(plot); 
     pdf.plotOn(plot,RooFit.LineColor(kBlue))
     pdf_eff.plotOn(plot,RooFit.LineColor(kRed))
     pdf_mom.plotOn(plot,RooFit.LineColor(kOrange))
     plot.Draw()
c2.Flush()

c3 = TCanvas()
c3.Divide(3,1)
for (i,var) in enumerate(angles) :
     c3.cd(1+i) # start argument of enumerate is python >= 2.6...
     origplot = var.frame()
     origdata.plotOn(origplot); 
     pdf.plotOn(origplot,RooFit.LineColor(kBlue))
     pdf_eff.plotOn(origplot,RooFit.LineColor(kRed))
     pdf_mom.plotOn(origplot,RooFit.LineColor(kOrange))
     origplot.Draw()
c3.Flush()
