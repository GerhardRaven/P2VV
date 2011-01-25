import RooFitDecorators
from itertools import product,count
from ModelBuilders import *
from ROOT import *
gSystem.Load("libp2vv.so")

class efficiency :
    def __init__(self,*args) : 
        (self.cpsi,self.ctheta,self.phi) = args[0] if len(args)==1 else args
    def accept(self) :
        x = self.cpsi.getVal()
        y = self.ctheta.getVal()
        from math import cos
        z = -1+2*cos( self.phi.getVal() )
        from random import random
        return random() > ( x*x*y ) # /( 1 if z<0 else 1-z ) )

fname="p2vv_9.root"
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

#replace input by inefficient data
eps = efficiency( angles )
inEffData = RooDataSet( "inEffData","inEffData", allObs )
for event in  data :
    allObs.assignValueOnly( event )
    if eps.accept() : inEffData.add( allObs )
data = inEffData;
          
# define the moments used to describe the efficiency
# for this, we need the PDF used to generate the data
marginalObs = pdf.getObservables( data.get() )
marginalObs.remove( angles )
# marginalize pdf over 'the rest' so we get the normalization of the moments right...
pdf_marginal = pdf.createProjection(marginalObs)

# compute the 'canonical' six moments
bnames = [ 'AzAz','AparApar','AperpAperp','AparAperp','AzAperp','AzApar' ]
sixmom = [ EffMoment( w['%s_basis'%n], 1., pdf_marginal, allObs ) for n in bnames ]
computeMoments(data,sixmom)
xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in sixmom ] )


# compute the Fourier series for the efficiency
moments = []
ab = abasis(w,angles)
for (i,l) in product(range(10),range(4)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, pdf_marginal, allObs ) for m in range(-l,l+1) ]

# loop over all data, determine moments
computeMoments(data,moments)


# compute the 'canonical' moments given the Fourier series
c = dict()
for m in moments : c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()

from math import sqrt,pi
xi_c = { 'AparApar_basis'   :   ( c[(0,0,0)]-  c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]-  c[(2,2,0)]/5) + sqrt(3./20)*(c[(0,2,2)]-  c[(2,2,2)]/5)  )
       , 'AzAz_basis'       :   ( c[(0,0,0)]+2*c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]+2*c[(2,2,0)]/5) - sqrt(3./20)*(c[(0,2,2)]+2*c[(2,2,2)]/5)  )
       , 'AperpAperp_basis' :   ( c[(0,0,0)]-  c[(2,0,0)]/5 - sqrt(1./ 5)*( c[(0,2,0)]-  c[(2,2,0)]/5 ) )
       , 'AparAperp_basis'  :   sqrt(3./5.)*( c[(0,2,-1)] - c[(2,2,-1)]/5 )
       , 'AzAperp_basis'    : - sqrt(6./5.)* 3*pi/32 *( c[(1,2, 1)] - c[(3,2, 1)]/4 - 5*c[(5,2, 1)]/128  - 7*c[(7,2, 1)]/512 - 105*c[(9,2, 1)]/16484) 
       , 'AzApar_basis'     :   sqrt(6./5.)* 3*pi/32 *( c[(1,2,-2)] - c[(3,2,-2)]/4 - 5*c[(5,2,-2)]/128  - 7*c[(7,2,-2)]/512 - 105*c[(9,2,-2)]/16484)
       }

# normalize moments and compare
def norm_xi( d ) :
    n =  (d['AparApar_basis'] + d['AzAz_basis'] + d['AperpAperp_basis'])/3
    for i in d.iterkeys() : d[i] = d[i]/n

norm_xi(xi_c)
norm_xi(xi_m)
for name in xi_c.iterkeys() :
    print '%s : direct moment: %s ;  moment computed from Fourier series: %s ; ratio = %s ' % ( name, xi_m[name], xi_c[name], xi_m[name]/xi_c[name])


## build PDF using the Fourier series efficiency...
pdf_eff = buildEff_x_PDF(w,'_eff',pdf,[ ( m.basis() , m.coefficient() ) for m in moments ] )

## build PDF using Fourier coefficients reverse engineered from the six moments...

### make some plots...
c = TCanvas()
c.Divide(3,1);
for (i,var) in enumerate(angles) :
     c.cd(1+i) # start argument of enumerate is python >= 2.6...
     plot = var.frame()
     data.plotOn(plot); 
     pdf.plotOn(plot,RooFit.LineColor(kBlue))
     pdf_eff.plotOn(plot,RooFit.LineColor(kRed))
     plot.Draw()
c.Flush()

