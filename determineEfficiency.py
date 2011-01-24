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
moments = []

ab = abasis(w,angles)
for (i,l) in product(range(6),range(4)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, pdf_marginal, allObs ) for m in range(-l,l+1) ]


# loop over all data, determine moments
pdf_eff = buildEffMomentsPDF(w, "_eff", pdf, data, moments)


# compare with the canonical six moments
bnames = [ 'AzAz','AparApar','AperpAperp','AparAperp','AzAperp','AzApar' ]
sixmom = [ EffMoment( w['%s_basis'%n], 1., pdf_marginal, allObs ) for n in bnames ]
computeMoments(data,sixmom)

# 
c = dict()
for m in moments :
   c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()

from math import sqrt,pi
xi = { 'AparApar_basis'   :  8*sqrt(pi)/9 * ( c[(0,0,0)]-  c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]-  c[(2,2,0)]/5) + sqrt(3./20)*(c[(0,2,2)]-  c[(2,2,2)]/5)  )
     , 'AzAz_basis'       :  8*sqrt(pi)/9 * ( c[(0,0,0)]+2*c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]+2*c[(2,2,0)]/5) - sqrt(3./20)*(c[(0,2,2)]+2*c[(2,2,2)]/5)  )
     , 'AperpAperp_basis' :  8*sqrt(pi)/9 * ( c[(0,0,0)]-  c[(2,0,0)]/5 - sqrt(1./ 5)*( c[(0,2,0)]-  c[(2,2,0)]/5 ) )
     , 'AparAperp_basis'  :  8*sqrt(pi)/9 * sqrt(3./5.)*( c[(0,2,-1)] - c[(2,2,-1)]/5 )
     , 'AzAperp_basis'    : -8*sqrt(pi)/9 * sqrt(6./5.)* 3*pi/32 *( c[(1,2,1)] - c[(3,2,1)]/4 - 5*c[(5,2,1)]/128 ) 
     , 'AzApar_basis'     :  8*sqrt(pi)/9 * sqrt(6./5.)* 3*pi/32 *( c[(1,2,-2)] - c[(3,2,-2)]/4  -5*c[(5,2,-2)]/128)
     }

for m in sixmom :
    name = m.basis().GetName()
    print '%s : moment: %s  computed: %s ; ratio = %s ' % ( name, m.coefficient(), xi[name], m.coefficient()/xi[name] )

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

#
#xi = { 'parpar' : 1
#     , '00'     : 1
#     , 'perpperp' : 1
#     , 'perp0'  : 0,
#     , 'par0'   : 0,
#     , 'parperp' : 0 ] # TODO: verify signconvention!!!
#
#coef = [ (0,0,0,   ( xi['parpar']+xi['00']+xi['perpperp'])/3 )
#       , (2,0,0,   ( xi['00']-xi['parpar'] )*float(5)/3      )
#       , (0,2,0,   ( xi['parpar']-xi['perpperp'] ) * sqrt(float(20)/9) )
#       , (0,2,-1,  ( xi['parperp'] ) * sqrt(float(5)/3) )
#       , (1,2,1,   ( xi['perp0'] ) * sqrt(float(5)/6) * float(32)/(3*pi) )
#       , (1,2,-1 , ( xi['par0'] ) * sqrt(float(5)/6)*float(32)/(3*pi) )
#       ]
#
# pdf = buildEff_x_PDF(w,'eff',pdf,[ ab.build('mom_eff',c[0],0,c[1],c[2],1.), c[3] for c in coef ] )
       
