########################################
### Author: Daan van Eijk
### Updated on: Jun 5 11
### Description: This script is an attempt to 'reverse engineer' a Fourier series from the given 6 angular efficiency moments.
###              This attempt never succeeded, and instead we used our own formalism of calculating the moments with respect to the Fourier bases
###              In the end the point estimates agree with analyses using the moments when using our own formalism
###              So this script is not so much needed anymore.
########################################
from itertools import count
from math import sqrt,pi
from ModelBuilders import *
from RooFitDecorators import *
from ROOT import *
from P2VVLoad import P2VVLibrary

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s", type="int", dest="seed")
(options, args) = parser.parse_args()

_seed = options.seed
print 'seed =', _seed

RooRandom.randomGenerator().SetSeed(_seed)

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

useTransversityAngles = True

read = False

if read:
    fname="p2vv_toydata.root"
    pdfName = "jpsiphipdf"
    dataName = "jpsiphipdfData"
    workspaceName = "ws"

    f = TFile(fname)
    ws = f.Get(workspaceName) 
    pdf = ws[pdfName]
    data = ws[dataName]    

else:
    ws = RooWorkspace('ws')
    declareObservables(ws,'Bs2Jpsiphi')
    definePolarAngularAmplitudes(ws)
    defineJPsiPhiPhysicsParams(ws)

    ws.factory("RooGaussModel::tres_sig(t,mu[0],sigma[0.05])")
    
    ws.factory("{wtag[0.0]}")


    pdf = buildJpsiphi(ws,'jpsiphipdf',useTransversityAngles)  ## for now we rely quite a bit on a naming convention -- 
                                     ## in future we should pass more information into the builder
                                     ## maybe a dictionary of what's what...

    obs = ws.set('transversityangles' if useTransversityAngles else 'helicityangles')
    obs.add( ws.argSet('t,tagdecision') )

    data = pdf.generate( obs, 10000)
    ws.put(data)

    file = TFile("p2vv_toydata.root","RECREATE")
    ws.Write("ws")
    file.Close()

def compute_moments( c ) :
    return { 'AparApar_basis'   :   ( c[(0,0,0)]-  c[(2,0,0)]/5. + sqrt(1./20.)*( c[(0,2,0)]-  c[(2,2,0)]/5.) + sqrt(3./20)*(c[(0,2,2)]-  c[(2,2,2)]/5.)  )
           , 'AzAz_basis'       :   ( c[(0,0,0)]+2.*c[(2,0,0)]/5. + sqrt(1./20.)*( c[(0,2,0)]+2.*c[(2,2,0)]/5.) - sqrt(3./20.)*(c[(0,2,2)]+2.*c[(2,2,2)]/5.)  )
           , 'AperpAperp_basis' :   ( c[(0,0,0)]-  c[(2,0,0)]/5. - sqrt(1./ 5.)*( c[(0,2,0)]-  c[(2,2,0)]/5. ) )
           , 'AparAperp_basis'  :   sqrt(3./5.)*( c[(0,2,-1)] - c[(2,2,-1)]/5. )
           , 'AzAperp_basis'    : - sqrt(6./5.)* 3.*pi/32. *( c[(1,2, 1)] - c[(3,2, 1)]/4.     )# - 5.*c[(5,2, 1)]/128.  - 7.*c[(7,2, 1)]/512. - 105.*c[(9,2, 1)]/16384. - 231.*c[(11,2, 1)]/65536. - 9009.*c[(13,2, 1)]/4194304. - 23595.*c[(15,2, 1)]/16777216. - 1042899.*c[(17,2, 1)]/1073741824.) 
           , 'AzApar_basis'     :   sqrt(6./5.)* 3.*pi/32. *( c[(1,2,-2)] - c[(3,2,-2)]/4.     )# - 5.*c[(5,2,-2)]/128.  - 7.*c[(7,2,-2)]/512. - 105.*c[(9,2,-2)]/16384. - 231.*c[(11,2,-2)]/65536. - 9009.*c[(13,2,-2)]/4194304. - 23595.*c[(15,2,-2)]/16777216. - 1042899.*c[(17,2,-2)]/1073741824.)
           }
def legendre(x) :
        return (1,x,0.5*(3*x*x-1),0.5*(5*x*x*x-3*x) )
def spharmonic(ct,phi) :
        st = sqrt(1-ct*ct)
        x=st*cos(phi)
        y=st*sin(phi)
        z=ct
        n0 = sqrt(1./(4*pi))
        n1 = sqrt(3./(4*pi))
        return ( (n0,), (-n1*y,n1*z,-n1*x),  )  # 0,1,2  -> -1,0,+1

#class efficiency :
#    def __init__(self,*args) : 
#        (self.cpsi,self.ctheta,self.phi) = args[0] if len(args)==1 else args
#    def accept(self) :
#        p = legendre( self.cpsi.getVal() )
#        y = spharmonic( self.ctheta.getVal(), self.phi.getVal() )
#        from random import random
#        #return random() < ( p[0]+0.2*p[1]+0.4*p[2] )/3 * ( y[0][0] + 0.1*y[1][-1+1] + 0.2*y[1][0+1] + 0.3*y[1][1+1] )
#        return random() < ( p[0]+p[2]) * ( y[0][0])# + 0.1*y[1][-1+1] + 0.2*y[1][0+1] + 0.3*y[1][1+1] ) 
    
class efficiency :
    def __init__(self,*args) : 
        (self.cpsi,self.ctheta,self.phi) = args[0] if len(args)==1 else args
    def accept(self) :
        x = self.cpsi.getVal()
        y = self.ctheta.getVal()
        from math import cos
        z = -1+2*cos( self.phi.getVal() )
        from random import random,seed
        return random() > ( x*x*y )/2 # /( 1 if z<0 else 1-z ) )

allObs = pdf.getObservables( data.get() )

angles = ws.set('transversityangles' if useTransversityAngles else 'helicityangles')
angles.Print()
angles.remove( ws['tagdecision'] )
angles.remove( ws['t'] )

angles.Print()

#replace input by inefficient data
eps = efficiency( angles )
inEffData = RooDataSet( "inEffData","inEffData", allObs )
for event in  data :
    allObs.assignValueOnly( event )
    if eps.accept() : inEffData.add( allObs )
origdata = data
data = inEffData;


# define the moments used to describe the efficiency
# for this, we need the PDF used to generate the data
marginalObs = pdf.getObservables( data.get() )
marginalObs.remove( angles )
# marginalize pdf over 'the rest' so we get the normalization of the moments right...
pdf_marginal = pdf.createProjection(marginalObs)

# compute the 'canonical' six moments
bnames = [ 'AzAz','AparApar','AperpAperp','AparAperp','AzAperp','AzApar' ]
sixmom = [ EffMoment( ws['%s_basis'%n], 1., pdf_marginal, allObs ) for n in bnames ]

xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in sixmom ] )
#print 'eerste keer: xi_m =', xi_m
computeMoments(data,sixmom)
xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in sixmom ] )
print 'Direct computed moments: ', xi_m

# compute the Fourier series for the efficiency
moments = []
ab = abasis(ws,angles)
for (i,l) in product(range(4),range(4)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, pdf_marginal, allObs ) for m in range(-l,l+1) ]

# loop over all data, determine moments
computeMoments(data,moments)

# compute the 'canonical' moments given the Fourier series
c = dict()
for m in moments : c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()
c_000 = c[(0,0,0)]

print 'Fourier coefficients are', c

signif = 5

print 'Using the following terms in the Fourier expansion (significance > %f)'%(signif)

for m in moments:
    if m.significance() > signif:
        print 'basis: ', m.basis().getTitle(), ', value = ', m.coefficient(), ', significance = ', m.significance()

#Correct the PDF with these Fourier coefficients 
#pdf_eff = buildEff_x_PDF(ws,'fourier_eff',pdf,[ ( m.basis() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] )
pdf_eff = buildEff_x_PDF(ws,'fourier_eff',pdf,[ ( m.basis() , m.coefficient() ) for m in moments if m.significance()>signif] )

#this calls the reverse engineering function compute_moments, that transforms the efficiency weights (Fourier coeff) into 6 xi's ('moments')
xi_c = compute_moments( c )

print 'These are the six moments computed from the Fourier series: xi_c =', xi_c
# normalize moments and compare
def norm_xi( d ) :
    n =  (d['AparApar_basis'] + d['AzAz_basis'] + d['AperpAperp_basis'])/3
    for i in d.iterkeys() : d[i] = d[i]/n

norm_xi(xi_c)
norm_xi(xi_m)

for name in xi_c.iterkeys() :
    print '%s : direct moment: %s ;  moment computed from Fourier series: %s ; ratio = %s ' % ( name, xi_m[name], xi_c[name], xi_m[name]/xi_c[name])

## build reverse engineered Fourier series using the direct moments (xi_m), to correct the pdf
#coef = [ (0,0,0,   ( xi_m['AparApar_basis']+xi_m['AzAz_basis']+xi_m['AperpAperp_basis'])/3 )
#       , (2,0,0,   ( xi_m['AzAz_basis']-xi_m['AparApar_basis'] )*float(5)/3      )
#       , (0,2,0,   ( xi_m['AparApar_basis']-xi_m['AperpAperp_basis'] ) * sqrt(float(20)/9) )  # 1/[sqrt(1/20) + sqrt(1/5) ] = 1/[ 3/2 sqrt(1/5) ] = 2*sqrt(5)/3
#       , (0,2,-1,  ( xi_m['AparAperp_basis'] ) * sqrt(float(5)/3) ) 
#       , (1,2,1,   ( xi_m['AzAperp_basis'] ) * sqrt(float(5)/6) * float(32)/(3*pi) )
#       , (1,2,-1 , ( xi_m['AzApar_basis'] ) * sqrt(float(5)/6)*float(32)/(3*pi) )
#       ]

## build reverse engineered Fourier series using the moments from the full Fourier series (xi_c), to correct the pdf
coef = [ (0,0,0,   ( xi_c['AparApar_basis']+xi_c['AzAz_basis']+xi_c['AperpAperp_basis'])/3 )
       , (2,0,0,   ( xi_c['AzAz_basis']-xi_c['AparApar_basis'] )*float(5)/3      )
       , (0,2,0,   ( xi_c['AparApar_basis']-xi_c['AperpAperp_basis'] ) * sqrt(float(20)/9) )  # 1/[sqrt(1/20) + sqrt(1/5) ] = 1/[ 3/2 sqrt(1/5) ] = 2*sqrt(5)/3
       , (0,2,-1,  ( xi_c['AparAperp_basis'] ) * sqrt(float(5)/3) ) 
       , (1,2,1,   ( xi_c['AzAperp_basis'] ) * sqrt(float(5)/6) * float(32)/(3*pi) )
       , (1,2,-2 , ( xi_c['AzApar_basis'] ) * sqrt(float(5)/6)*float(32)/(3*pi) )
       ]

print 'reverse engineerd c_ijk:'
print coef
pdf_mom = buildEff_x_PDF(ws,'reve_mom_eff',pdf,[ ( ab.build('reve_mom_eff',c[0],0,c[1],c[2],1.), c[3] ) for c in coef ] )

####################
### FitSet #########
####################

dG_OrigFit = RooRealVar("dG_OrigFit","dG_OrigFit",-1, 1)
dG_OrigFitErr = RooRealVar("dG_OrigFitErr","dG_OrigFitErr",-1, 1)

dG_WrongFit = RooRealVar("dG_WrongFit","dG_WrongFit",-1, 1)
dG_WrongFitErr = RooRealVar("dG_WrongFitErr","dG_WrongFitErr",-1, 1)

dG_AccFit = RooRealVar("dG_AccFit","dG_AccFit",-1, 1)
dG_AccFitErr = RooRealVar("dG_AccFitErr","dG_AccFitErr",-1, 1)

dG_FourAccFit = RooRealVar("dG_FourAccFit","dG_FourAccFit",-1, 1)
dG_FourAccFitErr = RooRealVar("dG_FourAccFitErr","dG_FourAccFitErr",-1, 1)


fitArgSet = RooArgSet(dG_OrigFit,
                      dG_OrigFitErr,
                      dG_WrongFit,
                      dG_WrongFitErr,
                      dG_AccFit,
                      dG_AccFitErr,
                      dG_FourAccFit,
                      dG_FourAccFitErr
                      )

fitDataSet = RooDataSet('fitDataSet','fitDataSet',fitArgSet)

print '*'*80
print '*'*20 + ' fitting original PDF on original data'
print '*'*80
ws['t_sig_dG'].setVal(0.05)
pdf.fitTo(origdata,RooFit.NumCPU(2))
dG_OrigFit.setVal(ws['t_sig_dG'].getVal())
dG_OrigFitErr.setVal(ws['t_sig_dG'].getError())
print '*'*80
print '*'*20 + ' fitting original PDF on inefficient data'
print '*'*80
ws['t_sig_dG'].setVal(0.05)
pdf.fitTo(data,RooFit.NumCPU(2))
dG_WrongFit.setVal(ws['t_sig_dG'].getVal())
dG_WrongFitErr.setVal(ws['t_sig_dG'].getError())
print '*'*80
print '*'*20 + ' fitting reverse engineerd moment PDF on inefficient data'
print '*'*80
ws['t_sig_dG'].setVal(0.05)
pdf_mom.fitTo(data,RooFit.NumCPU(2))
dG_AccFit.setVal(ws['t_sig_dG'].getVal())
dG_AccFitErr.setVal(ws['t_sig_dG'].getError())
print '*'*80
print '*'*20 + ' fitting Fourier expansion efficiency PDF on inefficient data'
print '*'*80
ws['t_sig_dG'].setVal(0.05)
pdf_eff.fitTo(data,RooFit.NumCPU(2))
dG_FourAccFit.setVal(ws['t_sig_dG'].getVal())
dG_FourAccFitErr.setVal(ws['t_sig_dG'].getError())

###################
### Fill FitSet ###
###################

fitDataSet.add(fitArgSet)

writefile = TFile.Open('fitDataSet.root',"RECREATE")
fitDataSet.Write('fitDataSet')
writefile.Close()


### make some plots...

canv = TCanvas()
canv.Divide(3,2)
for (i,var) in enumerate(angles) :
     canv.cd(1+i) # start argument of enumerate is python >= 2.6...
     plot = var.frame()
     origdata.plotOn(plot)
     pdf.plotOn(plot)
     plot.Draw()
for (i,var) in enumerate(angles) :
    canv.cd(4+i)
    plot2 = var.frame()
    data.plotOn(plot2)
    pdf_mom.plotOn(plot2,RooFit.LineColor(kRed))
    pdf_eff.plotOn(plot2,RooFit.LineColor(kGreen))
    plot2.Draw()
canv.Flush()

canv.SaveAs('canv.eps')
