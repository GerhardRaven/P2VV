from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi

from RooFitDecorators import *
import rootStyle
from ModelBuilders import _buildAngularFunction
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

##############################################################################
##########################   There we go!!!!! ################################
##############################################################################
signif = 1

from ModelBuilders import *

ws = RooWorkspace("ws")

useTransversityAngles = True

mcfilename = '/data/bfys/dveijk/MC/ReducedMCNTuple.root'

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))

ws.factory("{ helcosthetaK[-1,1], helcosthetaL[-1,1], helphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))

if useTransversityAngles: 
    ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
    angles = ws.set('transversityangles')
else:
    ws.defineSet("helicityangles","helcosthetaK,helcosthetaL,helphi")
    angles = ws.set('helicityangles')

ab = abasis(ws,angles)

ws.factory("{rz2[0.601,0.4,0.7],rperp2[0.16,0.1,0.5]")
ws.factory("RooFormulaVar::rpar2('1-@0-@1',{rz2,rperp2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")

ws.factory("{deltaz[0.],deltapar[2.5,%f,%f],deltaperp[-0.17,%f,%f]}"%(-2*pi,2*pi,-2*pi,2*pi))

##########################
### physics parameters ###
##########################
ws.factory("{#Gamma[0.68,0.4,0.9]}")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")

ws.factory("{t_sig_dG[0.060,-2,2]}")

ws.factory("{t_sig_dm[17.8,15,20]}")

ws.factory('{phis[-0.04,%f,%f]}'%(-pi,pi))

ws.factory("{expr::S('-1*sin(@0)',{phis}),expr::D('cos(@0)',{phis}),C[0]}")

###############################
### Experimental parameters ###
###############################
# For determination of efficiency from MC, put RooTruthModel, later replace this by realistic Resolution model. Does that imply building the JpsiPhi pdf again?

#ws.factory("RooTruthModel::tres_sig(t)")
ws.factory("RooGaussModel::tres_sig(t,mu[0],sigma[0.05])")

# For determination of efficiency from MC, put wtag to 0, later integrate out (or set to 0.5 as I used to do before). Does that imply rebuilding the JpsiPhi pdf?
ws.factory("{wtag[0.5,0.,0.501]}")

##################################### building the NEW pdf ###################################
# create observables
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

if useTransversityAngles:
    newpdf = buildJpsiphi(ws,'newpdf', True) 
else:
    newpdf = buildJpsiphi(ws,'newpdf', False)

ws.factory("{m[5200,5550]}")

if useTransversityAngles:
    ws.defineSet("observables","m,t,trcospsi,trcostheta,trphi,tagdecision")
else:
    ws.defineSet("observables","m,t,helcosthetaL,helcosthetaK,helphi,tagdecision")

MCdatafile = TFile(mcfilename)
NTupletree = MCdatafile.Get('MyTree')
# we take the subset of events with MC matching
tmpfile = TFile('tmp.root','RECREATE')

#Remove the non-matched candidates!!!
reducedtree = NTupletree.CopyTree('TRUEt>0')
reducedtree.SetDirectory(tmpfile)
reducedtree.Write()
MCdata = RooDataSet('MCdata','MCdata',reducedtree,ws.set('observables'),'t==t && m==m')
ws.put(MCdata)
print 'Number of truth matched MC events', MCdata.numEntries()

pdf = ws.pdf('newpdf')

allObs = pdf.getObservables( MCdata.get() )

angles = pdf.getObservables( MCdata.get() )
angles.remove(ws.var('t'))
angles.remove(ws.cat('tagdecision'))
angles.remove(ws.var('m'))
angles.remove(ws.var('wtag'))

print 'angles: ', [ i.GetName() for i in angles ]

print 'observables:', [ i.GetName() for i in allObs ]

# define the moments used to describe the efficiency
# for this, we need the PDF used to generate the data
marginalObs = pdf.getObservables( MCdata.get() )
marginalObs.remove( angles )

print 'marginal obs: ', [i.GetName() for i in marginalObs]
# marginalize pdf over 'the rest' so we get the normalization of the moments right...
pdf_marginal = pdf.createProjection(marginalObs)

# compute the 'canonical' six moments
bnames = [ 'AzAz','AparApar','AperpAperp','AparAperp','AzAperp','AzApar' ]
sixmom = [ EffMoment( ws['%s_basis'%n], 1., pdf_marginal, allObs ) for n in bnames ]
computeMoments(MCdata,sixmom)
xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in sixmom ] )

print 'Direct Moments xi_m =', xi_m

def norm_xi( d ) :
    n =  (d['AparApar_basis'] + d['AzAz_basis'] + d['AperpAperp_basis'])/3
    for i in d.iterkeys() : d[i] = d[i]/n

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
           , 'AzAperp_basis'    : - sqrt(6./5.)* 3*pi/32 *( c[(1,2, 1)] - c[(3,2, 1)]/4  )#- 5*c[(5,2, 1)]/128  - 7*c[(7,2, 1)]/512 - 105*c[(9,2, 1)]/16384) 
           , 'AzApar_basis'     :   sqrt(6./5.)* 3*pi/32 *( c[(1,2,-2)] - c[(3,2,-2)]/4  )#- 5*c[(5,2,-2)]/128  - 7*c[(7,2,-2)]/512 - 105*c[(9,2,-2)]/16384)
           }


norm_xi(xi_m)

print 'Normalized Direct Moments norm(xi_m) =', xi_m

# compute the Fourier series for the efficiency
moments = []
ab = abasis(ws,angles)
for (i,l) in product(range(4),range(4)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, pdf_marginal, allObs ) for m in range(-l,l+1) ]

# loop over all data, determine moments
computeMoments(MCdata,moments)

# compute the 'canonical' moments given the Fourier series
c = dict()
for m in moments : c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()
c_000 = c[(0,0,0)]

if useTransversityAngles:

    xi_c = compute_moments( c )
    # normalize moments and compare
   
    norm_xi(xi_c)
    norm_xi(xi_m)
    for name in xi_c.iterkeys() :
        print '%s : direct moment: %s ;  from Fourier series: %s ; ratio = %s ' % \
          ( name, xi_m[name], xi_c[name], xi_m[name]/xi_c[name])

# build PDF using the Fourier series efficiency...
# TODO: normalize relative to 0000 so that c_000 = 1
print 'using the following terms in Fourier expansion: '

for n,c in [ ( m.basis().GetName() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] :
    print '%s : %s ' % (n,c)

if False:
    pdf_eff = buildEff_x_PDF(ws,'fourier_eff',pdf,[ ( m.basis() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] )

    print 'effTimesPdfName: ', pdf_eff.GetName() 
    #ws.put( pdf_eff )

    # make some plots...
    c = TCanvas("c1","Use MC data to determine efficiency",900,300)
    c.Divide(3,1);
    for (i,var) in enumerate(angles) :
        c.cd(1+i) # start argument of enumerate is python >= 2.6...
        plot = var.frame()
        MCdata.plotOn(plot); 
        pdf.plotOn(plot,RooFit.LineColor(kBlue))
        pdf_eff.plotOn(plot,RooFit.LineColor(kRed))
        #pdf_mom.plotOn(plot,RooFit.LineColor(kOrange))
        plot.Draw()

    c.Flush()
    c.Update()
    c.Print("anglesJpsiKstarMC.eps")       

##############################
### Now read data and fit! ###
##############################

t = ws.var('t')

tmin = 0.3
tmax = 14.
t.setRange(tmin,tmax)

if useTransversityAngles:
    trcostheta = ws.var('trcostheta')
    trcospsi = ws.var('trcospsi')
    trphi = ws.var('trphi')
else:
    cthetaK = ws.var('helcosthetaK')
    cthetaL = ws.var('helcosthetaL')
    phi = ws.var('helphi')

tagdecision = ws.cat('tagdecision')

#wide
m = ws.var('m')

mwidemin = 5200
mwidemax = 5550
m.setRange(mwidemin,mwidemax)
#narrow
mnarrowmin = 5321.67
mnarrowmax = 5411.67
#m.setRange(mnarrowmin,mnarrowmax)

#################
### Load Data ###
#################

#Using My file with latest tagging for tagged fit
datafile = TFile('/data/bfys/dveijk/Data/Bs2JpsiPhiForTaggedFit.root')
NTupletree = datafile.Get('MyTree')

ws.set('observables').add(ws.var('wtag'))

if useTransversityAngles:
    data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m && trcospsi==trcospsi && trcostheta==trcostheta && trphi==trphi && tagdecision==tagdecision && wtag==wtag')
else:
    data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m && helcosthetaL==helcosthetaL && helcosthetaK==helcosthetaK && helphi==helphi && tagdecision==tagdecision && wtag ==wtag')

print 'Number of unbiased events', data.numEntries()

ws.put(data)
    
data.table(tagdecision).Print('v')

#######################
### Build the PDF's ###
#######################

ws.factory("Gaussian::m_sig(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")

#background B mass pdf
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

#background propertime 
# TODO: split resolution in tres_sig and tres_nonpsi!!
ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres_sig,SingleSided)")
ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres_sig,SingleSided)")
ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")

#background angles: 

#if useTransversityAngles:
#    ws.factory("Uniform::bkgang({trcospsi,trcostheta,trphi})")
#else:
#    ws.factory("Uniform::bkgang({helcosthetaL,helcosthetaK,helphi})")

#if useTransversityAngles:
#    ws.factory("Chebychev::bkg_trcospsi(trcospsi,{c0_trcospsi[-0.13,-1,1]})")
#    ws.factory("Chebychev::bkg_trcostheta(trcostheta,{c0_trcostheta[0.08,-1,1]})")
#    ws.factory("Chebychev::bkg_trphi(trphi,{c0_trphi[0.10,-1,1]})")
#    ws.factory("PROD::bkgang(bkg_trcospsi,bkg_trcostheta,bkg_trphi)")
#else:
#    ws.factory("Chebychev::bkg_helcosthetaK(helcosthetaK,{c0_helcosthetaK[-0.13,-1,1]})")
#    ws.factory("Chebychev::bkg_helcosthetaL(helcosthetaL,{c0_helcosthetaL[0.08,-1,1]})")
#    ws.factory("Chebychev::bkg_helphi(helphi,{c0_helphi[0.10,-1,1]})")
#    ws.factory("PROD::bkgang(bkg_helcosthetaL,bkg_helcosthetaK,bkg_helphi)")

#Build angular background with RooP2VVAngleBasis
_ba = lambda name,comp : _buildAngularFunction(ws,ab,name,comp)

_ba("Bkg_0000",  [ ( 0,0,0,0,1.) ] )

_ba("Bkg_0010",  [ ( 0,0,1,0,1.) ] )
_ba("Bkg_0011",  [ ( 0,0,1,1,1.) ] )
_ba("Bkg_001m1",  [ ( 0,0,1,-1,1.) ] )

_ba("Bkg_1000",  [ ( 1,0,0,0,1.) ] )
_ba("Bkg_2000",  [ ( 2,0,0,0,1.) ] )
_ba("Bkg_3000",  [ ( 3,0,0,0,1.) ] )

c0000 = RooRealVar('c0000','c0000',1.)

c0010 = RooRealVar('c0010','c0010',0.,-1.,1.)
c0011 = RooRealVar('c0011','c0011',0.,-1.,1.)
c001m1 = RooRealVar('c001m1','c001m1',0.,-1.,1.)

c1000 = RooRealVar('c1000','c1000',0.,-1.,1.)
c2000 = RooRealVar('c2000','c2000',0.,-1.,1.)
c3000 = RooRealVar('c3000','c3000',0.,-1.,1.)

#Flat background in all angles
bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis']),RooArgList(c0000))

#Only higher order in cospsi
#bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c1000,c2000,c3000))

#Add higher order in costheta
#bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c1000,c2000,c3000))

#Hiher roder in all angles
#bkg = RooRealSumPdf('bkg','bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_0011_basis'],ws['Bkg_001m1_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c0011,c001m1,c1000,c2000,c3000))

ws.put(bkg)

ws.factory("PROD::sig_pdf( m_sig, newpdf)")

#ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, bkgang)")
ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, bkg)")

ws.factory("SUM::pdf_ext(Nsig[543,0,1000]*sig_pdf,Nbkg[812,0,1000]*bkg_pdf)")
ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sig_pdf,bkg_pdf)")

pdf_ext = ws.pdf('pdf_ext')

print 'GOING TO BUILD THE ACCEPTANCE CORRECTED FULL PDF!'
angcorrpdf = buildEff_x_PDF(ws,'angcorrpdf',pdf_ext,[ ( m.basis() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] )
ws.put(angcorrpdf)

wsfile = TFile('TaggedWS.root','RECREATE')
ws.Write()
wsfile.Close()
