#Build 2011 Unbiased fit with Coeff parameterization
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

from ModelBuilders import *

############################
### Flags and rootfiles  ###
############################

signif = 1
swave = True

#mcfilename = '/data/bfys/dveijk/MC/ReducedMCNTuple.root'
mcfilename = '/data/bfys/dveijk/MC/2011/MC2011.root'
#datafilename = '/data/bfys/dveijk/DataJpsiPhi/2011/FromWouter.root'
#datafilename = '/data/bfys/dveijk/DataJpsiPhi/2011/Pass3.root'
#datafilename = '/data/bfys/dveijk/DataJpsiPhi/2011/Pass3WidePhiMass.root'
datafilename = '/data/bfys/dveijk/DataJpsiPhi/2011/Pass3Version2.root'
datasetname = 'MyTree'

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))
ws.factory("{m[5200,5550],mdau2[986,1050]}")
ws.factory("{biased[Biased=+1,NotBiased=0],unbiased[Unbiased=+1,NotUnbiased=0],fullybiased[FullyBiased =+1,NotFullyBiased=0]}")

ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
angles = ws.set('transversityangles')

ab = abasis(ws,angles)

_Az2         = 0.60
_deltaz      = 0.
_Aperp2      = 0.16
_deltaperp   = -0.17
_As2         = 0.02
_deltas      = 0.5
_Apar2       = 1. - _As2 - _Az2 - _Aperp2
_deltapar    = 2.5

#Normalize:
_rz2      = 1.
_rs2      = _As2/_Az2
_rperp2   = _Aperp2/_Az2
_rpar2    = _Apar2/_Az2 

def calculateRe(firstamp2,secondamp2,firstdelta,seconddelta):
   Re = sqrt(firstamp2*secondamp2)*cos(seconddelta-firstdelta)
   return Re

def calculateIm(firstamp2,secondamp2,firstdelta,seconddelta):
   Im = sqrt(firstamp2*secondamp2)*sin(seconddelta-firstdelta)
   return Im

ws.factory("{NAzAz[%f],NAparApar[%f,-1,1],NAperpAperp[%f,-1,1],NAsAs[%f,-1,1]}"%(_rz2,_rpar2,_rperp2,_rs2))

ws.factory("{ReAparAperp[%f,-1,1]}"%(calculateRe(_rpar2,_rperp2,_deltapar,_deltaperp)))
ws.factory("{ReAzAperp[%f,-1,1]}"%(calculateRe(_rz2,_rperp2,_deltaz,_deltaperp)))
ws.factory("{ReAzApar[%f,-1,1]}"%(calculateRe(_rz2,_rpar2,_deltaz,_deltapar)))

ws.factory("{ImAparAperp[%f,-1,1]}"%(calculateIm(_rpar2,_rperp2,_deltapar,_deltaperp)))
ws.factory("{ImAzAperp[%f,-1,1]}"%(calculateIm(_rz2,_rperp2,_deltaz,_deltaperp)))
ws.factory("{ImAsAperp[%f,-1,1]}"%(calculateIm(_rs2,_rperp2,_deltas,_deltaperp)))
  
if swave:
    ws.factory("{ReAsAz[%f,-1,1]}"%(calculateRe(_rs2,_rz2,_deltas,_deltaz)))
    ws.factory("{ReAsApar[%f,-1,1]}"%(calculateRe(_rs2,_rpar2,_deltas,_deltapar)))
    ws.factory("{ImAsAz[%f,-1,1]}"%(calculateIm(_rs2,_rz2,_deltas,_deltaz)))
    ws.factory("{ImAsApar[%f,-1,1]}"%(calculateIm(_rs2,_rpar2,_deltas,_deltapar)))
  
##########################
### physics parameters ###
##########################
ws.factory("{#Gamma[0.68,0.4,0.9]}")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")

ws.factory("{t_sig_dG[0.060,-2,2]}")

ws.factory("{t_sig_dm[17.8,15,20]}")

#Parametrize in terms of S,D,C:
#ws.factory("{S[-0.04,-2.,2.],D[1.,-2.,2.],C[0]}")

#Parametrize in terms of phis:
ws.factory('{phis[-0.04,%f,%f]}'%(-2*pi,2*pi))
ws.factory("{expr::S('-1*sin(@0)',{phis}),expr::D('cos(@0)',{phis}),C[0]}")

#Parametrize in terms of lambda
#ws.factory('{phis[-0.04,%f,%f]}'%(-2*pi,2*pi))
#ws.factory("{relambda[1.,-10,10],imlambda[-0.04,-10,10]}")
#ws.factory("{expr::relambda('cos(@0)',{phis}),expr::imlambda('-1*sin(@0)',{phis})}")
#ws.factory("{expr::lambda2('@0*@0+@1*@1',{relambda,imlambda}),expr::S('(2*@0)/(1+@1)',{imlambda,lambda2}),expr::D('(2*@0)/(1+@1)',{relambda,lambda2}),expr::C('(1-@0)/(1+@0)',{lambda2})}")
#ws.factory("{expr::lambda2('@0*@0+@1*@1',{relambda,imlambda}),expr::S('(2*@0)/(1+@1)',{imlambda,lambda2}),expr::D('(2*@0)/(1+@1)',{relambda,lambda2}),C[0]}")
#ws.factory("{lambda2[1],expr::S('(2*@0)/(1+@1)',{imlambda,lambda2}),expr::D('(2*@0)/(1+@1)',{relambda,lambda2}),expr::C('(1-@0)/(1+@0)',{lambda2})}")

###############################
### Experimental parameters ###
###############################
# For determination of efficiency from MC, put RooTruthModel, later replace this by realistic Resolution model. Does that imply building the JpsiPhi pdf again? Yes, you should, but you never did it yet. Apparantly it doesn't matter for the result, but fix it!!!

ws.factory("RooTruthModel::tres_MC(t)")

#2011
ws.factory("GaussModel::tres_3(t,tres_mu[0],tres_s3[0.513],tres_SF[1.00,0.9,1.1])")
ws.factory("GaussModel::tres_2(t,tres_mu,tres_s2[0.0853],tres_SF)")
ws.factory("GaussModel::tres_1(t,tres_mu,tres_s1[0.0434],tres_SF)")
ws.factory("AddModel::tres({tres_3,tres_2,tres_1},{tres_f3[0.0017],tres_f2[0.165]})")

#ws.factory("GaussModel::tres(t,tres_mean[0.0],tres_sigma[0.06])")

# For determination of efficiency from MC, put wtag to 0, later integrate out (or set to 0.5 as I used to do before). Does that imply rebuilding the JpsiPhi pdf?
ws.factory("tagomega[0.,0.,0.501]")
ws.factory("expr::wtag('tagomega',tagomega)")

#Build the MC PDF
if swave:
    MCpdf = buildJpsiphiSWave(ws,'MCpdf', True,'tres_MC')
else:
    MCpdf = buildJpsiphi(ws,'MCpdf', True,'tres_MC')

#Build the signal PDF
if swave:
    newpdf = buildJpsiphiSWave(ws,'newpdf', True,'tres')
else:
    newpdf = buildJpsiphi(ws,'newpdf', True,'tres')

ws.defineSet("MCobservables","t,trcospsi,trcostheta,trphi")

MCdatafile = TFile(mcfilename)
NTupletree = MCdatafile.Get('MyTree')
## # we take the subset of events with MC matching
## tmpfile = TFile('tmp.root','RECREATE')

## #Remove the non-matched candidates!!!
## reducedtree = NTupletree.CopyTree('TRUEt>0')
## reducedtree.SetDirectory(tmpfile)
## reducedtree.Write()
## MCdata = RooDataSet('MCdata','MCdata',reducedtree,ws.set('observables'),'t==t && m==m')
MCdata = RooDataSet('MCdata','MCdata',NTupletree,ws.set('MCobservables'),'t==t && trcospsi==trcospsi && trcostheta == trcostheta && trphi==trphi')
ws.put(MCdata)
print 'Number of MC events', MCdata.numEntries()

allObs = MCpdf.getObservables( MCdata.get() )

angles = MCpdf.getObservables( MCdata.get() )
angles.remove(ws.var('t'))
angles.remove(ws.cat('tagdecision'))

print 'angles: ', [ i.GetName() for i in angles ]

print 'MCobservables:', [ i.GetName() for i in allObs ]

# define the moments used to describe the efficiency
# for this, we need the PDF used to generate the data
marginalObs = MCpdf.getObservables( MCdata.get() )
marginalObs.remove( angles )

print 'marginal obs: ', [i.GetName() for i in marginalObs]

# marginalize pdf over 'the rest' so we get the normalization of the moments right...
MCpdf_marginal = MCpdf.createProjection(marginalObs)

# compute the 'canonical' six moments
bnames = [ 'AzAz','AparApar','AperpAperp','AparAperp','AzAperp','AzApar' ]
sixmom = [ EffMoment( ws['%s_basis'%n], 1., MCpdf_marginal, allObs ) for n in bnames ]
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
for (i,l) in product(range(4),range(3)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, MCpdf_marginal, allObs ) for m in range(-l,l+1) ]

# loop over all data, determine moments
computeMoments(MCdata,moments)

# compute the 'canonical' moments given the Fourier series
c = dict()
for m in moments : c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()
c_000 = c[(0,0,0)]



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
    MCpdf_eff = buildEff_x_PDF(ws,'fourier_eff',MCpdf,[ ( m.basis() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] )

    print 'effTimesPdfName: ', MCpdf_eff.GetName() 
    #ws.put( MCpdf_eff )

    # make some plots...
    c = TCanvas("c1","Use MC data to determine efficiency",900,300)
    c.Divide(3,1);
    for (i,var) in enumerate(angles) :
        c.cd(1+i) # start argument of enumerate is python >= 2.6...
        plot = var.frame()
        MCdata.plotOn(plot); 
        MCpdf.plotOn(plot,RooFit.LineColor(kBlue))
        MCpdf_eff.plotOn(plot,RooFit.LineColor(kRed))
        #MCpdf_mom.plotOn(plot,RooFit.LineColor(kOrange))
        plot.Draw()

    c.Flush()
    c.Update()
    c.Print("AngularAcceptanceCorrectionMC.eps")

#####################
### Now read data ###
#####################
tmin = 0.3
tmax = 14.
ws['t'].setRange(tmin,tmax)

#wide
mwidemin = 5200
mwidemax = 5550
ws['m'].setRange(mwidemin,mwidemax)
#narrow
mnarrowmin = 5321.67
mnarrowmax = 5411.67
#m.setRange(mnarrowmin,mnarrowmax)

#################
### Load Data ###
#################

datafile = TFile(datafilename)
NTupletree = datafile.Get(datasetname)
ws.defineSet("observables","t,trcospsi,trcostheta,trphi,m,tagdecision")

ws.set('observables').add(ws['tagomega'])
ws.set('observables').add(ws['biased'])
ws.set('observables').add(ws['unbiased'])
ws.set('observables').add(ws['mdau2'])

data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m && trcospsi==trcospsi && trcostheta==trcostheta && trphi==trphi && tagdecision==tagdecision && tagomega==tagomega && biased==biased && unbiased==unbiased && mdau2 == mdau2')

print 'Number of events', data.numEntries()

ws.put(data)
    
data.table(ws['tagdecision']).Print('v')
print 'TAGDECISION FOR UNBIASED EVENTS ONLY'
data.table(ws['tagdecision'],'unbiased == 1').Print('v')

print 'TRIGGER SUMMARY FOR ALL EVENTS'
data.table(ws['biased']).Print('v')
data.table(ws['unbiased']).Print('v')

# Make SuperCategory from (triggeredByUnbiasedHlt1AndHlt2,triggeredByBiasedHlt1AndHlt2)
TypeCat = RooSuperCategory('TypeCat','TypeCat',RooArgSet(ws['unbiased'],ws['biased']))
data.table(TypeCat).Print('v')
ws.put(TypeCat)
# Make MappedCategory from SuperCategory to split in unbiased and fullybiased
fitcat = RooMappedCategory('fitcat','fitcat',ws['TypeCat'],'00')#'00' means NotUnbiased && NotBiased
fitcat.map("{Unbiased;NotBiased}","AllUnbiased")
fitcat.map("{Unbiased;Biased}","AllUnbiased") 
fitcat.map("{NotUnbiased;Biased}","FullyBiased")

addcolumn = ws['data'].addColumn(fitcat)
addcolumn.SetName('fitcat')
ws.put(addcolumn)

ws['data'].table(ws['fitcat']).Print('v')

ws['fitcat'].setRange("unbiased","AllUnbiased")
ws['fitcat'].setRange("fullybiased","FullyBiased")

unbiaseddata = ws['data'].reduce(RooFit.CutRange('unbiased'))
fullybiaseddata = ws['data'].reduce(RooFit.CutRange('fullybiased'))

unbiaseddata.SetTitle('unbiaseddata')
unbiaseddata.SetName('unbiaseddata')
ws.put(unbiaseddata)
fullybiaseddata.SetTitle('fullybiaseddata')
fullybiaseddata.SetName('fullybiaseddata')
ws.put(fullybiaseddata)

#######################
### Build the PDF's ###
#######################
# Signal mass
ws.factory("Gaussian::m_sig_1(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")
ws.factory("expr::m_sig_sigma_2('2.14*@0',{m_sig_sigma_1})")
ws.factory("Gaussian::m_sig_2(m,m_sig_mean,m_sig_sigma_2)")
ws.factory("SUM::m_sig(m_bkg_f[0.83]*m_sig_1,m_sig_2)")

# Bkg mass
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

# Bkg time 
ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")
ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres,SingleSided)")
ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")

# Bkg angles
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
#ang_bkg = RooRealSumPdf('ang_bkg','ang_bkg',RooArgList(ws['Bkg_0000_basis']),RooArgList(c0000))
#ws.put(ang_bkg)

#Only higher order in cospsi
#ang_bkg = RooRealSumPdf('ang_bkg','ang_bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c1000,c2000,c3000))

#Add higher order in costheta
#ang_bkg = RooRealSumPdf('ang_bkg','ang_bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c1000,c2000,c3000))

#Hiher roder in all angles
#ang_bkg = RooRealSumPdf('ang_bkg','ang_bkg',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_0011_basis'],ws['Bkg_001m1_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c0011,c001m1,c1000,c2000,c3000))

ws.factory("Uniform::ang_bkg({trcostheta,trcospsi,trphi})")

###########################
### Putting it together ###
###########################
# Mass only PDF
ws.factory("SUM::m_pdf(Nsig[1000,0,16000]*m_sig,Nbkg[5000,0,16000]*m_bkg)")
#ws['m_pdf'].fitTo(ws['data'],RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False))

# Signal PDF
ws.factory("PROD::sig_pdf( m_sig, newpdf)")

# Bkg PDF
ws.factory("PROD::bkg_pdf( m_bkg, t_bkg, ang_bkg)")

ws.factory("SUM::pdf_ext(Nsig*sig_pdf,Nbkg*bkg_pdf)")
ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sig_pdf,bkg_pdf)")

print 'GOING TO BUILD THE ANGULAR ACCEPTANCE CORRECTED FULL PDF!'
angcorrpdf = buildEff_x_PDF(ws,'angcorr',ws['pdf_ext'],[ ( m.basis() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] )
ws.put(angcorrpdf)

#2011 values
ws.factory("Gaussian::p0(p0var[0.,0.5],p0mean[0.384],p0sigma[0.010])")
ws.factory("Gaussian::p1(p1var[-2.0,2.0],p1mean[1.037],p1sigma[0.081])")
ws.factory("expr:wtag_syst('p0var+p1var*(tagomega-etamean)',tagomega,etamean[0.379],p0var,p1var)")

###########################
### Tagging systematics ###
###########################

wtag = ws.function('wtag')
wtag_syst = ws.function('wtag_syst')

#pdf
customizer = RooCustomizer(ws['pdf_ext'],'inc_tag_syst')
customizer.replaceArg( wtag, wtag_syst )
tagpdf = customizer.build()
getattr(ws,'import')(tagpdf,RooFit.RecycleConflictNodes())

#angcorrpdf
customizer = RooCustomizer(angcorrpdf,'inc_tag_syst')
customizer.replaceArg( wtag, wtag_syst )
angcorrtagpdf = customizer.build()
getattr(ws,'import')(angcorrtagpdf,RooFit.RecycleConflictNodes())

#Constrain deltams
ws.factory("Gaussian::dmsconstraint(t_sig_dm,t_sig_dm_mean[17.63],t_sig_dm_sigma[0.11])")

#Constrain tres_SF
ws.factory("Gaussian::tres_SFconstraint(tres_SF,tres_SF_mean[1.00],tres_SF_sigma[0.02])")

wsfile = TFile('UnbiasedWS_Coeff.root','RECREATE')
ws.Write()
wsfile.Close()

