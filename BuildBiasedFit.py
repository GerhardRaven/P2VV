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

mcfilename =  '/data/bfys/dveijk/MC/2011/MC2011_ForAngAccComp.root'

datafilename = '/data/bfys/dveijk/DataJpsiPhi/2011/Pass3Version2.root'
datasetname = 'MyTree'

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,14.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))
ws.factory("{m[5200,5550],mdau1[3030,3150],mdau2[1007.46,1031.46]}")
ws.factory("{biased[Biased=+1,NotBiased=0],unbiased[Unbiased=+1,NotUnbiased=0],fullybiased[FullyBiased =+1,NotFullyBiased=0]}")

ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
angles = ws.set('transversityangles')

ab = abasis(ws,angles)

ws.factory("{rz2[0.601,0.4,0.7],rperp2[0.16,0.1,0.5],rs2[0.,0.,0.15]}")
ws.factory("RooFormulaVar::rpar2('1-@0-@1-@2',{rz2,rperp2,rs2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")
ws.factory("RooFormulaVar::rs('sqrt(@0)',{rs2})")

ws.factory("{deltaz[0.],deltapar[2.5,%f,%f],deltaperp[-0.17,%f,%f],deltas[0.,%f,%f]}"%(-2*pi,2*pi,-2*pi,2*pi,-2*pi,2*pi))

#Parametrize differently
#ws.factory("{NAzAz[1],NAparApar[0.44,-1,1],NAperpAperp[0.4,-1,1],ReAparAperp[-0.114,-1,1],ReAzAperp[0.08,-1,1],ReAzApar[-0.66,-1,1],ImAparAperp[0.410,-1,1],ImAzAperp[-0.63,-1,1]}")

ws.factory("expr::ReAz   ('@0    * cos(@1)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('@0    * sin(@1)',   {rz,deltaz})")
ws.factory("expr::ReApar ('@0    * cos(@1)', {rpar,deltapar})")
ws.factory("expr::ImApar ('@0    * sin(@1)', {rpar,deltapar})")
ws.factory("expr::ReAperp('@0    * cos(@1)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('@0    * sin(@1)',{rperp,deltaperp})")

ws.factory("expr::ReAs('@0 * cos(@1)',{rs,deltas})")
ws.factory("expr::ImAs('@0 * sin(@1)',{rs,deltas})")

# define the relevant combinations of strong amplitudes
ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,    ImAz                      })")  # |A_z|^2
ws.factory("expr::NAparApar  ('( @0 * @0 + @1 * @1 )',{ReApar,  ImApar                    })")  # |A_par|^2
ws.factory("expr::NAperpAperp('( @0 * @0 + @1 * @1 )',{ReAperp, ImAperp                   })")  # |A_perp|^2
ws.factory("expr::ReAparAperp('( @0 * @2 + @1 * @3 )',{ReApar,  ImApar,  ReAperp, ImAperp })")  # |A_par||A_perp| cos(delta_perp - delta_par)
ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,    ImAz,    ReAperp, ImAperp })")  # |A_z||A_perp|   cos(delta_perp - delta_z)
ws.factory("expr::ReAzApar   ('( @0 * @2 + @1 * @3 )',{ReAz,    ImAz,    ReApar,  ImApar  })")  # |A_z||A_par|    cos(delta_par  - delta_z)
ws.factory("expr::ImAparAperp('( @0 * @3 - @1 * @2 )',{ReApar,  ImApar,  ReAperp, ImAperp })")  # |A_par|A_perp|  sin(delta_perp - delta_par)
ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,    ImAz,    ReAperp, ImAperp })")  # |A_z||A_perp|   sin(delta_perp - delta_z)

ws.factory("expr::NAsAs      ('( @0 * @0 + @1 * @1 )',{ReAs,    ImAs                      })")  # |A_s|^2
ws.factory("expr::ReAsAz     ('( @0 * @2 + @1 * @3 )',{ReAs,    ImAs,    ReAz,    ImAz    })")  # |A_s||A_z| cos(delta_z - delta_s)
ws.factory("expr::ReAsApar   ('( @0 * @2 + @1 * @3 )',{ReAs,    ImAs,    ReApar,  ImApar  })")  # |A_s||A_par| cos(delta_par - delta_s)
ws.factory("expr::ImAsAperp  ('( @0 * @3 - @1 * @2 )',{ReAs,    ImAs,    ReAperp, ImAperp })")  # |A_s||A_perp| sin(delta_perp - delta_s)
ws.factory("expr::ImAsAz     ('( @0 * @3 - @1 * @2 )',{ReAs,    ImAs,    ReAz,    ImAz    })")  # |A_s||A_z| sin(delta_z - delta_s)
ws.factory("expr::ImAsApar   ('( @0 * @3 - @1 * @2 )',{ReAs,    ImAs,    ReApar,  ImApar  })")  # |A_s||A_par| sin(delta_par - delta_s)
    
##########################
### physics parameters ###
##########################
ws.factory("{#Gamma[0.68,0.4,0.9]}")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")

ws.factory("{t_sig_dG[0.06852,-2,2]}")

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

ws.factory("GaussModel::tres_3(t,tres_mu[-0.0027],tres_s3[0.513],tres_SF[1.00,0.5,1.5])")
ws.factory("GaussModel::tres_2(t,tres_mu,tres_s2[0.0853],tres_SF)")
ws.factory("GaussModel::tres_1(t,tres_mu,tres_s1[0.0434],tres_SF)")
ws.factory("AddModel::tres({tres_3,tres_2,tres_1},{tres_f3[0.0017],tres_f2[0.165]})")

#ws.factory("GaussModel::tres(t,tres_mean[0.0],tres_sigma[0.05])")

# For determination of efficiency from MC, put wtag to 0, later integrate out (or set to 0.5 as I used to do before). Does that imply rebuilding the JpsiPhi pdf?
ws.factory("tagomega[0.,0.,0.5]")
ws.factory("expr::wtag('tagomega',tagomega)")

#Build the MC PDF
MCpdf = buildJpsiphiSWave(ws,'MCpdf', True,'tres_MC')

#Build the signal PDF
newpdf = buildJpsiphiSWave(ws,'newpdf', True,'tres')

ws.defineSet("MCobservables","t,trcospsi,trcostheta,trphi,tagdecision,tagomega")

MCdatafile = TFile(mcfilename)
NTupletree = MCdatafile.Get('MyTree')

MCdata = RooDataSet('MCdata','MCdata',NTupletree,ws.set('MCobservables'),'t==t && trcospsi==trcospsi && trcostheta == trcostheta && trphi==trphi && tagomega == tagomega')
ws.put(MCdata)
print 'Number of MC events', MCdata.numEntries()

allObs = MCpdf.getObservables( MCdata.get() )
print 'MCobservables:', [ i.GetName() for i in allObs ]

angles = MCpdf.getObservables( MCdata.get() )
angles.remove(ws.var('t'))
angles.remove(ws.cat('tagdecision'))
angles.remove(ws.var('tagomega'))
print 'angles: ', [ i.GetName() for i in angles ]

# define the moments used to describe the efficiency
# for this, we need the PDF used to generate the data
#marginalObs = MCpdf.getObservables( MCdata.get() )
#marginalObs.remove(angles)
#marginalObs.remove(ws['t'])
#print 'marginal obs: ', [i.GetName() for i in marginalObs]
# marginalize pdf over 'the rest'
#MCpdf_marginal = MCpdf.createProjection(marginalObs)

# compute the 'canonical' ten moments
bnames = [ 'AzAz','AparApar','AperpAperp','AparAperp','AzAperp','AzApar','AsAs','AsAz','AsApar','AsAperp']
tenmom = [ EffMoment( ws['%s_basis'%n], 1., MCpdf, angles) for n in bnames ]

computeMoments(MCdata,tenmom)
xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in tenmom ] )
print 'Direct Moments xi_m =', xi_m

def norm_xi( d ) :
    n =  (1.)/(8.*sqrt(pi))
    for i in d.iterkeys() : d[i] = d[i]*n

norm_xi(xi_m)
print 'Normalized direct moments xi_m =', xi_m

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
    return { 'AzAz_basis'       :   4*( c[(0,0,0)]+2*c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]+2*c[(2,2,0)]/5) - sqrt(3./20)*(c[(0,2,2)]+2*c[(2,2,2)]/5)  )
             ,'AparApar_basis'   :   4*( c[(0,0,0)]-  c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]-  c[(2,2,0)]/5) + sqrt(3./20)*(c[(0,2,2)]-  c[(2,2,2)]/5)  )
             , 'AperpAperp_basis' :   4*( c[(0,0,0)]-  c[(2,0,0)]/5 - sqrt(1./ 5)*( c[(0,2,0)]-  c[(2,2,0)]/5 ) )
             , 'AparAperp_basis'  :   4*sqrt(3./5.)*( c[(0,2,-1)] - c[(2,2,-1)]/5 )
             , 'AzApar_basis'     :   4*sqrt(6./5.)* 3*pi/32 *( c[(1,2,-2)] - c[(3,2,-2)]/4)# - 5*c[(5,2,-2)]/128  - 7*c[(7,2,-2)]/512 - 105*c[(9,2,-2)]/16384)
             , 'AzAperp_basis'    : - 4*sqrt(6./5.)* 3*pi/32 *( c[(1,2, 1)] - c[(3,2, 1)]/4)# - 5*c[(5,2, 1)]/128  - 7*c[(7,2, 1)]/512 - 105*c[(9,2, 1)]/16384)
             , 'AsAs_basis'       :   2*(2*c[(0,0,0)]+sqrt(1./5)*c[(0,2,0)]-sqrt(3./5)*c[(0,2,2)])
             , 'AsApar_basis'     :   12*sqrt(2./5.)*pi/8 *( c[(0,2,-2)] - c[(2,2,-2)]/8)# - c[(4,2,-2)]/64 - 5*pi*c[(6,2,-2)]/1024 -35*pi*c[(8,2,-2)]/16384)
             , 'AsAperp_basis'     :   -12*sqrt(2./5.)*pi/8 *( c[(0,2,1)] - c[(2,2,1)]/8)# - c[(4,2,1)]/64 - 5*pi*c[(6,2,1)]/1024 -35*pi*c[(8,2,1)]/16384)
             , 'AsAz_basis'       :   (2./3)*(4*sqrt(3)*c[(1,0,0)]+2*sqrt(3./5)*c[(1,2,0)]-6*sqrt(1./5)*c[(1,2,2)])
           }

# compute the Fourier series for the efficiency
moments = []
ab = abasis(ws,angles)
for (i,l) in product(range(4),range(3)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, MCpdf, angles ) for m in range(-l,l+1) ]

# loop over all data, determine moments
computeMoments(MCdata,moments)

# compute the 'canonical' moments given the Fourier series
c = dict()
for m in moments :
    c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()

xi_c = compute_moments( c )
norm_xi(xi_c)

for name in xi_c.iterkeys() :
    print '%s : direct moment: %s ;  from Fourier series: %s ; ratio = %s ' % \
          ( name, xi_m[name], xi_c[name], xi_m[name]/xi_c[name])

print 'using the following terms in Fourier expansion: '

for n,coeff in [ ( m.basis().GetName() , m.coefficient() ) for m in moments if m.significance()>signif] :
    print '%s : %s ' % (n,coeff)


if False:
    MCpdf_eff = buildEff_x_PDF(ws,'fourier_eff',MCpdf,[ ( m.basis() , m.coefficient() ) for m in moments if m.significance()>signif] )
    ws.put( MCpdf_eff )

    effpdfs = RooArgList()
    effcoeffs = RooArgList()
    list = [ [ m.basis() , c[(m.basis().i(),m.basis().l(),m.basis().m())] ] for m in moments if m.significance()>signif]

    for i in list:
        ws.factory("const_%s[%s]"%(i[0].GetName(),i[1]))
        effpdfs.add(i[0])
        effcoeffs.add(ws['const_%s'%(i[0].GetName())])

    effshape = RooAddition('effshape','effshape',effpdfs,effcoeffs)

    # make some plots...
    c = TCanvas("c1","Use MC data to determine efficiency",900,600)
    c.Divide(3,2);
    for (i,var) in enumerate(angles) :
        c.cd(1+i)
        plot2 = var.frame()
        effshape.plotOn(plot2)
        plot2.Draw()
        c.cd(4+i) # start argument of enumerate is python >= 2.6...
        plot = var.frame()
        MCdata.plotOn(plot)
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

#################
### Load Data ###
#################

datafile = TFile(datafilename)
NTupletree = datafile.Get(datasetname)
ws.defineSet("observables","t,trcospsi,trcostheta,trphi,m,tagdecision")

ws.set('observables').add(ws['tagomega'])
ws.set('observables').add(ws['biased'])
ws.set('observables').add(ws['unbiased'])
ws.set('observables').add(ws['mdau1'])
ws.set('observables').add(ws['mdau2'])

data = RooDataSet('data','data',NTupletree,ws.set('observables'),'t==t && m ==m && trcospsi==trcospsi && trcostheta==trcostheta && trphi==trphi && tagdecision==tagdecision && tagomega==tagomega && biased==biased && unbiased==unbiased && mdau1 == mdau1 && mdau2 == mdau2')

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

##############################
### Proper Time Acceptance ###
##############################
ws.factory("expr::effshape('1/(1+(a*t)**(-c))',t,a[1.45],c[2.37])")
effhist = ws['effshape'].createHistogram('effhist',ws['t'],RooFit.Binning(10,ws['t'].getMin(),ws['t'].getMax()))
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effhist)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,effdatahist)")

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
#Single exponential for biased sample!!!!
ws.factory("RooDecay::t_bkg_B(t,t_bkg_B_tau[0.2,0.01,2.0],tres,SingleSided)")

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
#bkg_ang = RooRealSumPdf('bkg_ang','bkg_ang',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c1000,c2000,c3000))

#Add higher order in costheta
#bkg_ang = RooRealSumPdf('bkg_ang','bkg_ang',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c1000,c2000,c3000))

#Hiher roder in all angles
#bkg_ang = RooRealSumPdf('bkg_ang','bkg_ang',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_0011_basis'],ws['Bkg_001m1_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c0011,c001m1,c1000,c2000,c3000))

ws['m'].setRange("leftsideband",5200,5330)
ws['m'].setRange("rightsideband",5410,5550)

sidebanddata = ws['data'].reduce(RooFit.CutRange('leftsideband'))
rightsidebanddata = ws['data'].reduce(RooFit.CutRange('rightsideband'))
sidebanddata.append(rightsidebanddata)
sidebanddata.SetName('sidebanddata')
sidebanddata.SetTitle('sidebanddata')

#Set binning to famous 7x5x9
ws['trcospsi'].setBins(7)
ws['trcostheta'].setBins(5)
ws['trphi'].setBins(9)

sidebandhist = RooDataHist('sidebandhist','sidebandhist', RooArgSet(ws['trcospsi'],ws['trcostheta'],ws['trphi']), sidebanddata)
#ang_bkg =  RooHistPdf('ang_bkg','ang_bkg',RooArgSet(ws['trcospsi'],ws['trcostheta'],ws['trphi']), sidebandhist)
#ws.put(ang_bkg)

ws.factory("Uniform::ang_bkg({trcostheta,trcospsi,trphi})")


###########################
### Putting it together ###
###########################
# Mass only PDF
ws.factory("SUM::m_pdf(Nsig_B[1000,0,2000]*m_sig,Nbkg_B[1000,0,2000]*m_bkg)")
#ws['m_pdf'].fitTo(ws['data'],RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False))

ws.factory("PROD::sigpdf( m_sig, newpdf)")
ws.factory("PROD::bkgpdf( m_bkg, t_bkg_B, ang_bkg)")

ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sigpdf,bkgpdf)")
ws.factory("SUM::pdf_ext(Nsig_B*sigpdf,Nbkg_B*bkgpdf)")

ws.factory("EffHistProd::accpdf(pdf,effpdf)")
ws.factory("EffHistProd::accpdf_ext(pdf_ext,effpdf)")

print 'GOING TO BUILD THE ACCEPTANCE CORRECTED FULL PDF!'
timeangcorrpdf = buildEff_x_PDF(ws,'angcorr',ws['accpdf_ext'],[ ( m.basis() , m.coefficient() ) for m in moments if m.significance()>signif] )
ws.put(timeangcorrpdf)


###########################
### Tagging systematics ###
###########################
ws.factory("Gaussian::p0(p0var[0.,0.5],p0mean[0.384],p0sigma[0.010])")
ws.factory("Gaussian::p1(p1var[-2.,2.],p1mean[1.037],p1sigma[0.081])")
ws.factory("expr:wtag_syst('p0var+p1var*(tagomega-etamean)',tagomega,etamean[0.379],p0var,p1var)")

wtag = ws.function('wtag')
wtag_syst = ws.function('wtag_syst')

#timecorrpdf
customizer = RooCustomizer(ws['accpdf_ext'],'inc_tag_syst')
customizer.replaceArg( wtag, wtag_syst )
timecorrtagpdf = customizer.build()
getattr(ws,'import')(timecorrtagpdf,RooFit.RecycleConflictNodes())

#timeangcorrpdf
customizer = RooCustomizer(ws['accpdf_ext_angcorr'],'inc_tag_syst')
customizer.replaceArg( wtag, wtag_syst )
timeangcorrtagpdf = customizer.build()
getattr(ws,'import')(timeangcorrtagpdf,RooFit.RecycleConflictNodes())

#Constrain deltams
ws.factory("Gaussian::dmsconstraint(t_sig_dm,t_sig_dm_mean[17.63],t_sig_dm_sigma[0.11])")

#Constrain tres_SF
ws.factory("Gaussian::tres_SFconstraint(tres_SF,tres_SF_mean[1.00],tres_SF_sigma[0.04])")

wsfile = TFile('BiasedWS.root','RECREATE')
ws.Write()
wsfile.Close()
