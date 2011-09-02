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

mcfilename = '/data/bfys/dveijk/MC/2011/MC2011.root'
#mcfilename = '/data/bfys/dveijk/MC/ReducedMCNTuple.root'
acceptancecurvefilename = '/data/bfys/dveijk/DataJpsiPhi/2011/effhistoOverlap.root'

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,14.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))
ws.factory("{m[5200,5550],mdau2[986,1050]}")
ws.factory("{biased[Biased=+1,NotBiased=0],unbiased[Unbiased=+1,NotUnbiased=0],fullybiased[FullyBiased =+1,NotFullyBiased=0]}")

ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
angles = ws.set('transversityangles')

ab = abasis(ws,angles)

if swave:
    ws.factory("{rz2[0.601,0.4,0.7],rperp2[0.16,0.1,0.5],rs2[0.00001,0.,0.15]}")
    ws.factory("RooFormulaVar::rpar2('1-@0-@1-@2',{rz2,rperp2,rs2})")
    ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
    ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
    ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")
    ws.factory("RooFormulaVar::rs('sqrt(@0)',{rs2})")

    ws.factory("{deltaz[0.],deltapar[2.5,%f,%f],deltaperp[-0.17,%f,%f],deltas[0.5,%f,%f]}"%(-2*pi,2*pi,-2*pi,2*pi,-2*pi,2*pi))
else:
    ws.factory("{rz2[0.601,0.4,0.7],rperp2[0.16,0.1,0.5]}")
    ws.factory("RooFormulaVar::rpar2('1-@0-@1',{rz2,rperp2})")
    ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
    ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
    ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")
 
    ws.factory("{deltaz[0.],deltapar[2.5,%f,%f],deltaperp[-0.17,%f,%f]}"%(-2*pi,2*pi,-2*pi,2*pi))

#Parametrize differently
#ws.factory("{NAzAz[1],NAparApar[0.44,-1,1],NAperpAperp[0.4,-1,1],ReAparAperp[-0.114,-1,1],ReAzAperp[0.08,-1,1],ReAzApar[-0.66,-1,1],ImAparAperp[0.410,-1,1],ImAzAperp[-0.63,-1,1]}")

ws.factory("expr::ReAz   ('@0    * cos(@1)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('@0    * sin(@1)',   {rz,deltaz})")
ws.factory("expr::ReApar ('@0    * cos(@1)', {rpar,deltapar})")
ws.factory("expr::ImApar ('@0    * sin(@1)', {rpar,deltapar})")
ws.factory("expr::ReAperp('@0    * cos(@1)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('@0    * sin(@1)',{rperp,deltaperp})")

if swave:
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

if swave:
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

#ws.factory("GaussModel::tres_3(t,tres_mu_fit[-0.0016],tres_s3_fit[0.183])")
#ws.factory("GaussModel::tres_2(t,tres_mu_fit,tres_s2_fit[0.06464])")
#ws.factory("GaussModel::tres_1(t,tres_mu_fit,tres_s1_fit[0.0337])")
#ws.factory("AddModel::tres({tres_3,tres_2,tres_1},{tres_f3[0.017],tres_f2[0.456],tres_f1[0.527]})")

ws.factory("GaussModel::tres(t,tres_mean[0.0],tres_sigma[0.05])")

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
    c.Print("anglesJpsiKstarMC.eps")       

##############################
### Proper Time Acceptance ###
##############################
timeaccfile = TFile(acceptancecurvefilename)
effhisto = timeaccfile.Get('effhisto')

effdatahist = RooDataHist('eff','eff',RooArgList(ws['t']),effhisto)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,eff)")

#######################
### Build the PDF's ###
#######################
# Signal mass
ws.factory("Gaussian::m_sig(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")

# Bkg mass
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

# Bkg time 
#ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")
#ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres,SingleSided)")
#ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")
ws.factory("RooDecay::t_bkg(t,t_bkg_tau[0.2,0.01,2.0],tres,SingleSided)")

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
bkg_ang = RooRealSumPdf('bkg_ang','bkg_ang',RooArgList(ws['Bkg_0000_basis']),RooArgList(c0000))
ws.put(bkg_ang)

#Only higher order in cospsi
#bkg_ang = RooRealSumPdf('bkg_ang','bkg_ang',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c1000,c2000,c3000))

#Add higher order in costheta
#bkg_ang = RooRealSumPdf('bkg_ang','bkg_ang',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c1000,c2000,c3000))

#Hiher roder in all angles
#bkg_ang = RooRealSumPdf('bkg_ang','bkg_ang',RooArgList(ws['Bkg_0000_basis'],ws['Bkg_0010_basis'],ws['Bkg_0011_basis'],ws['Bkg_001m1_basis'],ws['Bkg_1000_basis'],ws['Bkg_2000_basis'],ws['Bkg_3000_basis']),RooArgList(c0000,c0010,c0011,c001m1,c1000,c2000,c3000))

###########################
### Putting it together ###
###########################
# Mass only PDF
ws.factory("SUM::m_pdf(Nsig[1000,0,16000]*m_sig,Nbkg[5000,0,16000]*m_bkg)")
#ws['m_pdf'].fitTo(ws['data'],RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False))

# Signal PDF
ws.factory("PROD::sig_pdf( m_sig, newpdf)")
accsig = RooEffHistProd('accsig','accsig',ws['sig_pdf'],ws['effpdf'])
ws.put(accsig)

# Bkg PDF
#Whoooo, this will become a nasty trick!
#RooEffHistProd can't deal with RooRealSumPdf apparantly. So instead of
#eff(t)[t*m*ang] = RooEffHistProd
#Do
#eff(t)[t*m] = RooEffHistProd and then multiply this object by bkg_ang
ws.factory("PROD::bkg_pdf_tm( m_bkg, t_bkg)")

accbkg_tm = RooEffHistProd('accbkg_tm','accbkg_tm',ws['bkg_pdf_tm'],ws['effpdf'])
ws.put(accbkg_tm)
ws.factory("PROD::accbkg(accbkg_tm, bkg_ang)")

ws.factory("SUM::accpdf_ext(Nsig*accsig,Nbkg*accbkg)")
ws.factory("SUM::accpdf(f_sig*accsig,accbkg)")

print 'GOING TO BUILD THE ACCEPTANCE CORRECTED FULL PDF!'

timeangcorrpdf = buildEff_x_PDF(ws,'angcorr',ws['accpdf_ext'],[ ( m.basis() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] )
ws.put(timeangcorrpdf)

ws.factory("Gaussian::p0(p0var[0.,0.5],p0mean[0.34],p0sigma[0.012])")
ws.factory("Gaussian::p1(p1var[-4.0,4.0],p1mean[1.01],p1sigma[0.12])")
ws.factory("expr:wtag_syst('p0var+p1var*(tagomega-etamean)',tagomega,etamean[0.339],p0var,p1var)")

###########################
### Tagging systematics ###
###########################

wtag = ws.function('wtag')
wtag_syst = ws.function('wtag_syst')

#timecorrpdf
customizer = RooCustomizer(ws['accpdf_ext'],'inc_tag_syst')
customizer.replaceArg( wtag, wtag_syst )
timecorrtagpdf = customizer.build()
getattr(ws,'import')(timecorrtagpdf,RooFit.RecycleConflictNodes())

#timeangcorrpdf
customizer = RooCustomizer(timeangcorrpdf,'inc_tag_syst')
customizer.replaceArg( wtag, wtag_syst )
timeangcorrtagpdf = customizer.build()
getattr(ws,'import')(timeangcorrtagpdf,RooFit.RecycleConflictNodes())

wsfile = TFile('AccWS.root','RECREATE')
ws.Write()
wsfile.Close()

