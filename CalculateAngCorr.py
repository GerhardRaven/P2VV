####################################################################################################
####################################################################################################
# This script builds a root file with the worskpace containing all information for the Tagged fit of Note6.
# Daan van Eijk, 29-09-2011

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

#2010
#mcfilename = '/data/bfys/dveijk/MC/ReducedMCNTuple.root'
#2011
#mcfilename = '/data/bfys/dveijk/MC/2011/MC2011.root'
#For fitter comparisons, from Greig
mcfilename = '/data/bfys/dveijk/MC/2011/MC2011_ForAngAccComp.root'

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[-2,20.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0],tagomega[0.,0.5]}"%(-pi,pi))
ws.factory("{biased[Biased=+1,NotBiased=0],unbiased[Unbiased=+1,NotUnbiased=0],fullybiased[FullyBiased =+1,NotFullyBiased=0]}")

ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
angles = ws.set('transversityangles')

ab = abasis(ws,angles)

#ws.factory("{rz2[0.601],rperp2[0.16],rs2[0.]}")
ws.factory("{rz2[0.60],rperp2[0.16],rs2[0.]}")
ws.factory("RooFormulaVar::rpar2('1-@0-@1-@2',{rz2,rperp2,rs2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")
ws.factory("RooFormulaVar::rs('sqrt(@0)',{rs2})")

ws.factory("{deltaz[0.],deltapar[2.5],deltaperp[-0.17],deltas[0.]}")

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
#ws.factory("{#Gamma[0.68]}")
ws.factory("{#Gamma[0.681]}")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")

#ws.factory("{t_sig_dG[0.06852]}")
ws.factory("{t_sig_dG[0.060]}")

ws.factory("{t_sig_dm[17.8]}")

#Parametrize in terms of phis:
ws.factory('{phis[-0.04]}')
ws.factory("{expr::S('-1*sin(@0)',{phis}),expr::D('cos(@0)',{phis}),C[0]}")

###############################
### Experimental parameters ###
###############################
# For determination of efficiency from MC, put RooTruthModel, later replace this by realistic Resolution model. Does that imply building the JpsiPhi pdf again? Yes, you should, but you never did it yet. Apparantly it doesn't matter for the result, but fix it!!!

ws.factory("RooTruthModel::tres_MC(t)")

#Build the MC PDF
MCpdf = buildJpsiphiSWave(ws,'MCpdf', True,'tres_MC')

ws.defineSet("MCobservables","t,trcospsi,trcostheta,trphi,tagdecision,tagomega")

MCdatafile = TFile(mcfilename)
NTupletree = MCdatafile.Get('MyTree')

#For ang acc corr cross check
tmin = 0.3
ws['t'].setRange(tmin,ws['t'].getMax())

MCdata = RooDataSet('MCdata','MCdata',NTupletree,ws.set('MCobservables'),'t==t && trcospsi==trcospsi && trcostheta == trcostheta && trphi==trphi && tagomega == tagomega')
ws.put(MCdata)
print 'Number of MC events', MCdata.numEntries()

allObs = MCpdf.getObservables( MCdata.get() )

angles = MCpdf.getObservables( MCdata.get() )
angles.remove(ws.var('t'))
angles.remove(ws.cat('tagdecision'))
angles.remove(ws.var('tagomega'))

print 'angles: ', [ i.GetName() for i in angles ]

print 'MCobservables:', [ i.GetName() for i in allObs ]

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

computeMoments(MCdata,MCpdf,tenmom)
xi_m = dict( [ (m.basis().GetName(),m.coefficient()) for m in tenmom ] )
print 'Direct Moments xi_m =', xi_m

def norm_xi( d ) :
    n =  (1.)/(8.*sqrt(pi))
    for i in d.iterkeys() : d[i] = d[i]*n

norm_xi(xi_m)
print 'Normalized direct moments xi_m =', xi_m
assert False
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
for (i,l) in product(range(10),range(3)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    #moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, MCpdf_marginal, allObs ) for m in range(-l,l+1) ]
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, MCpdf, angles ) for m in range(-l,l+1) ]

# loop over all data, determine moments
computeMoments(MCdata,MCpdf,moments)

# compute the 'canonical' moments given the Fourier series
c = dict()
for m in moments :
    c[ ( m.basis().i(),m.basis().l(),m.basis().m() ) ] = m.coefficient()

c_000 = c[(0,0,0)]

xi_c = compute_moments( c )
# normalize moments and compare
norm_xi(xi_c)

for name in xi_c.iterkeys() :
    print '%s : direct moment: %s ;  from Fourier series: %s ; ratio = %s ' % \
          ( name, xi_m[name], xi_c[name], xi_m[name]/xi_c[name])

# build PDF using the Fourier series efficiency...
# TODO: normalize relative to 0000 so that c_000 = 1
print 'using the following terms in Fourier expansion: '

for n,coeff in [ ( m.basis().GetName() , m.coefficient() ) for m in moments if m.significance()>signif] :
    print '%s : %s ' % (n,coeff)

MCpdf_eff = buildEff_x_PDF(ws,'fourier_eff',MCpdf,[ ( m.basis() , m.coefficient() ) for m in moments if m.significance()>signif] )

print 'effTimesPdfName: ', MCpdf_eff.GetName() 
ws.put( MCpdf_eff )

effpdfs = RooArgList()
effcoeffs = RooArgList()
list = [ [ m.basis() , c[(m.basis().i(),m.basis().l(),m.basis().m())] ] for m in moments if m.significance()>signif]

#ws.factory("const_0[%s]"%(list[0][1]))
#effpdfs.add(list[0][0])
#effcoeffs.add(ws['const_0'])

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

## testc = TCanvas("testc","Use MC data to determine efficiency",900,300)
## testc.Divide(3,2);

## testc.cd(1)
## trcospsiframe = ws['trcospsi'].frame()
## effshape.plotOn(trcospsiframe)
## trcospsiframe.Draw()

## testc.cd(2)
## trcosthetaframe = ws['trcostheta'].frame()
## effshape.plotOn(trcosthetaframe)
## trcosthetaframe.Draw()

## testc.cd(3)
## trphiframe = ws['trphi'].frame()
## effshape.plotOn(trphiframe)
## trphiframe.Draw()

## testc.cd(4)
## trcospsiframe = ws['trcospsi'].frame()
## effshape.plotOn(trcospsiframe,RooFit.Project(RooArgSet(ws['trcostheta'],ws['trphi'])))
## trcospsiframe.Draw()

## testc.cd(5)
## trcosthetaframe = ws['trcostheta'].frame()
## effshape.plotOn(trcosthetaframe,RooFit.Project(RooArgSet(ws['trcospsi'],ws['trphi'])))
## trcosthetaframe.Draw()

## testc.cd(6)
## trphiframe = ws['trphi'].frame()
## effshape.plotOn(trphiframe,RooFit.Project(RooArgSet(ws['trcostheta'],ws['trcospsi'])))
## trphiframe.Draw()
