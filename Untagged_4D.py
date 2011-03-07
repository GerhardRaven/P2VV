from ROOT import *
gSystem.Load("libp2vv")
from math import pi

import rootStyle
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

from ModelBuilders import *
from ModelBuilders import _buildAngularFunction
from RooFitDecorators import *

ncpu = RooCmdArg( RooFit.NumCPU(16) )

sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))
lw = RooCmdArg(RooFit.LineWidth(2))
ms = RooCmdArg(RooFit.MarkerSize(0.4))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))
stash = []

##############################################################################
##########################   There we go!!!!! ################################
##############################################################################

#####################################
### Angular acceptance correction ###
#####################################

signif = 1

blinded = False

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

ws.factory("{deltaz[0.],deltapar[2.5,%f,%f],deltaperp[-0.17]}"%(-2*pi,2*pi))

# create observables
ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")

##########################
### physics parameters ###
##########################
ws.factory("{#Gamma[0.68,0.3,0.9]}")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")

ws.factory("{t_sig_dG[0.060,-2,2]}")

ws.factory("{t_sig_dm[17.8]}")

ws.factory('{phis[-0.04]}')

ws.factory("{expr::S('-1*sin(phis)',{phis}),expr::D('cos(phis)',{phis}),C[0]}")
ws.factory("{expr::Sold('sin(phis)',{phis}),expr::Dold('cos(phis)',{phis}),Cold[0]}")

###############################
### Experimental parameters ###
###############################
# For determination of efficiency from MC, put RooTruthModel, later replace this by realistic Resolution model. Does that imply building the JpsiPhi pdf again?

#ws.factory("RooTruthModel::tres_sig(t)")
ws.factory("RooGaussModel::tres_sig(t,mu[0],sigma[0.05])")

# For determination of efficiency from MC, put wtag to 0, later integrate out (or set to 0.5 as I used to do before). Does that imply rebuilding the JpsiPhi pdf?
ws.factory("{wtag[0.0]}")

########################
### Building the PDF ###
########################

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

print 'angles: ', [ i.GetName() for i in angles ]

print 'observables:', [ i.GetName() for i in allObs ]

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
    c.Print("anglesJpsiPhiMC.eps")       

######################
### Build the PDFs ###
######################

#Build new workspace, since the PDF is going to look quite different from the signal PDF that we used to determine the angular acceptance correction

w = RooWorkspace("w")

###################
### Observables ###
###################
mode = 'Bs2Jpsiphi'

declareObservables(w,mode)

files = { 'Bu2JpsiK'  : 'Bu2JpsiKTuple.root'
          , 'Bs2Jpsiphi': '/data/bfys/dveijk/Data/Bs2JpsiPhiForTaggedFit.root'
          , 'Bd2JpsiKstar' : 'Bd2JpsiKstarTuple.root'
        }

tmin = -1
tmax = 14

w.var('t').setMin(tmin)
w.var('t').setMax(tmax)

mmin = 5200
mmax = 5550

w.var('m').setMin(mmin)
w.var('m').setMax(mmax)

obs = w.set('observables')

#angles = w.set('helicityangles')
angles = w.set('transversityangles')

obs.add( angles )

file = TFile(files[mode])
#data = RooDataSet('data','data',file.Get('dataset'),obs,'t==t && m==m')
data = RooDataSet('data','data',file.Get('MyTree'),obs,'t==t && m==m && sigmat==sigmat')

mb = MassPdfBuilder(w,w['m'],w['mdau1'],w['mdau2'] if mode != 'Bu2JpsiK' else None,mode)

####
massc = TCanvas()
massc.Divide(2,3)

### TODO: move these plots into the mass PDF builder...
if True :
    mdau1pdf = mb.dau1Pdf()
    #mdau1pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mdau1pdf.fitTo(data,ncpu)
    for i in mdau1pdf.getParameters(data): i.setConstant(True)
    comps = { 'mpsi_sig'       : ( sigcolor,dashed )
            , 'mpsi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( massc.cd(1), mb.dau1Obs(), data, mdau1pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True and mode != 'Bu2JpsiK' :
    massc.cd(2)
    mdau2pdf = mb.dau2Pdf()
    #mdau2pdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mdau2pdf.fitTo(data,ncpu)
    for i in mdau2pdf.getParameters(data) : i.setConstant(True)
    comps = { 'mphi_sig'       : ( sigcolor,dashed )
            , 'mphi_bkg'       : ( bkgcolor,dashed ) 
            }
    plot ( massc.cd(2), mb.dau2Obs(), data, mdau2pdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

if True :
    mpdf = mb.Pdf()
    #mpdf.getParameters(data).readFromFile('initialvalues_%s.txt'%mode,'ReadFromFile','Mass')
    mpdf.fitTo(data,ncpu)
    for i in mpdf.getParameters(data) : i.setConstant(True)
    ## TODO: at some point, we need to add the KK mass to the story...
    for i,obs in enumerate( [ mb.Obs(), mb.dau1Obs()])  : #  , mb.dau2Obs() ] ):
        # TODO: plot current in signal of other, sideband of other....
        comps = { 'm_sig'       : ( sigcolor,dashed )
                , 'm_psibkg'    : ( bkgcolor,dashed ) 
                , 'm_nonpsibkg' : ( nonpsicolor,dashed )
                }
        plot( massc.cd(3+i), obs, data, mpdf, comps, frameOpts = None, dataOpts = (ms,xes), pdfOpts = (lw,) )

massc.Print("InitialMassFits.eps")
#import sys
#sys.exit(0)

### Now that we have the masses fitted, let's make angle SPlots...
## TODO: move this into the mass PDF builder...
for i in mpdf.getParameters(data) : i.setConstant( i not in mb.yields() )
#mpdf.Print("T")
#mpdf.getParameters(data).Print("V")
splot = RooStats.SPlot("splotdata","splotdata",data,mpdf,mb.yields())
wdata = splot.GetSDataSet()

######################
### To show the sPlots for all the observables, per sample
######################

if False :
    c = TCanvas()
    c2 = TCanvas()
    stash.append(c)
    stash.append(c2)
    observables = RooArgSet(angles)
    observables.add(w['t'])
    observables.add(w['sigmat'])
    observables.add(w['wtag'])
    observables.add(w['mdau2'])
    c.Divide(len(observables),3)
    c2.Divide(1,3)
    for (i,sample) in enumerate([ j.GetName() for j in mb.yields() ]):
        dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","%s_sw"%sample) # need a dummy cut, as passing a (const char*)0 is kind of difficult...
        # should only make splots for observables not used in fit...
        for (j,observable) in enumerate( observables ) :
            f = observable.frame()
            dataw.plotOn(f, RooFit.DataError(RooAbsData.SumW2) )
            c.cd(i*len(observables)+j+1)
            f.Draw()
        c2.cd(1+i)
        f = w['t'].frame(-2.,2.,100)
        dataw.plotOn(f)
        f.Draw()


################
### If you want to describe background in terms of nth order RooP2VVAngleBasis, it calls the BkgAnglePdfBuilder
################
        
if False :
    # and now, build the angular distributions for psi and non-psi background, using the SWeights we
    # just got...
    ab = abasis(w,angles)
    x = BkgAnglePdfBuilder( w, ab, data
                          , { 'psibkg'    : { 'ranges' : (range(3),range(3),range(-3,4))  , 'weight' : 'N_psibkg_sw'    }
                            , 'nonpsibkg' : { 'ranges' : (range(3),range(3),range(-3,4)) , 'weight' : 'N_nonpsibkg_sw' } 
                            } )
    psibkg = x.psibkgPdf()
    nonpsibkg = x.nonpsibkgPdf()
    (c1,c2) = x.makeplots()

#################
### For now, use a simple flat background (also RooP2VVAngleBasis, such that it gets corrected with signla efficiency to compare with other fits....)
#################
ab = abasis(w,angles)
_ba = lambda name,comp : _buildAngularFunction(w,ab,name,comp)

_ba("psibkg_0000",  [ ( 0,0,0,0,1.) ] )
_ba("nonpsibkg_0000",  [ ( 0,0,0,0,1.) ] )
c0000 = RooRealVar('[c0000','c0000',1.)
psiangbkg = RooRealSumPdf('psiangbkg','psiangbkg',RooArgList(w['psibkg_0000_basis']),RooArgList(c0000))
nonpsiangbkg = RooRealSumPdf('nonpsiangbkg','nonpsiangbkg',RooArgList(w['nonpsibkg_0000_basis']),RooArgList(c0000))
w.put(psiangbkg)
w.put(nonpsiangbkg)

################
### next, we pick up the sigmat distributions for our three components...
################

sigmat = w['sigmat']
sigmat.setBins(40)
p= sigmat.frame()
stdata = {}
stpdf = {}
c = TCanvas()
for (f,sample) in enumerate([ 'sig','psibkg','nonpsibkg' ]):
      dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0","N_%s_sw"%sample) 
      stdata[sample] = RooDataHist("sigmat_%s_data"%sample,"hist Err Per Ev",RooArgSet(sigmat),dataw)
      stpdf[sample] = w.put(RooHistPdf("sigmat_%s"%sample,"sigmat_%s"%sample,RooArgSet(sigmat),stdata[sample])) # ,2)
      print stpdf[sample]
      dataw.plotOn(p) # python has trouble finding the right plotOn...
      stpdf[sample].plotOn(p,RooFit.LineColor( [kRed,kBlue,kBlack][f]))
p.Draw()

#################
### Build time resolution model
#################

y = TimeResolutionBuilder(w, w['t'],w['sigmat'])

#################
### Build time PDFs for the background species
#################

z = BkgTimePdfBuilder(w, y, stpdf)

w.factory('{phis[-0.04]}')
w.factory("{expr::S('-1*sin(phis)',{phis}),expr::D('cos(phis)',{phis}),C[0]}")
w.factory("{#Gamma[0.68,0.4,0.9]}")
w.factory("{dG[0.05,-0.3,0.3]}")

if blinded:
   #Building blinded parameters
   w.factory("RooUnblindUniform::#Gamma_unblind('BsCalvin',0.4,#Gamma)")
   w.factory("expr::t_sig_tau('1/@0',{#Gamma_unblind})")

   w.factory("RooUnblindUniform::t_sig_dG('BsHobbes',0.2,dG)")

w.factory("expr::t_sig_tau('1/@0',{#Gamma})")
w.factory("expr::t_sig_dG('@0',dG)")

w.factory("{t_sig_dm[17.7]}")

definePolarAngularAmplitudes(w)
print 'about to build J/psi phi signal'
sigjpsiphi = buildJpsiphi(w,'jpsiphisignal',False)
print 'about to createProjection of J/psi phi signal over tagdecision'
#t_sig = sigjpsiphi # w.put(sigjpsiphi.createProjection( w.argSet('tagdecision') ))
t_sig = w.put(sigjpsiphi.createProjection( w.argSet('tagdecision') ))

######################
### Multiplying... ###
######################

m_sig = w.pdf('m_sig')

sig = RooProdPdf('sig','sig', m_sig,t_sig)
w.put(sig)

#w.factory("PROD:psibkg(    m_psibkg,    t_psibkg,    %s )" % x.psibkgPdf().GetName())
#w.factory("PROD:nonpsibkg( m_nonpsibkg, t_nonpsibkg, %s )" % x.nonpsibkgPdf().GetName())

w.factory("PROD:psibkg(    m_psibkg,    t_psibkg,    psiangbkg )" )
w.factory("PROD:nonpsibkg( m_nonpsibkg, t_nonpsibkg, nonpsiangbkg )")

w.factory("SUM::pdf(f_sig[0.1,0,0.4]*sig, SUM(f_psi[0.5,0.1,0.9]*psibkg,nonpsibkg))")

uncorrpdf = w['pdf']
w['m_sig_sigma'].setConstant(False)
#w['m_sig_f1'].setConstant(False)
#w['m_sig_sigma2'].setConstant(False)

#Set back phis to zero for untagged fit
w['phis'].setVal(0)

#############
### Apply angular acceptance correction
##############
pdf = buildEff_x_PDF(w,'pdf',uncorrpdf,[ ( m.basis() , m.coefficient()/c_000 ) for m in moments if m.significance()>signif] )
w.put(pdf)

result = pdf.fitTo(data,ncpu,RooFit.Save(True),RooFit.Minos(false))
pdf.getParameters(data).writeToFile('fitresult_%s.txt'%mode)

paramlist = pdf.getParameters(data)
writeFitParamsLatex(paramlist)

dict = writeCorrMatrixLatex(result)

t = w['t']
m = w['m']
t.setRange('largeTime',0.3,t.getMax())

c = TCanvas()
c.Divide(5,2)

#===========================================================================================================
bkgcomps = { 'psibkg'    : ( bkgcolor,dashed ) 
           , 'nonpsibkg' : ( nonpsicolor,dashed )
           }
comps = { 'sig'       : ( sigcolor,dashed )
        , 'psibkg'    : ( bkgcolor,dashed ) 
        , 'nonpsibkg' : ( nonpsicolor,dashed )
        }

#===========================================================================================================
_c1 = plot( c.cd(1),w['mdau1'],data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), )
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c2 = plot( c.cd(2),w['mdau1'],data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m(#mu#mu) (t>0.3)') )
          , dataOpts = ( ms, xes, RooFit.CutRange('largeTime') )
          , pdfOpts = ( lw, RooFit.ProjectionRange('largeTime') ) 
          )
#===========================================================================================================
_c3 = plot( c.cd(3),m,data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m') )
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c4 = plot( c.cd(4),m,data,pdf,comps
          , frameOpts = ( RooFit.Bins(50), RooFit.Title('m (t>0.3)') )
          , dataOpts = ( ms, xes, RooFit.CutRange('largeTime') )
          , pdfOpts = ( lw, RooFit.ProjectionRange('largeTime') ) 
          )
#===========================================================================================================
_c5 = plot( c.cd(5),sigmat,data,pdf,comps
          , dataOpts = ( ms, xes )
          , pdfOpts = ( lw, ) 
          )
#===========================================================================================================
_c6 = plot( c.cd(6), t, data, pdf, comps
          , frameOpts = ( RooFit.Range(-0.4,0.4), RooFit.Bins(100), RooFit.Title('proper time, full mass range') )
          , dataOpts = ( ms,xes )
          , pdfOpts = ( lw, )
          )
#===========================================================================================================
_c7 = plot( c.cd(7), t, data, pdf, comps
          , frameOpts = ( RooFit.Title('proper time, signal region'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('sigRegion') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('sigRegion') )
          , logy = True
          )
#===========================================================================================================
_c8 = plot( c.cd(8), t, data, pdf, bkgcomps
          , frameOpts = ( RooFit.Title('proper time, lower sideband'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('leftSideband') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('leftSideband') )
          , logy = True
          )
#===========================================================================================================
_c9 = plot( c.cd(9), t, data, pdf, bkgcomps
          , frameOpts = ( RooFit.Title('proper time, upper sideband'), )
          , dataOpts = ( ms,xes,RooFit.CutRange('rightSideband') )
          , pdfOpts = ( lw,RooFit.ProjectionRange('rightSideband') )
          , logy = True
          )


c.Print("FinalPlots.eps")
