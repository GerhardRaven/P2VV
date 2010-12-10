from ROOT import *
import RooFitDecorators
gSystem.Load("libp2vv")

feelTheNeedForSpeed = True
if feelTheNeedForSpeed:
    ### experimental fast(er) toy generator...
    RooMultiCatGenerator.registerSampler( RooNumGenFactory.instance() )
    RooNumGenConfig.defaultConfig().methodND(False,True).setLabel( "RooMultiCatGenerator" )
    RooNumGenConfig.defaultConfig().methodND(False,True).Print()
    RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))

class apybasis : # TODO: can also implement this by returning a 'bound' function instead...
    def __init__(self,w,*args) :
       self.w = w
       (cpsi,cheta,phi) = args if len(args)==3 else args[0]

       def _f(x) :
            if type(x) is str : return w[x]
            if not self.w.function(x.GetName()) : x = w.put(x)
            return x
       self.cpsi   = _f(cpsi)
       self.ctheta = _f(cheta)
       self.phi    = _f(phi)

    def build(self,label,i,j,k,l,c) :
        name = "%s_%d_%d_%d_%d" % (label,i,j,k,l)
        name.replace("-","m")
        b = self.w.function(name) # workaround a bug in ROOT 5.26 -- if name not present, w.obj(name) will SEGV...
        if not b : 
            self.w.put( RooP2VVAngleBasis(name,name,self.cpsi,self.ctheta,self.phi,i,j,k,l,c) )
            b = self.w[name]
        return b


def buildAngularBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    def _ba(name,comp) :
        n = name + '_basis'
        s = RooArgSet()
        for c in comp : s.add( ab.build(name,c[0],c[1],c[2],c[3],c[4]) )
        return ws.put(RooAddition_( n, n, s ) )

    return ( _ba("AzAz",       [ ( 0,0,0, 0, 2.), ( 0,0,2,0,  sqrt(1./ 5.)), ( 0,0,2,2, -sqrt( 3./5.))
                               , ( 2,0,0, 0, 4.), ( 2,0,2,0,  sqrt(4./ 5.)), ( 2,0,2,2, -sqrt(12./5.)) ] )
           , _ba("AparApar",   [ ( 2,2,0, 0, 1.), ( 2,2,2,0,  sqrt(1./20.)), ( 2,2,2,2,  sqrt( 3./20.)) ] )
           , _ba("AperpAperp", [ ( 2,2,0, 0, 1.), ( 2,2,2,0, -sqrt(1./ 5.)) ] )
           , _ba("AparAperp",  [ ( 2,2,2,-1,  sqrt(3./5.)) ] )
           , _ba("AzAperp",    [ ( 2,1,2, 1,  sqrt(6./5.)) ] )
           , _ba("AzApar",     [ ( 2,1,2,-2, -sqrt(6./5.)) ] )
           )

def buildJpsiphi(ws, name) :
    basis = buildAngularBasis(ws, apybasis(ws,ws.set('transversityangles')))

    # define the relevant combinations of strong amplitudes
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")

    ws.put(RooFormulaVar("qtag_","@0*(1-2*@1)",RooArgList( ws['tagdecision'],ws['tagomega']) ) )

    ws.factory("expr::N('1-@0*@1',{qtag_,C})")
    ws.factory("Minus[-1]")
    ws.factory("$Alias(Addition_,sum_)")

    # TODO: move this bit into a derivative of RooBDecay, and do tagdecision explicitly
    #       -- at that point, FOAM will do the angles, and we avoid the max search
    # generate untagged, then do tag
    # for this we need to pass qtag into the pdf
    # this can be done generically if we pass 8 instead of 4 factors
    # into RooBDecay -- 4 for tag = +1 and 4 for tag = -1 (tag = 0 would take the sum)
    # then generate time according to the sum over tag
    # and do the tag conditionally given the time...
    # (i.e. we generate not the time distributions of tagged events,
    # but first the one for untagged events, and then we generate the
    # asymmetry, which is quick...)
    # Next, how to do Jpsi K* if we do tag,rec instead of (un)mix...?
    # in that case, we have three asymmetries (of which only one, mix/unmix,
    # is non-zero)
    # Note that we can use a RooCustomizer to automate the replacement of
    # fjpsiphi_sinh and fjpsiphi_sin, but the qtag in N is more tricky...
    ws.factory("sum_::fjpsiphi_cosh({ prod(N,NAzAz,                    AzAz_basis)"
                                   ", prod(N,NAparApar,                AparApar_basis)"
                                   ", prod(N,NAperpAperp,              AperpAperp_basis)"
                                   ", prod(N,ImAparAperp,            C,AparAperp_basis)"
                                   ", prod(N,ImAzAperp,              C,AzAperp_basis)"
                                   ", prod(N,ReAzApar,                 AzApar_basis)"
                                   "})")
    ws.factory("sum_::fjpsiphi_cos ({ prod(N,NAzAz,            qtag_,C,AzAz_basis)"
                                   ", prod(N,NAparApar,        qtag_,C,AparApar_basis)"
                                   ", prod(N,NAperpAperp,      qtag_,C,AperpAperp_basis)"
                                   ", prod(N,ImAparAperp,      qtag_,  AparAperp_basis)"
                                   ", prod(N,ImAzAperp,        qtag_,  AzAperp_basis)"
                                   ", prod(N,ReAzApar,         qtag_,C,AzApar_basis)"
                                   "})")
    ws.factory("sum_::fjpsiphi_sinh({ prod(N,NAzAz,      Minus,      D,AzAz_basis)"
                                   ", prod(N,NAparApar,  Minus,      D,AparApar_basis)"
                                   ", prod(N,NAperpAperp,            D,AperpAperp_basis)"
                                   ", prod(N,ReAparAperp,            S,AparAperp_basis)"
                                   ", prod(N,ReAzAperp,              S,AzAperp_basis)"
                                   ", prod(N,ReAzApar,   Minus,      D,AzApar_basis)"
                                   "})")
    ws.factory("sum_::fjpsiphi_sin ({ prod(N,NAzAz,      Minus,qtag_,S,AzAz_basis)"
                                   ", prod(N,NAparApar,  Minus,qtag_,S,AparApar_basis)"
                                   ", prod(N,NAperpAperp,      qtag_,S,AperpAperp_basis)"
                                   ", prod(N,ReAparAperp,Minus,qtag_,D,AparAperp_basis)"
                                   ", prod(N,ReAzAperp,  Minus,qtag_,D,AzAperp_basis)"
                                   ", prod(N,ReAzApar,   Minus,qtag_,S,AzApar_basis)"
                                   "})")
    ws.factory("BDecay::%s(t,tau,dG,fjpsiphi_cosh,fjpsiphi_sinh,fjpsiphi_cos,fjpsiphi_sin,dm,tres_sig,SingleSided)" % name)
    return ws.pdf(name)


def buildJpsikstar(ws, name) :
    buildAngularBasis(ws, apybasis(ws,ws.set('transversityangles')))
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")

    ws.put(RooFormulaVar("qmix","@0*@1*(1-2*@2)",RooArgList( ws["qrec"],ws["tagdecision"],ws["wmistag"] ) ) )

    ws.factory("{Minus[-1],Zero[0],One[1]}")
    # in J/psi K*, things factorize because deltaGamma = 0 -- use this !
    ws.factory("PROD::%s( RealSumPdf( { AzAz_basis , AparApar_basis, AperpAperp_basis, AparAperp_basis, AzAperp_basis, AzApar_basis} "
                       "            , { NAzAz,       NAparApar,      NAperpAperp,      ImAparAperp,     ImAzAperp,     ReAzApar } )"
                       ", BDecay(t,tau,Zero,One,Zero,qmix,Zero,dm,tres_sig,SingleSided))"%name)
    return ws.pdf(name)


## Looping over data in python is quite a bit slower than in C++
## Why?? How to improve this??
def computeMoments( data, moments ) :
    if not moments : return None
    allObs = moments[0].basis().getObservables(data)
    for event in data : 
        allObs.assignValueOnly( event )
        for m in moments : m.inc()

def buildMomentPDF(w,name,data,moments) :
    if not moments : return None
    computeMoments( data, moments ) 
    coef = RooArgList()
    fact = RooArgList()
    for m in moments :
        C = 'C_%f' % m.coefficient()
        w.factory( '%s[%f]'%(C,m.coefficient() ) )
        coef.add( w[C] )
        fact.add( m.basis() )
    w.put( RooRealSumPdf(name,name,fact,coef) )
    return w.pdf(name)

def buildMoment_x_PDF(w,name,pdf,moments) :
   if not moments : return pdf
   # now we need to multiply all relevant components (i.e. all RooP2VVAngleBasis ones) 
   # of "pdf" with their efficiency corrected versions, multiply them with the right moment basis & coefficient...
   customizer = RooCustomizer(pdf,name)
   for c in pdf.getComponents() :
        if c is not RooP2VVAngleBasis : continue
        name = "%s_eff" % c.GetName()
        s = RooArgSet()
        [ s.add( c.createProduct( m.basis() , m.coefficient()) ) for m in moments ]
        rep = w.put( RooAddition_( name, name, s, True ) )  # hand over ownership & put in workspace...
        customizer.replaceArg( c, rep )
   return customizer.build(True)

def buildEffMomentsPDF(w,name,pdf,data,moments) :
    computeMoments(data,moments)
    return buildMoment_x_PDF(w,name,pdf,moments)

def buildMassPDFs(ws):
    #signal B mass pdf
    ws.factory("Gaussian::m_sig(m,m_sig_mean[5280,5200,5400],m_sig_sigma[10,3,30])")
    
    # background B mass pdf
    #ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.1,0.1])")
    # signal phi mass pdf
    #ws.factory("Voigtian::mphi_phisig(mKK,mphi_phi_mean[1919.455],mphi_phi_width[4.26],mphi_phi_sigma[1,0.1,10])")
    # background phi mass pdf (would like to use a nice threshold function here ... RooDstD0BG seems a bit complicated
    #ws.factory("DstD0BG::mphi_combbkg(mKK,mphi_bkg_m0[987.4],mphi_bkg_C[10,1,30],mphi_bkg_B[1,0.1,10],const_zero)")
    
    # signal J/psi mass pdf
    ws.factory("CBShape::mpsi_sig(mdau1,mpsi_sig_mean[3097,3085,3110],mpsi_sig_sigma[13.1,5,20],"
               "mpsi_sig_alpha[1.36,0.5,3],mpsi_sig_n[3])")
    
    # background J/psi mass pdf
    ws.factory("Exponential::mpsi_bkg(mdau1,mpsi_bkg_exp[-0.0002,-0.01,0.01])")

def build3GaussianResoModel(ws,name):
    ws.factory("GaussModel::tres_%s_g1(t,tres_%s_m1[0,-0.2,0.2],tres_%s_s1[1.1,0.3,2], 1, sigmat)" % (name,name,name))
    ws.factory("GaussModel::tres_%s_g2(t,tres_%s_m1,            tres_%s_s2[2.0,1.5,10],1, sigmat)" % (name,name,name))
    ws.factory("GaussModel::tres_%s_g3(t,tres_%s_m1,            tres_%s_s3[0.54,0.1,3.0])" % (name,name,name))
    ws.factory("AddModel::tres_%s({tres_%s_g3,tres_%s_g2,tres_%s_g1},{tres_%s_f3[0.001,0.00,0.01],tres_%s_f2[0.2,0.01,1]})" % (name,name,name,name,name,name))
    
def buildResoModels(ws):
    build3GaussianResoModel(ws,'sig')
    build3GaussianResoModel(ws,'nonpsi')

def declareObservables( ws ):
    from math import pi

    # transvercity angles
    ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f]}"%(-pi,pi))
    ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
    # helicity angles. we can also compute these from the transversity angles and add them as columns
    ws.factory("{ helcosthetaK[-1,1], helcosthetaL[-1,1], helphi[%f,%f]}"%(-pi,pi))
    ws.defineSet("helicityangles","helcosthetaK,helcosthetaL,helphi")
    # tag 
    ws.factory("tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]")
    ws.factory("tagomega[0,0.5]")
    # B, jpsi, phi mass
    ws.factory("m[5200,5450]")
    ws.factory("mdau1[%f,%f]"%(3097-60,3097+60))
    ws.factory("mdau2[0,2000]")
    # time , time-error
    ws.factory("{t[-2,14],sigmat[0.005,0.1]}")

    # define a set for the observables (we can also define this from the pdf and data.
    ws.defineSet("observables","m,t,sigmat,mdau1,mdau2,trcostheta,trcospsi,trphi,tagomega,tagdecision") 

    # the next is something we may need to switch on or off, depending on whether we use a pdf for sigmat
    ws.defineSet("conditionalobservables","sigmat")


def buildABkgPdf( ws, name, resname, psimasspdfname ):
    from itertools import repeat

    # build the dependencies if needed
    if not ws.function(resname)     : buildResoModels(ws)
    if not ws.function("m_%s"%name) : buildMassPDFs(ws)

    #background B mass pdf
    ws.factory("Exponential::m_%s(m,m_%s_exp[-0.001,-0.01,-0.0001])"% tuple(repeat(name,2)) )

    #background propertime
    ws.factory("RooDecay::t_sl_%s(t,0,%s,SingleSided)"%(name,resname))
    ws.factory("RooDecay::t_ml_%s(t,t_ml_%s_tau[0.21,0.1,0.5],%s,SingleSided)"%(name,name,resname))
    ws.factory("RooDecay::t_ll_%s(t,t_ll_%s_tau[1.92,1.,2.5],%s,SingleSided)"%(name,name,resname))
    ws.factory("SUM::t_%s(t_%s_fll[0.004,0,1]*t_ll_%s,t_%s_fml[0.02,0,1]*t_ml_%s,t_sl_%s)"% tuple(repeat(name,6)))
    
    #background angles: 
    ws.factory("Chebychev::trcospsipdf_%s(trcospsi,{tcp_0_%s[-0.13,-1,1],tcp_1_%s[0.23,-1,1],tcp_2_%s[-0.057,-1,1],tcp_3_%s[-0.0058,-1,1],tcp_4_%s[-0.0154,-1,1]})"% tuple(repeat(name, 6)) )
    ws.factory("Chebychev::trcosthetapdf_%s(trcostheta,{tct_0_%s[0.08,-1,1],tct_1_%s[-0.22,-1,1],tct_2_%s[-0.022,-1,1],tct_3_%s[0.21,-1,1],tct_4_%s[0.0125,-1,1]})"% tuple(repeat(name, 6)) )
    ws.factory("Chebychev::trphipdf_%s(trphi,{tp_0_%s[0.10,-1,1],tp_1_%s[0.328,-1,1],tp_2_%s[0.081,-1,1],tp_3_%s[0.316,-1,1],tp_4_%s[0.044,-1,1]})"% tuple(repeat(name, 6)) )

    #now multiply
    ws.factory("PROD::%s_pdf(trcosthetapdf_%s,trcospsipdf_%s,trphipdf_%s, t_%s, m_%s, %s )"%(tuple(repeat(name,6)) + (psimasspdfname,)))


def buildBkgPdf( ws ):
    
    # assume that the resolution models and psi mass models have been built     
    buildABkgPdf(ws,'nonpsi','tres_nonpsi','mpsi_bkg')
    buildABkgPdf(ws,'psi','tres_sig','mpsi_sig')
    # add them
    ws.factory("SUM::bkg_pdf(f_psi[0.5,0.01,1]*psi_pdf,nonpsi_pdf)")
    
    return ws.pdf('bkg_pdf')

def definePolarAngularAmplitudes(ws):

    ##choice: either fit for the Re&Im of the 3 amplitudes (and then
    ##        constrain one phase and the sum of magnitudes)
    ##        or fit in terms of angles and relative magnitudes
    ##         Note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009
    # ws.factory("{rz[0.556],rpar[0.211],rperp[0.233]}")
    ws.factory("{rz[0.463,0.1,0.9],rpar[0.211],rperp[0.347,0.1,0.9]}")
    ws.factory("{deltaz[0],deltapar[-2.93],deltaperp[2.91]}")
    ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
    ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
    ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
    ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
    ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
    ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")
    
    # this doesn't belong here
    ws.factory("{gamma[0.68,0.4,0.9],dm[17.7],dG[0.05,-0.3,0.3]}")
    ws.factory("expr::tau('1/@0',{gamma})")
    ##choice: either fit for the three degrees of freedom independently
    ##        i.e. make S,D,C independent parameters
    ##ws.factory("{S[0.717,-1,1],D[0.696,-1,1],C[0,-1,1]}")
    ##        or write S,D,C in terms of phi_s
    #ws.factory("{expr::S('-sin(phis)',{phis[0.]}),expr::D('cos(phis)',{phis}),C[0]}")
    ws.factory("{expr::S('-sin(phis)',{phis[0.8]}),expr::D('cos(phis)',{phis}),C[0]}")
    ##        The no-CP violation case:
    ##ws.factory("{S[0],C[0],D[1]}")

def buildFullJpsiPhiPdf( ws ):
    # naming convention for pdfnames:
    # pdfname = observable_component
    # naming convention for pdf params:
    # parname = pdfname_parname

    declareObservables(ws)
    buildResoModels(ws)
    buildMassPDFs(ws)
    buildBkgPdf(ws)

    # this builds the angular pdf for J/psi phi
    definePolarAngularAmplitudes(ws)
    buildJpsiphi(ws,'jpsiphipdf')
    
    # still need to combine with the mass
    ws.factory("PROD::sig_pdf(jpsiphipdf,m_sig)")
    
    # make the final pdf
    ws.factory("SUM::pdf_ext(Nsig[1200,0,1000]*sig_pdf,Nbkg[81200,0,1000000]*bkg_pdf)")
    ws.factory("SUM::pdf(f_sig[0.01,0.,1.0]*sig_pdf,bkg_pdf)")


def readParameters( ws, filename, pdfname='pdf_ext'):
    pdf = ws[pdfname]
    pdf.getParameters(ws.set('observables')).readFromFile( filename )
    data = ws['data']
    if not data:
        print 'warning: no dataset in workspace. cannot initialize yields'
    else:
        from math import sqrt
        fsig = ws['Nsig'].getVal() /  ws['Nbkg'].getVal()
        N = data.numEntries()
        ws['Nsig'].setVal( N * fsig )
        ws['Nsig'].setError( sqrt( N * fsig ) )
        ws['Nbkg'].setVal( N * (1-fsig) )
        ws['Nbkg'].setError( sqrt( N * (1-fsig) ) )


### TODO: make a python version of buildEfficiencyPDF....
### TODO: sample efficiency from 3D angle histogram -- i.e. make a Fourier transform...
