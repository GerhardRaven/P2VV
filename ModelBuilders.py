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

class abasis : # TODO: can also implement this by returning a 'bound' function instead...
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
       print 'using %s,%s,%s' % (self.cpsi,self.ctheta,self.phi)

    def angles(self) : 
        return RooArgList(self.cpsi,self.ctheta,self.phi)
    def build(self,label,i,j,k,l,c) :
        name = "%s_%d_%d_%d_%d" % (label,i,j,k,l)
        name.replace("-","m")
        b = self.w.function(name) # workaround a bug in ROOT 5.26 -- if name not present, w.obj(name) will SEGV...
        if not b : 
            b = self.w.put( RooP2VVAngleBasis(name,name,self.cpsi,self.ctheta,self.phi,i,j,k,l,c) )
        return b

def _buildAngularFunction(ws,ab,name,comp) :
        n = name + '_basis'
        s = RooArgSet()
        for c in comp : s.add( ab.build(name,c[0],c[1],c[2],c[3],c[4]) )
        return ws.put(RooAddition_( n, n, s ) )

def _buildTransversityBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    # transversity amplitudes in terms of transversity angles
    _ba = lambda name,comp : _buildAngularFunction(ws,ab,name,comp)

    return ( _ba("AzAz",       [ ( 0,0,0, 0, 2.), ( 0,0,2,0,  sqrt(1./ 5.)), ( 0,0,2,2, -sqrt( 3./ 5.))
                               , ( 2,0,0, 0, 4.), ( 2,0,2,0,  sqrt(4./ 5.)), ( 2,0,2,2, -sqrt(12./ 5.)) ] )
           , _ba("AparApar",   [ ( 2,2,0, 0, 1.), ( 2,2,2,0,  sqrt(1./20.)), ( 2,2,2,2,  sqrt( 3./20.)) ] )
           , _ba("AperpAperp", [ ( 2,2,0, 0, 1.), ( 2,2,2,0, -sqrt(1./ 5.)) ] )
           , _ba("AparAperp",  [ ( 2,2,2,-1,  sqrt(3./5.)) ] )
           , _ba("AzAperp",    [ ( 2,1,2, 1,  sqrt(6./5.)) ] )
           , _ba("AzApar",     [ ( 2,1,2,-2, -sqrt(6./5.)) ] )
           )

def _buildHelicityBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    # transversity amplitudes in terms of helicity angles
    _ba = lambda name,comp : _buildAngularFunction(ws,ab,name,comp)

    return ( _ba("AzAz",       [ ( 2,2,0, 0, 2.), (2,2,2,0, -sqrt(4./5.))
                               , ( 2,0,0, 0, 2.), (2,0,2,0, -sqrt(4./5.)) ] )
           , _ba("AparApar",   [ ( 2,2,0, 0, 1.), (2,2,2, 0, sqrt(1./20.)), ( 2,2,2,2,  -sqrt(3./20.)) ] )
           , _ba("AperpAperp", [ ( 2,2,0, 0, 2.), ( 2,2,2,0, sqrt(1./ 5.)), (2,2,2,2,sqrt(3./5.)) ] )
           , _ba("AparAperp",  [ ( 2,2,2,-2,  sqrt(3./5.)) ] )
           , _ba("AzAperp",    [ ( 2,1,2,-1, -sqrt(12./5.)) ] )
           , _ba("AzApar",     [ ( 2,1,2, 1,  sqrt(6./5.)) ] )
           )
## TODO: replace hardwired use of 'tagomega' with a passable rule, which by default returns 'tagomega'
##       this is needed when fitting for mistag rate in tagging categories...
def buildJpsiphi(ws, name, transversity ) : # TODO: add tagsplit
    afb = AngleFunctionBuilder(ws, 'transversity' if transversity else 'helicity' )
    basis = afb.basis()

    #(tagcat,tagpdf) = buildTagging(ws,'sigtag',tagsplit)

    # define the relevant combinations of strong amplitudes
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,    ImAz                      })")  # |A_z|^2
    ws.factory("expr::NAparApar  ('( @0 * @0 + @1 * @1 )',{ReApar,  ImApar                    })")  # |A_par|^2
    ws.factory("expr::NAperpAperp('( @0 * @0 + @1 * @1 )',{ReAperp, ImAperp                   })")  # |A_perp|^2
    ws.factory("expr::ReAparAperp('( @0 * @2 + @1 * @3 )',{ReApar,  ImApar,  ReAperp, ImAperp })")  # |A_par||A_perp| cos(delta_perp - delta_par)
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,    ImAz,    ReAperp, ImAperp })")  # |A_z||A_perp|   cos(delta_perp - delta_z)
    ws.factory("expr::ReAzApar   ('( @0 * @2 + @1 * @3 )',{ReAz,    ImAz,    ReApar,  ImApar  })")  # |A_z||A_par|    cos(delta_par  - delta_z)
    ws.factory("expr::ImAparAperp('( @0 * @3 - @1 * @2 )',{ReApar,  ImApar,  ReAperp, ImAperp })")  # |A_par|A_perp|  sin(delta_perp - delta_par)
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,    ImAz,    ReAperp, ImAperp })")  # |A_z||A_perp|   sin(delta_perp - delta_z)

    ws.factory("expr::qtag_('@0*(1-2*@1)',{tagdecision,wtag})")
    ws.factory("expr::N('1.0/(1.0+@0*@1)',{tagdecision,C})")
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
    # is normally non-zero)
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

    ws.factory("BDecay::%s(t,t_sig_tau,t_sig_dG,fjpsiphi_cosh,fjpsiphi_sinh,fjpsiphi_cos,fjpsiphi_sin,t_sig_dm,tres_sig,SingleSided)" % name)
    return ws.pdf(name)


def buildJpsikstar(ws, name) :
    _buildTransversityBasis(ws, abasis(ws,ws.set('transversityangles')))
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
## So we adapt the arguments, and then defer to the C++ _computeMoments
def computeMoments( data, moments ) :
    if not moments : return None
    vecmom = std.vector('IMoment*')()
    for m in moments : vecmom.push_back(m)
    return _computeMoments( data, vecmom )

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
    return w.put( RooRealSumPdf(name,name,fact,coef) )

def buildMoment_x_PDF(w,name,pdf,moments) :
   if not moments : return pdf
   # now we need to multiply all relevant components (i.e. all RooP2VVAngleBasis ones) 
   # of "pdf" with their efficiency corrected versions, multiply them with the right moment basis & coefficient...
   customizer = RooCustomizer(pdf,name)
   for c in pdf.getComponents() :
        if type(c) is not RooP2VVAngleBasis : continue
        name = "%s_eff" % c.GetName()
        s = RooArgSet()
        [ s.add( c.createProduct( m.basis() , m.coefficient()) ) for m in moments ]
        rep = w.put( RooAddition_( name, name, s, True ) )  # hand over ownership & put in workspace...
        customizer.replaceArg( c, rep )
   return customizer.build(True)

def buildEffMomentsPDF(w,name,pdf,data,moments) :
    computeMoments(data,moments)
    return buildMoment_x_PDF(w,name,pdf,moments)

class TimeResolutionBuilder :
    # TODO: build a PDF for sigmat (eg. RooHistPdf, RooThresholdPdf... or the sum of two gamma functions....)
    # gamma = RooRealVar("gamma","gamma",15,10,30)
    # beta = RooRealVar("beta","beta",2.21456e-03,0,0.1)
    # mu = RooRealVar("mu","mu",0,0.001,0.01)
    # pdf = RooGammaPdf("pdf","pdf",st,gamma,beta,mu)
    def __init__(self,ws, t, sigmat) :
        if type(t)      is str : t      = ws[t]
        if type(sigmat) is str : sigmat = ws[sigmat]
        for name in [ 'sig','nonpsi' ] :
            ws.factory("GaussModel::tres_%s_g1(%s,tres_%s_m[0,-0.2,0.2],tres_%s_s1[1.1,0.3,2 ], 1, %s)" % (name,t.GetName(),name,name,sigmat.GetName()))
            ws.factory("GaussModel::tres_%s_g2(%s,tres_%s_m,            tres_%s_s2[20.,1.5,30], 1, %s)" % (name,t.GetName(),name,name,sigmat.GetName()))
            ws.factory("GaussModel::tres_%s_g3(%s,tres_%s_m,            tres_s3[0.54,0.1,3.0])" % (name,t.GetName(),name) )
            ws.factory("AddModel::tres_%s({tres_%s_g3,tres_%s_g2,tres_%s_g1},{tres_%s_f3[0.001,0.00,0.05],tres_%s_f2[0.2,0.01,1]})" % (name,name,name,name,name,name))
            ws['tres_%s_f3'%name].setVal(0)
            ws['tres_%s_f3'%name].setConstant(True)
            ws['tres_s3'].setVal(1)
            ws['tres_s3'].setConstant(True)
        self._sigres = ws['tres_sig']
        self._nonpsi = ws['tres_nonpsi']
    def signal(self) : return self._sigres
    def nonpsi(self) : return self._nonpsi

    def signalSigmaT(self) : return self._sigmaT['signal']
    def psibkgSigmaT(self) : return self._sigmaT['psibkg']
    def nonpsibkgSigmaT(self) : return self._sigmaT['nonpsibkg']

class BkgTimePdfBuilder : #background propertime
    def __init__(self, ws, resbuilder, sigmatpdf ) :
        for name,resname in { 'nonpsibkg': resbuilder.nonpsi().GetName() , 'psibkg' : resbuilder.signal().GetName() }.iteritems() :
            ws.factory("PROD::t_%s_sl(Decay(t,0,                        %s,SingleSided)|sigmat,%s)"%(name,     resname,sigmatpdf[name].GetName()))
            ws.factory("PROD::t_%s_ml(Decay(t,t_%s_ml_tau[0.21,0.1,0.5],%s,SingleSided)|sigmat,%s)"%(name,name,resname,sigmatpdf[name].GetName()))
            ws.factory("PROD::t_%s_ll(Decay(t,t_%s_ll_tau[1.92,1.0,2.5],%s,SingleSided)|sigmat,%s)"%(name,name,resname,sigmatpdf[name].GetName()))
            ws.factory("SUM::t_%s(t_%s_fll[0.004,0,1]*t_%s_ll,t_%s_fml[0.02,0,1]*t_%s_ml,t_%s_sl)"% (name,name,name,name,name,name) )
        # fix fraction of  ll and lifetime for nonpsi to zero...
        ws['t_nonpsibkg_fll'].setVal(0)
        ws['t_nonpsibkg_fll'].setConstant(True)
        ws['t_nonpsibkg_ll_tau'].setVal(0)
        ws['t_nonpsibkg_ll_tau'].setConstant(True)
        ws['t_psibkg_fll'].setVal(0)
        ws['t_psibkg_fll'].setConstant(True)
        ws['t_psibkg_ll_tau'].setVal(0)
        ws['t_psibkg_ll_tau'].setConstant(True)
        self._nonpsi = ws['t_nonpsibkg']
        self._psi = ws['t_psibkg']

    def psibkgPdf(self) : return self._psi
    def nonpsibkgPdf(self) : return self._nonpsi

   
class AngleFunctionBuilder :
    def __init__( self, ws, basis ) :
       # TODO: move knowledge of which set to use into the basis builder itself...
       lookup = { 'transversity' : ( _buildTransversityBasis, 'transversityangles' ) 
                , 'helicity'     : ( _buildHelicityBasis,     'helicityangles' ) 
                }
       (b,v) = lookup[basis]
       self._basis = b( ws, abasis(ws, ws.set(v) ) )
    def basis(self) : return self._basis

class BkgAnglePdfBuilder :
    def __init__(self,ws,basis,data, opt) : 
        self._angles = basis.angles()
        self._dataw = {}
        self._pdf = {}
        for name in [ 'psibkg', 'nonpsibkg'] :
            dataw = RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"1>0",opt[name]['weight']) # need a dummy cut, as passing a (const char*)0 is kind of difficult...
            self._dataw[name] = dataw
            moments = []
            from itertools import product
            ranges = opt[name]['ranges']
            for (i,l,m) in product(ranges[0],ranges[1],ranges[2]) :
                  if abs(m)>l : continue
                  #  Warning: the Y_lm are orthonormal, but the P_i are orthogonal, with dot product 2/(2*i+1)
                  moments.append( Moment( basis.build(name+'_mom',i,0,l,m,1.), float(2*i+1)/2 ) )
            self._pdf[name] = buildMomentPDF( ws, 'angles_%s'%name, dataw, moments )
    def psibkgPdf(self) : return self._pdf['psibkg']
    def nonpsibkgPdf(self) : return self._pdf['nonpsibkg']    

    def makeplots(self) : 
        from ROOT import gStyle
        gStyle.SetOptStat(0)

        canvas = []
        for i in ['psibkg','nonpsibkg'] : 
            c = TCanvas('angle_%s'%i)
            canvas.append(c)
            c.Divide(3,2)
            for (f,v) in enumerate( self._angles ) :
                c.cd(1+f)
                frame = v.frame()
                self._dataw[i].plotOn(frame)
                self._pdf[i].plotOn(frame)
                frame.Draw()

                c.cd(4+f)
                others = RooArgList( self._angles )
                others.remove( v )
                hist = self._pdf[i].createHistogram( others.names() )
                self._pdf[i].fillHistogram( hist,others,1., RooArgSet(v))
                hist.Draw('COLZ')
                # create residuals in 2D
                #datahist = data.createHistogram( others.name() )
                #self._dataw[i].fillHistogram( datahist )
        return (canvas[0],canvas[1])

class MassPdfBuilder :
    ## TODO: investigate use of per-event mass error... or add a 2nd Gaussian to the b-mass
    ## TODO: integrate SPLOT functionality into the MassPdfBuilder...
    def __init__(self,ws,m,m_dau1,m_dau2) : # assume B-mass, J/psi mass, phi mass
        #### define J/psi mass observable & corresponding PDF
        self._mdau1 = m_dau1
        # signal J/psi mass pdf
        ws.factory("CBShape::mpsi_sig(%s,mpsi_sig_mean[3094,3090,3105],mpsi_sig_sigma[13.2,8,18],mpsi_sig_alpha[1.39,0.8,2],mpsi_sig_n[3])"%m_dau1.GetName())
        self._mdau1_sig = ws['mpsi_sig']
        # background J/psi mass pdf
        ws.factory("Exponential::mpsi_bkg(%s,mpsi_bkg_exp[-0.0005,-0.001,0.0])"%m_dau1.GetName())
        self._mdau1_bkg = ws['mpsi_bkg']
        # overall J/psi mass pdf
        ws.factory("SUM::mpsi(mpsi_fjpsi[0.5,0.2,0.8]*mpsi_sig,mpsi_bkg)")
        self._mdau1_pdf = ws['mpsi']


        # signal phi mass pdf
        self._mdau2 = m_dau2
        ws.factory("Voigtian::mphi_phisig(%s,mphi_phi_mean[1019.455],mphi_phi_width[4.26],mphi_phi_sigma[1.2,0.1,5])"%m_dau2.GetName())
        self._mdau2_sig = ws['mphi_phisig']
        # background phi mass pdf (would like to use a nice threshold function here ... RooDstD0BG seems a bit complicated
        # On top of that, someone cut at phi-20 MeV/c^2, so we don't see the KK threshold anyway...
        # so we might as well just do a 2nd order polynomial or so...
        # Even more so, by only using a +- 10 MeV window, we kill half the background ;-)
        # So until we actually include the phi mass, we take a linear function...
        #ws.factory("DstD0BG::mphi_combbkg(%s,mphi_bkg_m0[987.4],mphi_bkg_C[6,1,10],mphi_bkg_B[16,8,30],zero[0])"%m_dau2.GetName())
        ws.factory("Chebychev::mphi_combbkg(%s,{mphi_bkg_1[0.2,-1,1],mphi_bkg_2[-0.01,-0.1,0.1]})"%m_dau2.GetName())
        self._mdau2_bkg = ws['mphi_combbkg']
        ws.factory("SUM::m_phi(mphi_fphi[0.2,0.05,0.8]*mphi_phisig,mphi_combbkg)")
        self._mdau2_pdf = ws['m_phi']

        #########################
        #signal B mass pdf
        # TODO: can we include sigmam without introducing a Punzi problem?
        #       note that the background PDF would not include sigmam...
        #       but both mean and sigma are different for t>0.3 and t<0.3...
        #       we could just take a sigmam distribution with t>0.3 and |m-m_bs|<25 as signal...
        #       and t<0.3 and |m-m_bs|>30 as background...
        self._m = m
        #ws.factory("PROD::m_psisig(SUM(m_sig_f1[0.9,0.1,0.99]*Gaussian(%s,m_sig_mean[5380,5200,5400],m_sig_sigma[10,3,30]),Gaussian(%s,m_sig_mean,m_sig_sigma2[15,10,35])),mpsi_sig)"%(m.GetName(),m.GetName()))
        #ws.factory("PROD::m_psisig(Gaussian(%s,m_sig_mean[5366,5350,5380],m_sig_sigma[10,3,30]),mpsi_sig)"%(m.GetName()))
        #ws.factory("PROD::m_nonpsisig(Gaussian(%s,m_sig_mean,m_sig_sigma),mpsi_bkg)"%m.GetName())
        #ws.factory("SUM::m_sig(m_sigfpsi[1.0]*m_psisig,m_nonpsisig)")
        ws.factory("PROD::m_sig(SUM(m_sig_f1[0.2,0.1,0.99]*Gaussian(%s,m_sig_mean[5366,5360,5375],m_sig_sigma[5,3,10]),Gaussian(%s,m_sig_mean,m_sig_sigma2[8,6,14])),SUM(m_sig_fpsi[0.95,0.8,1]*mpsi_sig,mpsi_bkg))"%(m.GetName(),m.GetName()))
        ws['m_sig_fpsi'].setVal(1)
        ws['m_sig_fpsi'].setConstant(True)
        self._m_sig = ws['m_sig']
        if False :
            ws.factory("PROD::m_sig(Gaussian(m,m_sig_mean,expr('@0*@1',{m_sig_s[0.5,5],sigmam}))|sigmam)")
        
        #background B mass pdf
        ws.factory("PROD::m_psibkg(Exponential(%s,m_psibkg_exp[-0.0003,-0.001,-0.0001]),mpsi_sig)"% m.GetName() )
        self._m_psibkg = ws['m_psibkg']
        ws.factory("PROD::m_nonpsibkg(Exponential(%s,m_nonpsibkg_exp[-0.0006,-0.001,-0.0001]),mpsi_bkg)"% m.GetName() )
        self._m_nonpsibkg = ws['m_nonpsibkg']
        ws.factory("SUM::m_bkg(m_bkgfpsi[0.1,0.01,0.99]*m_psibkg,m_nonpsibkg)")
        self._m_bkg = ws['m_bkg']

        ### now it becomes a bit tricky. 
        ## we can have various combinations of {sig,bkg} x {sig,bkg} x {sig,bkg} here...
        ## for now, we just split the b bkg into psi vs non-psi, but at some point we
        ## need to split b sig into phi vs kk...
        #ws.factory("SUM::m_b(m_fb[0.1,0.01,0.99]*m_sig,m_bkg)")
        # we make an extended PDF so that we can make SPlots!
        # but we cannot write the background in terms of fN and (1-f)N and still do splots...
        #ws.factory("expr::N_psibkg('@0*@1',{N_fpsi[0.5,0,1],N_bkg[20000,0,30000]})")
        #ws.factory("expr::N_nonpsibkg('(1-@0)*@1',{N_fpsi,N_bkg})")
        #ws.factory("SUM::m_b(N_sig[1000,0,10000]*m_sig,N_psibkg*m_psibkg,N_nonpsibkg*m_nonpsibkg)")
        ws.factory("SUM::m_b(N_sig[1000,0,10000]*m_sig,N_psibkg[15000,0,30000]*m_psibkg,N_nonpsibkg[15000,0,30000]*m_nonpsibkg)")
        self._m_pdf = ws['m_b']

        self._yields = RooArgList(ws.argSet('N_sig,N_psibkg,N_nonpsibkg'))

        sigmid = 5367.4
        sigwid = 15
        self._m.setRange('sigRegion',sigmid-sigwid,sigmid+sigwid)
        self._m.setRange('leftSideband',self._m.getMin(),sigmid-sigwid)
        self._m.setRange('rightSideband',sigmid+sigwid,self._m.getMax())
        
    def sigPdf(self)    : return self._m_sig
    def bkgPdf(self)    : return self._m_bkg
    def psibkgPdf(self)    : return self._m_psibkg
    def nonpsibkgPdf(self)    : return self._m_nonpsibkg
    def Pdf(self)       : return self._m_pdf
    def Obs(self)       : return self._m
    def yields(self)    : return self._yields

    def sigDau1Pdf(self) : return self._mdau1_sig
    def bkgDau1Pdf(self) : return self._mdau1_bkg
    def dau1Pdf(self)    : return self._mdau1_pdf
    def dau1Obs(self)    : return self._mdau1

    def sigDau2Pdf(self) : return self._mdau2_sig
    def bkgDau2Pdf(self) : return self._mdau2_bkg
    def dau2Pdf(self)    : return self._mdau2_pdf
    def dau2Obs(self)    : return self._mdau2




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
    ws.factory("tagomega[0,0.500001]")


    # B, jpsi, phi mass
    #ws.factory("m[5200,5450]")
    #ws.factory("m[5250,5450]")
    ws.factory("m[%f,%f]"%(5366-50,5366+50))
    ws.factory("mdau1[%f,%f]"%(3097-60,3097+40))
    #ws.factory("mdau1[%f,%f]"%(3097-60,3097+50))
    #ws.factory("mdau2[%f,%f]"%(1019.455-20,1019.455+20)) 
    ws.factory("mdau2[%f,%f]"%(1019.455-10,1019.455+10)) ### Note: +- 10 Mev/c^2 keeps all of the phi signal, and kills 1/2 of the background!
    # time , time-error
    ws.factory("{t[-4,10],sigmat[0.005,0.1]}")

    # define a set for the observables (we can also define this from the pdf and data.
    ws.defineSet("observables","m,t,sigmat,mdau1,mdau2,trcostheta,trcospsi,trphi,tagomega,tagdecision") 

    # the next is something we may need to switch on or off, depending on whether we use a pdf for sigmat
    #ws.defineSet("conditionalobservables","sigmat,tagomega")
    ws.defineSet("conditionalobservables","sigmat")

##
## call with buildTagging(ws, name, [ 0.25, 0.35, 0.45 ] ) 
## to define tagging categories corresponding to the intervals [0,0.25),[0.25,0.35),[0.35,0.45),[0.45,0.5]
## note that the final interval is defined 'by construction' and that the pdf returned gives 
## the efficiency to be in the i-th category, with the last category having an efficiency 1-sum_i eff_i
class TagPdfBuilder :
    def __init__(ws,tagcatdef,tagomega='tagomega')  :
        self._ws = ws
        ws.factory("ThresholdCategory::tagcat('tagomega','untagged',0)" )
        self._tagcat = ws['tagcat']
        self._sig = ws.put(RooThresholdPdf('sig_tagcat','sig_tagcat',ws[tagomega]))
        self._bkg = ws.put(RooThresholdPdf('bkg_tagcat','bkg_tagcat',ws[tagomega]))

        for (i,upper) in enumerate( tagcatdef ) :
            cname = 'tagcat_%d' % i
            tagcat.addThreshold(upper,cname)
            sname = 'sig_%s' % cname
            bname = 'bkg_%s' % cname
            sig_eff = ws.put(RooRealVar( sname, sname , 0.2, 0., 1.))
            bkg_eff = ws.put(RooRealVar( bname, bname , 0.2, 0., 1.))
            self._sig.addThreshold(upper, sig_eff )
            self._bkg.addThreshold(upper, bkg_eff )

    def sigPdf(self) : return self._sig

    def bkgPdf(self) : return self._bkg

    def tagCat(self) : return self._tagCat

def buildTagging( ws, name, tagcatdef ) :
    # either make PDF conditional on tagomega distribution
    # and use a fittable version RooHistPdf for tagOmega,
    # different for sig and bkg
    #
    # or split tagomega distribution in discrete categories,
    # and multiply by efficiency for each category, seperate 
    # for signal and background... -- or just make the fit
    # extended, and treat each bin as Poisson bkg + Poisson sig
    ws.factory("wtag[0,0,0.5]")
    tagcat = ws.factory("ThresholdCategory::%s('tagomega','untagged',0)" % name+"cat")
    pdf = ws.put( RooThresholdPdf(name+'effpdf',name+'effpdf',ws['tagomega']) )
    for (i,upper) in enumerate( tagcatdef ) :
        cname = '%s%d' % (name,i)
        ename = cname + '_eff'
        eff = ws.put(RooRealVar( ename, ename , 0.2, 0., 1.))
        tagcat.addThreshold(upper,cname)
        pdf.addThreshold(upper, eff )

    return (tagcat,pdf)
    # 
  



def buildABkgPdf( ws, name, resname, psimasspdfname ):
    from itertools import repeat

    # build the dependencies if needed
    if not ws.function(resname)     : buildResoModels(ws)
    if not ws.function('m_%s'%name) : buildMassPDFs(ws)

    #
    tpb = TimePdfBuilder(ws, name,resname)
    #background angles: 
    ws.factory("Chebychev::trcospsipdf_%s(trcospsi,{tcp_0_%s[-0.13,-1,1],tcp_1_%s[0.23,-1,1],tcp_2_%s[-0.057,-1,1],tcp_3_%s[-0.0058,-1,1],tcp_4_%s[-0.0154,-1,1]})"% tuple(repeat(name, 6)) )
    ws.factory("Chebychev::trcosthetapdf_%s(trcostheta,{tct_0_%s[0.08,-1,1],tct_1_%s[-0.22,-1,1],tct_2_%s[-0.022,-1,1],tct_3_%s[0.21,-1,1],tct_4_%s[0.0125,-1,1]})"% tuple(repeat(name, 6)) )
    ws.factory("Chebychev::trphipdf_%s(trphi,{tp_0_%s[0.10,-1,1],tp_1_%s[0.328,-1,1],tp_2_%s[0.081,-1,1],tp_3_%s[0.316,-1,1],tp_4_%s[0.044,-1,1]})"% tuple(repeat(name, 6)) )

    # apb = BkgAnglePdfBuilder(ws, 
    #now multiply
    ws.factory("PROD::%s_pdf(trcosthetapdf_%s,trcospsipdf_%s,trphipdf_%s, t_%s, m_%s, %s )"%(tuple(repeat(name,6)) + (psimasspdfname,)))
    return ws['%s_pdf'%name]


def buildBkgPdf( ws, name = 'bkg_pdf' ):
    
    # assume that the resolution models and psi mass models have been built     
    nonpsibkg = buildABkgPdf(ws,'nonpsi','tres_nonpsi','mpsi_bkg')
    psibkg    = buildABkgPdf(ws,'psi',   'tres_sig',   'mpsi_sig')
    # add them
    ws.factory("SUM::%s(f_psi[0.5,0.01,1]*%s,%s)"%(name,psibkg.GetName(),nonpsibkg.GetName()))
    return ws.pdf(name)

def definePolarAngularAmplitudes(ws):
    ##choice: either fit for the Re&Im of the 3 amplitudes (and then
    ##        constrain one phase and the sum of magnitudes)
    ##        or fit in terms of angles and relative magnitudes
    ##         Note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009
    # ws.factory("{rz[0.556],rpar[0.211],rperp[0.233]}")
    ws.factory("{rz[0.463,0.0,1.0],rpar[0.211],rperp[0.347,0.0,1.0]}")
    ws.factory("{deltaz[0],deltapar[-2.93],deltaperp[2.91]}")
    ws.factory("expr::ReAz   ('rz    * cos(deltaz)',   {rz,deltaz})")
    ws.factory("expr::ImAz   ('rz    * sin(deltaz)',   {rz,deltaz})")
    ws.factory("expr::ReApar ('rpar  * cos(deltapar)', {rpar,deltapar})")
    ws.factory("expr::ImApar ('rpar  * sin(deltapar)', {rpar,deltapar})")
    ws.factory("expr::ReAperp('rperp * cos(deltaperp)',{rperp,deltaperp})")
    ws.factory("expr::ImAperp('rperp * sin(deltaperp)',{rperp,deltaperp})")
    
def defineJPsiPhiPhysicsParams(ws):
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
    buildJpsiphi(ws,[0.3,0.35,0.45],'jpsiphipdf')
    
    # still need to combine with the mass
    ws.factory("PROD::sig_pdf(jpsiphipdf,m_sig)")
    
    # make the final pdf
    ws.factory("SUM::pdf_ext(N_sig[1200,0,10000]*sig_pdf,N_bkg[81200,0,1000000]*bkg_pdf)")
    ws.factory("SUM::pdf(f_sig[0.01,0.,1.0]*sig_pdf,bkg_pdf)")
    return ws['pdf_ext']


def readParameters( ws, filename, pdfname='pdf_ext'):
    pdf = ws[pdfname]
    pdf.getParameters(ws.set('observables')).readFromFile( filename )
    data = ws['data']
    if not data:
        print 'warning: no dataset in workspace. cannot initialize yields'
        return
    from math import sqrt
    fsig = ws['N_sig'].getVal() /  ws['N_bkg'].getVal()
    N = data.numEntries()
    ws['N_sig'].setVal( N * fsig )
    ws['N_sig'].setError( sqrt( N * fsig ) )
    ws['N_bkg'].setVal( N * (1-fsig) )
    ws['N_bkg'].setError( sqrt( N * (1-fsig) ) )

### backwards compatibility stub    
global _timeresbuilder
_timeresbuilder = None
def buildResoModels(ws):
    global _timeresbuilder
    if not _timeresbuilder : _timeresbuilder = TimeResolutionBuilder(ws)

##### backwards compatibility glue... making this a singleton...
global _MassPDFBuilder 
_MassPDFBuilder = None
def buildMassPDFs(ws) :
    global _MassPDFBuilder
    if not _MassPDFBuilder :
        m_B = ws['m']
        m_mumu = ws['mdau1']
        m_KK = ws['mdau2']
        _MassPDFBuilder = MassPdfBuilder(ws,m_B,m_mumu,m_KK) # todo: make this a J/psi phi builder, so that we can also have a J/psi K* one ;-]
    return _MassPDFBuilder
