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

    return ( _ba("AzAz",       [ ( 0,0,0, 0, 2.), (0,0,2,0, -sqrt(4. /5.))
                               , ( 2,0,0, 0, 4.), (2,0,2,0, -sqrt(16./5.)) ] )
           , _ba("AparApar",   [ ( 2,2,0, 0, 1.), (2,2,2, 0, sqrt(1. /20.)), ( 2,2,2,2,  -sqrt(3./20.)) ] )
           , _ba("AperpAperp", [ ( 2,2,0, 0, 1.), ( 2,2,2,0, sqrt(1. /20.)), (2,2,2,2,sqrt(3./20.)) ] )
           , _ba("AparAperp",  [ ( 2,2,2,-2,  sqrt(3./5.)) ] )
           , _ba("AzAperp",    [ ( 2,1,2,-1, -sqrt(6./5.)) ] )
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

def buildEff_x_PDF(w,name,pdf,eff) :
   if not eff : return pdf
   # now we need to multiply all relevant components (i.e. all RooP2VVAngleBasis ones) 
   # of "pdf" with their efficiency corrected versions, multiply them with the right basis fcn & coefficent
   # those are assumed to be in eff....
   customizer = RooCustomizer(pdf,name)
   for c in pdf.getComponents() :
        if type(c) is not RooP2VVAngleBasis : continue
        name = "%s_eff" % c.GetName()
        s = RooArgSet()
        [ s.add( c.createProduct( fijk, cijk ) )  for (fijk,cijk) in eff ]
        rep = w.put( RooAddition_( name, name, s, True ) )  # hand over ownership & put in workspace...
        customizer.replaceArg( c, rep )
   return customizer.build(True)

def buildEffMomentsPDF(w,name,pdf,data,moments) :
   computeMoments(data,moments)
   return buildEff_x_PDF(w,name,pdf,[ ( m.basis() , m.coefficient() ) for m in moments ] )

class TimeResolutionBuilder :
    # TODO: build a PDF for sigmat (eg. RooHistPdf, RooThresholdPdf... or the sum of two gamma functions....)
    # gamma = RooRealVar("gamma","gamma",15,10,30)
    # beta = RooRealVar("beta","beta",2.21456e-03,0,0.1)
    # mu = RooRealVar("mu","mu",0,0.001,0.01)
    # pdf = RooGammaPdf("pdf","pdf",st,gamma,beta,mu)
    def __init__(self,ws, t, sigmat) :
        if type(t)      is str : t      = ws[t]
        if type(sigmat) is str : sigmat = ws[sigmat]
        # define outlier catcher
        if false :
            ws.factory("GaussModel::tres_3(%s,zero[0],tres_3_s[3,1,5])" % (t.GetName()) )
        else :
           ws.factory("{tres_3_l[1.7,0.9,3.0],tres_3_s[1,0.5,2]}")
           ws.factory("GExpModel::tres_3_gexpr(%s,tres_3_s,tres_3_l,kFALSE,Normal)"%(t.GetName()))
           ws.factory("GExpModel::tres_3_gexpl(%s,tres_3_s,tres_3_l,kFALSE,Flipped)"%(t.GetName()))
           ws.factory(" AddModel::tres_3({tres_3_gexpl,tres_3_gexpr},{half[0.5]})")

        for name in [ 'sig','nonpsi' ] :
            ws.factory("GaussModel::tres_%s_1(%s,tres_%s_m[0,-0.2,0.2],tres_%s_s1[1.1,0.3,2 ], 1, %s)" % (name,t.GetName(),name,name,sigmat.GetName()))
            if False :
               ws.factory("GaussModel::tres_%s_2(%s,tres_%s_m,            tres_%s_s2[20.,1.5,30], 1, %s)" % (name,t.GetName(),name,name,sigmat.GetName()))
            else :
               # try GExp instead of G2  -- so that the limit GExp -> G (i.e. 'lifetime'->0) it returns to the original
               # for now, the mean is forced to zero by the GExpModel code...
               ws.factory("{tres_%s_s2[1.5,0.9,3.0],tres_%s_l[1.,20.0]}"%(name,name))
               # choice: either scale lifetime with error or not... let's first try an absolute lifetime...
               # try to use the same width as in the first Gaussian!
               ws.factory("GExpModel::tres_%s_2_gexpr(%s,tres_%s_s2,tres_%s_l,%s,%s,kFALSE,Normal)"%(name,t.GetName(),name,name,sigmat.GetName(),sigmat.GetName()))
               ws.factory("GExpModel::tres_%s_2_gexpl(%s,tres_%s_s2,tres_%s_l,%s,%s,kFALSE,Flipped)"%(name,t.GetName(),name,name,sigmat.GetName(),sigmat.GetName()))
               ws.factory(" AddModel::tres_%s_2({tres_%s_2_gexpl,tres_%s_2_gexpr},{half[0.5]})"%(name,name,name))

            ws.factory("AddModel::tres_%s({tres_3,tres_%s_2,tres_%s_1},{tres_%s_f3[0.001,0.00,0.02],tres_%s_f2[0.2,0.01,1]})" % (name,name,name,name,name))
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
            if False :
                # this results in horrible wrong plots....
                ws.factory("Decay::t_%s_sl(t,0,                        %s,SingleSided)"%(name,     resname))
                ws.factory("Decay::t_%s_ml(t,t_%s_ml_tau[0.21,0.1,0.5],%s,SingleSided)"%(name,name,resname))
                ws.factory("Decay::t_%s_ll(t,t_%s_ll_tau[1.92,1.0,2.5],%s,SingleSided)"%(name,name,resname))
                ws.factory("PROD::t_%s(SUM(t_%s_fll[0.004,0,1]*t_%s_ll,t_%s_fml[0.02,0,1]*t_%s_ml,t_%s_sl)|sigmat,%s)"% (name,name,name,name,name,name,sigmatpdf[name].GetName()) )
            else :
                ws.factory("PROD::t_%s_sl(Decay(t,0,                        %s,SingleSided)|sigmat,%s)"%(name,     resname,sigmatpdf[name].GetName()))
                ws.factory("PROD::t_%s_ml(Decay(t,t_%s_ml_tau[0.21,0.1,0.5],%s,SingleSided)|sigmat,%s)"%(name,name,resname,sigmatpdf[name].GetName()))
                ws.factory("PROD::t_%s_ll(Decay(t,t_%s_ll_tau[1.92,1.0,2.5],%s,SingleSided)|sigmat,%s)"%(name,name,resname,sigmatpdf[name].GetName()))
                ws.factory("SUM::t_%s(t_%s_fll[0.004,0,1]*t_%s_ll,t_%s_fml[0.02,0,1]*t_%s_ml,t_%s_sl)"% (name,name,name,name,name,name) )
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
    def __init__(self,ws,m,m_dau1,m_dau2,mode) : # assume B-mass, J/psi mass, phi mass
        #### define J/psi mass observable & corresponding PDF
        self._mdau1 = m_dau1
        # signal J/psi mass pdf
        ws.factory("CBShape::mpsi_sig(%s,mpsi_sig_mean[3094,3090,3105],mpsi_sig_sigma[13.2,8,18],mpsi_sig_alpha[1.39,0.8,2],mpsi_sig_n[3])"%m_dau1.GetName())
        self._mdau1_sig = ws['mpsi_sig']
        # background J/psi mass pdf
        # given the narrow window, might as well take a 1st (2nd?) order polynomial...
        ws.factory("Exponential::mpsi_bkg(%s,mpsi_bkg_exp[-0.0005,-0.001,0.0])"%m_dau1.GetName())
        #  ws.factory("Chebychev::mpsi_bkg(%s,{mpsi_bkg_p1[0.2,-1,1],mpsi_bkg_p2[-0.01,-0.1,0.1]})"%m_dau1.GetName())
        self._mdau1_bkg = ws['mpsi_bkg']
        # overall J/psi mass pdf
        ws.factory("SUM::mpsi(mpsi_fjpsi[0.5,0.2,0.8]*mpsi_sig,mpsi_bkg)")
        self._mdau1_pdf = ws['mpsi']


        if mode == 'Bs2Jpsiphi':
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
            ws.factory("Chebychev::mphi_bkg(%s,{mphi_bkg_p1[0.2,-1,1],mphi_bkg_p2[-0.01,-0.1,0.1]})"%m_dau2.GetName())
            self._mdau2_bkg = ws['mphi_bkg']
            ws.factory("SUM::m_phi(mphi_fphi[0.2,0.05,0.8]*mphi_phisig,mphi_bkg)")
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
        if 'Bs' in mode :
            sigmid = 5367.4
            sigwid = 15  # in J/psi phi we have a 7 MeV resolution
        if 'Bu' in mode or 'Bd' in mode:
            sigmid = 5279.17
            sigwid = 20  # in J/psi K+ we have a 9 MeV resolution

        ws.factory('m_sig_mean[%s,%s,%s]'%(sigmid,sigmid-0.5*sigwid,sigmid+0.5*sigwid))
        ws.factory("PROD::m_sig(SUM(m_sig_f1[1]*Gaussian(%s,m_sig_mean,m_sig_sigma[7,4,12]),Gaussian(%s,m_sig_mean,m_sig_sigma2[14])),SUM(m_sig_fpsi[1]*mpsi_sig,mpsi_bkg))"%(m.GetName(),m.GetName()))
        self._m_sig = ws['m_sig']
        #if False :
        #    ws.factory("PROD::m_sig(Gaussian(m,m_sig_mean,expr('@0*@1',{m_sig_s[0.5,5],sigmam}))|sigmam)")
        
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
        if mode == 'Bs2Jpsiphi':
            ws.factory('{N_sig[1000,0,10000],N_psibkg[15000,0,30000],N_nonpsibkg[15000,0,30000]}')
        if mode == 'Bu2JpsiK' or mode == 'Bd2JpsiKstar' :
            ws.factory('{N_sig[13000,0,100000],N_psibkg[150000,0,300000],N_nonpsibkg[150000,0,300000]}')
        
        ws.factory("SUM::m_b(N_sig*m_sig,N_psibkg*m_psibkg,N_nonpsibkg*m_nonpsibkg)")
        self._m_pdf = ws['m_b']

        self._yields = RooArgList(ws.argSet('N_sig,N_psibkg,N_nonpsibkg'))

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




def declareObservables( ws, mode ):
    from math import pi

    if mode in [ 'Bs2Jpsiphi', 'Bd2JpsiKstar' ] :
        # transvercity angles
        ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f]}"%(-pi,pi))
        ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
        # helicity angles. we can also compute these from the transversity angles and add them as columns
        ws.factory("{ helcosthetaK[-1,1], helcosthetaL[-1,1], helphi[%f,%f]}"%(-pi,pi))
        ws.defineSet("helicityangles","helcosthetaK,helcosthetaL,helphi")

    # trigger
    ws.factory("unbiased[yes=+1,no=0]")
    # tag 
    ws.factory("tagdecision[bbar=+1,b=-1,untagged=0]")
    ws.factory("tagomega[0,0.500001]")

    ws.factory("decaytype[JpsiKplus=10,JpsiKmin=11,JpsiKstar0=20,JpsiKstarbar0=21,Jpsiphi=40]")


    # B, jpsi, phi mass
    #ws.factory("m[5200,5450]")
    #ws.factory("m[5250,5450]")
    ws.factory("mdau1[%f,%f]"%(3097-60,3097+40))
    #ws.factory("mdau1[%f,%f]"%(3097-60,3097+50))
    if mode == 'Bs2Jpsiphi' :
        ws.factory("m[%f,%f]"%(5366-50,5366+50))
        ws.factory("mdau2[%f,%f]"%(1019.455-12,1019.455+12)) ### Note: +- 10 Mev/c^2 keeps all of the phi signal, and kills 1/2 of the background -- but the roadmap uses +- 12 instead...
    if mode == 'Bd2JpsiKstar' :
        ws.factory("mdau2[%f,%f]"%(892-50,892+50))
    if mode in [ 'Bu2JpsiK','Bd2JpsiKstar' ]:
        ws.factory("m[%f,%f]"%(5279-50,5279+50))
    # time , time-error
    ws.factory("{t[-5,14],sigmat[0.005,0.1]}")  # Note 2 uses [0.3,14] -- so maybe switch to [-4,14] instead...

    # define a set for the observables (we can also define this from the pdf and data.
    if mode in [ 'Bs2Jpsiphi', 'Bd2JpsiKstar' ] :
        ws.defineSet("observables","m,t,sigmat,mdau1,mdau2,trcostheta,trcospsi,trphi,tagomega,tagdecision,decaytype,unbiased") 
    else :
        ws.defineSet("observables","m,t,sigmat,mdau1,tagomega,tagdecision,decaytype,unbiased") 

##
## call with buildTagging(ws, name, [ 0.25, 0.35, 0.45 ] ) 
## to define tagging categories corresponding to the intervals [0,0.25),[0.25,0.35),[0.35,0.45),[0.45,0.5]
## note that the final interval is defined 'by construction' and that the pdf returned gives 
## the efficiency to be in the i-th category, with the last category having an efficiency 1-sum_i eff_i
class TagPdfBuilder :
    def __init__(self,ws,tagcatdef,tagomega='tagomega')  :
        self._ws = ws
        if type(tagomega) == str : 
            tagomega = ws[tagomega]
        elif ws[tagomega.GetName()] != tagomega :
            raise LogicError('tagomega in ws does not match given RooAbsReal')

        ws.factory("ThresholdCategory::tagcat(%s,'untagged',0)"%tagomega.GetName() ) 
        self._tagcat = ws['tagcat']
        self._sig = ws.put(RooThresholdPdf('tagcat_sig','tagcat_sig',tagomega))# should we worry about the fact that the value is eff/binwidth, and the binwidth is not constant???
        self._psibkg = ws.put(RooThresholdPdf('tagcat_psibkg','tagcat_psibkg',tagomega))
        self._nonpsibkg = ws.put(RooThresholdPdf('tagcat_nonpsibkg','tagcat_nonpsibkg',tagomega))

        for (i,upper) in enumerate( tagcatdef ) :
            self._tagcat.addThreshold(upper,'tagcat_%d' % i)
            for (n,pdf) in [('sig',self._sig),('psibkg',self._psibkg),('nonpsibkg',self._nonpsibkg) ] :
                name = 'tagcat_%s_eff%s' % (n,i)
                eff = ws.put(RooRealVar( name, name , float(1)/(1+len(tagcatdef)), 0., 1.))
                pdf.addThreshold(upper, eff )

    def sigPdf(self)       : return self._sig
    def psibkgPdf(self)    : return self._bkg
    def nonpsibkgPdf(self) : return self._bkg

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
    ws.factory("{gamma[0.68,0.4,0.9],t_sig_dm[17.7],t_sig_dG[0.05,-0.3,0.3]}")
    ws.factory("expr::t_sig_tau('1/@0',{gamma})")
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
