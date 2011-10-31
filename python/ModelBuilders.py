from ROOT import *
import RooFitDecorators
gSystem.Load("libP2VV")
from math import pi

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

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
        return ws.put(RooAddition( n, n, s ) )

def _buildTransversityBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    # transversity amplitudes in terms of transversity angles
    _ba = lambda name,comp : _buildAngularFunction(ws,ab,name,comp)

    return ( _ba("AzAz",       [ ( 0,0,0, 0, 2.), ( 0,0,2,0,  sqrt(1./ 5.)), ( 0,0,2,2, -sqrt( 3./ 5.))
                               , ( 2,0,0, 0, 4.), ( 2,0,2,0,  2.*sqrt(1./ 5.)), ( 2,0,2,2, -2.*sqrt(3./ 5.)) ] )
             , _ba("AparApar",   [ ( 2,2,0, 0, 1.), ( 2,2,2,0,  (1./2.)*sqrt(1./5.)), ( 2,2,2,2,  (1./2.)*sqrt( 3./5.)) ] )
             , _ba("AperpAperp", [ ( 2,2,0, 0, 1.), ( 2,2,2,0, -sqrt(1./ 5.)) ] )
             , _ba("AparAperp",  [ ( 2,2,2,-1,  sqrt(3./5.)) ] )
             , _ba("AzApar",     [ ( 2,1,2,-2, -sqrt(6./5.)) ] )
             , _ba("AzAperp",    [ ( 2,1,2, 1,  sqrt(6./5.)) ] )
             , _ba("AsAs",       [ ( 0,0,0, 0, 2.), ( 0,0,2,0,  sqrt(1./ 5.)), ( 0,0,2,2, -sqrt( 3./ 5.)) ] )
             , _ba("AsApar",     [ (1,1,2,-2, -3.*sqrt(2./5.))] )
             , _ba("AsAperp",    [ (1,1,2,1, 3.*sqrt(2./5.))] )
             , _ba("AsAz",       [ (1,0,0,0, 4.*sqrt(3)),(1,0,2,0, 2.*sqrt(3./5.)),(1,0,2,2, -6.*sqrt(1./5.))] )
             )

def _buildHelicityBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    # transversity amplitudes in terms of helicity angles
    _ba = lambda name,comp : _buildAngularFunction(ws,ab,name,comp)

    return ( _ba("AzAz",       [ ( 0,0,0, 0, 2.), (0,0,2,0, -sqrt(4./5.))
                               , ( 2,0,0, 0, 4.), (2,0,2,0, -sqrt(16./5.)) ] )
             , _ba("AparApar",   [ ( 2,2,0, 0, 1.), (2,2,2, 0, sqrt(1./20.)), ( 2,2,2,2,  -sqrt(3./20.)) ] )
             , _ba("AperpAperp", [ ( 2,2,0, 0, 1.), ( 2,2,2,0, sqrt(1./ 20.)), (2,2,2,2,sqrt(3./20.)) ] )
             , _ba("AparAperp",  [ ( 2,2,2,-2,  sqrt(3./5.)) ] )
             , _ba("AzApar",     [ ( 2,1,2, 1,  sqrt(6./5.)) ] )
             , _ba("AzAperp",    [ ( 2,1,2,-1, -sqrt(6./5.)) ] )
             , _ba("AsAs",       [ ( 0,0,0, 0, 2.), ( 0,0,2,0,  -2.*sqrt(1./ 5.)) ] )
             , _ba("AsApar",     [ (1,1,2,1, 3.*sqrt(2./5.)) ] )
             , _ba("AsAperp",    [ (1,1,2,-1, -3.*sqrt(2./5.))] )
             , _ba("AsAz",       [ (1,0,0,0, 4.*sqrt(3)), (1,0,2,0, -4.*sqrt(3./5.)) ] )
             )
## TODO: replace hardwired use of 'tagomega' with a passable rule, which by default returns 'tagomega'
##       this is needed when fitting for mistag rate in tagging categories...
def buildJpsiphi(ws, name, transversity,resoname) : # TODO: add tagsplit
    afb = AngleFunctionBuilder(ws, 'transversity' if transversity else 'helicity' )
    basis = afb.basis()

    ws.factory("expr::qtag_('@0*(1-2*@1)',{tagdecision,wtag})")
    ws.factory("expr::N('1.0/(1.0+@0*@1)',{tagdecision,C})")
    ws.factory("Minus[-1]")

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

    ws.factory("sum::fjpsiphi_cosh( prod(N,NAzAz,                    AzAz_basis)"
                                 ", prod(N,NAparApar,                AparApar_basis)"
                                 ", prod(N,NAperpAperp,              AperpAperp_basis)"
                                 ", prod(N,ImAparAperp,            C,AparAperp_basis)"
                                 ", prod(N,ImAzAperp,              C,AzAperp_basis)"
                                 ", prod(N,ReAzApar,                 AzApar_basis)"
                                 ")")
    ws.factory("sum::fjpsiphi_cos ( prod(N,NAzAz,            qtag_,C,AzAz_basis)"
                                 ", prod(N,NAparApar,        qtag_,C,AparApar_basis)"
                                 ", prod(N,NAperpAperp,      qtag_,C,AperpAperp_basis)"
                                 ", prod(N,ImAparAperp,      qtag_,  AparAperp_basis)"
                                 ", prod(N,ImAzAperp,        qtag_,  AzAperp_basis)"
                                 ", prod(N,ReAzApar,         qtag_,C,AzApar_basis)"
                                 ")")
    ws.factory("sum::fjpsiphi_sinh( prod(N,NAzAz,      Minus,      D,AzAz_basis)"
                                 ", prod(N,NAparApar,  Minus,      D,AparApar_basis)"
                                 ", prod(N,NAperpAperp,            D,AperpAperp_basis)"
                                 ", prod(N,ReAparAperp,            S,AparAperp_basis)"
                                 ", prod(N,ReAzAperp,              S,AzAperp_basis)"
                                 ", prod(N,ReAzApar,   Minus,      D,AzApar_basis)"
                                 ")")
    ws.factory("sum::fjpsiphi_sin ( prod(N,NAzAz,      Minus,qtag_,S,AzAz_basis)"
                                 ", prod(N,NAparApar,  Minus,qtag_,S,AparApar_basis)"
                                 ", prod(N,NAperpAperp,      qtag_,S,AperpAperp_basis)"
                                 ", prod(N,ReAparAperp,Minus,qtag_,D,AparAperp_basis)"
                                 ", prod(N,ReAzAperp,  Minus,qtag_,D,AzAperp_basis)"
                                 ", prod(N,ReAzApar,   Minus,qtag_,S,AzApar_basis)"
                                 ")")

    ws.factory("BDecay::%s(t,t_sig_tau,t_sig_dG,fjpsiphi_cosh,fjpsiphi_sinh,fjpsiphi_cos,fjpsiphi_sin,t_sig_dm,%s,SingleSided)" %(name,resoname))
    return ws.pdf(name)

#############
### SWAVE ###
#############

def buildJpsiphiSWave(ws, name, transversity,resoname) : # TODO: add tagsplit
    afb = AngleFunctionBuilder(ws, 'transversity' if transversity else 'helicity' )
    basis = afb.basis()
    
    ws.factory("expr::qtag_('@0*(1-2*@1)',{tagdecision,wtag})")
    ws.factory("expr::N('1.0/(1.0+@0*@1)',{tagdecision,C})")
    ws.factory("Minus[-1]")

    ws.factory("sum::fjpsiphi_cosh( prod(N,NAzAz,                    AzAz_basis)"
                                 ", prod(N,NAparApar,                AparApar_basis)"
                                 ", prod(N,NAperpAperp,              AperpAperp_basis)"
                                 ", prod(N,ImAparAperp,            C,AparAperp_basis)"
                                 ", prod(N,ImAzAperp,              C,AzAperp_basis)"
                                 ", prod(N,ReAzApar,                 AzApar_basis)"

                                 ", prod(N,NAsAs,                    AsAs_basis)"
                                 ", prod(N,ImAsAperp,                AsAperp_basis)"
                                 ", prod(N,ReAsAz,                 C,AsAz_basis)"
                                 ", prod(N,ReAsApar,               C,AsApar_basis)"
                                 ")")
    ws.factory("sum::fjpsiphi_cos ( prod(N,NAzAz,            qtag_,C,AzAz_basis)"
                                 ", prod(N,NAparApar,        qtag_,C,AparApar_basis)"
                                 ", prod(N,NAperpAperp,      qtag_,C,AperpAperp_basis)"
                                 ", prod(N,ImAparAperp,      qtag_,  AparAperp_basis)"
                                 ", prod(N,ImAzAperp,        qtag_,  AzAperp_basis)"
                                 ", prod(N,ReAzApar,         qtag_,C,AzApar_basis)"

                                 ", prod(N,NAsAs,            qtag_,C,AsAs_basis)"
                                 ", prod(N,ImAsAperp,        qtag_,C,AsAperp_basis)"
                                 ", prod(N,ReAsAz,           qtag_,  AsAz_basis)"
                                 ", prod(N,ReAsApar,         qtag_,  AsApar_basis)"
                                 ")")
    ws.factory("sum::fjpsiphi_sinh( prod(N,NAzAz,      Minus,      D,AzAz_basis)"
                                 ", prod(N,NAparApar,  Minus,      D,AparApar_basis)"
                                 ", prod(N,NAperpAperp,            D,AperpAperp_basis)"
                                 ", prod(N,ReAparAperp,            S,AparAperp_basis)"
                                 ", prod(N,ReAzAperp,              S,AzAperp_basis)"
                                 ", prod(N,ReAzApar,   Minus,      D,AzApar_basis)"

                                 ", prod(N,NAsAs,                  D,AsAs_basis)"
                                 ", prod(N,ImAsAperp,              D,AsAperp_basis)"
                                 ", prod(N,ImAsAz,                 S,AsAz_basis)"
                                 ", prod(N,ImAsApar,               S,AsApar_basis)"
                                 ")")
    ws.factory("sum::fjpsiphi_sin ( prod(N,NAzAz,      Minus,qtag_,S,AzAz_basis)"
                                 ", prod(N,NAparApar,  Minus,qtag_,S,AparApar_basis)"
                                 ", prod(N,NAperpAperp,      qtag_,S,AperpAperp_basis)"
                                 ", prod(N,ReAparAperp,Minus,qtag_,D,AparAperp_basis)"
                                 ", prod(N,ReAzAperp,  Minus,qtag_,D,AzAperp_basis)"
                                 ", prod(N,ReAzApar,   Minus,qtag_,S,AzApar_basis)"

                                 ", prod(N,NAsAs,            qtag_,S,AsAs_basis)"
                                 ", prod(N,ImAsAperp,        qtag_,S,AsAperp_basis)"
                                 ", prod(N,ImAsAz,     Minus,qtag_,D,AsAz_basis)"
                                 ", prod(N,ImAsApar,   Minus,qtag_,D,AsApar_basis)"
                                 ")")
    
    
    ws.factory("BDecay::%s(t,t_sig_tau,t_sig_dG,fjpsiphi_cosh,fjpsiphi_sinh,fjpsiphi_cos,fjpsiphi_sin,t_sig_dm,%s,SingleSided)"%(name,resoname))
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

def buildJpsikstarAnglePdf(ws, name, transversity) :
    afb = AngleFunctionBuilder(ws, 'transversity' if transversity else 'helicity' )
    basis = afb.basis()

    #_buildTransversityBasis(ws, abasis(ws,ws.set('transversityangles')))
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")

        # in J/psi K*, things factorize because deltaGamma = 0 -- use this !
    ws.factory("RealSumPdf::%s( { AzAz_basis , AparApar_basis, AperpAperp_basis, AparAperp_basis, AzAperp_basis, AzApar_basis} "
               "            , { NAzAz,       NAparApar,      NAperpAperp,      ImAparAperp,     ImAzAperp,     ReAzApar } )"%name)
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
        n = "%s_%s_eff" % (name,c.GetName())
        s = RooArgSet()
        [ s.add( c.createProduct( fijk, cijk ) )  for (fijk,cijk) in eff ]
        rep = w.put( RooAddition( n, n, s, True ) )  # hand over ownership & put in workspace...
        customizer.replaceArg( c, rep )
   return customizer.build(True)

def buildEffMomentsPDF(w,name,pdf,data,moments) :
   computeMoments(data,moments)
   return buildEff_x_PDF(w,name,pdf,[ ( m.basis() , m.coefficient() ) for m in moments ] )

class AngleFunctionBuilder :
    def __init__( self, ws, basis ) :
       # TODO: move knowledge of which set to use into the basis builder itself...
       lookup = { 'transversity' : ( _buildTransversityBasis, 'transversityangles' ) 
                  , 'helicity'     : ( _buildHelicityBasis,     'helicityangles' )
                  }
       (b,v) = lookup[basis]
       self._basis = b( ws, abasis(ws, ws.set(v) ) )
    def basis(self) : return self._basis

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
