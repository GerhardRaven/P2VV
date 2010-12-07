from ROOT import *
gSystem.Load("libp2vv")

RooWorkspace.put = getattr(RooWorkspace,'import')

def buildAngularBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    def _ba(name,comp) :
        n = name + '_basis'
        s = RooArgSet()
        for c in comp : s.add( ab(name,c[0],c[1],c[2],c[3],c[4]) )
        ws.put(RooAddition_( n, n, s ) )
        return ws.function(n)

    return ( _ba("AzAz",       [ ( 0,0,0, 0, 2.), ( 0,0,2,0, sqrt(1./ 5.)), ( 0,0,2,2, -sqrt( 3./ 5.)) 
                               , ( 2,0,0, 0, 4.), ( 2,0,2,0, sqrt(4./ 5.)), ( 2,0,2,2, -sqrt(12./ 5.)) ] )
           , _ba("AparApar",   [ ( 2,2,0, 0, 1.), ( 2,2,2,0, sqrt(1./20.)), ( 2,2,2,2,  sqrt( 3./20.)) ] )
           , _ba("AperpAperp", [ ( 2,2,0, 0, 1.), ( 2,2,2,0,-sqrt(1./ 5.)) ] )
           , _ba("AparAperp",  [ ( 2,2,2,-1, sqrt(3./5.)) ] )
           , _ba("AzAperp",    [ ( 2,1,2, 1,sqrt(6./5.)) ] )
           , _ba("AzApar",     [ ( 2,1,2,-2,-sqrt(6./5.)) ] )
           )

def buildJpsiphi(ws, name) :
    basis = buildAngularBasis(ws, abasis(ws,'trcospsi','trcostheta','trphi') ) 

    # define the relevant combinations of strong amplitudes
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")

    ws.put(RooFormulaVar("qtag_","@0*(1-2*@1)",RooArgList( ws.cat('tagdecision'),ws.var('wmistag')) ) )

    ws.factory("expr::N('1-@0*@1',{qtag_,C})")
    ws.factory("Minus[-1]")
    ws.factory("$Alias(Addition_,sum_)") 

    ws.factory("sum_::fjpsiphi_cosh({ prod(N,NAzAz,                    AzAz_basis)"
                                   ", prod(N,NAparApar,                AparApar_basis)"
                                   ", prod(N,NAperpAperp,              AperpAperp_basis)"
                                   ", prod(N,ImAparAperp,      qtag_,C,AparAperp_basis)"
                                   ", prod(N,ImAzAperp,        qtag_,C,AzAperp_basis)"
                                   ", prod(N,ReAzApar,                 AzApar_basis)"
                                   "})")
    ws.factory("sum_::fjpsiphi_cos ({ prod(N,NAzAz,            qtag_,C,AzAz_basis)"
                                   ", prod(N,NAparApar,        qtag_,C,AparApar_basis)"
                                   ", prod(N,NAperpAperp,      qtag_,C,AperpAperp_basis)"
                                   ", prod(N,ImAparAperp,              AparAperp_basis)"
                                   ", prod(N,ImAzAperp,                AzAperp_basis)"
                                   ", prod(N,ReAzApar,         qtag_,C,AzApar_basis)"
                                   "})") 
    ws.factory("sum_::fjpsiphi_sinh({ prod(N,NAzAz,      Minus,      D,AzAz_basis)"
                                   ", prod(N,NAparApar,  Minus,      D,AparApar_basis)"
                                   ", prod(N,NAperpAperp,            D,AperpAperp_basis)"
                                   ", prod(N,ReAparAperp,      qtag_,S,AparAperp_basis)"
                                   ", prod(N,ReAzAperp,        qtag_,S,AzAperp_basis)"
                                   ", prod(N,ReAzApar,   Minus,      D,AzApar_basis)"
                                   "})")
    ws.factory("sum_::fjpsiphi_sin ({ prod(N,NAzAz,      Minus,qtag_,S,AzAz_basis)"
                                   ", prod(N,NAparApar,  Minus,qtag_,S,AparApar_basis)"
                                   ", prod(N,NAperpAperp,      qtag_,S,AperpAperp_basis)"
                                   ", prod(N,ReAparAperp,Minus,      D,AparAperp_basis)"
                                   ", prod(N,ReAzAperp,  Minus,      D,AzAperp_basis)"
                                   ", prod(N,ReAzApar,   Minus,qtag_,S,AzApar_basis)"
                                   "})")
    ws.factory("BDecay::%s(t,tau,dG,fjpsiphi_cosh,fjpsiphi_sinh,fjpsiphi_cos,fjpsiphi_sin,dm,res,SingleSided)" % name)
    return ws.pdf(name)


def buildJpsikstar(ws, name) :
    buildAngularBasis(ws, abasis(ws,'trcospsi','trcostheta','trphi') )
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")

    ws.put(RooFormulaVar("qmix","@0*@1*(1-2*@2)",RooArgList( ws.cat("qrec"),ws.cat("tagdecision"),ws.var("wmistag") ) ) )

    ws.factory("{Minus[-1],Zero[0],One[1]}")
    # in J/psi K*, things factorize because deltaGamma = 0 -- use this !
    ws.factory("PROD::%s( RealSumPdf( { AzAz_basis , AparApar_basis, AperpAperp_basis, AparAperp_basis, AzAperp_basis, AzApar_basis} "
                       "            , { NAzAz,       NAparApar,      NAperpAperp,      ImAparAperp,     ImAzAperp,     ReAzApar } )"
                       ", BDecay(t,tau,Zero,One,Zero,qmix,Zero,dm,res,SingleSided))"%name)
    return ws.pdf(name)

def buildMomentPDF(w,name,data,moments) :
    if not moments : return None
    allObs = moments[0].basis().getObservables(data)
    for i in range( data.numEntries() ) :
        allObs.assignValueOnly( data.get(i) )
        for m in moments : m.inc()
    coef = RooArgList()
    fact = RooArgList()
    for m in moments :
        C = 'C_%f' % m.coefficient()
        w.factory( '%s[%f]'%(C,m.coefficient() ) )
        coef.add( w.var( C ) )
        fact.add( m.basis() )
    w.put(  RooRealSumPdf(name,name,fact,coef) )
    return w.pdf(name)
