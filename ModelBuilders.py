from ROOT import *
gSystem.Load("libp2vv")

def buildAngularBasis(ws, ab) :
    #definition of the angular part of the PDF in terms of basis functions...
    ws.put(RooAddition_("AzAz_basis",            "AzAz_basis",        RooArgSet( ab("AzAz",       0,0,0, 0, 2.), ab("AzAz",       0,0,2,0, sqrt(1./ 5.)), ab("AzAz",      0,0,2,2, -sqrt( 3./ 5.))
                                                                               , ab("AzAz",       2,0,0, 0, 4.), ab("AzAz",       2,0,2,0, sqrt(4./ 5.)), ab("AzAz",      2,0,2,2, -sqrt(12./ 5.)) )))
    ws.put(RooAddition_("AparApar_basis",        "AparApar_basis",    RooArgSet( ab("AparApar",   2,2,0, 0, 1.), ab("AparApar",   2,2,2,0, sqrt(1./20.)), ab("AparApar",  2,2,2,2,  sqrt( 3./20.)) )))
    ws.put(RooAddition_("AperpAperp_basis",      "AperpAperp_basis",  RooArgSet( ab("AperpAperp", 2,2,0, 0, 1.), ab("AperpAperp", 2,2,2,0,-sqrt(1./ 5.)))))
    ws.put(RooAddition_("AparAperp_basis",       "AparAperp_basis",   RooArgSet( ab("AparAperp",  2,2,2,-1, sqrt( 9./15.)))))
    ws.put(RooAddition_("AzAperp_basis",         "AzAperp_basis",     RooArgSet( ab("AzAperp",    2,1,2, 1,-sqrt(18./15.)))))
    ws.put(RooAddition_("AzApar_basis",          "AzApar_basis",      RooArgSet( ab("AzApar",     2,1,2,-2, sqrt(18./15.)))))

def buildJpsiphi(ws, name) :
    buildAngularBasis(ws, abasis(ws,'trcospsi','trcostheta','trphi') )

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

    #ws.factory("expr::N('1-@0*@1',{qtag_,C})") #in J/psi K*, we need to drop this factor...
    ws.factory("expr::N('1/(1+@0*@1)',{qtag_,C})") #in J/psi K*, we need to drop this factor...

    ws.factory("Minus[-1]");
    ws.factory("$Alias(Addition_,sum_)") ;

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
    ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAparApar  ('( @4 * @4 + @5 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::NAperpAperp('( @2 * @2 + @3 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAparAperp('( @4 * @2 + @5 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ReAzApar   ('( @0 * @4 + @1 * @5 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAparAperp('( @4 * @3 - @5 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")
    ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,ImAz,ReAperp,ImAperp,ReApar,ImApar})")

    ws.put(RooFormulaVar("qmix","@0*@1*(1-2*@2)",RooArgList( ws.cat("qrec"),ws.cat("qtag"),ws.var("wmistag") ) ) )

    ws.factory("{Minus[-1],Zero[0],One[1]}");
    # in J/psi K*, things factorize because deltaGamma = 0 -- use this !
    ws.factory("PROD::%s( RealSumPdf( { AzAz_basis , AparApar_basis, AperpAperp_basis, AparAperp_basis, AzAperp_basis, AzApar_basis} "
                       "           , { NAzAz,       NAparApar,      NAperpAperp,      ImAparAperp,     ImAzAperp,     ReAzApar } )"
                       ", BDecay(t,tau,Zero,One,Zero,qmix,Zero,dm,res,SingleSided))"%name)
    return ws.pdf(name)
