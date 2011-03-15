###############################################################################
## function that returns a P2VV PDF                           ##
## mode and options are specified in the configuration object ##
################################################################
def getP2VVPDF(config) :



#################
## signal PDFs ##
###############################################################################
## buildJpsiV: function to build B->J/psiV signal PDFs               ##
##   assumes observables, parameters and resolution model (tres_sig) ##
##   exist in workspace                                              ##
###############################################################################
def buildJpsiV(config) :
  # get general settings
  mode = config.value('P2VVMode')
  BDecClass = config.value('BDecayClass')
  pdfName = config.value('sigPDFName')

  # get workspace
  ws = config.workspace()
  ws.factory("$Alias(Addition_, sum_)")

  # get amplitudes
  ReA0    = config['ReA0'].name()
  ImA0    = config['ImA0'].name()
  ReApar  = config['ReApar'].name()
  ImApar  = config['ImApar'].name()
  ReAperp = config['ReAperp'].name()
  ImAperp = config['ImAperp'].name()
  ReAS    = config['ReAS'].name()
  ImAS    = config['ImAS'].name()

  # get observables
  BLifetime = config['BLifetime'].name()
  iTag      = config['iTag'].name()
  misTag    = config['misTag'].name()
  if mode is 'Bd2JpsiKstar' : 
    fTag = config['fTag'].name()

  # get lifetime and mixing parameters
  BMeanLife = config['BMeanLife'].name()
  dGamma    = config['dGamma'].name()
  dm        = config['dm'].name()

  # get CP violation paramters
  if BDecClass is 'RooBTagDecay' :
    CCP = config['CCP'].name()
  if mode is 'Bs2Jpsiphi' :
    DCP = config['DCP'].name()
    SCP = config['SCP'].name()

  # get tagging parameters
  if BDecClass is 'RooBTagDecay' :
    dilution   = config['tagDilution'].name()
    ADilMisTag = config['ADilMisTag'].name()
    avgCEven   = config['avgCEven'].name()
    avgCOdd    = config['avgCOdd'].name()
    if mode is 'Bd2JpsiKstar' :
      ANorm = config['ANorm'].name()

  # build angular functions
  AngleFunctionBuilder(config)

  # define the relevant combinations of amplitudes
  ws.factory("expr::A0Sq('(@0 * @0 + @1 * @1)', {%s, %s})"\
      % (ReA0, ImA0))                                  # |A_0|^2
  ws.factory("expr::AparSq('(@0 * @0 + @1 * @1)', {%s, %s})"\
      % (ReApar, ImApar))                              # |A_par|^2
  ws.factory("expr::AperpSq('(@0 * @0 + @1 * @1)', {%s, %s})"\
      % (ReAperp, ImAperp))                            # |A_perp|^2
  ws.factory("expr::ASSq('(@0 * @0 + @1 * @1)', {%s, %s})"\
      % (ReAS, ImAS))                                  # |A_S|^2
  ws.factory("expr::ReA0Apar('(@0 * @2 + @1 * @3)', {%s, %s, %s, %s})"\
      % (ReA0, ImA0, ReApar, ImApar))                  # Re(A_0* A_par)
  ws.factory("expr::ReA0Aperp('(@0 * @2 + @1 * @3 )', {%s, %s, %s, %s})"\
      % (ReA0, ImA0, ReAperp, ImAperp))                # Re(A_0* A_perp)
  ws.factory("expr::ImA0Aperp('(@0 * @3 - @1 * @2)', {%s, %s, %s, %s})"\
      % (ReA0, ImA0, ReAperp, ImAperp))                # Im(A_0* A_perp)
  ws.factory("expr::ReA0AS('(@0 * @2 + @1 * @3)', {%s, %s, %s, %s})"\
      % (ReA0, ImA0, ReAS, ImAS))                      # Re(A_0* A_S)
  ws.factory("expr::ReAparAperp('(@0 * @2 + @1 * @3)', {%s, %s, %s, %s})"\
      % (ReApar, ImApar, ReAperp, ImAperp))            # Re(A_par* A_perp)
  ws.factory("expr::ImAparAperp('(@0 * @3 - @1 * @2)', {%s, %s, %s, %s})"\
      % (ReApar, ImApar, ReAperp, ImAperp))            # Im(A_par* A_perp)
  ws.factory("expr::ReAparAS('(@0 * @2 + @1 * @3)', {%s, %s, %s, %s})"\
      % (ReApar, ImApar, ReAS, ImAS))                  # Re(A_par* A_S)
  ws.factory("expr::ReAperpAS('(@0 * @2 + @1 * @3)', {%s, %s, %s, %s})"\
      % (ReAperp, ImAperp, ReAS, ImAS))                # Re(A_perp* A_S)
  ws.factory("expr::ImAperpAS('(@0 * @3 - @1 * @2)', {%s, %s, %s, %s})"\
      % (ReAperp, ImAperp, ReAS, ImAS))                # Im(A_perp* A_S)

  # build the signal PDF
  if mode is 'Bd2JpsiKstar' :
    # B0->J/psiK*

    # build time function coefficient
    ws.factory("RealSumPdf::fjpsikstar("
        "{A0Sq_basis,        AparSq_basis,    AperpSq_basis, ASSq_basis,"
        " ReA0Apar_basis,    ImA0Aperp_basis, ReA0AS_basis,             "
        " ImAparAperp_basis, ReAparAS_basis,  ImAperpAS_basis},         "
        "{A0Sq,              AparSq,          AperpSq,       ASSq,      "
        " ReA0Apar,          ImA0Aperp,       ReA0AS,                   "
        " ImAparAperp,       ReAparAS,        ImAperpAS})               ")

    if BDecClass is 'RooBDecay' :
      # use RooBDecay

      # define factors that depend on the flavour tags
      ws.factory("expr::cTag('@0 * @1 * (1. - 2. * @2)', {%s, %s, %s})"\
          % (iTag, fTag, misTag))

      # build PDF
      ws.factory("PROD::%s(fjpsikstar,\
         BDecay(%s, %s, %s, one, zero, cTag, zero, %s, tres_sig, SingleSided)"\
         % (pdfName, BLifetime, BMeanLife, dGamma, dm))

    else :
      # use RooBTagDecay
      ws.factory("PROD::%s(fjpsikstar,\
          BTagDecay(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, one,\
          tres_sig, SingleSided)"\
          % (pdfName, BLifetime, iTag, fTag, BMeanLife, dGamma, dm, dilution,
             ADilMisTag, ANorm, avgCEven, avgCOdd))

  elif mode is 'Bs2Jpsiphi' :
    # B_s0->J/psiphi

    # set strings to build time function coefficients
    coshCStr =\
        "prod(A0Sq,                          A0Sq_basis),"
        "prod(AparSq,                        AparSq_basis),"
        "prod(AperpSq,                       AperpSq_basis),"
        "prod(ASSq,                          ASSq_basis),"
        "prod(ReA0Apar,                      ReA0Apar_basis),"
        "prod(ImA0Aperp,                 {C} ImA0Aperp_basis),"
        "prod(ReA0AS,                        ReA0AS_basis),"
        "prod(ImAparAperp,               {C} ImAparAperp_basis),"
        "prod(ReAparAS,                      ReAparAS_basis),"
        "prod(ImAperpAS,                 {C} ImAperpAS_basis)"
    cosCStr =\
        "prod(A0Sq,               {cTag} {C} A0Sq_basis),"
        "prod(AparSq,             {cTag} {C} AparSq_basis),"
        "prod(AperpSq,            {cTag} {C} AperpSq_basis),"
        "prod(ASSq,               {cTag} {C} ASSq_basis),"
        "prod(ReA0Apar,           {cTag} {C} ReA0Apar_basis),"
        "prod(ImA0Aperp,          {cTag}     ImA0Aperp_basis),"
        "prod(ReA0AS,             {cTag} {C} ReA0AS_basis),"
        "prod(ImAparAperp,        {cTag}     ImAparAperp_basis),"
        "prod(ReAparAS,           {cTag} {C} ReAparAS_basis),"
        "prod(ImAperpAS,          {cTag}     ImAperpAS_basis)"
    sinhCStr =\
        "prod(A0Sq,        minus,        {D} A0Sq_basis),"
        "prod(AparSq,      minus,        {D} AparSq_basis),"
        "prod(AperpSq,                   {D} AperpSq_basis),"
        "prod(ASSq,        minus,        {D} ASSq_basis),"
        "prod(ReA0Apar,    minus,        {D} ReA0Apar_basis),"
        "prod(ReA0Aperp,                 {S} ImA0Aperp_basis),"
        "prod(ReA0AS,      minus,        {D} ReA0AS_basis),"
        "prod(ReAparAperp,               {S} ImAparAperp_basis),"
        "prod(ReAparAS,    minus,        {D} ReAparAS_basis),"
        "prod(ReAperpAS,   minus,        {S} ImAperpAS_basis)"
    sinCStr =\
        "prod(A0Sq,        minus, {cTag} {S} A0Sq_basis),"
        "prod(AparSq,      minus, {cTag} {S} AparSq_basis),"
        "prod(AperpSq,            {cTag} {S} AperpSq_basis),"
        "prod(ASSq,        minus, {cTag} {S} ASSq_basis),"
        "prod(ReA0Apar,    minus, {cTag} {S} ReA0Apar_basis),"
        "prod(ReA0Aperp,   minus, {cTag} {D} ImA0Aperp_basis),"
        "prod(ReA0AS,      minus, {cTag} {S} ReA0AS_basis),"
        "prod(ReAparAperp, minus, {cTag} {D} ImAparAperp_basis),"
        "prod(ReAparAS,    minus, {cTag} {S} ReAparAS_basis),"
        "prod(ReAperpAS,          {cTag} {D} ImAperpAS_basis)"

    if BDecClass is 'RooBDecay' :
      # use RooBDecay

      # define factors that depend on initial state flavour tag
      ws.factory("expr::cTag('@0 * (1. - 2. * @1)', {%s, %s})"\
          % (iTag, misTag))

      # format time function strings
      coshCStr = coshCoefStr.format(C='')
      cosCStr  =  cosCoefStr.format(C='', cTag='cTag,')
      sinhCStr = sinhCoefStr.format(D=DCP+',', S=SCP+',')
      sinCStr  =  sinCoefStr.format(D=DCP+',', S=SCP+',', cTag='cTag,')

      # build time function coefficients
      ws.factory('sum_::fjpsiphi_cosh({' + coshCStr + '})')
      ws.factory('sum_::fjpsiphi_cos({'  + cosCStr  + '})')
      ws.factory('sum_::fjpsiphi_sinh({' + sinhCStr + '})')
      ws.factory('sum_::fjpsiphi_sin({'  + sinCStr  + '})')

      # build PDF
      ws.factory("BDecay::%s(%s, %s, %s, fjpsiphi_cosh, fjpsiphi_sinh,\
                 fjpsiphi_cos, fjpsiphi_sin, %s, tres_sig, SingleSided))"\
          % (pdfName, BLifetime, BMeanLife, dGamma, dm))

    else :
      # use RooBTagDecay

      # format time function strings
      coshCStr = coshCoefStr.format(C=CCP+',')
      cosCStr  =  cosCoefStr.format(C=CCP+',', cTag='')
      sinhCStr = sinhCoefStr.format(D=DCP+',', S=SCP+',')
      sinCStr  =  sinCoefStr.format(D=DCP+',', S=SCP+',', cTag='')

      # build time function coefficients
      ws.factory('sum_::fjpsiphi_cosh({' + coshCStr + '})')
      ws.factory('sum_::fjpsiphi_cos({'  + cosCStr  + '})')
      ws.factory('sum_::fjpsiphi_sinh({' + sinhCStr + '})')
      ws.factory('sum_::fjpsiphi_sin({'  + sinCStr  + '})')

      # build PDF
      ws.factory("BTagDecay::%s(%s, %s, %s, %s, %s, %s, %s, %s, %s,\
          fjpsiphi_cosh, fjpsiphi_sinh, fjpsiphi_cos, fjpsiphi_sin,\
          tres_sig, SingleSided)"\
          % (pdfName, BLifetime, iTag, BMeanLife, dGamma, dm, dilution,
             ADilMisTag, avgCEven, avgCOdd))


  return ws.pdf(pdfName)


#######################
## angular functions ##
###############################################################################
class AngleFunctionBuilder :
  def __init__(self, config) :
    import RooFitDecorators
    from ROOT import RooArgSet, RooAddition_, RooP2VVAngleBasis

    # set configuration member
    self._config = config

    # get type of angles
    angType = self._config.value('anglesType')
    if angType and angType[0] is 'trans' :
      self._anglesType = 'trans'
    else :
      self._anglesType = 'hel'

    # get angles from workspace
    workspace = self._config.workspace()
    if self._anglesType is 'trans' :
      self._cpsi   = workspace[self._config['trcpsi'].name()]
      self._ctheta = workspace[self._config['trctheta'].name()]
      self._phi    = workspace[self._config['trphi'].name()]
    else :
      self._cpsi   = workspace[self._config['helcthetaK'].name()]
      self._ctheta = workspace[self._config['helcthetal'].name()]
      self._phi    = workspace[self._config['helphi'].name()]

    print 'INFO: AngleFunctionBuilder.__init__(): using %s, %s, %s'\
        % (self._cpsi, self._ctheta, self._phi)

    # specify components of angular functions
    angTerms = []
    if self._anglesType is 'trans' :
      # using transversity angles
      angTerms.append(('A0Sq',        [(0, 0, 0,  0,  2.             ),
                                       (0, 0, 2,  0,  sqrt( 1. /  5.)),
                                       (0, 0, 2,  2, -sqrt( 3. /  5.)),
                                       (2, 0, 0,  0,  4.             ),
                                       (2, 0, 2,  0,  sqrt( 4. /  5.)),
                                       (2, 0, 2,  2, -sqrt(12. /  5.))]))
      angTerms.append(('AparSq',      [(2, 2, 0,  0,  1.             ),
                                       (2, 2, 2,  0,  sqrt( 1. / 20.)),
                                       (2, 2, 2,  2,  sqrt( 3. / 20.))]))
      angTerms.append(('AperpSq',     [(2, 2, 0,  0,  1.             ),
                                       (2, 2, 2,  0, -sqrt( 1. /  5.))]))
      angTerms.append(('ReA0Apar',    [(2, 1, 2, -2, -sqrt( 6. /  5.))]))
      angTerms.append(('ImA0Aperp',   [(2, 1, 2,  1,  sqrt( 6. /  5.))]))
      angTerms.append(('ImAparAperp', [(2, 2, 2, -1,  sqrt( 3. /  5.))]))

      if self._config.value('incKSWave') :
        angTerms.append(('ASSq',        [(0, 0, 0,  0,  2.             ),
                                         (0, 0, 2,  0,  sqrt( 1. /  5.)),
                                         (0, 0, 2,  2, -sqrt( 3. /  5.))]))
        angTerms.append(('ReA0AS',      [(1, 0, 0,  0,  sqrt(48.      )),
                                         (1, 0, 2,  0,  sqrt(12. /  5.)),
                                         (1, 0, 2,  2, -sqrt(36. /  5.))]))
        angTerms.append(('ReAparAS',    [(1, 1, 2, -2, -sqrt(18. /  5.))]))
        angTerms.append(('ImAperpAS',   [(1, 1, 2,  1, -sqrt(18. /  5.))]))

    else :
      # using helicity angles
      angTerms.append(('A0Sq',        [(0, 0, 0,  0,  2.             ),
                                       (0, 0, 2,  0, -sqrt( 4. /  5.)),
                                       (2, 0, 0,  0,  4.             ),
                                       (2, 0, 2,  0, -sqrt(16. /  5.))]))
      angTerms.append(('AparSq',      [(2, 2, 0,  0,  1.             ),
                                       (2, 2, 2,  0,  sqrt( 1. / 20.)),
                                       (2, 2, 2,  2, -sqrt( 3. / 20.))]))
      angTerms.append(('AperpSq',     [(2, 2, 0,  0,  1.             ),
                                       (2, 2, 2,  0,  sqrt( 1. / 20.)),
                                       (2, 2, 2,  2,  sqrt( 3. / 20.))]))
      angTerms.append(('ReA0Apar',    [(2, 1, 2,  1,  sqrt( 6. /  5.))]))
      angTerms.append(('ImA0Aperp',   [(2, 1, 2, -1, -sqrt( 6. /  5.))]))
      angTerms.append(('ImAparAperp', [(2, 2, 2, -2,  sqrt( 3. /  5.))]))

      if self._config.value('incKSWave') :
        angTerms.append(('ASSq',        [(0, 0, 0,  0,  2.             ),
                                         (0, 0, 2,  0, -sqrt( 4. /  5.))]))
        angTerms.append(('ReA0AS',      [(1, 0, 0,  0,  sqrt(48.      )),
                                         (1, 0, 2,  0, -sqrt(48. /  5.))]))
        angTerms.append(('ReAparAS',    [(1, 1, 2,  1,  sqrt(18. /  5.))]))
        angTerms.append(('ImAperpAS',   [(1, 1, 2, -1,  sqrt(18. /  5.))]))

    # build angular function for each term in the signal PDF
    self._basis = []
    for angTerm in angTerms :
      name = angTerm[0] + '_basis'
      basesSet = RooArgSet()
      for comps in angTerm[1] :
        basesSet.add(self.buildBasisFunc(name, comps[0], comps[1], comps[2],
            comps[3], comps[4]))

      self._basis.append(self._config.workspace().put(RooAddition_(name, name,
          basesSet))

  def angles(self) : 
    from ROOT import RooArgList
    return RooArgList(self._cpsi, self._ctheta, self._phi)

  def basis(self) :
    return self._basis

  def buildBasisFunc(self, name, i, j, k, l, c) :
    import RooFitDecorators
    from ROOT import RooP2VVAngleBasis

    name = "%s_%d_%d_%d_%d" % (name, i, j, k, l)
    name.replace("-", "m")

    basisFunc = self._config.workspace().function(name)
    if not basisFunc : 
      basisFunc = self._config.workspace().put(RooP2VVAngleBasis(name, name,
          self._cpsi, self._ctheta, self._phi, i, j, k, l, c))

    return basisFunc


#############
## moments ##
###############################################################################
## Looping over data in python is quite a bit slower than in C++
## So we adapt the arguments, and then defer to the C++ _computeMoments
def computeMoments( data, moments ) :
    if not moments : return None
    vecmom = std.vector('IMoment*')()
    for m in moments : vecmom.push_back(m)
    return _computeMoments( data, vecmom )


###############################################################################
def buildMomentPDF(w,name,data,moments) :
    import RooFitDecorators
    from ROOT import RooArgList, RooRealSumPdf

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


##################
## efficiencies ##
###############################################################################
def buildEff_x_PDF(w,name,pdf,eff) :
   import RooFitDecorators
   from ROOT import RooArgSet, RooCustomizer, RooP2VVAngleBasis, RooAddition_

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
        rep = w.put( RooAddition_( n, n, s, True ) )  # hand over ownership & put in workspace...
        customizer.replaceArg( c, rep )
   return customizer.build(True)


###############################################################################
def buildEffMomentsPDF(w,name,pdf,data,moments) :
   computeMoments(data,moments)
   return buildEff_x_PDF(w,name,pdf,[ ( m.basis() , m.coefficient() ) for m in moments ] )


#####################
## time resolution ##
###############################################################################
class TimeResolutionBuilder :
  # TODO: build a PDF for sigmat (eg. RooHistPdf, RooThresholdPdf... or the sum of two gamma functions....)
  # gamma = RooRealVar("gamma","gamma",15,10,30)
  # beta = RooRealVar("beta","beta",2.21456e-03,0,0.1)
  # mu = RooRealVar("mu","mu",0,0.001,0.01)
  # pdf = RooGammaPdf("pdf","pdf",st,gamma,beta,mu)
  def __init__(self,ws, t, sigmat) :
    import RooFitDecorators

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


###############################################################################
### backwards compatibility stub    
global _timeresbuilder
_timeresbuilder = None
def buildResoModels(ws):
    global _timeresbuilder
    if not _timeresbuilder : _timeresbuilder = TimeResolutionBuilder(ws)


###############
## mass PDFs ##
###############################################################################
class MassPdfBuilder :
    ## TODO: investigate use of per-event mass error... or add a 2nd Gaussian to the b-mass
    ## TODO: integrate SPLOT functionality into the MassPdfBuilder...
    def __init__(self,ws,m,m_dau1,m_dau2,mode) : # assume B-mass, J/psi mass, phi mass
        import RooFitDecorators
        from ROOT import RooArgList

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


###############################################################################
##### backwards compatibility glue... making this a singleton...
global _MassPDFBuilder 
_MassPDFBuilder = None
def buildMassPDFs(ws) :
    import RooFitDecorators

    global _MassPDFBuilder
    if not _MassPDFBuilder :
        m_B = ws['m']
        m_mumu = ws['mdau1']
        m_KK = ws['mdau2']
        _MassPDFBuilder = MassPdfBuilder(ws,m_B,m_mumu,m_KK) # todo: make this a J/psi phi builder, so that we can also have a J/psi K* one ;-]
    return _MassPDFBuilder


#####################
## background PDFs ##
###############################################################################
class BkgTimePdfBuilder : #background propertime
    def __init__(self, ws, resbuilder, sigmatpdf ) :
        import RooFitDecorators

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


###############################################################################
class BkgAnglePdfBuilder :
    def __init__(self,ws,basis,data, opt) : 
        from ROOT import RooDataSet

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
        from ROOT import gStyle, RooArgList, RooArgSet
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


###############################################################################
def buildABkgPdf( ws, name, resname, psimasspdfname ):
    import RooFitDecorators
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


###############################################################################
def buildBkgPdf( ws, name = 'bkg_pdf' ):
    
    # assume that the resolution models and psi mass models have been built     
    nonpsibkg = buildABkgPdf(ws,'nonpsi','tres_nonpsi','mpsi_bkg')
    psibkg    = buildABkgPdf(ws,'psi',   'tres_sig',   'mpsi_sig')
    # add them
    ws.factory("SUM::%s(f_psi[0.5,0.01,1]*%s,%s)"%(name,psibkg.GetName(),nonpsibkg.GetName()))
    return ws.pdf(name)


#############
## tagging ##
###############################################################################
## call with buildTagging(ws, name, [ 0.25, 0.35, 0.45 ] ) 
## to define tagging categories corresponding to the intervals [0,0.25),[0.25,0.35),[0.35,0.45),[0.45,0.5]
## note that the final interval is defined 'by construction' and that the pdf returned gives 
## the efficiency to be in the i-th category, with the last category having an efficiency 1-sum_i eff_i
class TagPdfBuilder :
    def __init__(self,ws,tagcatdef,tagomega='tagomega')  :
        import RooFitDecorators
        from ROOT import RooRealVar, RooThresholdPdf

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


###############################################################################
def buildTagging( ws, name, tagcatdef ) :
    # either make PDF conditional on tagomega distribution
    # and use a fittable version RooHistPdf for tagOmega,
    # different for sig and bkg
    #
    # or split tagomega distribution in discrete categories,
    # and multiply by efficiency for each category, seperate 
    # for signal and background... -- or just make the fit
    # extended, and treat each bin as Poisson bkg + Poisson sig
    import RooFitDecorators
    from ROOT import RooRealVar, RooThresholdPdf

    ws.factory("misTag[0,0,0.5]")
    tagcat = ws.factory("ThresholdCategory::%s('tagomega','untagged',0)" % name+"cat")
    pdf = ws.put( RooThresholdPdf(name+'effpdf',name+'effpdf',ws['tagomega']) )
    for (i,upper) in enumerate( tagcatdef ) :
        cname = '%s%d' % (name,i)
        ename = cname + '_eff'
        eff = ws.put(RooRealVar( ename, ename , 0.2, 0., 1.))
        tagcat.addThreshold(upper,cname)
        pdf.addThreshold(upper, eff )

    return (tagcat,pdf)

