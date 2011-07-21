###############################################################################
## P2VVConfiguration: configuration of P2VV settings                         ##
##                                                                           ##
## authors:                                                                  ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                ##
##                                                                           ##
###############################################################################

def getP2VVConfig(mode = '', options = [], createWS = True) :
  """function that returns default P2VV configuration for a given mode
  """

  from math import pi, sqrt, cos, sin
  from ROOT import RooFit

  # create configuration instance
  if createWS :
    config = P2VVConfiguration('', '')
  else :
    config = P2VVConfiguration(False, '')

  # check arguments
  if type(mode) is not str :
    print "P2VV - ERROR: getP2VVConfig: argument 'mode' is not a string"
    return config
  if type(options) is not list :
    print "P2VV - ERROR: getP2VVConfig: argument 'options' is not a list"
    return config

  # set mode in configuration
  config.addSetting('P2VVMode', P2VVSetting('mode',
      'P2VV mode', mode))

  # get options
  optDict = {}
  for opt in options :
    if type(opt) is str :
      optSplit = opt.split('=')
      optDict[optSplit[0]] = optSplit[1] if len(optSplit) > 1 else ''

  config.addSetting('P2VVOptions', P2VVSetting('options',
      'P2VV options', optDict))


  ######################
  ## general settings ##
  ######################
  config.addSetting('sigPDFName', P2VVSetting('sigPDFName',
      'ROOT name of the signal PDF', mode + 'PDF'))

  onlySignal = False
  if 'onlySignal' in options : onlySignal = True
  config.addSetting('onlySignal', P2VVSetting('onlySignal',
      'build only the signal PDF?', onlySignal))

  if 'noFactFlavSpec' in options :
    config.addSetting('noFactFlavSpec', P2VVSetting('noFactFlavSpec',
        'factorize PDF for flavour-specific states?', True))

  config.addSetting('P2VVFrameDrawOpts', P2VVSetting('P2VVFrameDrawOpts',
      'options for drawing plot frames', []))

  config.addSetting('P2VVDataPlotOpts', P2VVSetting('P2VVDataPlotOpts',
      'options for plotting data', [RooFit.MarkerStyle(8),
      RooFit.MarkerSize(0.4), RooFit.DrawOption("P")]))

  config.addSetting('P2VVPDFPlotOpts', P2VVSetting('P2VVPDFPlotOpts',
      'options for plotting PDFs', [RooFit.LineWidth(2)]))


  ###############
  ## constants ##
  ###############
  config.addSetting('zeroConst', RooRealSetting('zero',
      'zero constant', 'par', 0., '', ''))
  config.addSetting('halfConst', RooRealSetting('half',
      'one half constant', 'par', 0.5, '', ''))
  config.addSetting('oneConst', RooRealSetting('one',
      'one constant', 'par', 1., '', ''))
  config.addSetting('minusConst', RooRealSetting('minus',
      'minus constant', 'par', -1., '', ''))


  ##################
  ## decay angles ##
  ##################
  anglesType = []
  if 'anglesType' in optDict and optDict['anglesType'] != '' :
    angType = optDict['anglesType'].lower()
    if 'hel'   in angType : anglesType.append('hel')
    if 'trans' in angType : anglesType.append('trans')

  elif mode in ['Bd2mumuKstar', 'Bd2JpsiKstar', 'Bs2Jpsiphi'] :
    anglesType = ['hel']

  # helicity angles
  if 'hel' in anglesType :
    config.addSetting('cpsiAng', RooRealSetting('helcthetaK',
        'cosine of kaon polarization angle', 'obs', 0., -1., 1.))
    config.addSetting('cthetaAng', RooRealSetting('helcthetaL',
        'cosine of lepton polarization angle', 'obs', 0., -1., 1.))
    config.addSetting('phiAng', RooRealSetting('helphi',
        'angle between decay planes', 'obs', 0., -pi, pi))

  # transversity angles
  if 'trans' in anglesType :
    config.addSetting('cpsiAng', RooRealSetting('trcpsi',
        'cosine of kaon polarization angle', 'obs', 0., -1., 1.))
    config.addSetting('cthetaAng', RooRealSetting('trctheta',
        'cosine of transversity polar angle', 'obs', 0., -1., 1.))
    config.addSetting('phiAng', RooRealSetting('trphi',
        'transversity azimuthal angle', 'obs', 0., -pi, pi))

  config.addSetting('anglesType', P2VVSetting('anglesType',
      'type of angles: helicity or transversity', anglesType))


  #########################
  ## B0->mumuK* settings ##
  #########################
  if mode == 'Bd2mumuKstar' :
    # amplitudes
    config.addSetting('incMuMass', P2VVSetting('incMuMass',
      'include finite muon mass?', True))
    config.addSetting('incMuMuT', P2VVSetting('incMuMuT',
      'include timelike dimuon contribution?', True))
    config.addSetting('incMuMuS', P2VVSetting('incMuMuS',
      'include scalar dimuon contribution?', True))


  #########################################################
  ## B0->J/psiK, B0->J/psiK* and B_s0->J/psiphi settings ##
  #########################################################
  elif mode in ['Bu2JpsiK', 'Bd2JpsiKstar', 'Bs2Jpsiphi'] :
    ### general ###

    # RooFit class for time PDF
    if 'RooBDecay' in optDict : BDecClass = 'RooBDecay'
    else : BDecClass = 'RooBTagDecay'
    config.addSetting('BDecayClass', P2VVSetting('BDecayClass',
        'RooFit class for time PDF', BDecClass))

    #if not onlySignal :
      # decay type
      #config.addSetting('decayType', RooCatSetting('decaytype',
      #    'J/psiX decay type', True, 'Jpsiphi', {10:'JpsiKplus',
      #    11:'JpsiKmin', 20:'JpsiKstar0', 21:'JpsiKstarbar0', 40:'Jpsiphi'}))

      # trigger
      #config.addSetting('unbiased', RooCatSetting('unbiased',
      #    'unbiased trigger?', True, 'yes', {1 : 'yes', 0 : 'no'}))

    ### B lifetime and mixing ###
    tResModel = 'Gauss'
    if 'tResModel' in optDict and optDict['tResModel'] != '' :
      tResModel = optDict['tResModel']

    config.addSetting('tResModel', P2VVSetting('tResModel',
      'time resolution model', tResModel))

    config.addSetting('BLifetime', RooRealSetting('t',
        'B lifetime (ps)', 'obs', 0., -0.5, 5.))
    if tResModel == 'Gauss' :
      config.addSetting('BLifetimeError', RooRealSetting('sigmat',
          'B lifetime error (ps)', 'par', 0.05, '', ''))
    elif tResModel[:6] == '3Gauss' :
      config.addSetting('BLifetimeError', RooRealSetting('sigmat',
          'B lifetime error (ps)', 'par', 0.05, 0.005, 0.1))

    if mode == 'Bd2JpsiKstar' :
      config.addSetting('Gamma', RooRealSetting('Gamma',
          'Gamma (ps^-1)', 'par', 0.65, 0.4, 0.9))
      config.addSetting('dGamma', RooRealSetting('dGamma',
          'delta Gamma_s (ps^-1)', 'par', 0., '', ''))
      config.addSetting('dm', RooRealSetting('dm',
          'delta m_s (ps^-1)', 'par', 0.51, -1., 2.))
    elif mode == 'Bs2Jpsiphi' :
      config.addSetting('Gamma', RooRealSetting('gamma',
          'Gamma_s (ps^-1)', 'par', 0.68, 0.4, 0.9))
      config.addSetting('dGamma', RooRealSetting('dGamma',
          'delta Gamma_s (ps^-1)', 'par', 0.05, -0.3, 0.3))
      config.addSetting('dm', RooRealSetting('dm',
          'delta m_s (ps^-1)', 'par', 17.8, 15., 20.))

    config.addSetting('BMeanLife', RooFormSetting('BMeanLife',
        'B mean lifetime', '1. / @0', ['Gamma']))

    if not onlySignal :
      ### masses ###
      #config.addSetting('mJpsi', RooRealSetting('mdau1',
      #    'J/psi mass (MeV)', 'obs', 3097., 3097. - 60., 3097. + 40.))
      if mode in ['Bu2JpsiK', 'Bd2JpsiKstar'] :
        config.addSetting('BMass', RooRealSetting('BMass',
            'B0 mass (MeV)', 'obs', 5279., 5279. - 50., 5279. + 50.))
        #if mode == 'Bd2JpsiKstar' :
        #  config.addSetting('mKstar', RooRealSetting('mdau2',
        #      'K*0 mass (MeV)', 'obs', 892., 892. - 50., 892. + 50.))
      elif mode == 'Bs2Jpsiphi' :
        config.addSetting('BMass', RooRealSetting('BMass',
            'B_s0 mass (MeV)', 'obs', 5366., 5366. - 50., 5366. + 50.))
        #config.addSetting('mphi', RooRealSetting('mdau2',
        #    'phi mass (MeV)', 'obs', 1019.455, 1019.455 - 12., 1019.455 + 12.))
        # Note: +-10 Mev/c^2 keeps all of the phi signal, and kills 1/2 of the
        # background -- but the roadmap uses +-12 instead...


    if mode in ['Bd2JpsiKstar', 'Bs2Jpsiphi'] :
      ### CP violation ###

      # initial CP violation parameter values
      lambdaSq = 1.
      phi_s    = -0.2

      ReLambda = sqrt(lambdaSq) * cos(-phi_s)
      ImLambda = sqrt(lambdaSq) * sin(-phi_s)

      CCP = (1. - lambdaSq) / (1. + lambdaSq)
      DCP = 2. * ReLambda   / (1. + lambdaSq)
      SCP = 2. * ImLambda   / (1. + lambdaSq)

      lambdaCPType = 'polar'
      if 'lambdaCPType' in optDict and optDict['lambdaCPType'] != '' :
        lambdaCPType = optDict['lambdaCPType']
      config.addSetting('lambdaCPType', P2VVSetting('lambdaCPType',
        'type of lambda (CP violation)', lambdaCPType))

      if lambdaCPType == 'polar' :
        config.addSetting('lambdaCPSq', RooRealSetting('lambda^2',
            'CP violation param. |lambda|^2', 'par', lambdaSq, 0., 5.))
        config.addSetting('CCP', RooFormSetting('C',
            'lambda param. C', '(1. - @0) / (1. + @0)', ['lambdaCPSq']))
        if mode == 'Bs2Jpsiphi' :
          config.addSetting('phiCP', RooRealSetting('phis',
              'CP violation param. phi_s', 'par', phi_s, -2 * pi, 2 * pi))
          config.addSetting('DCP', RooFormSetting('D',
              'B0_s lambda param. D', '2. * sqrt(@0) * cos(@1) / (1. + @0)',
              ['lambdaCPSq', 'phiCP']))
          config.addSetting('SCP', RooFormSetting('S',
              'B0_s lambda param. S', '-2. * sqrt(@0) * sin(@1) / (1. + @0)',
              ['lambdaCPSq', 'phiCP']))

      elif lambdaCPType == 'cartesian' :
        config.addSetting('ReLambdaCP', RooRealSetting('ReLambda',
            'CP violation param. Re(lambda)', 'par', ReLambda, -2., 2.))
        config.addSetting('ImLambdaCP', RooRealSetting('ImLambda',
            'CP violation param. Im(lambda)', 'par', ImLambda, -2., 2.))
        config.addSetting('lambdaCPSq', RooFormSetting('lambdaSq',
            'CP violation param. |lambda|^2', '@0 * @0 + @1 * @1',
            ['ReLambdaCP', 'ImLambdaCP']))
        config.addSetting('CCP', RooFormSetting('C',
            'B0_s lambda param. C', '(1. - @0) / (1. + @0)', ['lambdaCPSq']))
        if mode == 'Bs2Jpsiphi' :
          config.addSetting('DCP', RooFormSetting('D',
              'B0_s lambda param. D', '2. * @1 / (1. + @0)',
              ['lambdaCPSq', 'ReLambdaCP']))
          config.addSetting('SCP', RooFormSetting('S',
              'B0_s lambda param. S', '2. * @1 / (1. + @0)',
              ['lambdaCPSq', 'ImLambdaCP']))

      else :
        config.addSetting('CCP', RooRealSetting('C',
            'B0_s lambda param. C', 'par', CCP, -3., 3.))
        if mode == 'Bs2Jpsiphi' :
          config.addSetting('DCP', RooRealSetting('D',
              'B0_s lambda param. D', 'par', DCP, -3., 3.))
          config.addSetting('SCP', RooRealSetting('S',
              'B0_s lambda param. S', 'par', SCP, -3., 3.))

      ### flavour tags ###
      if 'allowITagZero' in optDict : allowITagZero = True
      else : allowITagZero = False
      config.addSetting('allowITagZero', P2VVSetting('allowITagZero',
          'allow zero value for init. state tag (untagged)', allowITagZero))

      asymType = 'direct'
      if 'asymType' in optDict and optDict['asymType'] != '' :
        asymType = optDict['asymType']
      config.addSetting('asymType', P2VVSetting('asymType',
        'type of CP asymmetries', asymType))

      tagType = 'singleCat'
      numTagCats = 1
      if 'tagType' in optDict and optDict['tagType'] != '' :
        tagType = optDict['tagType']

        if BDecClass == 'RooBDecay' and tagType != 'singleCat' :
          print "P2VV - ERROR: getP2VVConfig: tagging categories are not implemented for RooBDecay"
          tagType = 'singleCat'

        if tagType[:10] == 'categories' :
          try : numTagCats = int(tagType[10:])
          except : numTagCats = 2
          tagType = tagType[:10]
          if numTagCats < 2 : numTagCats = 2

      config.addSetting('tagType', P2VVSetting('tagType',
          'flavour tagging implementation', tagType))
      config.addSetting('numTagCats', P2VVSetting('numTagCats',
          'number of flavour tagging categories', numTagCats))

      # initial state flavour tag
      if allowITagZero :
        config.addSetting('iTag', RooCatSetting('iTag',
            'initial state flavour tag', True, 'B',
            {+1 : 'B', -1 : 'Bbar', 0 : 'untagged'}))
      else :
        config.addSetting('iTag', RooCatSetting('iTag',
            'initial state flavour tag', True, 'B', {+1 : 'B', -1 : 'Bbar'}))

      # final state flavour tag
      if mode == 'Bd2JpsiKstar' :
        config.addSetting('fTag', RooCatSetting('fTag',
            'final state flavour tag', True, 'B', {+1 : 'B', -1 : 'Bbar'}))

      # production and normalization asymmetries
      if not 'coef' in asymType.lower()  :
        config.addSetting('AProd', RooRealSetting('AProd',
            'production asymmetry', 'par', 0., -1., 1.))
      if mode == 'Bd2JpsiKstar' or not 'coef' in asymType.lower() :
        config.addSetting('ANorm', RooFormSetting('ANorm',
            'normalization asymmetry', '-@0', ['CCP']))

      # tagging parameters
      if tagType == 'categories' :
        tagCats = {}
        for tagCatIter in range(numTagCats) :
          tagCat = '{0:02d}'.format(tagCatIter)
          tagCats[tagCatIter] = 'cat' + tagCat

          if tagCatIter == 0 :
            # wrong tag parameters (untagged)
            config.addSetting('tagDilution' + tagCat,
                RooRealSetting('tagDilution' + tagCat,
                'tagging dilution', 'par', 0., '', ''))
            config.addSetting('ADilWTag' + tagCat,
                RooRealSetting('ADilWTag' + tagCat,
                'dilution/wrong tag asymmetry', 'par', 0., '', ''))

            # average even and odd coefficients
            if asymType in ['coefficients', 'equalCoefs'] :
              config.addSetting('avgCEvenSum', RooRealSetting('avgCEvenSum',
                  'sum of CP average even coefficients', 'par', 1., '', ''))
              config.addSetting('avgCOddSum', RooRealSetting('avgCOddSum',
                  'sum of CP average odd coefficients', 'par', 0., -2., 2.))
            else :
              config.addSetting('avgCEvenSum', RooFormSetting('avgCEvenSum',
                  'sum of CP average even coefficients', '1. + @0*@1',
                  ['AProd', 'ANorm']))
              config.addSetting('avgCOddSum', RooFormSetting('avgCOddSum',
                  'sum of CP average odd coefficients', '@0 + @1',
                  ['AProd', 'ANorm']))

          else :
            # tagging category coefficients
            config.addSetting('tagCatCoef' + tagCat,
                RooRealSetting('tagCatCoef' + tagCat,
                'tagging category coefficient', 'par',
                (1. - float(tagCatIter) / float(numTagCats)) / float(numTagCats),
                0., 1.))

            # wrong tag parameters
            config.addSetting('wTag' + tagCat,
                RooRealSetting('wTag' + tagCat,
                'wrong tag fraction B', 'par',
                0.5 * (1. - float(tagCatIter) / float(numTagCats)), '', ''))
            config.addSetting('wTagBar' + tagCat,
                RooRealSetting('wTagBar' + tagCat,
                'wrong tag fraction anti-B', 'par',
                0.5 * (1. - float(tagCatIter) / float(numTagCats)), '', ''))
            config.addSetting('tagDilution' + tagCat,
                RooFormSetting('tagDilution' + tagCat,
                'tagging dilution', '1. - @0 - @1',
                ['wTag' + tagCat, 'wTagBar' + tagCat]))
            config.addSetting('ADilWTag' + tagCat,
                RooFormSetting('ADilWTag' + tagCat,
                'dilution/wrong tag asymmetry', '(@0 - @1) / @2',
                ['wTag' + tagCat, 'wTagBar' + tagCat, 'tagDilution' + tagCat]))

            # average even and odd coefficients
            if asymType == 'coefficients' :
              config.addSetting('avgCEven' + tagCat,
                  RooRealSetting('avgCEven' + tagCat,
                  'CP average even coefficients', 'par', 1., 0., 2.))
              config.addSetting('avgCOdd' + tagCat,
                  RooRealSetting('avgCOdd' + tagCat,
                  'CP average odd coefficients', 'par', 0., -2., 2.))
            elif asymType != 'equalCoefs' :
              config.addSetting('ATagEff' + tagCat,
                  RooRealSetting('ATagEff' + tagCat,
                  'tagging efficiency asymmetry', 'par', 0., -1., 1.))
              config.addSetting('avgCEven' + tagCat,
                  RooFormSetting('avgCEven' + tagCat,
                  'CP average even coefficients', '1. + @0*@1 + @0*@2 + @1*@2',
                  ['AProd', 'ANorm', 'ATagEff' + tagCat]))
              config.addSetting('avgCOdd' + tagCat,
                  RooFormSetting('avgCOdd' + tagCat,
                  'CP average odd coefficients', '@0 + @1 + @2 + @0*@1*@2',
                  ['AProd', 'ANorm', 'ATagEff' + tagCat]))


        # tagging category
        config.addSetting('tagCat', RooCatSetting('tagCat',
            'tagging category', True, 'cat01', tagCats))

      else :
        # wrong tag parameters
        config.addSetting('wTag', RooRealSetting('wTag',
            'wrong tag fraction B', 'par', 0., '', ''))
        config.addSetting('wTagBar', RooRealSetting('wTagBar',
            'wrong tag fraction anti-B', 'par', 0., '', ''))
        config.addSetting('tagDilution',
            RooFormSetting('tagDilution', 'tagging dilution',
            '1. - @0 - @1', ['wTag', 'wTagBar']))
        config.addSetting('ADilWTag',
            RooFormSetting('ADilWTag', 'dilution/wrong tag asymmetry',
            '(@0 - @1) / @2', ['wTag', 'wTagBar', 'tagDilution']))

        # average even and odd coefficients
        if asymType == 'coefficients' :
          config.addSetting('avgCEven', RooRealSetting('avgCEven',
              'CP average even coefficients', 'par', 1., '', ''))
          config.addSetting('avgCOdd', RooRealSetting('avgCOdd',
              'CP average odd coefficients', 'par', 0., -2., 2.))
        elif asymType != 'equalCoefs' :
          config.addSetting('avgCEven', RooFormSetting('avgCEven',
              'CP average even coefficients', '1. + @0*@1',
              ['AProd', 'ANorm']))
          config.addSetting('avgCOdd', RooFormSetting('avgCOdd',
              'CP average odd coefficients', '@0 + @1',
              ['AProd', 'ANorm']))


      ### amplitudes ###

      # initial values (J/psiK*) from Phys. Rev. D76 031102 (2007)
      # (arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009)
      # S-wave amplitude from first LHCb results
      A0Sq    = 0.556
      AparSq  = 0.211
      AparPh  = -(-2.93)
      AperpSq = 0.233
      AperpPh = pi - 2.91
      ASSq    = 0.05
      ASPh    = 2.2

      ReA0    = sqrt(A0Sq)
      ReApar  = sqrt(AparSq)  * cos(AparPh)
      ImApar  = sqrt(AparSq)  * sin(AparPh)
      ReAperp = sqrt(AperpSq) * cos(AperpPh)
      ImAperp = sqrt(AperpSq) * sin(AperpPh)
      ReAS    = sqrt(ASSq)    * cos(ASPh)
      ImAS    = sqrt(ASSq)    * sin(ASPh)

      KSWave = 'none'
      if 'KSWave' in optDict and type(optDict['KSWave']) is str\
          and optDict['KSWave'] != '' :
        KSWave = optDict['KSWave']
        if mode == 'Bs2Jpsiphi' and KSWave == 'include' : KSWave = 'includeOdd'
      config.addSetting('KSWave', P2VVSetting('KSWave',
        'include KK or Kpi S-wave?', KSWave))

      ampsType = 'transCartesian'
      if 'ampsType' in optDict and optDict['ampsType'] != '' :
        ampsType = optDict['ampsType']
      config.addSetting('ampsType', P2VVSetting('ampsType',
        'type of amplitudes', ampsType))

      if ampsType == 'transCartesian' or ampsType == 'transPolar':
        if ampsType == 'transPolar' :
          # polar transversity amplitudes
          config.addSetting('A0Mag2', RooRealSetting('A0Mag2',
              '|A0|^2', 'par', A0Sq, 0., 1.))
          config.addSetting('A0Ph', RooRealSetting('delta0',
              'delta_0', 'par', 0., '', ''))
          config.addSetting('AparMag2', RooFormSetting('AparMag2',
              '|A_par|^2', '1. - @0 - @1', ['A0Mag2', 'AperpMag2']))
          config.addSetting('AparPh', RooRealSetting('deltaPar',
              'delta_par', 'par', AparPh, -2. * pi, 2. * pi))
          config.addSetting('AperpMag2', RooRealSetting('AperpMag2',
              '|A_perp|^2', 'par', AperpSq, 0., 1.))
          config.addSetting('AperpPh', RooRealSetting('deltaPerp',
              'delta_perp', 'par', AperpPh, -2. * pi, 2. * pi))
          if KSWave[:7] == 'include' :
            config.addSetting('ASMag2', RooRealSetting('ASMag2',
                '|A_S|^2', 'par', ASSq, 0., 1.))
            config.addSetting('ASPh', RooRealSetting('deltaS',
                'delta_S', 'par', ASPh, -2. * pi, 2. * pi))

          # construct cartesian amplitudes with polar parameters
          config.addSetting('ReA0', RooFormSetting('ReA0',
              'Re(A_0)',    'sqrt(@0) * cos(@1)', ['A0Mag2', 'A0Ph']))
          config.addSetting('ImA0', RooFormSetting('ImA0',
              'Im(A_0)',    'sqrt(@0) * sin(@1)', ['A0Mag2', 'A0Ph']))
          config.addSetting('ReApar', RooFormSetting('ReApar',
              'Re(A_par)',  'sqrt(@0) * cos(@1)', ['AparMag2', 'AparPh']))
          config.addSetting('ImApar', RooFormSetting('ImApar',
              'Im(A_par)',  'sqrt(@0) * sin(@1)', ['AparMag2', 'AparPh']))
          config.addSetting('ReAperp', RooFormSetting('ReAperp',
              'Re(A_perp)', 'sqrt(@0) * cos(@1)', ['AperpMag2', 'AperpPh']))
          config.addSetting('ImAperp', RooFormSetting('ImAperp',
              'Im(A_perp)', 'sqrt(@0) * sin(@1)', ['AperpMag2', 'AperpPh']))
          if KSWave[:7] == 'include' :
            config.addSetting('ReAS', RooFormSetting('ReAS',
                'Re(A_S)', 'sqrt(@0) * cos(@1)', ['ASMag2', 'ASPh']))
            config.addSetting('ImAS', RooFormSetting('ImAS',
                'Im(A_S)', 'sqrt(@0) * sin(@1)', ['ASMag2', 'ASPh']))

        elif ampsType == 'transCartesian' :
          # cartesian transversity amplitudes
          config.addSetting('ReA0', RooRealSetting('ReA0',
              'Re(A_0)', 'par', 1., '', ''))
          config.addSetting('ImA0', RooRealSetting('ImA0',
              'Im(A_0)', 'par', 0., '', ''))
          config.addSetting('ReApar', RooRealSetting('ReApar',
              'Re(A_par)', 'par', ReApar / ReA0, -1., 1.))
          config.addSetting('ImApar', RooRealSetting('ImApar',
              'Im(A_par)', 'par', ImApar / ReA0, -1., 1.))
          config.addSetting('ReAperp', RooRealSetting('ReAperp',
              'Re(A_perp)', 'par', ReAperp / ReA0, -1., 1.))
          config.addSetting('ImAperp', RooRealSetting('ImAperp',
              'Im(A_perp)', 'par', ImAperp / ReA0, -1., 1.))
          if KSWave[:7] == 'include' :
            config.addSetting('ReAS', RooRealSetting('ReAS',
                'Re(A_S)', 'par', ReAS / ReA0, -1., 1.))
            config.addSetting('ImAS', RooRealSetting('ImAS',
                'Im(A_S)', 'par', ImAS / ReA0, -1., 1.))

        # define angular coefficients with transversity amplitudes
        if mode == 'Bd2JpsiKstar' :
          # Bd0->J/psiK* coefficients
          config.addSetting('J_0020x0020_0', RooFormSetting('J_0020x0020_0',
              'J_0020x0020_0 coefficient', '@0 * @0 + @1 * @1',
              ['ReA0', 'ImA0']))                            # |A_0|^2
          config.addSetting('J_22x002022_0', RooFormSetting('J_22x002022_0',
              'J_22x002022_0 coefficient', '@0 * @0 + @1 * @1',
              ['ReApar', 'ImApar']))                        # |A_par|^2
          config.addSetting('J_22x002022_1', RooFormSetting('J_22x002022_1',
              'J_22x002022_1 coefficient', '@0 * @0 + @1 * @1',
              ['ReAperp', 'ImAperp']))                      # |A_perp|^2
          config.addSetting('J_21x21_0', RooFormSetting('J_21x21_0',
              'J_21x21_0 coefficient', '@0 * @2 + @1 * @3',
              ['ReA0', 'ImA0', 'ReApar', 'ImApar']))        # Re(A_0* A_par)
          config.addSetting('J_21x2m1_0', RooFormSetting('J_21x2m1_0',
              'J_21x2m1_0 coefficient', '@0 * @3 - @1 * @2',
              ['ReA0', 'ImA0', 'ReAperp', 'ImAperp']))      # Im(A_0* A_perp)
          config.addSetting('J_22x2m2_0', RooFormSetting('J_22x2m2_0',
              'J_22x2m2_0 coefficient', '@0 * @3 - @1 * @2',
              ['ReApar', 'ImApar', 'ReAperp', 'ImAperp']))  # Im(A_par* A_perp)
          if KSWave[:7] == 'include' :
            config.addSetting('J_00x0020_0', RooFormSetting('J_00x0020_0',
                'J_00x0020_0 coefficient', '@0 * @0 + @1 * @1',
                ['ReAS', 'ImAS']))                          # |A_S|^2
            config.addSetting('J_10x0020_0', RooFormSetting('J_10x0020_0',
                'J_10x0020_0 coefficient', '@0 * @2 + @1 * @3',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS']))          # Re(A_0* A_S)
            config.addSetting('J_11x21_0', RooFormSetting('J_11x21_0',
                'J_11x21_0 coefficient', '@0 * @2 + @1 * @3',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS']))      # Re(A_par* A_S)
            config.addSetting('J_11x2m1_0', RooFormSetting('J_11x2m1_0',
                'J_11x2m1_0 coefficient', '@0 * @3 - @1 * @2',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS']))    # Im(A_perp* A_S)

        elif mode == 'Bs2Jpsiphi' :
          # Bs0->J/psiphi coefficients
          config.addSetting('J_0020x0020_0_cosh',     # +|A_0|^2 * 1
              RooFormSetting('J_0020x0020_0_cosh',
              'J_0020x0020_0 cosh coefficient',
              '@0 * @0 + @1 * @1',
              ['ReA0', 'ImA0']))
          config.addSetting('J_0020x0020_0_cos',      # +|A_0|^2 * C
              RooFormSetting('J_0020x0020_0_cos',
              'J_0020x0020_0 cos coefficient',
              '(@0 * @0 + @1 * @1) * @2',
              ['ReA0', 'ImA0', 'CCP']))
          config.addSetting('J_0020x0020_0_sinh',     # -|A_0|^2 * D
              RooFormSetting('J_0020x0020_0_sinh',
              'J_0020x0020_0 sinh coefficient',
              '-(@0 * @0 + @1 * @1) * @2',
              ['ReA0', 'ImA0', 'DCP']))
          config.addSetting('J_0020x0020_0_sin',      # -|A_0|^2 * S
              RooFormSetting('J_0020x0020_0_sin',
              'J_0020x0020_0 sin coefficient',
              '-(@0 * @0 + @1 * @1) * @2',
              ['ReA0', 'ImA0', 'SCP']))

          config.addSetting('J_22x002022_0_cosh',     # +|A_par|^2 * 1
              RooFormSetting('J_22x002022_0_cosh',
              'J_22x002022_0 cosh coefficient',
              '@0 * @0 + @1 * @1',
              ['ReApar', 'ImApar']))
          config.addSetting('J_22x002022_0_cos',      # +|A_par|^2 * C
              RooFormSetting('J_22x002022_0_cos',
              'J_22x002022_0 cos coefficient',
              '(@0 * @0 + @1 * @1) * @2',
              ['ReApar', 'ImApar', 'CCP']))
          config.addSetting('J_22x002022_0_sinh',     # -|A_par|^2 * D
              RooFormSetting('J_22x002022_0_sinh',
              'J_22x002022_0 sinh coefficient',
              '-(@0 * @0 + @1 * @1) * @2',
              ['ReApar', 'ImApar', 'DCP']))
          config.addSetting('J_22x002022_0_sin',      # -|A_par|^2 * S
              RooFormSetting('J_22x002022_0_sin',
              'J_22x002022_0 sin coefficient',
              '-(@0 * @0 + @1 * @1) * @2',
              ['ReApar', 'ImApar', 'SCP']))

          config.addSetting('J_22x002022_1_cosh',     # +|A_perp|^2 * 1
              RooFormSetting('J_22x002022_1_cosh',
              'J_22x002022_1 cosh coefficient',
              '@0 * @0 + @1 * @1',
              ['ReAperp', 'ImAperp']))
          config.addSetting('J_22x002022_1_cos',      # +|A_perp|^2 * C
              RooFormSetting('J_22x002022_1_cos',
              'J_22x002022_1 cos coefficient',
              '(@0 * @0 + @1 * @1) * @2',
              ['ReAperp', 'ImAperp', 'CCP']))
          config.addSetting('J_22x002022_1_sinh',     # +|A_perp|^2 * D
              RooFormSetting('J_22x002022_1_sinh',
              'J_22x002022_1 sinh coefficient',
              '(@0 * @0 + @1 * @1) * @2',
              ['ReAperp', 'ImAperp', 'DCP']))
          config.addSetting('J_22x002022_1_sin',      # +|A_perp|^2 * S
              RooFormSetting('J_22x002022_1_sin',
              'J_22x002022_1 sin coefficient',
              '(@0 * @0 + @1 * @1) * @2',
              ['ReAperp', 'ImAperp', 'SCP']))

          config.addSetting('J_21x21_0_cosh',         # +Re(A_0* A_par) * 1
              RooFormSetting('J_21x21_0_cosh',
              'J_21x21_0 cosh coefficient',
              '@0 * @2 + @1 * @3',
              ['ReA0', 'ImA0', 'ReApar', 'ImApar']))
          config.addSetting('J_21x21_0_cos',          # +Re(A_0* A_par) * C
              RooFormSetting('J_21x21_0_cos',
              'J_21x21_0 cos coefficient',
              '(@0 * @2 + @1 * @3) * @4',
              ['ReA0', 'ImA0', 'ReApar', 'ImApar', 'CCP']))
          config.addSetting('J_21x21_0_sinh',         # -Re(A_0* A_par) * D
              RooFormSetting('J_21x21_0_sinh',
              'J_21x21_0 sinh coefficient',
              '-(@0 * @2 + @1 * @3) * @4',
              ['ReA0', 'ImA0', 'ReApar', 'ImApar', 'DCP']))
          config.addSetting('J_21x21_0_sin',          # -Re(A_0* A_par) * S
              RooFormSetting('J_21x21_0_sin',
              'J_21x21_0 sin coefficient',
              '-(@0 * @2 + @1 * @3) * @4',
              ['ReA0', 'ImA0', 'ReApar', 'ImApar', 'SCP']))

          config.addSetting('J_21x2m1_0_cosh',        # +Im(A_0* A_perp) * C
              RooFormSetting('J_21x2m1_0_cosh',
              'J_21x2m1_0 cosh coefficient',
              '(@0 * @3 - @1 * @2) * @4',
              ['ReA0', 'ImA0', 'ReAperp', 'ImAperp', 'CCP']))
          config.addSetting('J_21x2m1_0_cos',         # +Im(A_0* A_perp) * 1
              RooFormSetting('J_21x2m1_0_cos',
              'J_21x2m1_0 cos coefficient',
              '@0 * @3 - @1 * @2',
              ['ReA0', 'ImA0', 'ReAperp', 'ImAperp']))
          config.addSetting('J_21x2m1_0_sinh',        # +Re(A_0* A_perp) * S
              RooFormSetting('J_21x2m1_0_sinh',
              'J_21x2m1_0 sinh coefficient',
              '(@0 * @2 + @1 * @3) * @4',
              ['ReA0', 'ImA0', 'ReAperp', 'ImAperp', 'SCP']))
          config.addSetting('J_21x2m1_0_sin',         # -Re(A_0* A_perp) * D
              RooFormSetting('J_21x2m1_0_sin',
              'J_21x2m1_0 sin coefficient',
              '-(@0 * @2 + @1 * @3) * @4',
              ['ReA0', 'ImA0', 'ReAperp', 'ImAperp', 'DCP']))

          config.addSetting('J_22x2m2_0_cosh',        # +Im(A_par* A_perp) * C
              RooFormSetting('J_22x2m2_0_cosh',
              'J_22x2m2_0 cosh coefficient',
              '(@0 * @3 - @1 * @2) * @4',
              ['ReApar', 'ImApar', 'ReAperp', 'ImAperp', 'CCP']))
          config.addSetting('J_22x2m2_0_cos',         # +Im(A_par* A_perp) * 1
              RooFormSetting('J_22x2m2_0_cos',
              'J_22x2m2_0 cos coefficient',
              '@0 * @3 - @1 * @2',
              ['ReApar', 'ImApar', 'ReAperp', 'ImAperp']))
          config.addSetting('J_22x2m2_0_sinh',        # +Re(A_par* A_perp) * S
              RooFormSetting('J_22x2m2_0_sinh',
              'J_22x2m2_0 sinh coefficient',
              '(@0 * @2 + @1 * @3) * @4',
              ['ReApar', 'ImApar', 'ReAperp', 'ImAperp', 'SCP']))
          config.addSetting('J_22x2m2_0_sin',         # -Re(A_par* A_perp) * D
              RooFormSetting('J_22x2m2_0_sin',
              'J_22x2m2_0 sin coefficient',
              '-(@0 * @2 + @1 * @3) * @4',
              ['ReApar', 'ImApar', 'ReAperp', 'ImAperp', 'DCP']))

          if KSWave == 'includeEven' :
            config.addSetting('J_00x0020_0_cosh',     # +|A_S|^2 * 1
                RooFormSetting('J_00x0020_0_cosh',
                'J_00x0020_0 cosh coefficient',
                '@0 * @0 + @1 * @1',
                ['ReAS', 'ImAS']))
            config.addSetting('J_00x0020_0_cos',      # +|A_S|^2 * C
                RooFormSetting('J_00x0020_0_cos',
                'J_00x0020_0 cos coefficient',
                '(@0 * @0 + @1 * @1) * @2',
                ['ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_00x0020_0_sinh',     # -|A_S|^2 * D
                RooFormSetting('J_00x0020_0_sinh',
                'J_00x0020_0 sinh coefficient',
                '-(@0 * @0 + @1 * @1) * @2',
                ['ReAS', 'ImAS', 'DCP']))
            config.addSetting('J_00x0020_0_sin',      # -|A_S|^2 * S
                RooFormSetting('J_00x0020_0_sin',
                'J_00x0020_0 sin coefficient',
                '-(@0 * @0 + @1 * @1) * @2',
                ['ReAS', 'ImAS', 'SCP']))

            config.addSetting('J_10x0020_0_cosh',     # +Re(A_0* A_S) * 1
                RooFormSetting('J_10x0020_0_cosh',
                'J_10x0020_0 cosh coefficient',
                '@0 * @2 + @1 * @3',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS']))
            config.addSetting('J_10x0020_0_cos',      # +Re(A_0* A_S) * C
                RooFormSetting('J_10x0020_0_cos',
                'J_10x0020_0 cos coefficient',
                '(@0 * @2 + @1 * @3) * @4',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_10x0020_0_sinh',     # -Re(A_0* A_S) * D
                RooFormSetting('J_10x0020_0_sinh',
                'J_10x0020_0 sinh coefficient',
                '-(@0 * @2 + @1 * @3) * @4',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS', 'DCP']))
            config.addSetting('J_10x0020_0_sin',      # -Re(A_0* A_S) * S
                RooFormSetting('J_10x0020_0_sin',
                'J_10x0020_0 sin coefficient',
                '-(@0 * @2 + @1 * @3) * @4',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS', 'SCP']))

            config.addSetting('J_11x21_0_cosh',       # +Re(A_par* A_S) * 1
                RooFormSetting('J_11x21_0_cosh',
                'J_11x21_0 cosh coefficient',
                '@0 * @2 + @1 * @3',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS']))
            config.addSetting('J_11x21_0_cos',        # +Re(A_par* A_S) * C
                RooFormSetting('J_11x21_0_cos',
                'J_11x21_0 cos coefficient',
                '(@0 * @2 + @1 * @3) * @4',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_11x21_0_sinh',       # -Re(A_par* A_S) * D
                RooFormSetting('J_11x21_0_sinh',
                'J_11x21_0 sinh coefficient',
                '-(@0 * @2 + @1 * @3) * @4',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS', 'DCP']))
            config.addSetting('J_11x21_0_sin',        # -Re(A_par* A_S) * S
                RooFormSetting('J_11x21_0_sin',
                'J_11x21_0 sin coefficient',
                '-(@0 * @2 + @1 * @3) * @4',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS', 'SCP']))

            config.addSetting('J_11x2m1_0_cosh',      # +Im(A_perp* A_S) * C
                RooFormSetting('J_11x2m1_0_cosh',
                'J_11x2m1_0 cosh coefficient',
                '(@0 * @3 - @1 * @2) * @4',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_11x2m1_0_cos',       # +Im(A_perp* A_S) * 1
                RooFormSetting('J_11x2m1_0_cos',
                'J_11x2m1_0 cos coefficient',
                '@0 * @3 - @1 * @2',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS']))
            config.addSetting('J_11x2m1_0_sinh',      # -Re(A_perp* A_S) * S
                RooFormSetting('J_11x2m1_0_sinh',
                'J_11x2m1_0 sinh coefficient',
                '-(@0 * @2 + @1 * @3) * @4',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS', 'SCP']))
            config.addSetting('J_11x2m1_0_sin',       # +Re(A_perp* A_S) * D
                RooFormSetting('J_11x2m1_0_sin',
                'J_11x2m1_0 sin coefficient',
                '(@0 * @2 + @1 * @3) * @4',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS', 'DCP']))

          elif KSWave == 'includeOdd' :
            config.addSetting('J_00x0020_0_cosh',     # +|A_S|^2 * 1
                RooFormSetting('J_00x0020_0_cosh',
                'J_00x0020_0 cosh coefficient',
                '@0 * @0 + @1 * @1',
                ['ReAS', 'ImAS']))
            config.addSetting('J_00x0020_0_cos',      # +|A_S|^2 * C
                RooFormSetting('J_00x0020_0_cos',
                'J_00x0020_0 cos coefficient',
                '(@0 * @0 + @1 * @1) * @2',
                ['ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_00x0020_0_sinh',     # +|A_S|^2 * D
                RooFormSetting('J_00x0020_0_sinh',
                'J_00x0020_0 sinh coefficient',
                '(@0 * @0 + @1 * @1) * @2',
                ['ReAS', 'ImAS', 'DCP']))
            config.addSetting('J_00x0020_0_sin',      # +|A_S|^2 * S
                RooFormSetting('J_00x0020_0_sin',
                'J_00x0020_0 sin coefficient',
                '(@0 * @0 + @1 * @1) * @2',
                ['ReAS', 'ImAS', 'SCP']))

            config.addSetting('J_10x0020_0_cosh',
                RooFormSetting('J_10x0020_0_cosh',    # +Re(A_0* A_S) * C
                'J_10x0020_0 cosh coefficient',
                '(@0 * @2 + @1 * @3) * @4',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_10x0020_0_cos',      # +Re(A_0* A_S) * 1
                RooFormSetting('J_10x0020_0_cos',
                'J_10x0020_0 cos coefficient',
                '@0 * @2 + @1 * @3',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS']))
            config.addSetting('J_10x0020_0_sinh',     # -Im(A_0* A_S) * S
                RooFormSetting('J_10x0020_0_sinh',
                'J_10x0020_0 sinh coefficient',
                '-(@0 * @3 - @1 * @2) * @4',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS', 'SCP']))
            config.addSetting('J_10x0020_0_sin',      # +Im(A_0* A_S) * D
                RooFormSetting('J_10x0020_0_sin',
                'J_10x0020_0 sin coefficient',
                '(@0 * @3 - @1 * @2) * @4',
                ['ReA0', 'ImA0', 'ReAS', 'ImAS', 'DCP']))

            config.addSetting('J_11x21_0_cosh',       # +Re(A_par* A_S) * C
                RooFormSetting('J_11x21_0_cosh',
                'J_11x21_0 cosh coefficient',
                '(@0 * @2 + @1 * @3) * @4',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_11x21_0_cos',        # +Re(A_par* A_S) * 1
                RooFormSetting('J_11x21_0_cos',
                'J_11x21_0 cos coefficient',
                '@0 * @2 + @1 * @3',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS']))
            config.addSetting('J_11x21_0_sinh',       # -Im(A_par* A_S) * S
                RooFormSetting('J_11x21_0_sinh',
                'J_11x21_0 sinh coefficient',
                '-(@0 * @3 - @1 * @2) * @4',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS', 'SCP']))
            config.addSetting('J_11x21_0_sin',        # +Im(A_par* A_S) * D
                RooFormSetting('J_11x21_0_sin',
                'J_11x21_0 sin coefficient',
                '(@0 * @3 - @1 * @2) * @4',
                ['ReApar', 'ImApar', 'ReAS', 'ImAS', 'DCP']))

            config.addSetting('J_11x2m1_0_cosh',      # +Im(A_perp* A_S) * 1
                RooFormSetting('J_11x2m1_0_cosh',
                'J_11x2m1_0 cosh coefficient',
                '@0 * @3 - @1 * @2',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS']))
            config.addSetting('J_11x2m1_0_cos',       # +Im(A_perp* A_S) * C
                RooFormSetting('J_11x2m1_0_cos',
                'J_11x2m1_0 cos coefficient',
                '(@0 * @3 - @1 * @2) * @4',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS', 'CCP']))
            config.addSetting('J_11x2m1_0_sinh',      # +Im(A_perp* A_S) * D
                RooFormSetting('J_11x2m1_0_sinh',
                'J_11x2m1_0 sinh coefficient',
                '(@0 * @3 - @1 * @2) * @4',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS', 'DCP']))
            config.addSetting('J_11x2m1_0_sin',       # +Im(A_perp* A_S) * S
                RooFormSetting('J_11x2m1_0_sin',
                'J_11x2m1_0 sin coefficient',
                '(@0 * @3 - @1 * @2) * @4',
                ['ReAperp', 'ImAperp', 'ReAS', 'ImAS', 'SCP']))

      else :
        # angular coefficients
        if mode == 'Bd2JpsiKstar' :
          # Bd0->J/psiK* coefficients
          config.addSetting('J_0020x0020_0', RooRealSetting('J_0020x0020_0',
              'J_0020x0020_0 coefficient', 'par',           # |A_0|^2
              1., '', ''))
          config.addSetting('J_22x002022_0', RooRealSetting('J_22x002022_0',
              'J_22x002022_0 coefficient', 'par',           # |A_par|^2
              AparSq / A0Sq, -5., 5.))
          config.addSetting('J_22x002022_1', RooRealSetting('J_22x002022_1',
              'J_22x002022_1 coefficient', 'par',           # |A_perp|^2
              AperpSq / A0Sq, -5., 5.))
          config.addSetting('J_21x21_0', RooRealSetting('J_21x21_0',
              'J_21x21_0 coefficient', 'par',               # Re(A_0* A_par)
              ReApar / ReA0, -5., 5.))
          config.addSetting('J_21x2m1_0', RooRealSetting('J_21x2m1_0',
              'J_21x2m1_0 coefficient', 'par',              # Im(A_0* A_perp)
              ImAperp / ReA0, -5., 5.))
          config.addSetting('J_22x2m2_0', RooRealSetting('J_22x2m2_0',
              'J_22x2m2_0 coefficient', 'par',              # Im(A_par* A_perp)
              (ReApar * ImAperp - ImApar * ReAperp) / A0Sq, -5., 5.))
          if KSWave[:7] == 'include' :
            config.addSetting('J_00x0020_0', RooRealSetting('J_00x0020_0',
                'J_00x0020_0 coefficient', 'par',           # |A_S|^2
                ASSq / A0Sq, -5., 5.))
            config.addSetting('J_10x0020_0', RooRealSetting('J_10x0020_0',
                'J_10x0020_0 coefficient', 'par',           # Re(A_0* A_S)
                ReAS / ReA0, -5., 5.))
            config.addSetting('J_11x21_0', RooRealSetting('J_11x21_0',
                'J_11x21_0 coefficient', 'par',             # Re(A_par* A_S)
                (ReApar * ReAS + ImApar * ImAS) / A0Sq, -5., 5.))
            config.addSetting('J_11x2m1_0', RooRealSetting('J_11x2m1_0',
                'J_11x2m1_0 coefficient', 'par',            # Im(A_perp* A_S)
                (ReAperp * ImAS - ImAperp * ReAS) / A0Sq, -5., 5.))

        elif mode == 'Bs2Jpsiphi' :
          # Bs0->J/psiphi coefficients
          config.addSetting('J_0020x0020_0_cosh',     # +|A_0|^2 * 1
              RooRealSetting('J_0020x0020_0_cosh',
              'J_0020x0020_0 cosh coefficient', 'par',
              1., '', ''))
          config.addSetting('J_0020x0020_0_cos',      # +|A_0|^2 * C
              RooRealSetting('J_0020x0020_0_cos',
              'J_0020x0020_0 cos coefficient', 'par',
              CCP, -5., 5.))
          config.addSetting('J_0020x0020_0_sinh',     # -|A_0|^2 * D
              RooRealSetting('J_0020x0020_0_sinh',
              'J_0020x0020_0 sinh coefficient', 'par',
              -DCP, -5., 5.))
          config.addSetting('J_0020x0020_0_sin',      # -|A_0|^2 * S
              RooRealSetting('J_0020x0020_0_sin',
              'J_0020x0020_0 sin coefficient', 'par',
              -SCP, -5., 5.))

          config.addSetting('J_22x002022_0_cosh',     # +|A_par|^2 * 1
              RooRealSetting('J_22x002022_0_cosh',
              'J_22x002022_0 cosh coefficient', 'par',
              AparSq / A0Sq, -5., 5.))
          config.addSetting('J_22x002022_0_cos',      # +|A_par|^2 * C
              RooRealSetting('J_22x002022_0_cos',
              'J_22x002022_0 cos coefficient', 'par',
              AparSq / A0Sq * CCP, -5., 5.))
          config.addSetting('J_22x002022_0_sinh',     # -|A_par|^2 * D
              RooRealSetting('J_22x002022_0_sinh',
              'J_22x002022_0 sinh coefficient', 'par',
              -AparSq / A0Sq * DCP, -5., 5.))
          config.addSetting('J_22x002022_0_sin',      # -|A_par|^2 * S
              RooRealSetting('J_22x002022_0_sin',
              'J_22x002022_0 sin coefficient', 'par',
              -AparSq / A0Sq * SCP, -5., 5.))

          config.addSetting('J_22x002022_1_cosh',     # +|A_perp|^2 * 1
              RooRealSetting('J_22x002022_1_cosh',
              'J_22x002022_1 cosh coefficient', 'par',
              AperpSq / A0Sq, -5., 5.))
          config.addSetting('J_22x002022_1_cos',      # +|A_perp|^2 * C
              RooRealSetting('J_22x002022_1_cos',
              'J_22x002022_1 cos coefficient', 'par',
              AperpSq / A0Sq * CCP, -5., 5.))
          config.addSetting('J_22x002022_1_sinh',     # +|A_perp|^2 * D
              RooRealSetting('J_22x002022_1_sinh',
              'J_22x002022_1 sinh coefficient', 'par',
              AperpSq / A0Sq * DCP, -5., 5.))
          config.addSetting('J_22x002022_1_sin',      # +|A_perp|^2 * S
              RooRealSetting('J_22x002022_1_sin',
              'J_22x002022_1 sin coefficient', 'par',
              AperpSq / A0Sq * SCP, -5., 5.))

          config.addSetting('J_21x21_0_cosh',         # +Re(A_0* A_par) * 1
              RooRealSetting('J_21x21_0_cosh',
              'J_21x21_0 cosh coefficient', 'par',
              ReApar / ReA0, -5., 5.))
          config.addSetting('J_21x21_0_cos',          # +Re(A_0* A_par) * C
              RooRealSetting('J_21x21_0_cos',
              'J_21x21_0 cos coefficient', 'par',
              ReApar / ReA0 * CCP, -5., 5.))
          config.addSetting('J_21x21_0_sinh',         # -Re(A_0* A_par) * D
              RooRealSetting('J_21x21_0_sinh',
              'J_21x21_0 sinh coefficient', 'par',
              -ReApar / ReA0 * DCP, -5., 5.))
          config.addSetting('J_21x21_0_sin',          # -Re(A_0* A_par) * S
              RooRealSetting('J_21x21_0_sin',
              'J_21x21_0 sin coefficient', 'par',
              -ReApar / ReA0 * SCP, -5., 5.))

          config.addSetting('J_21x2m1_0_cosh',        # +Im(A_0* A_perp) * C
              RooRealSetting('J_21x2m1_0_cosh',
              'J_21x2m1_0 cosh coefficient', 'par',
              ImAperp / ReA0 * CCP, -5., 5.))
          config.addSetting('J_21x2m1_0_cos',         # +Im(A_0* A_perp) * 1
              RooRealSetting('J_21x2m1_0_cos',
              'J_21x2m1_0 cos coefficient', 'par',
              ImAperp / ReA0, -5., 5.))
          config.addSetting('J_21x2m1_0_sinh',        # +Re(A_0* A_perp) * S
              RooRealSetting('J_21x2m1_0_sinh',
              'J_21x2m1_0 sinh coefficient', 'par',
              ReAperp / ReA0 * SCP, -5., 5.))
          config.addSetting('J_21x2m1_0_sin',         # -Re(A_0* A_perp) * D
              RooRealSetting('J_21x2m1_0_sin',
              'J_21x2m1_0 sin coefficient', 'par',
              -ReAperp / ReA0 * DCP, -5., 5.))

          config.addSetting('J_22x2m2_0_cosh',        # +Im(A_par* A_perp) * C
              RooRealSetting('J_22x2m2_0_cosh',
              'J_22x2m2_0 cosh coefficient', 'par',
              (ReApar * ImAperp - ImApar * ReAperp) / A0Sq * CCP, -5., 5.))
          config.addSetting('J_22x2m2_0_cos',         # +Im(A_par* A_perp) * 1
              RooRealSetting('J_22x2m2_0_cos',
              'J_22x2m2_0 cos coefficient', 'par',
              (ReApar * ImAperp - ImApar * ReAperp) / A0Sq, -5., 5.))
          config.addSetting('J_22x2m2_0_sinh',        # +Re(A_par* A_perp) * S
              RooRealSetting('J_22x2m2_0_sinh',
              'J_22x2m2_0 sinh coefficient', 'par',
              (ReApar * ReAperp + ImApar * ImAperp) / A0Sq * SCP, -5., 5.))
          config.addSetting('J_22x2m2_0_sin',         # -Re(A_par* A_perp) * D
              RooRealSetting('J_22x2m2_0_sin',
              'J_22x2m2_0 sin coefficient', 'par',
              -(ReApar * ReAperp + ImApar * ImAperp) / A0Sq * DCP, -5., 5.))

          if KSWave == 'includeEven' :
            config.addSetting('J_00x0020_0_cosh',     # +|A_S|^2 * 1
                RooRealSetting('J_00x0020_0_cosh',
                'J_00x0020_0 cosh coefficient', 'par',
                ASSq / A0Sq, -5., 5.))
            config.addSetting('J_00x0020_0_cos',      # +|A_S|^2 * C
                RooRealSetting('J_00x0020_0_cos',
                'J_00x0020_0 cos coefficient', 'par',
                ASSq / A0Sq * CCP, -5., 5.))
            config.addSetting('J_00x0020_0_sinh',     # -|A_S|^2 * D
                RooRealSetting('J_00x0020_0_sinh',
                'J_00x0020_0 sinh coefficient', 'par',
                -ASSq / A0Sq * DCP, -5., 5.))
            config.addSetting('J_00x0020_0_sin',      # -|A_S|^2 * S
                RooRealSetting('J_00x0020_0_sin',
                'J_00x0020_0 sin coefficient', 'par',
                -ASSq / A0Sq * SCP, -5., 5.))

            config.addSetting('J_10x0020_0_cosh',     # +Re(A_0* A_S) * 1
                RooRealSetting('J_10x0020_0_cosh',
                'J_10x0020_0 cosh coefficient', 'par',
                ReAS / ReA0, -5., 5.))
            config.addSetting('J_10x0020_0_cos',      # +Re(A_0* A_S) * C
                RooRealSetting('J_10x0020_0_cos',
                'J_10x0020_0 cos coefficient', 'par',
                ReAS / ReA0 * CCP, -5., 5.))
            config.addSetting('J_10x0020_0_sinh',     # -Re(A_0* A_S) * D
                RooRealSetting('J_10x0020_0_sinh',
                'J_10x0020_0 sinh coefficient', 'par',
                -ReAS / ReA0 * DCP, -5., 5.))
            config.addSetting('J_10x0020_0_sin',      # -Re(A_0* A_S) * S
                RooRealSetting('J_10x0020_0_sin',
                'J_10x0020_0 sin coefficient', 'par',
                -ReAS / ReA0 * SCP, -5., 5.))

            config.addSetting('J_11x21_0_cosh',       # +Re(A_par* A_S) * 1
                RooRealSetting('J_11x21_0_cosh',
                'J_11x21_0 cosh coefficient', 'par',
                (ReApar * ReAS + ImApar * ImAS) / A0Sq, -5., 5.))
            config.addSetting('J_11x21_0_cos',        # +Re(A_par* A_S) * C
                RooRealSetting('J_11x21_0_cos',
                'J_11x21_0 cos coefficient', 'par',
                (ReApar * ReAS + ImApar * ImAS) / A0Sq * CCP, -5., 5.))
            config.addSetting('J_11x21_0_sinh',       # -Re(A_par* A_S) * D
                RooRealSetting('J_11x21_0_sinh',
                'J_11x21_0 sinh coefficient', 'par',
                -(ReApar * ReAS + ImApar * ImAS) / A0Sq * DCP, -5., 5.))
            config.addSetting('J_11x21_0_sin',        # -Re(A_par* A_S) * S
                RooRealSetting('J_11x21_0_sin',
                'J_11x21_0 sin coefficient', 'par',
                -(ReApar * ReAS + ImApar * ImAS) / A0Sq * SCP, -5., 5.))

            config.addSetting('J_11x2m1_0_cosh',      # +Im(A_perp* A_S) * C
                RooRealSetting('J_11x2m1_0_cosh',
                'J_11x2m1_0 cosh coefficient', 'par',
                (ReAperp * ImAS - ImAperp * ReAS) / A0Sq * CCP, -5., 5.))
            config.addSetting('J_11x2m1_0_cos',       # +Im(A_perp* A_S) * 1
                RooRealSetting('J_11x2m1_0_cos',
                'J_11x2m1_0 cos coefficient', 'par',
                (ReAperp * ImAS - ImAperp * ReAS) / A0Sq, -5., 5.))
            config.addSetting('J_11x2m1_0_sinh',      # -Re(A_perp* A_S) * S
                RooRealSetting('J_11x2m1_0_sinh',
                'J_11x2m1_0 sinh coefficient', 'par',
                -(ReAperp * ReAS + ImAperp * ImAS) / A0Sq *SCP, -5., 5.))
            config.addSetting('J_11x2m1_0_sin',       # +Re(A_perp* A_S) * D
                RooRealSetting('J_11x2m1_0_sin',
                'J_11x2m1_0 sin coefficient', 'par',
                (ReAperp * ReAS + ImAperp * ImAS) / A0Sq * DCP, -5., 5.))

          elif KSWave == 'includeOdd' :
            config.addSetting('J_00x0020_0_cosh',     # +|A_S|^2 * 1
                RooRealSetting('J_00x0020_0_cosh',
                'J_00x0020_0 cosh coefficient', 'par',
                ASSq / A0Sq, -5., 5.))
            config.addSetting('J_00x0020_0_cos',      # +|A_S|^2 * C
                RooRealSetting('J_00x0020_0_cos',
                'J_00x0020_0 cos coefficient', 'par',
                ASSq / A0Sq * CCP, -5., 5.))
            config.addSetting('J_00x0020_0_sinh',     # +|A_S|^2 * D
                RooRealSetting('J_00x0020_0_sinh',
                'J_00x0020_0 sinh coefficient', 'par',
                ASSq / A0Sq * DCP, -5., 5.))
            config.addSetting('J_00x0020_0_sin',      # +|A_S|^2 * S
                RooRealSetting('J_00x0020_0_sin',
                'J_00x0020_0 sin coefficient', 'par',
                ASSq / A0Sq * SCP, -5., 5.))

            config.addSetting('J_10x0020_0_cosh',     # +Re(A_0* A_S) * C
                RooRealSetting('J_10x0020_0_cosh',
                'J_10x0020_0 cosh coefficient', 'par',
                ReAS / ReA0 * CCP, -5., 5.))
            config.addSetting('J_10x0020_0_cos',      # +Re(A_0* A_S) * 1
                RooRealSetting('J_10x0020_0_cos',
                'J_10x0020_0 cos coefficient', 'par',
                ReAS / ReA0, -5., 5.))
            config.addSetting('J_10x0020_0_sinh',     # -Im(A_0* A_S) * S
                RooRealSetting('J_10x0020_0_sinh',
                'J_10x0020_0 sinh coefficient', 'par',
                -ImAS / ReA0 * SCP, -5., 5.))
            config.addSetting('J_10x0020_0_sin',      # +Im(A_0* A_S) * D
                RooRealSetting('J_10x0020_0_sin',
                'J_10x0020_0 sin coefficient', 'par',
                ImAS / ReA0 * DCP, -5., 5.))

            config.addSetting('J_11x21_0_cosh',       # +Re(A_par* A_S) * C
                RooRealSetting('J_11x21_0_cosh',
                'J_11x21_0 cosh coefficient', 'par',
                (ReApar * ReAS + ImApar * ImAS) / A0Sq * CCP, -5., 5.))
            config.addSetting('J_11x21_0_cos',        # +Re(A_par* A_S) * 1
                RooRealSetting('J_11x21_0_cos',
                'J_11x21_0 cos coefficient', 'par',
                (ReApar * ReAS + ImApar * ImAS) / A0Sq, -5., 5.))
            config.addSetting('J_11x21_0_sinh',       # -Im(A_par* A_S) * S
                RooRealSetting('J_11x21_0_sinh',
                'J_11x21_0 sinh coefficient', 'par',
                -(ReApar * ImAS - ImApar * ReAS) / A0Sq * SCP, -5., 5.))
            config.addSetting('J_11x21_0_sin',        # +Im(A_par* A_S) * D
                RooRealSetting('J_11x21_0_sin',
                'J_11x21_0 sin coefficient', 'par',
                (ReApar * ImAS - ImApar * ReAS) / A0Sq * DCP, -5., 5.))

            config.addSetting('J_11x2m1_0_cosh',      # +Im(A_perp* A_S) * 1
                RooRealSetting('J_11x2m1_0_cosh',
                'J_11x2m1_0 cosh coefficient', 'par',
                (ReAperp * ImAS - ImAperp * ReAS) / A0Sq, -5., 5.))
            config.addSetting('J_11x2m1_0_cos',       # +Im(A_perp* A_S) * C
                RooRealSetting('J_11x2m1_0_cos',
                'J_11x2m1_0 cos coefficient', 'par',
                (ReAperp * ImAS - ImAperp * ReAS) / A0Sq * CCP, -5., 5.))
            config.addSetting('J_11x2m1_0_sinh',      # +Im(A_perp* A_S) * D
                RooRealSetting('J_11x2m1_0_sinh',
                'J_11x2m1_0 sinh coefficient', 'par',
                (ReAperp * ImAS - ImAperp * ReAS) / A0Sq * DCP, -5., 5.))
            config.addSetting('J_11x2m1_0_sin',       # +Im(A_perp* A_S) * S
                RooRealSetting('J_11x2m1_0_sin',
                'J_11x2m1_0 sin coefficient', 'par',
                (ReAperp * ImAS - ImAperp * ReAS) / A0Sq * SCP, -5., 5.))


      ### efficiencies ###

      effType = 'noEff'
      if 'effType' in optDict and optDict['effType'] != '' :
        effType = optDict['effType']
      config.addSetting('effType', P2VVSetting('effType',
        'type of efficiency', effType))

      if effType == 'noEff' :
        config.addSetting('angEffBasisFuncs', P2VVSetting('angEffBasisFuncs',
          'angular efficiency basis functions', ''))
      elif effType == 'angular' :
        config.addSetting('angEffBasisFuncs', P2VVSetting('angEffBasisFuncs',
          'angular efficiency basis functions', (4, 4)))


  return config


###############################################################################

class P2VVConfiguration :
  """class for configuration of the P2VV framework

  contains a dictionary with P2VV settings and a RooFit workspace
  """

  def __init__(self, workspace = '', WSFilePath = '') :
    from ROOT import RooWorkspace

    # settings dictionary
    self._settingsDict = {}

    # RooFit workspace
    if type(workspace) is str or type(workspace) is RooWorkspace :
      self.setWorkspace(workspace)
    else :
      self._workspace = ''

    # workspace file path
    if type(WSFilePath) is str and len(WSFilePath) > 0 :
      self.setWSFilePath(WSFilePath)

    # model builders
    self._modelBuilders = {}

    # plots stash
    self._plotsStash = []


  ## settings dictionary methods ##

  def __contains__(self, key) :
    return key in self._settingsDict

  def __getitem__(self, key) :
    if key in self._settingsDict :
      return self._settingsDict[key]
    else :
      return None

  def __iter__(self) :
    for key in self._settingsDict :
      yield key

  def iteritems(self) :
    for keyValue in self._settingsDict.iteritems() :
      yield keyValue

  def iterkeys(self) :
    for key in self._settingsDict.iterkeys() :
      yield key

  def itervalues(self) :
    for value in self._settingsDict.itervalues() :
      yield value

  def __len__(self) :
    return len(self._settingsDict)


  ## P2VV configuration methods ##

  def value(self, key) :
    """get the value of a setting
    """

    # check type of the setting key
    if type(key) is not str :
      print "P2VV - ERROR: P2VVConfiguration.value: setting key is not of type 'str'"
      return None

    # get value of setting
    setting = self[key]
    if setting : return setting.value()
    else : return None    

  def set(self, key, setting) :
    """set a setting
    """
    self.addSetting(key, setting, True)

  def addSetting(self, key, setting, overWriteOld = True) :
    """add a setting to the settings dictionary or overwrite an old setting
    """

    # check type of the setting key
    if type(key) is not str :
      print "P2VV - ERROR: P2VVConfiguration.addSetting: setting key is not a string"
      return

    # check type of the setting
    if not issubclass(setting.__class__, P2VVSetting) :
      print "P2VV - ERROR: P2VVConfiguration.addSetting: class of 'setting' does not inherit from 'P2VVSetting'"
      return

    # add new setting to settings dictionary or overwrite old setting
    if not key in self._settingsDict or overWriteOld:
      self._settingsDict[key] = setting

  def removeSetting(self, key) :
    """remove a setting from the settings dictionary
    """

    # check type of the setting key
    if type(key) is not str :
      print "P2VV - ERROR: P2VVConfiguration.removeSetting: setting key is not a string"
      return

    # check if setting is in dictionary
    if key not in self._settingsDict :
      print "P2VV - ERROR: P2VVConfiguration.removeSetting: '%s': no sucth setting"\
          % key
      return

    del self._settingsDict[key]


  ## RooFit workspace methods ##

  def workspace(self) :
    """get workspace
    """
    return self._workspace

  def setWorkspace(self, workspace = '') :
    """set workspace
    """

    from ROOT import RooWorkspace

    if type(workspace) is RooWorkspace :
      self._workspace = workspace
    elif type(workspace) is str :
      if workspace == '' :
        self._workspace = RooWorkspace('P2VVRooWS')
      else :
        self._workspace = RooWorkspace(workspace)
    else :
      print "P2VV - WARNING: P2VVConfiguration.setWorkspace: argument is not a RooWorkspace: no workspace set"
      self._workspace = ''

  def writeWorkspace(self, WSPath = '') :
    """write workspace to ROOT file
    """

    if self._workspace == '' :
      print "P2VV - ERROR: P2VVConfiguration.writeWorkspace: no workspace set"
      return

    if type(WSPath) is str and len(WSPath) > 0 :
      self._workspace.writeToFile(WSPath)
    elif 'WSPath' in self._settingsDict\
        and type(self._settingsDict['WSPath'].value()) is str\
        and len(self._settingsDict['WSPath'].value()) > 0 :
      self._workspace.writeToFile(self._settingsDict['WSPath'].value())
    else :
      print "P2VV - ERROR: P2VVConfiguration.writeWorkspace: no workspace file path set"

  def WSFilePath(self) :
    """get workspace file path
    """

    if 'WSPath' in self._settingsDict :
      return self._settingsDict['WSPath'].value()
    else :
      return ''

  def setWSFilePath(self, WSPath) :
    """set workspace file path
    """

    if type(WSPath) is str and len(WSPath) > 0 :
      self.addSetting('WSPath', P2VVSetting('WSPath', 'workspace file path',
          WSPath), True)
    else :
      print "P2VV - ERROR: P2VVConfiguration.setWSFilePath: argument is not a string or an empty string: no workspace file path set"

  def declareRooVars(self, varType = 'var') :
    """declares RooFit variables

    declares the RooFit variables that are in the settings dictionary and puts
    them in workspace
    """

    if self._workspace == '' :
      print "P2VV - ERROR: P2VVConfiguration.declareRooVars: no workspace set: can't declare any variables"
      return

    # loop over settings and find RooSettings to put in workspace
    formSettings   = []
    obsSettingKeys = []
    for key, setting in self._settingsDict.iteritems() :
      if not issubclass(setting.__class__, RooSetting) : continue

      if setting.type() is 'RooFormSetting' and varType in ['var', 'form'] :
        # put RooFormSettings in a list for later declaration
        formSettings.append(setting)

      elif (setting.observable() and varType in ['var', 'obs'])\
          or (not setting.observable() and varType in ['var', 'par']) :
        # declare variable
        setting.declare(self._workspace)
        if setting.observable() : obsSettingKeys.append(key)

    # define sets of observables
    self.defineRooSet('observables', obsSettingKeys)
    self.defineRooSet('angles', ['cpsiAng', 'cthetaAng', 'phiAng'])

    # declare RooFormulaVars
    if varType in ['var', 'form'] :
      while True :
        nDeclared = 0
        for setting in formSettings :
          if setting.declared() : continue

          # check if all the needed variables have been declared
          declare = True
          for var in setting.variables() :
            if not var in self._settingsDict\
                or not self._settingsDict[var].declared() :
              declare = False
              break

          # declare formula
          if declare :
            nDeclared += 1
            setting.declare(self._workspace, self)

        # exit loop if nothing was declared any more
        if nDeclared == 0 : break

  def defineRooSet(self, name, settingKeysList) :
    """defines a set of RooFit variables in the workspace
    """

    if self._workspace == '' :
      print "P2VV - ERROR: P2VVConfiguration.defineRooSet: no workspace set: can't define set"
      return

    if type(name) is not str or len(name) < 1 :
      print "P2VV - ERROR: P2VVConfiguration.defineRooSet: 'name' is not a string or an empty string"
      return

    if type(settingKeysList) is not list :
      print "P2VV - ERROR: P2VVConfiguration.defineRooSet: 'settingsList' is not a list"
      return

    # define set of RooVars in workspace
    varString = ''
    for key in settingKeysList :
      if key in self._settingsDict :
        varString += self._settingsDict[key].name() + ','

    self._workspace.defineSet(name, varString[:-1])


  ## model builders methods ##

  def modelBuilder(self, MBType) :
    """returns the model builder of type MBType
    """

    if MBType in self._modelBuilders :
      return self._modelBuilders[MBType]
    else :
      from P2VVModelBuilders import createModelBuilder

      modelBuilder = createModelBuilder(self, MBType)
      if modelBuilder : self._modelBuilders[MBType] = modelBuilder
      return modelBuilder


  ## plot methods ##

  def addPlotObj(self, obj) :
    if obj: self._plotsStash.append(obj)

  def plotsStash(self) :
    return self._plotsStash[:]


###############################################################################
class P2VVSetting :
  """general P2VV setting
  """

  def __init__(self, name, description = '', value = '') :
    self.setName(name)
    self.setDescription(description)
    self.setValue(value)

  def name(self) :
    return self._name

  def description(self) :
    return self._description

  def value(self) :
    return self._value

  def setName(self, name) :
    if type(name) is str and len(name) > 0:
      self._name = name
    else :
      print "P2VV - ERROR: P2VVSetting.setName: argument 'name' is not a string or an empty string: setting name 'noName'"
      self._name = 'noName'

  def setDescription(self, description) :
    if type(description) is str :
      self._description = description
    else :
      print "P2VV - ERROR: P2VVSetting.setDescription(%s): argument 'description' is not a string"\
        % self._name
      self._description = ''

  def setValue(self, value) :
    self._value = value

  def type(self) :
    return self.__class__.__name__


###############################################################################

class RooSetting(P2VVSetting) :
  """general RooFit variable setting
  """

  def __init__(self, name, description = '', observable = False) :
    P2VVSetting(name, description, '')
    self.setObservable(observable)
    self._declared = False

  def name(self, nameSel = '') :
    return self._name

  def observable(self) :
    return self._observable

  def setName(self, name) :
    if self._declared :
      print "P2VV - ERROR: RooSetting.setName(%s): variable has already been declared"\
          % self._name
      return

    P2VVSetting.setName(self, name)

  def setObservable(self, observable) :
    if observable :
      self._observable = True
    else :
      self._observable = False

  def declared(self) :
    return self._declared

  def _declare(self, workspace, factoryString) :
    import RooFitDecorators

    # check if variable has not been declared yet
    if self._declared :
      print "P2VV - ERROR: RooSetting.declare(%s): variable has already been declared"\
          % self._name
      return

    # declare variable if it's not already in workspace
    if self.name('unblind') not in workspace :
      workspace.factory(factoryString)

    self._declared = True


###############################################################################

class RooRealSetting(RooSetting) :
  """RooRealVar setting
  """

  def __init__(self, name, description = '', realType = '', value = 0.,
      minValue = '', maxValue = '',
      blindParams = ('RooUnblindUniform', 'blindString', 1.)) :
    self._declared = False
    self.setName(name)
    self.setDescription(description)
    self._observable  = False
    self._blinded     = False
    self._blindParams = ()
    self._value       = 0.
    self._minValue    = ''
    self._maxValue    = ''

    RooSetting.__init__(self, name, description, False)
    self.setValue(value)
    self.setMinMax(minValue, maxValue)
    self.setRealType(realType)
    self.setBlindParams(blindParams)

  def name(self, nameSel = 'blind') :
    if self._blinded and nameSel == 'blind': return self._name + '_blind'
    elif self._blinded and nameSel == 'unblind': return self._name + '_unblind'
    else : return self._name

  def blinded(self) :
    return self._blinded

  def minValue(self) :
    return self._minValue

  def maxValue(self) :
    return self._maxValue

  def set(self, name = '', descr = '', realType = 'type', val = '',
       min = 'min', max = 'max', blindParams = '') :
    if self._declared :
      print "P2VV - ERROR: RooRealSetting.set(%s): variable has already been declared"\
          % self._name
      return

    if name  != ''        : self.setName(name)
    if descr != ''        : self.setDescription(descr)
    if realType != 'type' : self.setRealType(realType)
    if val   != ''        : self.setValue(val)

    if min != 'min' and max != 'max' : self.setMinMax(min, max)
    elif min != 'min' : self.setMinMax(min, self.maxValue())
    elif max != 'max' : self.setMinMax(self.minValue(), max)

    if blindParams != '' : self.setBlindParams(blindParams)

  def setRealType(self, realType) :
    if realType == 'obs' :
      self.setObservable(True)
    else :
      self.setObservable(False)
      if realType == 'blind' : self.setBlinded(True)
      else : self.setBlinded(False)

  def setBlindParams(self, blindParams) :
    if type(blindParams) is tuple and len(blindParams) == 3 :
      self._blindParams = blindParams
    else :
      print "P2VV - ERROR: RooRealSetting.setBlindParams(%s): parameters do not have the required format"\
          % self._name

  def setValue(self, value) :
    if self._declared :
      print "P2VV - ERROR: RooRealSetting.setValue(%s): variable has already been declared"\
          % self._name
      return

    if type(value) is float or type(value) is int :
      self._value = float(value)
      self._forceMinMax()
    else :
      print "P2VV - ERROR: RooRealSetting.setValue(%s): value is not a float"\
          % self._name
      if type(self._minValue) is float :
        self._value = self._minValue
      else :
        self._value = 0.

  def setMinMax(self, minValue, maxValue) :
    if self._declared :
      print "P2VV - ERROR: RooRealSetting.setMinMax(%s): variable has already been declared"\
          % self._name
      return

    # set minimum and maximum values
    if (type(minValue) is float or type(minValue) is int)\
        and (type(maxValue) is float or type(maxValue) is int)\
        and float(maxValue) > float(minValue) :
      self._minValue = float(minValue)
      self._maxValue = float(maxValue)

      # force value between new minimum and new maximum
      self._forceMinMax()
    else :
      self._minValue = ''
      self._maxValue = ''

      # variable can't be observable without a range
      if self._observable :
        print "P2VV - WARNING: RooRealSetting.setMinMax(%s): a constant variable can't be observable: removing observable flag"\
            % self._name
        self.setObservable(False)

  def setObservable(self, observable) :
    if observable :
      # check range for a RooRealVar
      if type(self._minValue) is float\
          and type(self._maxValue is float) :
        self._observable = True
        self.setBlinded(False)
      else :
        print "P2VV - ERROR: RooRealSetting.setObservable(%s): a constant variable can't be observable"\
            % self._name
        self._observable = False
    else :
      self._observable = False

  def _forceMinMax(self) :
    if self._declared :
      print "P2VV - ERROR: RooRealSetting._forceMinMax(%s): variable has already been declared"\
          % self._name
      return

    if type(self._minValue) is float and type(self._maxValue) is float:
      if self._value < self._minValue :
        print "P2VV - WARNING: RooRealSetting._forceMinMax(%s): value is less than minimum value: forcing value into range"\
            % self._name
        self._value = self._minValue
      elif self._value > self._maxValue :
        print "P2VV - WARNING: RooRealSetting._forceMinMax(%s): value is greater than maximum value: forcing value into range"\
            % self._name
        self._value = self._maxValue

  def setBlinded(self, blind) :
    if blind :
      if not self._observable :
        self._blinded = True
      else :
        print "P2VV - ERROR: RooRealSetting.setBlinded(%s): an observable can't be blinded"\
            % self._name

    else :
      self._blinded = False

  def declare(self, workspace) :
    if self._minValue != '' and self._maxValue != '' :
      factStr = "%s[%f, %f, %f]"\
          % (self.name('blind'), self._value, self._minValue, self._maxValue)
    else :
      factStr = "%s[%f]"\
          % (self.name('blind'), self._value)

    if self._blinded :
      factStr =  'UnblindUniform::%s(%s, %f, %s)' % (self.name('unblind'),
          self._blindParams[1], self._blindParams[2], factStr)

    self._declare(workspace, factStr)


###############################################################################

class RooCatSetting(RooSetting) :
  """RooCategory setting
  """

  def __init__(self, name, description = '', observable = False, value = '',
      catTypesDict = {}) :
    self._declared = False
    self.setName(name)
    self.setDescription(description)
    self._observable   = False
    self._value        = ''
    self._catTypesDict = {}

    RooSetting.__init__(self, name, description, False)
    self.setCatTypesDict(catTypesDict)
    self.setValue(value)
    self.setObservable(observable)

  def catTypesDict(self) :
    return self._catTypesDict.copy()

  def set(self, name = '', descr = '', obs = '', val = '', dict = '') :
    if self._declared :
      print "P2VV - ERROR: RooCatSetting.set(%s): variable has already been declared"\
          % self._name
      return

    if name  != '' : self.setName(name)
    if descr != '' : self.setDescription(descr)
    if obs   != '' : self.setObservable(obs)
    if val   != '' : self.setValue(val)
    if dict  != '' : self.setCatTypesDict(dict)

  def setObservable(self, observable) :
    if observable :
      # check number of types for a RooCategory
      if len(self._catTypesDict) > 1 :
        self._observable = True
      else :
        print "P2VV - ERROR: RooCatSetting.setObservable(%s): a constant variable can't be observable"\
            % self._name
        self._observable = False
    else :
      self._observable = False

  def setCatTypesDict(self, catTypesDict) :
    if self._declared :
      print "P2VV - ERROR: RooCatSetting.setCatTypesDict(%s): variable has already been declared"\
          % self._name
      return

    # set dictionary of category types and indices
    self._catTypesDict = {}

    if type(catTypesDict) is not dict :
      # don't define category types if catTypesDict is not a dictionary
      print "P2VV - ERROR: RooCatSetting.setCatTypesDict(%s): argument 'catTypesDict' is not a dictionary: not defining any types"\
          % self._name
    else :
      for catTypeIndex, catType in catTypesDict.iteritems() :
        # check category type
        if type(catType) is not str or len(catType) <= 0 :
          print "P2VV - ERROR: RooCatSetting.setCatTypesDict(%s): category type is not a string or an empty string: not defining type"\
              % self._name
          continue

        # check category index
        if type(catTypeIndex) is not int :
          print "P2VV - ERROR: RooCatSetting.setCatTypesDict(%s): index of category type '%s' is not an integer: not defining type"\
              % (self._name, catType)
          continue

        # set category type
        self._catTypesDict[catTypeIndex] = catType

    if len(self._catTypesDict) < 2 and self._observable:
      # variable can't be observable with less than two types
      print "P2VV - WARNING: RooCatSetting.setCatTypesDict(%s): a constant variable can't be observable: removing observable flag"\
          % catType
      self.setObservable(False)

    # set value
    self.setValue(self._value)

  def setValue(self, value) :
    if self._declared :
      print "P2VV - ERROR: RooCatSetting.setValue(%s): variable has already been declared"\
          % self._name
      return

    if len(self._catTypesDict) == 0 :
      self._value = ''
    elif type(value) is not str\
        or value not in self._catTypesDict.values() :
      self._value = sorted(self._catTypesDict.keys())[0]
    else :
      self._value = value

  def declare(self, workspace) :
    factoryString = "%s[" % self._name
    for catTypeIndex in sorted(self._catTypesDict.keys()) :
      factoryString += "%s=%d, "\
          % (self._catTypesDict[catTypeIndex], catTypeIndex)
    factoryString = factoryString[:-2] + ']'

    self._declare(workspace, factoryString)


###############################################################################

class RooFormSetting(RooSetting) :
  """RooFormulaVar setting
  """

  def __init__(self, name, description = '', value = '1.', variables = []) :
    self._declared = False
    self.setName(name)
    self.setDescription(description)
    self._observable = False

    RooSetting.__init__(self, name, description, False)
    self.setValue(value)
    self.setVariables(variables)

  def variables(self) :
    return self._variables

  def set(self, name = '', descr = '', val = '', vars = '') :
    if self._declared :
      print "P2VV - ERROR: RooCatSetting.set(%s): variable has already been declared"\
          % self._name
      return

    if name  != '' : self.setName(name)
    if descr != '' : self.setDescription(descr)
    if val   != '' : self.setValue(val)
    if vars  != '' : self.setVariables(vars)

  def setValue(self, value) :
    if self._declared :
      print "P2VV - ERROR: RooFormSetting.setValue(%s): formula has already been declared"\
          % self._name
      return

    if type(value) is str :
      self._value = value
    else :
      print "P2VV - ERROR: RooFormSetting.setValue(%s): value is not a string"\
          % self._name
      self._value = ''

  def setObservable(self, observable) :
    if observable :
      print "P2VV - ERROR: RooFormSetting.setObservable(%s): a formula cannot be 'observable'"\
          % self._name

  def setVariables(self, variables) :
    if self._declared :
      print "P2VV - ERROR: RooFormSetting.setVariables(%s): formula has already been declared"\
          % self._name
      return

    self._variables = []

    if type(variables) is list :
      for var in variables :
        if type(var) is str and len(var) > 0 :
          self._variables.append(var)
        else :
          print "P2VV - ERROR: RooFormSetting.setVariables(%s): 'variable' is not a string or an empty string"\
              % self._name

    else :
      print "P2VV - ERROR: RooFormSetting.setVariables(%s): variables are not in a list"\
          % self._name

  def declare(self, workspace, settings) :
    factoryString = "expr::%s('%s', {" % (self._name, self._value)
    for var in self._variables :
      factoryString += "%s, " % settings[var].name('unblind')
    factoryString = factoryString[:-2] + '})'

    self._declare(workspace, factoryString)

