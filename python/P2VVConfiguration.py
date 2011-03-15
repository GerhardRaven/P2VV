###############################################################################
## function that returns default P2VV configuration for a given mode         ##
###############################################################################
def getP2VVConfig(mode = '', options = [], createWS = True) :
  from math import pi

  # create configuration instance
  if createWS is True :
    config = P2VVConfiguration('', '')
  else :
    config = P2VVConfiguration(False, '')

  # check arguments
  if type(mode) is not str :
    print "ERROR: getP2VVConfig: argument 'mode' is not a string"
    return config
  if type(options) is not list :
    print "ERROR: getP2VVConfig: argument 'options' is not a list"
    return config

  # set mode in configuration
  config.addSetting('P2VVMode', P2VVSetting('mode',
      'P2VV mode', mode))

  # get options
  optList = []
  for opt in options :
    if type(opt) is str :
      optList.append(opt.strip())

  config.addSetting('P2VVOptions', P2VVSetting('options',
      'P2VV options', optList))


  ######################
  ## general settings ##
  ######################
  config.addSetting('sigPDFName', P2VVSetting('sigPDFName',
      'ROOT name of the signal PDF', mode + 'PDF'))


  ###############
  ## constants ##
  ###############
  config.addSetting('zeroConst', RooRealSetting('zero',
      'zero constant', False, 0., '', ''))
  config.addSetting('oneConst', RooRealSetting('one',
      'one constant', False, 1., '', ''))
  config.addSetting('minusConst', RooRealSetting('minus',
      'minus constant', False, -1., '', ''))


  ##################
  ## decay angles ##
  ##################
  if 'noAngles' in optList\
      or not mode in ['Bd2mumuKstar', 'Bd2JpsiKstar', 'Bs2Jpsiphi']:
    anglesType = []
  else :
    if 'transAngles' in optList :
      if 'helAngles' in optList :
        anglesType = ['hel', 'trans']
      else :
        anglesType = ['trans']
    else :
      anglesType = ['hel']

    # helicity angles
    if 'hel' in anglesType :
      config.addSetting('helcthetaK', RooRealSetting('helcthetaK',
          'cosine of kaon polarization angle', True, 0., -1., 1.))
      config.addSetting('helcthetal', RooRealSetting('helcthetaL',
          'cosine of lepton polarization angle', True, 0., -1., 1.))
      config.addSetting('helphi', RooRealSetting('helphi',
          'angle between decay planes', True, 0., -pi, pi))

    # transversity angles
    if 'trans' in anglesType :
      config.addSetting('trcpsi', RooRealSetting('trcpsi',
          'cosine of kaon polarization angle', True, 0., -1., 1.))
      config.addSetting('trctheta', RooRealSetting('trctheta',
          'cosine of transversity polar angle', True, 0., -1., 1.))
      config.addSetting('trphi', RooRealSetting('trphi',
          'transversity azimuthal angle', True, 0., -pi, pi))

  config.addSetting('anglesType', P2VVSetting('anglesType',
      'type of angles: helicity or transversity', anglesType))


  #########################
  ## B0->mumuK* settings ##
  #########################
  if mode is 'Bd2mumuKstar' :
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
    BDecClass = 'RooBTagDecay'
    config.addSetting('BDecayClass', P2VVSetting('BDecayClass',
        'RooFit class for time PDF', BDecClass))


    ### observables ###

    # decay type
    config.addSetting('decayType', RooCatSetting('decaytype',
        'J/psiX decay type', True, 'Jpsiphi', {10:'JpsiKplus',
        11:'JpsiKmin', 20:'JpsiKstar0', 21:'JpsiKstarbar0', 40:'Jpsiphi'}))

    # trigger
    config.addSetting('unbiased', RooCatSetting('unbiased',
        'unbiased trigger?', True, 'yes', {1 : 'yes', 0 : 'no'}))

    # flavour tags
    config.addSetting('iTag', RooCatSetting('iTag',
        'initial state flavour tag', True, 'bbar', {-1 : 'b', +1 : 'bbar'}))
    config.addSetting('misTag', RooRealSetting('misTag',
        'mis-tag fraction', False, 0., 0., 0.5))
    if mode is 'Bd2JpsiKstar' :
      config.addSetting('fTag', RooCatSetting('fTag',
          'final state flavour tag', True, 'bbar', {-1 : 'b', +1 : 'bbar'}))

    # B lifetime
    config.addSetting('BLifetime', RooRealSetting('t',
        'B lifetime (ps)', True, 0., -5., 14.))
    # Note 2 uses [0.3,14] -- so maybe switch to [-4,14] instead...
    config.addSetting('BLifetimeError', RooRealSetting('sigmat',
        'B lifetime error (ps)', True, 0.005, 0.005, 0.1))

    # masses
    config.addSetting('mJpsi', RooRealSetting('mdau1',
        'J/psi mass (MeV)', True, 3097., 3097. - 60., 3097. + 40.))
    if mode in ['Bu2JpsiK', 'Bd2JpsiKstar'] :
      config.addSetting('mB', RooRealSetting('m',
          'B0 mass (MeV)', True, 5279., 5279. - 50., 5279. + 50.))
      if mode is 'Bd2JpsiKstar' :
        config.addSetting('mKstar', RooRealSetting('mdau2',
            'K*0 mass (MeV)', True, 892., 892. - 50., 892. + 50.))
    elif mode is 'Bs2Jpsiphi' :
      config.addSetting('m', RooRealSetting('m',
          'B_s0 mass (MeV)', True, 5366., 5366. - 50., 5366. + 50.))
      config.addSetting('mphi', RooRealSetting('mdau2',
          'phi mass (MeV)', True, 1019.455, 1019.455 - 12., 1019.455 + 12.))
      # Note: +-10 Mev/c^2 keeps all of the phi signal, and kills 1/2 of the
      # background -- but the roadmap uses +-12 instead...


    ### physics parameters ###
    if mode in ['Bd2JpsiKstar', 'Bs2Jpsiphi'] :
      # amplitudes
      # note: initial values from arXiv:0704.0522v2 [hep-ex] BaBar PUB-07-009
      config.addSetting('incKSWave', P2VVSetting('incKSWave',
        'include KK or Kpi S-wave?', True))

      config.addSetting('ReA0', RooRealSetting('ReA0',
          'Re(A_0)', False, 1., '', ''))
      config.addSetting('ImA0', RooRealSetting('ImA0',
          'Im(A_0)', False, 0., '', ''))
      config.addSetting('ReApar', RooRealSetting('ReApar',
          'Re(A_par)', False, -0.602, -1., 1.))
      config.addSetting('ImApar', RooRealSetting('ImApar',
          'Im(A_par)', False, -0.129, -1., 1.))
      config.addSetting('ReAperp', RooRealSetting('ReAperp',
          'Re(A_perp)', False, -0.630, -1., 1.))
      config.addSetting('ImAperp', RooRealSetting('ImAperp',
          'Im(A_perp)', False, 0.149, -1., 1.))
      config.addSetting('ReAS', RooRealSetting('ReAS',
          'Re(A_S)', False, -0.168, -1., 1.))
      config.addSetting('ImAS', RooRealSetting('ImAS',
          'Im(A_S)', False, 0.252, -1., 1.))

      #config.addSetting('A0Mag2', RooRealSetting('A0Mag2',
      #    '|A0|^2', False, 1., '', ''))
      #config.addSetting('AparMag2', RooRealSetting('AparMag2',
      #    '|A_par|^2', False, 0.379, 0., 1.2))
      #config.addSetting('AperpMag2', RooRealSetting('AperpMag2',
      #    '|A_perp|^2', False, 0.419, 0., 1.2))
      #config.addSetting('ASMag2', RooRealSetting('ASMag2',
      #    '|A_S|^2', False, 0.092, 0., 0.6))

      #config.addSetting('A0Ph', RooRealSetting('deltaz',
      #    'delta_0', False, 0., '', ''))
      #config.addSetting('AparPh', RooRealSetting('deltapar',
      #    'delta_par', False, -2.93, -2. * pi, 2. * pi))
      #config.addSetting('AperpPh', RooRealSetting('deltaperp',
      #    'delta_perp', False, 2.91, -2. * pi, 2. * pi))
      #config.addSetting('ASPh', RooRealSetting('deltas',
      #    'delta_S', False, 2.2, -2. * pi, 2. * pi))

      #config.addSetting('ReA0', RooFormSetting('ReA0',
      #    'Re(A_0)', 'sqrt(@0) * cos(@1)', ['A0Mag2', 'A0Ph']))
      #config.addSetting('ImA0', RooFormSetting('ImA0',
      #    'Im(A_0)', 'sqrt(@0) * sin(@1)', ['A0Mag2', 'A0Ph']))
      #config.addSetting('ReApar', RooFormSetting('ReApar',
      #    'Re(A_par)', 'sqrt(@0) * cos(@1)', ['AparMag2', 'AparPh']))
      #config.addSetting('ImApar', RooFormSetting('ImApar',
      #    'Im(A_par)', 'sqrt(@0) * sin(@1)', ['AparMag2', 'AparPh']))
      #config.addSetting('ReAperp', RooFormSetting('ReAperp',
      #    'Re(A_perp)', 'sqrt(@0) * cos(@1)', ['AperpMag2', 'AperpPh']))
      #config.addSetting('ImAperp', RooFormSetting('ImAperp',
      #    'Im(A_perp)', 'sqrt(@0) * sin(@1)', ['AperpMag2', 'AperpPh']))
      #config.addSetting('ReAS', RooFormSetting('ReAS',
      #    'Re(A_S)', 'sqrt(@0) * cos(@1)', ['ASMag2', 'ASPh']))
      #config.addSetting('ImAS', RooFormSetting('ImAS',
      #    'Im(A_S)', 'sqrt(@0) * sin(@1)', ['ASMag2', 'ASPh']))

      # lifetime and mixing
      if mode is 'Bd2JpsiKstar' :
        config.addSetting('Gamma', RooRealSetting('gamma',
            'Gamma (ps^-1)', False, 0.65, 0.4, 0.9))
        config.addSetting('dGamma', RooRealSetting('t_sig_dG',
            'delta Gamma_s (ps^-1)', False, 0., '', ''))
        config.addSetting('dm', RooRealSetting('t_sig_dm',
            'delta m_s (ps^-1)', False, 0.51, '', ''))

      elif mode is 'Bs2Jpsiphi' :
        config.addSetting('Gamma', RooRealSetting('gamma',
            'Gamma_s (ps^-1)', False, 0.68, 0.4, 0.9))
        config.addSetting('dGamma', RooRealSetting('t_sig_dG',
            'delta Gamma_s (ps^-1)', False, 0.05, -0.3, 0.3))
        config.addSetting('dm', RooRealSetting('t_sig_dm',
            'delta m_s (ps^-1)', False, 17.8, '', ''))

      config.addSetting('BMeanLife', RooFormSetting('t_sig_tau',
          'B mean lifetime', '1. / @0', ['Gamma']))

      # CP violation
      if mode is 'Bd2JpsiKstar' and BDecClass is 'RooBTagDecay' :
        config.addSetting('CCP', RooRealSetting('C',
            'B0 lambda parameter C', False, 0., '', ''))
      elif mode is 'Bs2Jpsiphi' :
        config.addSetting('phiCP', RooRealSetting('phis',
            'CP violation parameter phi_s', False, 0.8, -pi, pi))
        if BDecClass is 'RooBTagDecay' :
          config.addSetting('CCP', RooRealSetting('C',
              'B0_s lambda parameter C', False, 0., '', ''))
        config.addSetting('DCP', RooFormSetting('D',
            'B0_s lambda parameter D', 'cos(@0)', ['phis']))
        config.addSetting('SCP', RooFormSetting('S',
            'B0_s lambda parameter S', '-sin(@0)', ['phis']))

      # tagging parameters
      if BDecClass is 'RooBTagDecay' :
        config.addSetting('tagDilution', RooFormSetting('tagDilution',
            'mis-tag dilution', '1 - 2. * @0', ['misTag']))
        config.addSetting('ADilMisTag', RooRealSetting('ADilMisTag',
            'dilution/mis-tag asymmetry', False, 0., '', ''))
        if mode is 'Bd2JpsiKstar' :
          config.addSetting('ANorm', RooRealSetting('ANorm',
            'normalization asymmetry', False, 0., '', ''))
        config.addSetting('avgCEven', RooRealSetting('avgCEven',
            'CP average even coefficients', False, 1., '', ''))
        config.addSetting('avgCOdd', RooFormSetting('avgCOdd',
            'CP average odd coefficients', '-@0', ['C']))


  return config


###############################################################################
## class for configuration of the P2VV framework                   ##
## contains a dictionary with P2VV settings and a RooFit workspace ##
#####################################################################
class P2VVConfiguration :
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

  # get the value of a setting
  def value(self, key) :
    # check type of the setting key
    if type(key) is not str :
      print "ERROR: P2VVConfiguration.value: setting key is not of type 'str'"
      return None

    # get value of setting
    setting = self[key]
    if setting : return setting.value()
    else : return None    

  # add a setting to the settings dictionary or overwrite an old setting
  def addSetting(self, key, setting, overWriteOld = True) :
    # check type of the setting key
    if type(key) is not str :
      print "ERROR: P2VVConfiguration.addSetting: setting key is not of type 'str'"
      return

    # check type of the setting
    if not issubclass(setting.__class__, P2VVSetting) :
      print "ERROR: P2VVConfiguration.addSetting: class of 'setting' does not inherit from 'P2VVSetting'"
      return

    # add new setting to settings dictionary or overwrite old setting
    if not key in self._settingsDict or overWriteOld:
      self._settingsDict[key] = setting


  ## RooFit workspace methods ##

  # get workspace
  def workspace(self) :
    return self._workspace

  # set workspace
  def setWorkspace(self, workspace = '') :
    from ROOT import RooWorkspace

    if type(workspace) is RooWorkspace :
      self._workspace = workspace
    elif type(workspace) is str :
      if workspace is '' :
        self._workspace = RooWorkspace('P2VVRooWS')
      else :
        self._workspace = RooWorkspace(workspace)
    else :
      print "WARNING: P2VVConfiguration.setWorkspace: argument is not a RooWorkspace: no workspace set"
      self._workspace = ''

  # write workspace to ROOT file
  def writeWorkspace(self, WSPath = '') :
    if self._workspace is '' :
      print "ERROR: P2VVConfiguration.writeWorkspace: no workspace set"
      return

    if type(WSPath) is str and len(WSPath) > 0 :
      self._workspace.writeToFile(WSPath)
    elif 'WSPath' in self._settingsDict\
        and type(self._settingsDict['WSPath'].value()) is str\
        and len(self._settingsDict['WSPath'].value()) > 0 :
      self._workspace.writeToFile(self._settingsDict['WSPath'].value())
    else :
      print "ERROR: P2VVConfiguration.writeWorkspace: no workspace file path set"

  # get workspace file path
  def WSFilePath(self) :
    if 'WSPath' in self._settingsDict :
      return self._settingsDict['WSPath'].value()
    else :
      return ''

  # set workspace file path
  def setWSFilePath(self, WSPath) :
    if type(WSPath) is str and len(WSPath) > 0 :
      self.addSetting('WSPath', P2VVSetting('WSPath', 'workspace file path',
          WSPath), True)
    else :
      print "ERROR: P2VVConfiguration.setWSFilePath: argument is not a string or an empty string: no workspace file path set"

  # declare RooFit variables in settings dictionary and put them in workspace
  def declareRooVars(self, varType = 'var') :
    if self._workspace is '' :
      print "ERROR: P2VVConfiguration.declareRooVars: no workspace set: can't declare any variables"
      return

    # loop over settings and find RooSettings to put in workspace
    formSettings     = []
    obsSettingKeys  = []
    for key, setting in self._settingsDict.iteritems() :
      if not issubclass(setting.__class__, RooSetting) : continue

      if setting.type() is 'RooFormSetting' and varType in ['var', 'form'] :
        # put RooFormSettings in a list for later declaration
        formSettings.append(setting)
        for var in setting.variables() :
          if var not in self._settingsDict :
            formSettings.remove(setting)
            break

      elif (setting.observable() and varType in ['var', 'obs'])\
          or (not setting.observable() and varType in ['var', 'par']) :
        # declare variable
        setting.declare(self._workspace)
        if setting.observable() : obsSettingKeys.append(key)

    # define sets of observables
    self.defineRooSet('observables', obsSettingKeys)
    if 'anglesType' in self._settingsDict :
      if 'hel' in self._settingsDict['anglesType'].value() :
        self.defineRooSet('helAngles',
            ['helcthetaK', 'helcthetal', 'helphi'])
      if 'trans' in self._settingsDict['anglesType'].value() :
        self.defineRooSet('transAngles',
            ['trcpsi', 'trctheta', 'trphi'])

    # declare RooFormulaVars
    if varType in ['var', 'form'] :
      while True :
        nDeclared = 0
        for setting in formSettings :
          if setting.declared() : continue

          # check if all the needed variables have been declared
          declare = True
          for var in setting.variables() :
            if not self._settingsDict[var].declared() :
              declare = False
              break

          # declare formula
          if declare :
            nDeclared += 1
            setting.declare(self._workspace, self)

        # exit loop if nothing was declared any more
        if nDeclared == 0 : break

  def defineRooSet(self, name, settingKeysList) :
    if self._workspace is '' :
      print "ERROR: P2VVConfiguration.defineRooSet: no workspace set: can't define set"
      return

    if type(name) is not str or len(name) < 1 :
      print "ERROR: P2VVConfiguration.defineRooSet: 'name' is not a string or an empty string"
      return

    if type(settingKeysList) is not list :
      print "ERROR: P2VVConfiguration.defineRooSet: 'settingsList' is not a list"
      return

    varString = ''
    for key in settingKeysList :
      if key in self._settingsDict :
        varString += self._settingsDict[key].name() + ','

    self._workspace.defineSet(name, varString[:-1])


###############################################################################
## general P2VV setting ##
##########################
class P2VVSetting :
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
      print "ERROR: P2VVSetting.setName: argument 'name' is not a string or an empty string: setting name 'noName'"
      self._name = 'noName'

  def setDescription(self, description) :
    if type(description) is str :
      self._description = description
    else :
      print "ERROR: P2VVSetting.setDescription(%s): argument 'description' is not a string"\
        % self._name
      self._description = ''

  def setValue(self, value) :
    self._value = value

  def type(self) :
    return self.__class__.__name__


###############################################################################
## general RooFit variable setting ##
#####################################
class RooSetting(P2VVSetting) :
  def __init__(self, name, description = '', observable = False) :
    P2VVSetting(name, description, '')
    self.setObservable(observable)
    self._declared = False

  def observable(self) :
    return self._observable

  def setName(self, name) :
    if self._declared :
      print "ERROR: RooSetting.setName(%s): variable has already been declared"\
          % self._name
      return

    P2VVSetting.setName(self, name)

  def setObservable(self, observable) :
    if observable is True :
      self._observable = True
    else :
      self._observable = False

  def declared(self) :
    return self._declared

  def _declare(self, workspace, factoryString) :
    if self._declared :
      print "ERROR: RooSetting.declare(%s): variable has already been declared"\
          % self._name
      return

    workspace.factory(factoryString)
    self._declared = True


###############################################################################
## RooRealVar setting ##
########################
class RooRealSetting(RooSetting) :
  def __init__(self, name, description = '', observable = False, value = 0.,
      minValue = '', maxValue = '') :
    self._declared = False
    self.setName(name)
    self.setDescription(description)
    self._observable = False
    self._value      = 0.
    self._minValue   = ''
    self._maxValue   = ''

    RooSetting.__init__(self, name, description, False)
    self.setValue(value)
    self.setRange(minValue, maxValue)
    self.setObservable(observable)

  def minValue(self, minValue) :
    return self._minValue

  def maxValue(self, maxValue) :
    return self._maxValue

  def setValue(self, value) :
    if self._declared :
      print "ERROR: RooRealSetting.setValue(%s): variable has already been declared"\
          % self._name
      return

    if type(value) is float or type(value) is int :
      self._value = float(value)
      self._forceIntoRange()
    else :
      print "ERROR: RooRealSetting.setValue(%s): value is not a float"\
          % self._name
      if type(self._minValue) is float :
        self._value = self._minValue
      else :
        self._value = 0.

  def setRange(self, minValue, maxValue) :
    if self._declared :
      print "ERROR: RooRealSetting.setRange(%s): variable has already been declared"\
          % self._name
      return

    # set minimum and maximum values
    if (type(minValue) is float or type(minValue) is int)\
        and (type(maxValue) is float or type(maxValue) is int)\
        and float(maxValue) > float(minValue) :
      self._minValue = float(minValue)
      self._maxValue = float(maxValue)

      # force value into new range
      self._forceIntoRange()
    else :
      self._minValue = ''
      self._maxValue = ''

      # variable can't be observable without a range
      if self._observable :
        print "WARNING: RooRealSetting.setRange(%s): a constant variable can't be observable: removing observable flag"\
            % self._name
        self.setObservable(False)

  def setObservable(self, observable) :
    if observable is True :
      # check range for a RooRealVar
      if type(self._minValue) is float\
          and type(self._maxValue is float) :
        self._observable = True
      else :
        print "ERROR: RooRealSetting.setObservable(%s): a constant variable can't be observable"\
            % self._name
        self._observable = False
    else :
      self._observable = False

  def _forceIntoRange(self) :
    if self._declared :
      print "ERROR: RooRealSetting._forceIntoRange(%s): variable has already been declared"\
          % self._name
      return

    if type(self._minValue) is float and type(self._maxValue) is float:
      if self._value < self._minValue :
        print "WARNING: RooRealSetting._forceIntoRange(%s): value is less than minimum value: forcing value into range"\
            % self._name
        self._value = self._minValue
      elif self._value > self._maxValue :
        print "WARNING: RooRealSetting._forceIntoRange(%s): value is greater than maximum value: forcing value into range"\
            % self._name
        self._value = self._maxValue

  def declare(self, workspace) :
    if self._minValue is not '' and self._maxValue is not '' :
      factoryString = "%s[%f, %f, %f]"\
          % (self._name, self._value, self._minValue, self._maxValue)
    else :
      factoryString = "%s[%f]" % (self._name, self._value)

    self._declare(workspace, factoryString)


###############################################################################
## RooCategory setting ##
#########################
class RooCatSetting(RooSetting) :
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

  def setObservable(self, observable) :
    if observable is True :
      # check number of types for a RooCategory
      if len(self._catTypesDict) > 1 :
        self._observable = True
      else :
        print "ERROR: RooCatSetting.setObservable(%s): a constant variable can't be observable"\
            % self._name
        self._observable = False
    else :
      self._observable = False

  def setCatTypesDict(self, catTypesDict) :
    # set dictionary of category types and indices
    self._catTypesDict = {}

    if self._declared :
      print "ERROR: RooCatSetting.setCatTypesDict(%s): variable has already been declared"\
          % self._name
      return

    if type(catTypesDict) is not dict :
      # don't define category types if catTypesDict is not a dictionary
      print "ERROR: RooCatSetting.setCatTypesDict(%s): argument 'catTypesDict' is not a dictionary: not defining any types"\
          % self._name
    else :
      for catTypeIndex, catType in catTypesDict.iteritems() :
        # check category type
        if type(catType) is not str or len(catType) <= 0 :
          print "ERROR: RooCatSetting.setCatTypesDict(%s): category type is not a string or an empty string: not defining type"\
              % self._name
          continue

        # check category index
        if type(catTypeIndex) is not int :
          print "ERROR: RooCatSetting.setCatTypesDict(%s): index of category type '%s' is not an integer: not defining type"\
              % (self._name, catType)
          continue

        # set category type
        self._catTypesDict[catTypeIndex] = catType

    if len(self._catTypesDict) < 2 and self._observable:
    # variable can't be observable without a range
      print "WARNING: RooCatSetting.setCatTypesDict(%s): a constant variable can't be observable: removing observable flag"\
          % catType
      self.setObservable(False)

    # set value
    self.setValue(self._value)

  def setValue(self, value) :
    if self._declared :
      print "ERROR: RooCatSetting.setValue(%s): variable has already been declared"\
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
    for catTypeIndex, catType in self._catTypesDict.iteritems() :
      factoryString += "%s=%d, " % (catType, catTypeIndex)
    factoryString = factoryString[:-2] + ']'

    self._declare(workspace, factoryString)


###############################################################################
## RooFormulaVar setting ##
###########################
class RooFormSetting(RooSetting) :
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

  def setValue(self, value) :
    if self._declared :
      print "ERROR: RooFormSetting.setValue(%s): formula has already been declared"\
          % self._name
      return

    if type(value) is str :
      self._value = value
    else :
      print "ERROR: RooFormSetting.setValue(%s): value is not a string"\
          % self._name
      self._value = ''

  def setObservable(self, observable) :
    if observable is True :
      print "ERROR: RooFormSetting.setObservable(%s): a formula cannot be 'observable'"\
          % self._name

  def setVariables(self, variables) :
    if self._declared :
      print "ERROR: RooFormSetting.setVariables(%s): formula has already been declared"\
          % self._name
      return

    self._variables = []

    if type(variables) is list :
      for var in variables :
        if type(var) is str and len(var) > 0 :
          self._variables.append(var)
        else :
          print "ERROR: RooFormSetting.setVariables(%s): 'variable' is not a string or an empty string"\
              % self._name

    else :
      print "ERROR: RooFormSetting.setVariables(%s): variables are not in a list"\
          % self._name

  def declare(self, workspace, settings) :
    factoryString = "expr::%s('%s', {" % (self._name, self._value)
    for var in self._variables :
      factoryString += "%s, " % settings[var].name()
    factoryString = factoryString[:-2] + '})'

    self._declare(workspace, factoryString)


###############################################################################
## test ##
##########
config = getP2VVConfig(mode = 'Bs2Jpsiphi', options = 'transAngles', createWS = True)
config.declareRooVars()
config.workspace().Print()

