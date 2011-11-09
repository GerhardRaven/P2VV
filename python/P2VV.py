###############################################################################
## P2VV: common tools                                                        ##
##                                                                           ##
## authors:                                                                  ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           ##
##                                                                           ##
###############################################################################

def loadP2VVLib() :
  """function that loads the P2VV library

  assumes $P2VVPATH/lib is in $LD_LIBRARYPATH
  """

  from ROOT import gSystem
  gSystem.Load("libP2VV")


###############################################################################

def setRooFitOutput() :
  """controls the RooFit output messages
  """

  from ROOT import RooFit, RooMsgService

  # get message service instance
  msgServ = RooMsgService.instance()

  # remove all output streams
  for stream in range(msgServ.numStreams()) :
    msgServ.deleteStream(stream)

  # add default streams for P2VV
  msgServ.addStream(RooFit.PROGRESS)
  msgServ.addStream(RooFit.INFO, RooFit.Topic(RooFit.Minimization))
  #msgServ.getStream(1).addTopic(RooFit.Plotting)
  msgServ.getStream(1).addTopic(RooFit.Fitting)
  msgServ.getStream(1).addTopic(RooFit.Eval)
  msgServ.getStream(1).addTopic(RooFit.Caching)
  msgServ.getStream(1).addTopic(RooFit.ObjectHandling)
  msgServ.getStream(1).addTopic(RooFit.InputArguments)
  msgServ.getStream(1).addTopic(RooFit.DataHandling)
  msgServ.getStream(1).addTopic(RooFit.NumIntegration)


###############################################################################

def registerMultiCatGen() :
  """registers an experimental fast(er) toy generator
  """

  from ROOT import RooMultiCatGenerator, RooNumGenFactory, RooNumGenConfig,\
      RooMsgService, RooFit

  RooMultiCatGenerator.registerSampler(RooNumGenFactory.instance())
  RooNumGenConfig.defaultConfig().methodND(False,True)\
      .setLabel("RooMultiCatGenerator")
  RooNumGenConfig.defaultConfig().methodND(False,True).Print()


###############################################################################

def readData(filePath, dataSetName, NTuple = False,
    observables = None, tagCat = '', initTag = '') :
  """reads data from file (RooDataSet or TTree(s))
  """

  if NTuple :
    from ROOT import RooDataSet, TChain

    # create data set from NTuple file(s)
    print "P2VV - INFO: readData: reading NTuple(s) '%s' from file(s) '%s'"\
        % (dataSetName, filePath)
    files = TChain(dataSetName)
    files.Add(filePath)
    data = RooDataSet(dataSetName, dataSetName, files, observables)

  else :
    from ROOT import TFile

    # get data set from file
    print "P2VV - INFO: readData: reading RooDataset '%s' from file '%s'"\
        % (dataSetName, filePath)
    file = TFile.Open(filePath, 'READ')
    data = file.Get(dataSetName)
    file.Close()

  return data


###############################################################################

def writeData(filePath, dataSetName, data, NTuple = False) :
  """writes data to file (RooDataSet or TTree)
  """

  from ROOT import TFile

  print "P2VV - INFO: writeData: writing RooDataSet '%s' to file '%s'"\
      % (dataSetName, filePath)

  file = TFile.Open(filePath, 'RECREATE')
  if NTuple : data.tree().Write(dataSetName)
  else : data.Write(dataSetName)
  file.Close()


###############################################################################

def convertLambda(config, fitResult, printToScreen = True) :
  """converts the CP violation parameter lambda from a fit result with
     'convertComplexPars'
  """

  import RooFitDecorators

  if config.value('lambdaCPType') == 'polar' :
    # convert from polar lambda to cartesian lambda
    params = [config['lambdaCPSq'].name(), config['phiCP'].name()]
    result = fitResult.result(params)

    polVals = [[result[1][0], -result[1][1]]]
    polCovs = [[result[2][0][0], -result[2][0][1]],
               [-result[2][1][0], result[2][1][1]]]

    convert = convertComplexPars(polVals, polCovs, 'polar, magSq', 'cart', 0)

    cartVals = convert[0]
    cartCovs = convert[1]

    polVals[0][1] = -polVals[0][1]
    polCovs[0][1] = -polCovs[0][1]
    polCovs[1][0] = -polCovs[1][0]

  else :
    # convert from cartesian lambda to polar lambda
    params = [config['ReLambdaCP'].name(), config['ImLambdaCP'].name()]
    result = fitResult.result(params)

    cartVals = [[result[1][0], -result[1][1]]]
    cartCovs = [[result[2][0][0], -result[2][0][1]],
               [-result[2][1][0], result[2][1][1]]]

    convert = convertComplexPars(cartVals, cartCovs, 'cart', 'polar, magSq', 0)

    cartVals[0][1] = -cartVals[0][1]
    cartCovs[0][1] = -cartCovs[0][1]
    cartCovs[1][0] = -cartCovs[1][0]

    polVals = convert[0]
    polCovs = convert[1]

  if printToScreen :
    # print cartesian values and covariance matrix to screen
    params = [('Re(lambda)', 'Im(lambda)')]

    print 'P2VV - INFO: convertLambda: cartesian lambda:'
    printComplexParams(params, cartVals, cartCovs)

    # print polar values and covariance matrix to screen
    params = [('|lambda|^2', 'phi')]

    print 'P2VV - INFO: convertLambda: polar lambda:'
    printComplexParams(params, polVals, polCovs)

  return convert


###############################################################################

def convertAmplitudes(config, fitResult, printToScreen = True) :
  """converts amplitudes from a fit result with 'convertComplexPars'
  """

  import RooFitDecorators

  KSWave = False
  if config.value('KSWave') and type(config.value('KSWave')) is str\
      and config.value('KSWave')[:7] == 'include' :
    KSWave = True

  if config.value('ampsType') == 'transPolar' :
    # convert from polar amplitudes to cartesian amplitudes
    params = []
    params.append(config['A0Mag2'].name())
    params.append(config['AparPh'].name())
    params.append(config['AperpMag2'].name())
    params.append(config['AperpPh'].name())
    if KSWave :
      params.append(config['ASMag2'].name())
      params.append(config['ASPh'].name())

    result = fitResult.result(params)

    polVals = [(result[1][0], 0.), (1. - result[1][0] - result[1][2],
        result[1][1]), (result[1][2], result[1][3])]
    if KSWave :
      polVals.append((result[1][4], result[1][5]))

    polCovs = [list(row) for row in result[2]]
    polCovs.insert(1, [0. for i in range(len(result[2]))])  # A0Ph
    polCovs.insert(2, [-result[2][0][i] - result[2][2][i]\
        for i in range(len(result[2]))])                    # AparMag2
    for rowIter, row in enumerate(polCovs) :
      row.insert(1, 0.)
      if rowIter == 1 :
        row.insert(2, 0.)
      elif rowIter == 2 :
        row.insert(2, result[2][0][0] + result[2][2][2] + 2. * result[2][0][2])
      else :
        row.insert(2, polCovs[2][rowIter])

    convert = convertComplexPars(polVals, polCovs, 'polar, magSq', 'cart', 1)
    cartVals = convert[0]
    cartCovs = convert[1]

  else :
    # convert from cartesian amplitudes to polar amplitudes
    params = []
    params.append(config['ReApar'].name())
    params.append(config['ImApar'].name())
    params.append(config['ReAperp'].name())
    params.append(config['ImAperp'].name())
    if KSWave :
      params.append(config['ReAS'].name())
      params.append(config['ImAS'].name())

    result = fitResult.result(params)

    cartVals = [(1., 0.), (result[1][0], result[1][1]), (result[1][2],
        result[1][3])]
    if KSWave :
      cartVals.append((result[1][4], result[1][5]))
  
    cartCovs = [list(row) for row in result[2]]
    cartCovs.insert(0, [0. for i in range(len(cartCovs[0]))])  # ReA0
    cartCovs.insert(1, [0. for i in range(len(cartCovs[0]))])  # ImA0
    for row in cartCovs :
      row.insert(0, 0.)
      row.insert(1, 0.)

    convert = convertComplexPars(cartVals, cartCovs, 'cart', 'polar, magSq', 3)
    polVals = convert[0]
    polCovs = convert[1]

  if printToScreen :
    # print cartesian values and covariance matrix to screen
    params = [('Re(A_0)', 'Im(A_0)'), ('Re(A_par)', 'Im(A_par)'),
        ('Re(A_perp)', 'Im(A_perp)')]
    if KSWave :
      params.append(('Re(A_S)', 'Im(A_S)'))

    print 'P2VV - INFO: convertAmplitudes: cartesian amplitudes:'
    printComplexParams(params, cartVals, cartCovs)

    # print polar values and covariance matrix to screen
    params = [('|A_0|^2', 'arg(A_0)'), ('|A_par|^2', 'arg(A_par)'),
        ('|A_perp|^2', 'arg(A_perp)')]
    if KSWave :
      params.append(('|A_S|^2', 'arg(A_S)'))

    print 'P2VV - INFO: convertAmplitudes: polar amplitudes:'
    printComplexParams(params, polVals, polCovs)

  return convert


###############################################################################

def printComplexParams(parameters, values, covariances = None) :
  """print values and covariances of complex parameters to screen
  """

  from math import sqrt

  # print values and their errors
  for parIter, par in enumerate(parameters) :
    print '  {0:<11} = {1:10.4g} +/- {2:9.4g}'.format(par[0],
        values[parIter][0], sqrt(covariances[2*parIter][2*parIter]))
    print '  {0:<11} = {1:10.4g} +/- {2:9.4g}'.format(par[1],
        values[parIter][1], sqrt(covariances[2*parIter+1][2*parIter+1]))

  if covariances :
    # print covariance matrix
    print '\n  covariance matrix:\n              ',
    for par in parameters :
      print '{0:>11}  {1:>11} '.format(par[0], par[1]),
    print

    for parIter0, par in enumerate(parameters) :
      print '  {0:<11} '.format(par[0]),
      for parIter1 in range(len(parameters)) :
        print '{0:11.3e}  {1:11.3e} '.format(\
            covariances[2*parIter0][2*parIter1],\
            covariances[2*parIter0][2*parIter1+1]),
      print

      print '  {0:<11} '.format(par[1]),
      for parIter1 in range(len(parameters)) :
        print '{0:11.3e}  {1:11.3e} '.format(\
            covariances[2*parIter0+1][2*parIter1],\
            covariances[2*parIter0+1][2*parIter1+1]),
      print

  print


###############################################################################

def convertComplexPars(values, covariances = None, inType = 'trans,cart',
    outType = 'trans,polar,magSq', normalize = 0) :
  """converts complex parameters (cartesian, polar, helicity and transversity)
  """

  # check input/output arguments
  if type(inType) is str : inType = inType.lower()
  else : inType = ''

  if type(outType) is str : outType = outType.lower()
  else : outType = ''

  if type(normalize) is not int : normalize = 0

  # determine operation
  inHel    = 'hel'   in inType
  inPolar  = 'polar' in inType
  inMagSq  =  inPolar  and 'magsq' in inType

  outHel   = 'hel'   in outType
  outPolar = 'polar' in outType
  outMagSq =  outPolar and 'magsq' in outType

  # check if we need to do something
  if inHel == outHel and inPolar == outPolar and inMagSq == outMagSq\
      and normalize < 1 :
    return (values, covariances)

  # check values
  nPars = len(values)
  if nPars < 1 :
    print "P2VV - ERROR: convertAmplitudes: no parameters specified"
    return None
  if nPars < 3 and inHel != outHel :
    print "P2VV - ERROR: convertAmplitudes: less than three parameters specified"
    return None
  if covariances and len(covariances) != 2 * nPars :
    print "P2VV - ERROR: convertAmplitudes: size of covariance matrix does not correspond to number of parameters"
    return None
  if normalize > nPars :
    print "P2VV - ERROR: convertAmplitudes: numbber of parameters for normalization is larger than total number of parameters"
    return None

  from array import array
  from math import sqrt, sin, cos, atan2 
  from ROOT import TMatrixD

  # define array of output values
  valsOut = list([list(par) for par in values])

  # get covariances
  if covariances :
    Jacobians = []
    # TODO: make this a TMatrixTSym
    covMatrixIn = TMatrixD(2 * nPars, 2 * nPars,
        array('d', [c_ij for c in covariances for c_ij in c[:2 * nPars]]))

  # convert to cartesian values
  if inPolar :
    if covariances :
      JacPolToCart = TMatrixD(2 * nPars, 2 * nPars)
      Jacobians.append(JacPolToCart)

    valsTemp = [par[:] for par in valsOut]
    for parIter, par in enumerate(valsTemp) :
      if inMagSq :
        ReA = sqrt(par[0]) * cos(par[1])
        ImA = sqrt(par[0]) * sin(par[1])
      else :
        ReA = par[0] * cos(par[1])
        ImA = par[0] * sin(par[1])

      valsOut[parIter][0] = ReA
      valsOut[parIter][1] = ImA

      if covariances :
        JacPolToCart[2 * parIter][2 * parIter]         = ReA / par[0]
        JacPolToCart[2 * parIter][2 * parIter + 1]     = -ImA
        JacPolToCart[2 * parIter + 1][2 * parIter]     = ImA / par[0]
        JacPolToCart[2 * parIter + 1][2 * parIter + 1] = ReA
        if inMagSq :
          JacPolToCart[2 * parIter][2 * parIter]     /= 2.
          JacPolToCart[2 * parIter + 1][2 * parIter] /= 2.

  # normalize
  if normalize > 0 :
    if covariances :
      JacNorm = TMatrixD(2 * nPars, 2 * nPars)
      Jacobians.append(JacNorm)

    valsTemp = [par[:] for par in valsOut]

    norm = 0.
    for parIter in range(normalize) :
      norm += valsTemp[parIter][0] * valsTemp[parIter][0]\
          + valsTemp[parIter][1] * valsTemp[parIter][1]
    norm = sqrt(norm)

    for parIter, par in enumerate(valsTemp) :
      valsOut[parIter][0] /= norm
      valsOut[parIter][1] /= norm

      if covariances :
        for newValIter in range(2 * parIter, 2 * parIter + 2) :
          for oldValIter in range(2 * nPars) :
            if oldValIter != newValIter :
              if oldValIter < 2 * normalize :
                oldVal1 = par[newValIter - 2 * parIter]
                oldVal2 = valsTemp[int(oldValIter / 2)][int(oldValIter % 2)]
                oldValNormSq = oldVal1 * oldVal2 / norm / norm
                JacNorm[newValIter][oldValIter] = -oldValNormSq / norm

              else :
                JacNorm[newValIter][oldValIter] = 0.

            else :
              if oldValIter < 2 * normalize :
                oldVal = par[newValIter - 2 * parIter]
                oldValNormSq = oldVal * oldVal / norm / norm
                JacNorm[newValIter][oldValIter] = (1. - oldValNormSq) / norm

              else :
                JacNorm[newValIter][oldValIter] = 1. / norm

  # convert to polar values
  if outPolar :
    if covariances :
      JacCartToPol = TMatrixD(2 * nPars, 2 * nPars)
      Jacobians.append(JacCartToPol)

    valsTemp = [par[:] for par in valsOut]
    for parIter, par in enumerate(valsTemp) :
      ASq = par[0] * par[0] + par[1] * par[1]
      Aph = atan2(par[1], par[0])

      if outMagSq :
        valsOut[parIter][0] = ASq
      else :
        valsOut[parIter][0] = sqrt(ASq)

      valsOut[parIter][1] = Aph

      if covariances :
        JacCartToPol[2 * parIter][2 * parIter]         = par[0]
        JacCartToPol[2 * parIter][2 * parIter + 1]     = par[1]
        JacCartToPol[2 * parIter + 1][2 * parIter]     = -par[1] / ASq
        JacCartToPol[2 * parIter + 1][2 * parIter + 1] =  par[0] / ASq
        if outMagSq :
          JacCartToPol[2 * parIter][2 * parIter]     *= 2.
          JacCartToPol[2 * parIter][2 * parIter + 1] *= 2.
        else :
          JacCartToPol[2 * parIter][2 * parIter]     /= sqrt(ASq)
          JacCartToPol[2 * parIter][2 * parIter + 1] /= sqrt(ASq)

  covariancesOut = None
  if covariances :
    # construct total transformation matrix
    fullJacobian = TMatrixD(2 * nPars, 2 * nPars, array('d',
       [1. if i == j else 0. for j in range(2*nPars) for i in range(2*nPars)]))
    for Jac in Jacobians :
      fullJacobianTemp = TMatrixD(fullJacobian)
      fullJacobian.Mult(Jac, fullJacobianTemp)

    # construct output covariance matrix
    # TODO: once covMatrixIn is a TMatrixTSym, use Similarity
    covMatrixTemp = TMatrixD(2 * nPars, 2 * nPars)
    covMatrixOut  = TMatrixD(2 * nPars, 2 * nPars)
    covMatrixTemp.Mult(fullJacobian, covMatrixIn)
    covMatrixOut.MultT(covMatrixTemp, fullJacobian)

    # get covariances from matrix
    covariancesOut = tuple([tuple([covMatrixOut[i][j]\
        for j in range(2 * nPars)]) for i in range(2 * nPars)])

  return (tuple([tuple(par) for par in valsOut]), covariancesOut)

