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

def registerMultiCatGen() :
  """registers an experimental fast(er) toy generator
  """

  from ROOT import RooMultiCatGenerator, RooNumGenFactory, RooNumGenConfig,\
      RooMsgService, RooFit

  RooMultiCatGenerator.registerSampler(RooNumGenFactory.instance())
  RooNumGenConfig.defaultConfig().methodND(False,True)\
      .setLabel("RooMultiCatGenerator")
  RooNumGenConfig.defaultConfig().methodND(False,True).Print()
  #RooMsgService.instance().addStream(RooFit.DEBUG,\
  #    RooFit.Topic(RooFit.Generation))


###############################################################################

def convertComplexPars(values, covariances = None, inType = 'trans,cart',
    outType = 'trans,polar', normalize = 0, magnSquared = True) :
  """converts complex paramters (cartesian, polar, helicity and transversity)
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
  outHel   = 'hel'   in outType
  outPolar = 'polar' in outType

  # check if we need to do something
  if inHel == outHel and inPolar == outPolar and normalize < 1 :
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
    covMatrixIn = TMatrixD(2 * nPars, 2 * nPars,
        array('d', [c_ij for c in covariances for c_ij in c[:2 * nPars]]))

  # convert to cartesian values
  if inPolar :
    if covariances :
      JacPolToCart = TMatrixD(2 * nPars, 2 * nPars)
      Jacobians.append(JacPolToCart)

    valsTemp = [par[:] for par in valsOut]
    for parIter, par in enumerate(valsTemp) :
      if magnSquared :
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
        if magnSquared :
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

      if magnSquared :
        valsOut[parIter][0] = ASq
      else :
        valsOut[parIter][0] = sqrt(ASq)

      valsOut[parIter][1] = Aph

      if covariances :
        JacCartToPol[2 * parIter][2 * parIter]         = par[0]
        JacCartToPol[2 * parIter][2 * parIter + 1]     = par[1]
        JacCartToPol[2 * parIter + 1][2 * parIter]     = -par[1] / ASq
        JacCartToPol[2 * parIter + 1][2 * parIter + 1] =  par[0] / ASq
        if magnSquared :
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
    covMatrixTemp = TMatrixD(2 * nPars, 2 * nPars)
    covMatrixOut  = TMatrixD(2 * nPars, 2 * nPars)
    covMatrixTemp.Mult(fullJacobian, covMatrixIn)
    covMatrixOut.MultT(covMatrixTemp, fullJacobian)

    # get covariances from matrix
    covariancesOut = tuple([tuple([covMatrixOut[i][j]\
        for j in range(2 * nPars)]) for i in range(2 * nPars)])

  return (tuple([tuple(par) for par in valsOut]), covariancesOut)

