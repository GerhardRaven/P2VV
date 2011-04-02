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

def convertAmplitudes(amplitudes, covariances = None, inType = 'trans,cart',
    outType = 'trans,polar', normalize = 3, magnSquared = True) :
  """converts decay amplitudes (cartesian, polar, helicity and transversity)
  """

  # check input/output arguments
  if type(inType) is str : inType = inType.lower()
  else : inType = ''

  if type(outType) is str : outType = outType.lower()
  else : outType = ''

  if type(normalize) is not int : normalize = 0

  # determine operation
  inHel   = 'hel'  in inType
  inCart  = 'cart' in inType
  outHel  = 'hel'  in outType
  outCart = 'cart' in outType

  # check if we need to do something
  if inHel == outHel and inCart == outCart and normalize < 1 :
    return (amplitudes, covariances)

  # check amplitudes
  nAmps = len(amplitudes)
  if nAmps < 1 :
    print "P2VV - ERROR: convertAmplitudes: no amplitudes specified"
    return None
  if nAmps < 3 and inHel != outHel :
    print "P2VV - ERROR: convertAmplitudes: less than three amplitudes specified"
    return None
  if covariances and len(covariances) != 2 * nAmps :
    print "P2VV - ERROR: convertAmplitudes: size of covariance matrix does not correspond to number of amplitudes"
    return None

  from array import array
  from math import sqrt, sin, cos, atan2 
  from ROOT import TMatrixD

  # define array of output amplitudes
  ampsOut = list([list(amp) for amp in amplitudes])

  # get covariances
  if covariances :
    transMatrices = []
    covMatrixIn = TMatrixD(2 * nAmps, 2 * nAmps,
        array('d', [c_ij for c in covariances for c_ij in c[:2 * nAmps]]))

  # convert to cartesian amplitudes
  if not inCart :
    if covariances :
      TPol2Cart = TMatrixD(2 * nAmps, 2 * nAmps)

    ampsTemp = [amp[:] for amp in ampsOut]
    for ampIter, amp in enumerate(ampsTemp) :
      if magnSquared :
        ReA = sqrt(amp[0]) * cos(amp[1])
        ImA = sqrt(amp[0]) * sin(amp[1])
      else :
        ReA = amp[0] * cos(amp[1])
        ImA = amp[0] * sin(amp[1])

      ampsOut[ampIter][0] = ReA
      ampsOut[ampIter][1] = ImA

      if covariances :
        TPol2Cart[2 * ampIter][2 * ampIter]         = 0.5 * ReA / amp[0]
        TPol2Cart[2 * ampIter][2 * ampIter + 1]     = -ImA
        TPol2Cart[2 * ampIter + 1][2 * ampIter]     = 0.5 * ImA / amp[0]
        TPol2Cart[2 * ampIter + 1][2 * ampIter + 1] = ReA
        if not magnSquared :
          TPol2Cart[2 * ampIter][2 * ampIter]     /= amp[0]
          TPol2Cart[2 * ampIter + 1][2 * ampIter] /= amp[0]

        transMatrices.append(TPol2Cart)

  # convert to polar amplitudes
  if not outCart :
    if covariances :
      TCart2Pol = TMatrixD(2 * nAmps, 2 * nAmps)

    ampsTemp = [amp[:] for amp in ampsOut]
    for ampIter, amp in enumerate(ampsTemp) :
      ASq = amp[0] * amp[0] + amp[1] * amp[1]
      Aph = atan2(amp[1], amp[0])

      if magnSquared :
        ampsOut[ampIter][0] = ASq
      else :
        ampsOut[ampIter][0] = sqrt(ASq)

      ampsOut[ampIter][1] = Aph

      if covariances :
        TCart2Pol[2 * ampIter][2 * ampIter]         = 2. * amp[0]
        TCart2Pol[2 * ampIter][2 * ampIter + 1]     = 2. * amp[1]
        TCart2Pol[2 * ampIter + 1][2 * ampIter]     = -amp[1] / ASq
        TCart2Pol[2 * ampIter + 1][2 * ampIter + 1] =  amp[0] / ASq

        transMatrices.append(TCart2Pol)

  covariancesOut = None
  if covariances :
    # construct total transformation matrix
    TMatrix = TMatrixD(2 * nAmps, 2 * nAmps, array('d',
       [1. if i == j else 0. for j in range(2*nAmps) for i in range(2*nAmps)]))
    for matrix in transMatrices :
      TMatrixTemp = TMatrixD(TMatrix)
      TMatrix.Mult(matrix, TMatrixTemp)

    # construct output covariance matrix
    covMatrixTemp = TMatrixD(2 * nAmps, 2 * nAmps)
    covMatrixOut  = TMatrixD(2 * nAmps, 2 * nAmps)
    covMatrixTemp.Mult(TMatrix, covMatrixIn)
    covMatrixOut.MultT(covMatrixTemp, TMatrix)

    # get covariances from matrix
    covariancesOut = tuple([tuple([covMatrixOut[i][j]\
        for j in range(2 * nAmps)]) for i in range(2 * nAmps)])

  return (tuple([tuple(amp) for amp in ampsOut]), covariancesOut)

