###############################################################################
## P2VV: common tools                                                        ##
##                                                                           ##
## authors:                                                                  ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                ##
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

