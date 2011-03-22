###############################################################################
## P2VV: common tools                                                        ##
##                                                                           ##
## authors:                                                                  ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                ##
##                                                                           ##
###############################################################################


###################################################
## function that loads the P2VV library          ##
## (assumes $P2VVPATH/lib is in $LD_LIBRARYPATH) ##
###############################################################################
def loadP2VVLib() :
  from ROOT import gSystem
  gSystem.Load("libP2VV")


##################################################
## register experimental fast(er) toy generator ##
###############################################################################
def registerMultiCatGen() :
    from ROOT import RooMultiCatGenerator, RooNumGenFactory, RooNumGenConfig,\
        RooMsgService, RooFit

    RooMultiCatGenerator.registerSampler(RooNumGenFactory.instance())
    RooNumGenConfig.defaultConfig().methodND(False,True)\
        .setLabel("RooMultiCatGenerator")
    RooNumGenConfig.defaultConfig().methodND(False,True).Print()
    #RooMsgService.instance().addStream(RooFit.DEBUG,\
    #    RooFit.Topic(RooFit.Generation))

