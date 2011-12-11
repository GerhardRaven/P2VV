"""registers MultiCatGen

Registers MultiCatGen, an experimental fast(er) toy generator for RooFit
"""

from P2VVLoad import P2VVLibrary
print "P2VV - INFO: MultiCatGen: registering generator"
from ROOT import RooMultiCatGenerator, RooNumGenFactory, RooNumGenConfig, RooMsgService, RooFit
RooMultiCatGenerator.registerSampler(RooNumGenFactory.instance())
RooNumGenConfig.defaultConfig().methodND( False, True ).setLabel('RooMultiCatGenerator')
RooNumGenConfig.defaultConfig().methodND( False, True ).Print()


