"""load the P2VV library

Load the P2VV library, assuming that $P2VVPATH/lib is in $LD_LIBRARYPATH
"""

## Add this in case a reflex dictionary is used.
from ROOT import gSystem
gSystem.Load('libCintex')
from ROOT import Cintex
Cintex.Enable()

print "P2VV - INFO: loading P2VV library"
from ROOT import gSystem
gSystem.Load('libRooFit')
gSystem.Load('libP2VV')


