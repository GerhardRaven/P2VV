"""load the P2VV library

Load the P2VV library, assuming that the library directory is in $LD_LIBRARYPATH
"""

from ROOT import gSystem
print "P2VV - INFO: loading P2VV library"

# enable support for reflex dictionaries
gSystem.Load('libCintex')
from ROOT import Cintex
Cintex.Enable()

# load required libraries
gSystem.Load('libRooFit')
gSystem.Load('libP2VVDict')
gSystem.Load('libP2VVLib')
