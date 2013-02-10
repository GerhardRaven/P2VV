"""load the P2VV library

Load the P2VV library, assuming that $P2VVPATH/lib is in $LD_LIBRARYPATH
"""

print "P2VV - INFO: loading P2VV library"
from ROOT import gSystem
gSystem.Load('libRooFit')
gSystem.Load('libP2VV')


