"""load the P2VV library

Load the P2VV library, assuming that the library directory is in $LD_LIBRARYPATH
"""

from ROOT import gSystem

#if not 'libRooFit' in gSystem.GetLibraries() :
#    print "P2VV - INFO: loading Cintex library"
#    gSystem.Load('libCintex')
#
#from ROOT import Cintex
#Cintex.Enable()

if not 'libRooFit' in gSystem.GetLibraries() :
    print "P2VV - INFO: loading RooFit library"
    gSystem.Load('libRooFit')

if not 'libP2VV' in gSystem.GetLibraries() :
    print "P2VV - INFO: loading P2VV library"
    try:
        from ROOT import gCling
        gSystem.Load("libGui")
        import os
        gCling.AddIncludePath(os.path.join(os.environ['P2VVROOT'], 'include'))
    except ImportError:
        pass
    gSystem.Load('libP2VV')
