"""load the P2VV library

Load the P2VV library, assuming that the library directory is in $LD_LIBRARYPATH
"""

print "P2VV - INFO: loading Cintex library"
import ROOT
ROOT.gSystem.Load("libCintex.so")
ROOT.Cintex.Enable()
ROOT.gSystem.Load("libP2VVDict.so")
