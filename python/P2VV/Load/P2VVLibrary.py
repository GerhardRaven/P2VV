"""load the P2VV library

Load the P2VV library, assuming that the library directory is in $LD_LIBRARYPATH
"""

print "P2VV - INFO: loading Cintex library"
import PyCintex
PyCintex.Cintex.Enable()
PyCintex.loadDictionary("P2VVDict")
