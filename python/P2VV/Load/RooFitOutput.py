"""controls RooFit output messages

Defines RooFit output streams and topics for P2VV
"""

import P2VV.RooFitDecorators
print "P2VV - INFO: RooFitOutput: setting RooFit output streams"
from ROOT import RooFit, RooMsgService

# get message service instance
msgServ = RooMsgService.instance()

# remove plotting streams with "INFO" as a minimum level
for stream in msgServ :
#    if stream.minLevel == RooFit.INFO : stream -= RooFit.Minimization
    if stream.minLevel == RooFit.INFO :
        stream -= RooFit.Plotting
        stream -= RooFit.Caching

