"""controls RooFit output messages

Defines RooFit output streams and topics for P2VV
"""

print "P2VV - INFO: RooFitOutput: setting RooFit output streams"
from ROOT import RooFit, RooMsgService

# get message service instance
msgServ = RooMsgService.instance()

# remove plotting streams with "INFO" as a minimum level
for stream in msgServ :
    if stream.minLevel == RooFit.INFO : stream -= RooFit.Plotting

## remove all output streams
#for stream in range(msgServ.numStreams()) : msgServ.deleteStream(stream)
#
## add default streams for P2VV
#msgServ.addStream(RooFit.PROGRESS)
#msgServ.addStream(RooFit.INFO, RooFit.Topic(RooFit.Minimization))
#msgServ.getStream(1) += ( RooFit.Plotting
#                        , RooFit.Fitting
#                        , RooFit.Eval
#                        , RooFit.Caching
#                        , RooFit.ObjectHandling
#                        , RooFit.InputArguments
#                        , RooFit.DataHandling
#                        , RooFit.NumIntegration )

