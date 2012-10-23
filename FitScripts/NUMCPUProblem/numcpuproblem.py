from ROOT import TFile,TCanvas, RooFit, RooArgSet

file = TFile('sigws.root')
sigworkspace= file.Get('sigworkspace')
file.Close()

sigpdf = sigworkspace.pdf('sig_t')
sigdata = sigworkspace.data('DecayTree')
t =  sigworkspace.var('time')
st = sigworkspace.var('sigmat')
SFconstraint = sigworkspace.pdf('timeResSF_constraint')

from ROOT import RooMsgService
RooMsgService.instance().getStream(1).removeTopic(RooFit.Caching)
RooMsgService.instance().getStream(1).removeTopic(RooFit.Eval)
                                                              
cond = RooArgSet(st)
cons = RooArgSet(SFconstraint)
siglifetimeresult = sigpdf.fitTo(sigdata
                                 , RooFit.ConditionalObservables(cond)
                                 , RooFit.ExternalConstraints(cons)
                                 , RooFit.NumCPU(16)
#                                 , RooFit.Optimize(2)
                                 )
