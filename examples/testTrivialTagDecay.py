from ROOT import *
 
gSystem.Load('libP2VV')
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))

w = RooWorkspace()
tag = w.factory('tag[unmixed=1,untagged=0,mixed=-1]')
t   = w.factory('t[0,20]')
tau = w.factory('tau[1.5,1.0,2.0]')
dm  = w.factory('dm[0.5,0.4,0.6]')
eps = w.factory('eps[0.5,0,1]')
one = w.factory('ConstVar::one(1)')
zero = w.factory('ConstVar::zero(0)')
pdf = w.factory('TrivialTagDecay::pdf(t,tag,tau,zero,dm,eps,one,zero,one,zero,TruthModel(t),SingleSided)')
pdf.Print("T")
data = pdf.generate(RooArgSet(tag,t),100)
data.table(tag).Print('V')
pdf.fitTo(data)
