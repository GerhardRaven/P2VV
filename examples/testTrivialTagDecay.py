from ROOT import *
 
gSystem.Load('libP2VV')
#RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))

w = RooWorkspace()
tag  = w.factory('tag[unmixed=1,untagged=0,mixed=-1]')
t    = w.factory('t[0,10]')
tau  = w.factory('tau[1.5,1.0,2.0]')
dm   = w.factory('dm[18.,16.,20.]')
eps  = w.factory('eps[0.8,0,1]')
one  = w.factory('ConstVar::one(1)')
zero = w.factory('ConstVar::zero(0)')
dgamma = w.factory('dgamma[-0.1,-0.8,0.8]')
Dilution = w.factory('Dilution[0.3,0,1]')
pdf  = w.factory('TrivialTagDecay::pdf(t,tag,tau,dgamma,dm,eps,one,zero,Dilution,zero,TruthModel(t),SingleSided)')
# opdf  = w.factory('BDecay::pdf(t,tau,dgamma,dm,one,zero,Dilution,zero,TruthModel(t),SingleSided)')

data = pdf.generate(RooArgSet(tag,t),100000)
data.table(tag).Print('V')

pdf.fitTo(data)

d = TCanvas()
d.Divide(1,5)

for (i,(tv,label)) in enumerate( [ (1,'unmixed'),(-1,'mixed'),(0,'untagged') ] ) :
   d.cd(2+i)
   f = t.frame(RooFit.Title(label))
   data.plotOn(f, RooFit.Cut('tag==%s'%tv))
   pdf.plotOn(f, RooFit.Slice(tag,label))
   f.Draw()

d.cd(5)
f = t.frame(RooFit.Title('Asymmetry'))
data.plotOn(f,RooFit.Asymmetry(tag))
pdf.plotOn(f,RooFit.Asymmetry(tag))
f.Draw()

d.cd(1)
f = t.frame(RooFit.Title('All'))
data.plotOn(f)
pdf.plotOn(f)
f.Draw()
