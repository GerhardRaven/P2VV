from ROOT import *
 
gSystem.Load('libP2VV')
RooMsgService.instance().addStream(RooFit.DEBUG,RooFit.Topic(RooFit.Generation))

w = RooWorkspace()
tag = w.factory('tag[unmixed=1,untagged=0,mixed=-1]')
t   = w.factory('t[0,10]')
tau = w.factory('tau[1.5,1.0,2.0]')
dm  = w.factory('dm[0.5,0.4,0.6]')
eps = w.factory('eps[0.333,0,1]')
one = w.factory('ConstVar::one(1)')
zero = w.factory('ConstVar::zero(0)')
pdf = w.factory('TrivialTagDecay::pdf(t,tag,tau,zero,dm,eps,one,zero,one,zero,TruthModel(t),SingleSided)')
data = pdf.generate(RooArgSet(tag,t),100000)
data.table(tag).Print('V')

I = pdf.createIntegral(RooArgSet(t))
Id = dict()
for i in [0,1,-1] :
    tag.setIndex(i)
    Id[ tag.getLabel() ]  = I.getVal()
print Id
norm = sum( Id.itervalues() )
for k,v in Id.iteritems() : Id[k] *= float(1)/norm
print Id['untagged']
print Id['mixed']/(Id['mixed']+Id['unmixed'])

c = TCanvas()
c.Divide(1,3)

for (i,(tv,label)) in enumerate( [ (1,'unmixed'),(-1,'mixed'),(0,'untagged') ] ) :
   c.cd(1+i)
   f = t.frame(RooFit.Title(label))
   data.plotOn(f, RooFit.Cut('tag==%s'%tv))
   pdf.plotOn(f, RooFit.Slice(tag,label))
   f.Draw()

pdf.fitTo(data)

d = TCanvas()
d.Divide(1,3)

for (i,(tv,label)) in enumerate( [ (1,'unmixed'),(-1,'mixed'),(0,'untagged') ] ) :
   d.cd(1+i)
   f = t.frame(RooFit.Title(label))
   data.plotOn(f, RooFit.Cut('tag==%s'%tv))
   pdf.plotOn(f, RooFit.Slice(tag,label))
   f.Draw()
