from ROOT import *
 
gSystem.Load('libP2VV')
w = RooWorkspace()
tag = w.factory('tag[unmixed=1,untagged=0,mixed=-1]')
t   = w.factory('t[0,20]')
tau = w.factory('tau[1.5]')
dm  = w.factory('dm[0.5]')
eps = w.factory('eps[0.5,0,1]')
one = w.factory('ConstVar::one(1)')
zero = w.factory('ConstVar::zero(0)')
pdf = w.factory('TrivialTagDecay::pdf(t,tag,tau,zero,dm,eps,one,zero,one,zero,TruthModel(t),SingleSided)')

