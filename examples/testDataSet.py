from P2VV.RooFitWrappers import *
ws = RooObject ( workspace = 'myws' )

x = RealVar('x',Value = 1 )
y = RealVar('y',Value = 2 )
z = RealVar('z',Value = 3 )

l = RooArgList( [ x._var,y._var,z._var ] )
l.Print()
s = RooArgSet( [ x._var,y._var,z._var ] )
s.Print()

a = RooDataSet('foo','foo',[ i._var for i in [ x,y,z] ] )
a.Print()

b = RooDataSet('bar','bar',l )
b.Print()

c = RooDataSet('bar','bar',s )
c.Print()
