


from RooFitWrappers import *
ws = RooObject ( workspace = 'myws' )

x = RealVar('x',Value = 1 )
y = RealVar('y',Value = 2 )
z = RealVar('z',Value = 3 )



a = RooDataSet('foo','foo',[ i._var for i in [ x,y,z] ] )
print a
a.Print()

