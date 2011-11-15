
from RooFitWrappers import *
ws = RooObject ()
ws.setWorkspace( RooWorkspace('ws') )
x = Category('mystate', States = [ 'aap','noot','mies' ] )
y = Category('mystate2', States = [ 'x','y' ] )
z = SuperCategory('sc',[x,y])
z.Print("V")
for i in z.states() : print i
a = MappedCategory('mc',x,{'one':('aap','noot')
                          ,'two':('mies', ) } )

for i in a.states() : print i
a.Print("V")
