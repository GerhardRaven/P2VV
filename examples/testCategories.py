from RooFitWrappers import *
ws = RooObject ( workspace = 'myws' )
x = Category('mystate', States = [ 'aap','noot','mies','vuur' ] )
y = Category('mystate2', States = [ 'x','y' ] )
z = SuperCategory('sc',[x,y])
#z.Print("V")
for i in x.states() : print i
a = MappedCategory('mc',x,{'one':('aap','noot')
                          ,'two':('mies', ) } )

for i in a.states() : print i
#a.Print("V")


for i in x.states() :
    x.setLabel(i)
    print x.getLabel(),a.getLabel()
