from RooFitWrappers import RooObject
import RooFitDecorators 

class Mapping(object):
    def __init__(self, mapping, dataset):
        self._mapping = dict([(k.GetName(), v) for k, v in mapping.iteritems()])
        self._dataset = dataset
        RooObject._ws.put( dataset )

        rootVars = dict( (var.GetName(),var) for var in dataset.get() )
        for v, n in mapping.iteritems():
            rv = rootVars[n] 
            # Test if they are the same
            # if getattr fails, add a corresponding method for the missing type....
            __test = getattr( self, '__test%s' % type(rv).__name__ )
            __test(v,rv)

            RooObject._ws._mappings[v.GetName()] = n

            # Delete currently associated RooRealVar and put the one from the data in
            RooObject._ws._objects.pop(v.GetName())
            RooObject._ws._objects[n] = rv
            # Variables from a DataSet must be observable??
            rv._observable = True
            v._var = rv

    def __testRooCategory(self, cat, rooCat):
        states = dict( (cat.GetName(),cat.GetVal()) for cat in self )
        assert cat.observable() == True
        assert cat['Name']
        assert cat['States'] == states

    def __testRooRealVar(self, var, rooVar):
        assert var.observable() == True
        assert var['MinMax'] == (rooVar.getMin(), rooVar.getMax())
        assert var['Unit'] == rooVar.getUnit()

    def __getitem__(self, k):
        if type(k) != str: k = k.GetName()
        if k == 'DataSet':
            return self._dataset
        elif k in self._mapping:
            return self._mapping[k]
        else:
            raise KeyError('This mapping does not contain %k' % k)

    def __contains__(self, k):
        return k == 'DataSet' or k.GetName() in self._mapping
