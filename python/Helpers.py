from RooFitWrappers import RooObject

class Mapping(object):
    def __init__(self, mapping, dataset):
        self._mapping = dict([(k.GetName(), v) for k, v in mapping.iteritems()])
        self._dataset = dataset

        imp = getattr(RooObject._ws, 'import')
        imp(dataset)

        dataVars = dataset.get()
        it = dataVars.iterator()
        rootVars = {}
        while True:
            var = it.Next()
            if not var:
                break
            rootVars[var.GetName()] = var
            
        for v, n in mapping.iteritems():
            rv = rootVars[n] 
            # Test if they are the same
            if rv.IsA().GetName().find('RooRealVar') != -1:
                self.__testRealVar(v, rv)
            if rv.IsA().GetName().find('RooCategory') != -1:
                self.__testCategory(v, rv)

            RooObject._ws._mappings[v.GetName()] = n

            # Delete currently associated RooRealVar and put the one from the data in
            RooObject._ws._objects.pop(v.GetName())
            RooObject._ws._objects[n] = rv
            # Variables from a DataSet must be observable??
            rv._observable = True
            v._var = rv

    def __testCategory(self, cat, rooCat):
        states = {}
        it = self.typeIterator()
        while True:
            cat = it.Next()
            if not cat:
                break
            states[cat.GetName()] = cat.getVal()
        assert cat.observable() == True
        assert cat['Name']
        assert cat['States'] == states

    def __testRealVar(self, var, rooVar):
        assert var.observable() == True
        assert var['MinMax'] == (rooVar.getMin(), rooVar.getMax())
        assert var['Unit'] == rooVar.getUnit()

    def __getitem__(self, k):
        if type(k) != str:
            k = k.GetName()
        if k == 'DataSet':
            return self._dataset
        elif k in self._mapping:
            return self._mapping[k]
        else:
            raise KeyError('This mapping does not contain %k' % k)

    def __contains__(self, k):
        if k == 'DataSet':
            return True
        else:
            return k.GetName() in self._mapping
