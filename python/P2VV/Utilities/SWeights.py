###########################################################################################################################################
## Utilities.SWeights: P2VV utilities for computing sWeights                                                                             ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

class SData( object ) :
    def __init__( self, **kwargs ) :
        # get input arguments
        def getKwArg( keyword, member, kwargs ) :
            if keyword in kwargs : setattr( self, '_' + member, kwargs.pop(keyword) )
            else : raise KeyError, 'P2VV - ERROR: SData.__init__(): key %s not found in input arguments' % keyword
        getKwArg( 'Name', 'name',      kwargs )
        getKwArg( 'Data', 'inputData', kwargs )
        getKwArg( 'Pdf',  'pdf',       kwargs )

        # initialize dictionary for weighted data sets per specie
        self._data = dict()

        # get yields and observables
        self._yields = [ par for par in self._pdf.Parameters() if par.getAttribute('Yield') ]
        self._observables = self._pdf.Observables()

        # calculate sWeights
        from ROOT import RooStats, RooArgList, RooSimultaneous
        if isinstance( self._pdf._var, RooSimultaneous ) and kwargs.pop( 'Simultaneous', True ) :
            # split data set in categories of the simultaneous PDF
            splitCat        = self._pdf.indexCat()
            splitData       = self._inputData.split(splitCat)
            self._sPlots    = [ ]
            self._sDataSets = [ ]
            sDataVars       = None
            from ROOT import RooFormulaVar
            for splitCatState in splitCat:
                # calculate sWeights per category
                cat = splitCatState.GetName()
                data = splitData.FindObject(cat)

                origYieldVals = [ ( par.GetName(), par.getVal(), par.getError() ) for par in self._yields if par.GetName().endswith(cat) ]
                self._sPlots.append(  RooStats.SPlot( self._name + '_sData_' + cat, self._name + '_sData_' + cat
                                                     , data, self._pdf.getPdf(cat)
                                                     , RooArgList( par._var for par in self._yields if par.GetName().endswith(cat) ) )
                                   )
                self._sDataSets.append( self._sPlots[-1].GetSDataSet() )

                print 'P2VV - INFO: SData.__init__(): yields category %s:' % cat
                print '    original:',
                for vals in origYieldVals : print '%s = %.2f +/- %.2f  ' % vals,
                print '\n    new:     ',
                for par in self._yields :
                    if par.GetName().endswith(cat) : print '%s = %.2f +/- %.2f  ' % ( par.GetName(), par.getVal(), par.getError() ),
                print

                # add column for splitting category/categories (was removed when data set was split)
                # FIXME: How can we do this more generally? These are special cases and it might go wrong here...
                from ROOT import RooSuperCategory
                splitCat.setLabel(cat)
                __dref = lambda o : o._target_() if hasattr(o,'_target_') else o
                if isinstance( __dref(splitCat), RooSuperCategory ) :
                    for fundCat in splitCat.inputCatList() :
                        if not self._sDataSets[-1].get().find( fundCat.GetName() ) : self._sDataSets[-1].addColumn(fundCat)
                elif splitCat.isFundamental() and not self._sDataSets[-1].get().find( splitCat.GetName() ) :
                    self._sDataSets[-1].addColumn(splitCat)

                # add general sWeight and PDF value columns (it must be possible to simplify this...)
                # FIXME: in some cases "par.GetName().strip( '_' + cat )" goes wrong:
                # use "par.GetName()[ : par.GetName().find(cat) - 1 ]" instead
                # (case: 'N_bkgMass_notExclBiased'.strip('_notExclBiased') --> 'N_bkgM' ?!!!!!!)
                weightVars = [ (  RooFormulaVar( par.GetName()[ : par.GetName().find(cat) - 1 ] + '_sw', '', '@0'
                                                , RooArgList( self._sDataSets[-1].get().find( par.GetName() + '_sw' ) ) )
                                , RooFormulaVar( 'L_' + par.GetName()[ : par.GetName().find(cat) - 1 ], '', '@0'
                                                , RooArgList( self._sDataSets[-1].get().find( 'L_' + par.GetName() ) ) )
                               ) for par in self._yields if par.GetName().endswith(cat)
                             ]

                for weight, pdfVal in weightVars :
                    self._sDataSets[-1].addColumn(weight)
                    self._sDataSets[-1].addColumn(pdfVal)

                if not sDataVars :
                    # get set of variables in data
                    sDataVars = self._sDataSets[-1].get()
                    for par in self._yields :
                        if cat in par.GetName() :
                            sDataVars.remove( sDataVars.find( par.GetName() + '_sw' ) )
                            sDataVars.remove( sDataVars.find( 'L_' + par.GetName()  ) )

            # merge data sets from categories
            from ROOT import RooDataSet
            self._sData = RooDataSet( self._name + '_splotdata', self._name + '_splotdata', sDataVars )
            for data in self._sDataSets : self._sData.append(data)

        else :
            # calculate sWeights with full data set
            if isinstance( self._pdf._var, RooSimultaneous ) :
                print 'P2VV - WARNING: SData.__init__(): computing sWeights with a simultaneous PDF'
            self._sPlot = RooStats.SPlot( self._name + '_splotdata', self._name + '_splotdata', self._inputData, self._pdf._var
                                         , RooArgList( par._var for par in self._yields ) )
            self._sData = self._sPlot.GetSDataSet()

        # check keyword arguments
        if kwargs : raise KeyError, 'P2VV - ERROR: SData.__init__(): got unknown keywords %s for %s' % ( kwargs, type(self) )

    def usedObservables( self ) : return self._observables
    def components( self )      : return [ y.GetName()[2:] for y in self._yields ]
    def Pdf( self )             : return self._pdf

    def Yield( self, Component ) :
        yName = 'N_%s' % Component
        for y in self._yields :
            if y.GetName() == yName : return y.getVal()

        raise KeyError, 'P2VV - ERROR: SData.__init__(): unknown component %s' % Component

    def data( self, Component = None ) :
        if not Component : return self._sData

        if Component not in self._data :
            # check if component exists
            yName = 'N_%s' % Component
            if not any( yName in y.GetName() for y in self._yields ) :
                raise KeyError, 'P2VV - ERROR: SData.__init__(): unknown component: %s' % Component
            wName = '%s_sw' % yName
            if wName not in [ w.GetName() for w in self._sData.get() ] :
                raise KeyError, 'no weight in dataset for component %s' % Component

            # create weighted data set
            dName = '%s_weighted_%s' % ( self._sData.GetName(), Component )
            from ROOT import RooDataSet
            from P2VV.ROOTDecorators import ROOTversion
            if ROOTversion[0] <= 5 and ROOTversion[1] <= 34 and ROOTversion[2] < 2 :
                self._data[Component] = RooDataSet( dName, dName, self._sData.get(), Import = self._sData, WeightVar = ( wName ) )
            else :
                self._data[Component] = RooDataSet( dName, dName, self._sData.get(), Import = self._sData, WeightVar = ( wName, True ) )

        return self._data[Component]


def createSData( **kwargs ) :
    # make sweighted dataset using J/psi phi mass
    Observables = kwargs.pop('Observables')
    Data        = kwargs.pop('Data')
    FitOpts     = kwargs.pop('FitOpts')
    Components  = kwargs.pop('Components')
    Name        = kwargs.pop('Name')
    from P2VV.RooFitWrappers import buildPdf
    pdf = buildPdf( Components, Observables = Observables, Name= Name + '_splot_pdf')
    c_m = pdf.fitTo(Data,**FitOpts)
    c_state = dict( ( p, p.isConstant() ) for p in pdf.Parameters() )
    for p in pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
    sdata =  SData(Pdf = pdf, Data = Data, Name = Name + '_sdata')
    for p,c in c_state.iteritems() : p.setConstant(c)
    return sdata
