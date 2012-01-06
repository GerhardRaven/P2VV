###########################################################################################################################################
## P2VVParameterizations.FlavourTagging: Flavour tagging parameters                                                                      ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin


class TaggingParams ( _util_parse_mixin ):
    def __init__( self, **kwargs ) :
        self._numTagCats = kwargs.pop( 'NumTagCats', 1 )
        self._dilutions  = kwargs.pop('Dilutions')
        self._ADilWTags  = kwargs.pop('ADilWTags')
        self._CEvenOdds  = kwargs.pop('CEvenOdds')
        if self._numTagCats > 1 : self._tagCatCoefs = kwargs.pop('TagCatCoefs')

    def __getitem__( self, kw ) :
        def raiseError(kw) : raise RuntimeError( 'TaggingParams.__getitem__(\'%s\'): need to specify tagging category' % kw )
        if kw == 'dilution'                  : return self._dilutions[0]     if self._numTagCats == 1 else raiseError(kw)
        if kw == 'ADilWTag'                  : return self._ADilWTags[0]     if self._numTagCats == 1 else raiseError(kw)
        if kw == 'CEvenOdd'                  : return self._CEvenOdds[0]     if self._numTagCats == 1 else raiseError(kw)
        if kw in [ 'avgCEven',  'avgCOdd' ]  : return self._CEvenOdds[0][kw] if self._numTagCats == 1 else raiseError(kw)
        if kw in [ 'avgCEvens', 'avgCOdds' ] : return [ CEvenOdd[ kw[:-1] ] for CEvenOdd in self._CEvenOdds ]

        return getattr( self, '_' + kw )

class Trivial_TaggingParams( TaggingParams ) :
    def __init__( self ) :
        from RooFitWrappers import ConstVar
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__( self
                              , Dilutions = [ ConstVar( Name = 'Tagging dilution', Value = 0. ) ]
                              , ADilWTags = [ ConstVar( Name = 'zero', Value = 0. ) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class WTag_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar,ConstVar
        self._parseArg( 'wTag', kwargs, Title = 'B wrong tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__( self
                              , Dilutions = [ FormulaVar(  'tagDilution', '1. - 2*@0 '
                                                         , [ self._wTag ], Title = 'Tagging dilution' ) ]
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class LinearEstWTag_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar,ConstVar
        self._parseArg( 'estWTag', kwargs, Title = 'Estimated wrong tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
        self._parseArg( 'p0', kwargs, Title = 'p0  tagging parameter',    Value = 0.384, Constant = True, MinMax = ( 0., 0.5 ) )
        self._parseArg( 'p1', kwargs, Title = 'p1  tagging parameter',    Value = 1.037, Constant = True, MinMax = ( 0.8, 1.2 ) )
        self._parseArg( 'avgEstWTag', kwargs, Title = 'Average estimated wrong tag probability',    Value = 0.379, Constant = True, MinMax = ( 0., 0.5 ) )
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__( self
                              , Dilutions = [ FormulaVar(  'tagDilution', '1. - 2. * ( @2 + @3 * ( @0 - @1 ) ) '
                                                         , [ self._estWTag, self._avgEstWTag, self._p0, self._p1 ], Title = 'Calibrated Per-Event Dilution' ) ]
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class Dilution_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'dilution', kwargs, Title = 'Tagging dilution', Value = 0., MinMax = ( 0., 1. ) )
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        from RooFitWrappers import ConstVar
        TaggingParams.__init__( self
                              , Dilutions = [ self._dilution ]
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )
        self._check_extraneous_kw( kwargs )

class WTagsCoefAsyms_TaggingParams( TaggingParams ) :
    """flavour tagging parameters with wrong-tag parameters and asymmetry coefficients
    """

    def __init__( self, **kwargs ) :
        # get wrong-tag parameters
        if 'wTag' in kwargs and 'wTagBar' in kwargs :
            # wrong tag + wrong tag bar
            self._parseArg( 'wTag',    kwargs, Title = 'B wrong tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
            self._parseArg( 'wTagBar', kwargs, Title = 'Bbar wrong tag probability', Value = 0.25, MinMax = ( 0., 0.5 ) )
        else :
            # average wrong tag + wrong tag asymmetry
            self._parseArg( 'WTag',  kwargs, Title = 'Average wrong tag probability',   Value = 0.25, MinMax = (  0., 0.5 ) )
            self._parseArg( 'AWTag', kwargs, Title = 'Wrong tag probability asymmetry', Value = 0.,   MinMax = ( -1., 1.  ) )

        # get average even and odd coefficients
        if 'CEvenOdd' in kwargs :
            # a coefficients object is specified
            CEvenOdd = kwargs.pop('CEvenOdd')

        else :
            from RooFitWrappers import RooObject
            if 'AvgCEven' in kwargs and 'AvgCOdd' in kwargs :
                # (values for the) coefficients are specified
                avgCEven = kwargs.pop('AvgCEven')
                avgCOdd  = kwargs.pop('AvgCOdd')
                if not isinstance( avgCEven, RooObject ) and not isinstance( avgCOdd, RooObject ) :
                    avgCOdd  /= avgCEven
                    avgCEven  = 1.

            elif 'AProd' in kwargs and 'ANorm' in kwargs :
                # values for the production and normalization asymmetries are specified
                self._AProdVal = kwargs.pop('AProd')
                self._ANormVal = kwargs.pop('ANorm')
                if isinstance( self._AProdVal, RooObject ) : self._AProdVal = self._AProdVal.getVal()
                if isinstance( self._ANormVal, RooObject ) : self._ANormVal = self._ANormVal.getVal()

                avgCEven = 1.
                avgCOdd  = ( self._AProdVal + self._ANormVal ) / ( 1. + self._AProdVal * self._ANormVal )

            else :
                # use default values
                avgCEven = 1.
                avgCOdd  = 0.

            # create coefficients object
            from P2VVParameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
            CEvenOdd = Coefficients_CEvenOdd( avgCEven = avgCEven, avgCOdd = avgCOdd )

        # check for remaining arguments and initialize
        self._check_extraneous_kw( kwargs )
        from RooFitWrappers import FormulaVar
        if hasattr( self, 'wTag' ) and hasattr( self, 'wTagBar' ) :
            TaggingParams.__init__( self
                                   , Dilutions = [ FormulaVar(  'tagDilution', '1. - @0 - @1'
                                                   , [ self._wTag, self._wTagBar ], Title = 'Average tagging dilution' ) ]
                                   , ADilWTags = [ FormulaVar(  'ADilWTag',    '(@0 - @1) / (1. - @0 - @1)'
                                                   , [ self._wTag, self._wTagBar ], Title = 'Dilution/wrong tag asymmetry' ) ]
                                   , CEvenOdds = [ CEvenOdd ]
                                  )
        else :
            TaggingParams.__init__( self
                                   , Dilutions = [ FormulaVar(  'tagDilution', '1. - 2. * @0'
                                                   , [ self._WTag ],              Title = 'Average tagging dilution' ) ]
                                   , ADilWTags = [ FormulaVar(  'ADilWTag',    '2. * @0 * @1 / (1. - 2. * @0)'
                                                   , [ self._WTag, self._AWTag ], Title = 'Dilution/wrong tag asymmetry' ) ]
                                   , CEvenOdds = [ CEvenOdd ]
                                  )


class WTagCatsCoefAsyms_TaggingParams( TaggingParams ) :
    """flavour tagging parameters with tagging categories, wrong-tag parameters and asymmetry coefficients
    """

    def __init__( self, **kwargs ) :
        # get number of tagging categories
        numTagCats = kwargs.pop('NumTagCats')
        if numTagCats < 1 :
            raise KeyError('WTagsCoefAsyms_TaggingParams: number of tagging categories must be greater than or equal to one')

        # initialize parameter lists
        tagCatCoefs       = [ ]
        dilutions         = [ ]
        ADilWTags         = [ ]
        CEvenOdds         = [ ]
        self._ATagEffVals = [ ]
        self._wTags       = [ ]
        self._wTagBars    = [ ]
        self._WTags       = [ ]
        self._AWTags      = [ ]

        # get wrong tag parameters category 0
        if ( 'wTag0' in kwargs and 'wTagBar0' in kwargs ) or ( 'WTag0' in kwargs and 'AWTag0' in kwargs ) :
            startTagCat = 0
        else :
            startTagCat = 1
            self._wTags.append(None)
            self._wTagBars.append(None)
            self._WTags.append(None)
            self._AWTag.append(None)

            from RooFitWrappers import ConstVar
            dilutions.append( ConstVar( Name = 'tagDilution0', Value = 0. ) )
            ADilWTags.append( ConstVar( Name = 'ADilWTag0',    Value = 0. ) )

        # get average even and average odd coefficients of category 0
        if 'CEvenOddSum' in kwargs :
            # a coefficients object is specified
            CEvenOddSum = kwargs.pop('CEvenOddSum')

        else :
            from RooFitWrappers import RooObject
            if 'AvgCEvenSum' in kwargs and 'AvgCOddSum' in kwargs :
                # (values for the) coefficients are specified
                avgCEvenSum = kwargs.pop('AvgCEvenSum')
                avgCOddSum  = kwargs.pop('AvgCOddSum')

                if not isinstance( avgCEvenSum, RooObject ) and not isinstance( avgCOddSum, RooObject ) :
                    avgCOddSum  /= avgCEvenSum
                    avgCEvenSum  = 1.

            elif 'AProd' in kwargs and 'ANorm' in kwargs :
                # values for the production and normalization asymmetries are specified
                self._AProdVal = kwargs.pop('AProd')
                self._ANormVal = kwargs.pop('ANorm')
                if isinstance( self._AProdVal, RooObject ) : self._AProdVal = self._AProdVal.getVal()
                if isinstance( self._ANormVal, RooObject ) : self._ANormVal = self._ANormVal.getVal()

                avgCEvenSum = 1.
                avgCOddSum  = ( self._AProdVal + self._ANormVal ) / ( 1. + self._AProdVal * self._ANormVal )

            else :
                # use default values
                avgCEvenSum = 1.
                avgCOddSum  = 0.

            from P2VVParameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
            CEvenOddSum = Coefficients_CEvenOdd(  avgCEven = avgCEvenSum if isinstance( avgCEvenSum, RooObject ) \
                                                             else { 'Name' : 'avgCEvenSum', 'Value' : avgCEven }
                                                , avgCOdd  = avgCOddSum if isinstance( avgCOddSum, RooObject )   \
                                                             else { 'Name' : 'avgCOddSum', 'Value' : avgCOdd }
                                               )

        CEvenOdds.append(CEvenOddSum)

        # loop over tagging categories
        for index in range( startTagCat, numTagCats ) :
            if index > 0 :
                # get tagging category coefficient
                self._parseArg(  'tagCatCoef%d' % index, kwargs, ContainerList = tagCatCoefs
                               , Title = 'Tagging category coefficient %d' % index
                               , Value = (1. - float(index) / float(numTagCats)) / float(numTagCats)
                               , MinMax = ( 0., 1. )
                              )

                # get average even and average odd coefficients
                if 'CEvenOdd%d' % index in kwargs :
                    # a coefficients object is specified
                    CEvenOdd = kwargs.pop('CEvenOdd%d' % index)

                else :
                    if 'AvgCEven%d' % index in kwargs and 'AvgCOdd%d' % index in kwargs :
                        # (values for the) coefficients are specified
                        avgCEven = kwargs.pop('AvgCEven%d' % index)
                        avgCOdd  = kwargs.pop('AvgCOdd%d'  % index)

                        if not isinstance( avgCEven, RooObject ) and not isinstance( avgCOdd, RooObject ) :
                            avgCOdd  /= 1. + self._AProdVal * self._ANormVal
                            avgCEven /= 1. + self._AProdVal * self._ANormVal

                    elif 'ATagEff%d' % index in kwargs and hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                        # values for the asymmetries are specified
                        ATagEffVal = kwargs.pop('ATagEff%d' % index)
                        if isinstance( ATagEffVal, RooObject ) : ATagEffVal = ATagEffVal.getVal()
                        self._ATagEffVals.append( ATagEffVal )

                        avgCEven = 1. + self._AProdVal * self._ANormVal + self._AProdVal * ATagEffVal + self._ANormVal * ATagEffVal
                        avgCOdd  = self._AProdVal + self._ANormVal + ATagEffVal + self._AProdVal * self._ANormVal * ATagEffVal
                        avgCEven /= 1. + self._AProdVal * self._ANormVal
                        avgCOdd  /= 1. + self._AProdVal * self._ANormVal

                    else :
                        # use values for tagging efficiency asymmetry = 0
                        avgCEven = CEvenOddSum['avgCEven'].getVal()
                        avgCEven = CEvenOddSum['avgCOdd'].getVal()

                    CEvenOdd = Coefficients_CEvenOdd( avgCEven = avgCEven, avgCOdd = avgCOdd )
                    CEvenOdd = Coefficients_CEvenOdd(  avgCEven = avgCEven if isinstance( avgCEven, RooObject ) \
                                                                  else { 'Name' : 'avgCEven%d' % index, 'Value' : avgCEven }
                                                     , avgCOdd  = avgCOdd if isinstance( avgCOdd, RooObject )   \
                                                                  else { 'Name' : 'avgCOdd%d' % index, 'Value' : avgCOdd }
                                                    )

                CEvenOdds.append(CEvenOdd)

            # get wrong tag parameters
            from RooFitWrappers import FormulaVar
            if 'wTag%d' % index in kwargs and 'wTagBar%d' % index in kwargs :
                self._WTags.append(None)
                self._AWTags.append(None)
                self._parseArg(  'wTag%d' % index, kwargs, ContainerList = self._wTags
                               , Title = 'B wrong tag probability %d' % index
                               , Value = 0.5 * (1. - float(index) / float(numTagCats))
                               , MinMax = ( 0., 0.5 )
                              )
                self._parseArg(  'wTagBar%d' % index, kwargs, ContainerList = self._wTagBars
                               , Title = 'Bbar wrong tag probability %d' % index
                               , Value = 0.5 * (1. - float(index) / float(numTagCats))
                               , MinMax = ( 0., 0.5 )
                              )
                dilutions.append( FormulaVar(  'tagDilution%d' % index, '1. - @0 - @1'
                                             , [ self._wTags[index], self._wTagBars[index] ]
                                             , Title = 'Average tagging dilution %d' % index
                                            )
                                )
                ADilWTags.append( FormulaVar(  'ADilWTag%d' % index, '(@0 - @1) / (1. - @0 - @1)'
                                             , [ self._wTags[index], self._wTagBars[index] ]
                                             , Title = 'Dilution/wrong tag asymmetry %d' % index
                                            )
                                )

            else :
                self._wTags.append(None)
                self._wTagBars.append(None)
                self._parseArg(  'WTag%d' % index, kwargs, ContainerList = self._WTags
                               , Title = 'Average wrong tag probability %d' % index
                               , Value = 0.5 * (1. - float(index) / float(numTagCats))
                               , MinMax = ( 0., 0.5 )
                              )
                self._parseArg(  'AWTag%d' % index, kwargs, ContainerList = self._AWTags
                               , Title = 'Wrong tag probability asymmetry %d' % index
                               , Value = 0.
                               , MinMax = ( -1., 1. )
                              )
                dilutions.append( FormulaVar(  'tagDilution%d' % index, '1. - 2. * @0'
                                             , [ self._WTags[index], self._AWTags[index] ]
                                             , Title = 'Average tagging dilution %d' % index
                                            )
                                )
                ADilWTags.append( FormulaVar(  'ADilWTag%d' % index, '2. * @0 * @1 / (1. - 2. * @0)'
                                             , [ self._WTags[index], self._AWTags[index] ]
                                             , Title = 'Dilution/wrong tag asymmetry %d' % index
                                            )
                                )

        # check for remaining keyword arguments and initialize
        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__(  self, NumTagCats = numTagCats
                               , TagCatCoefs = tagCatCoefs
                               , Dilutions   = dilutions
                               , ADilWTags   = ADilWTags
                               , CEvenOdds   = CEvenOdds
                              )

class TaggingCategories( _util_parse_mixin ) :
    def __init__( self, **kwargs ) :
        # get tagging category variable
        self._numTagCats = kwargs.pop('NumTagCats')
        self._parseArg(  'tagCat', kwargs, ObjectType = 'Category', SingleArgKey = 'Name'
                       , Title = 'Tagging Category', Observable = True
                       , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, self._numTagCats ) ]
                      )

        self._tagCatCoefs = kwargs.pop('TagCatCoefs')
        self._ATagEffs    = kwargs.pop('ATagEffs'   )
        self._WTags       = kwargs.pop('WTags'      )
        self._AWTags      = kwargs.pop('AWTags'     )

        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )

    def tagCatsDict( self ) :
        tagCatsDict = dict(  [ ( 'NumTagCats', self._numTagCats ) ]
                           + [ ( 'tagCatCoef%d' % ( cat + 1 ), coef  ) for cat, coef  in enumerate( self._tagCatCoefs[ 1 : ] ) ]
                           + [ ( 'ATagEff%d'    % ( cat + 1 ), asym  ) for cat, asym  in enumerate( self._ATagEffs[ 1 : ]    ) ]
                           + [ ( 'WTag%d'       % ( cat + 1 ), WTag  ) for cat, WTag  in enumerate( self._WTags[ 1 : ]       ) ]
                           + [ ( 'AWTag%d'      % ( cat + 1 ), AWTag ) for cat, AWTag in enumerate( self._AWTags[ 1 : ]      ) ]
                          )


class Independent_TaggingCategories( TaggingCategories ) :
    def __init__( self, **kwargs ) :
        # get number of tagging categories
        if 'NumTagCats' not in kwargs : 
            raise KeyError('Independent_TaggingCategories: did not find "NumTagCats" argument')
        else :
            numTagCats = kwargs.pop('NumTagCats')

        # get tagging category variable (or its name)
        tagCat = kwargs.pop( 'tagCat', 'tagCat' )

        # get category parameters
        def getCatParam(name) :
            params = kwargs.pop( name, [ ] )
            if params and len(params) != numTagCats - 1 :
                raise AssertionError(  'Independent_TaggingCategories: length of %s list (%d) ' % ( name, len(params) ) \
                                     + 'is not equal to number of categories - 1 (%d - 1)' % numTagCats
                                    )
            return params

        catParams = [ getCatParam(params) for params in [ 'TagCatCoefs', 'ATagEffs', 'WTags', 'AWTags' ] ]

        if numTagCats == 6 :
            # set default parameters for the six (standard) categories
            if not catParams[0] : catParams[0] = [ 0.15, 0.07, 0.03, 0.01, 0.003 ]
            if not catParams[1] : catParams[1] = 5 * [ 0. ]
            if not catParams[2] : catParams[2] = [ 0.40, 0.35, 0.27, 0.24, 0.12  ]
            if not catParams[3] : catParams[3] = 5 * [ 0. ]
        else :
            # loop over tagging categories and set default parameters
            for cat in range( numTagCats - 1 ) :
                if len(catParams[0]) == cat : pass
                if len(catParams[1]) == cat : catParams[1].append(0.)
                if len(catParams[2]) == cat : pass
                if len(catParams[3]) == cat : catParams[3].append(0.)

        # check for remaining arguments and initialize
        self._check_extraneous_kw( kwargs )
        TaggingCategories.__init__(  self, NumTagCats = numTagCats, tagCat = tagCat, TagCatCoefs = catParams[0], ATagEffs = catParams[1]
                                   , WTags = catParams[2], AWTags = catParams[3]
                                  )


class Linear_TaggingCategories( Independent_TaggingCategories ) :
    def __init__( self, **kwargs ) :
        # set tagging category binning in estimated wrong-tag probability (eta)
        if 'TagCatBins' in kwargs :
            self._tagCatBins = kwargs.pop('TagCatBins')
        else :
            self._tagCatBins = [  ( 'Untagged', 0, 0.500001, 0.50 )
                                , ( 'tagCat1',  1, 0.499999, 0.43 )
                                , ( 'tagCat2',  2, 0.38,     0.35 )
                                , ( 'tagCat3',  3, 0.31,     0.28 )
                                , ( 'tagCat4',  4, 0.24,     0.21 )
                                , ( 'tagCat5',  5, 0.17,     0.14 )
                               ]

        TaggingCategories.__init__( self, **dict( list(kwargs.items()) + list(addArgs.items()) ) )


class Trivial_Background_Tag( _util_parse_mixin ) :
    def pdf(self) :
        return self._pdf
    def __init__( self, tagdecision, **kwargs ) :
        triple =  set([-1,0,+1])==set([ i for i in tagdecision.states().itervalues()] )
        double =  set([-1,+1])==set([ i for i in tagdecision.states().itervalues()] )
        assert triple or double
        if triple : self._parseArg('bkg_tag_eps',   kwargs, Title = 'background tagging efficiency ', Value = 0.25, MinMax = (0.0,1.0) )
        self._parseArg('bkg_tag_delta', kwargs, Title = 'background tagging asymmetry ', Value = 0.0, MinMax = (-0.5,0.5) )
        from RooFitWrappers import GenericPdf
        name = kwargs.pop('Name','Trivial_Background_TagPdf')
        if triple :
            self._pdf = GenericPdf( name, Formula = '(@0==0)*(1-@1)+(@0!=0)*@1*0.5*(1+@0*@2)'
                                        , Arguments = [ tagdecision,self._bkg_tag_eps,self._bkg_tag_delta ] )
        else :
            self._pdf = GenericPdf( name, Formula = '0.5*(1+@0*@1)'
                                        , Arguments = [ tagdecision,self._bkg_tag_delta ] )

        for (k,v) in kwargs.iteritems() :
            setattr(self,'_'+k,v)

class Uniform_Background_Tag( _util_parse_mixin ) :
    def pdf(self) :
        return self._pdf
    def __init__( self, tagdecision, **kwargs ) :
        from RooFitWrappers import GenericPdf
        name = kwargs.pop('Name','Uniform_Background_TagPdf')
        self._pdf = GenericPdf( name, Formula = '1'
                                , Arguments = [ tagdecision ] )
        for (k,v) in kwargs.iteritems() :
            setattr(self,'_'+k,v)
