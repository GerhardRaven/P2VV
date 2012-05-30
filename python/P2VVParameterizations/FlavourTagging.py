###########################################################################################################################################
## P2VVParameterizations.FlavourTagging: Flavour tagging parameters                                                                      ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin


###########################################################################################################################################
## Tagging parameters classes ##
################################

class TaggingParams ( _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin ):
    def __init__( self, **kwargs ) :
        self._numTagCats  = kwargs.pop( 'NumTagCats', [ 0, 1 ] )
        if type(self._numTagCats) == int : self._numTagCats = ( 0, self._numTagCats )
        assert self._numTagCats[0] >= 0 and self._numTagCats[1] > 0,\
               'TaggingParams.__init__(): invalid numbers of tagging categories: %d and %d' % ( self._numTagCats[0], self._numTagCats[1] )

        self._dilutions = kwargs.pop('Dilutions')
        self._ADilWTags = kwargs.pop('ADilWTags')
        self._CEvenOdds = kwargs.pop('CEvenOdds')
        if self._numTagCats[0] == 0 : self._dilutions = ( [ ], self._dilutions )
        if self._numTagCats[0] == 0 : self._ADilWTags = ( [ ], self._ADilWTags )
        if self._numTagCats[0] == 0 : self._CEvenOdds = [ self._CEvenOdds ]

        if self._numTagCats[0] > 1 or self._numTagCats[1] > 1 :
            self._tagCatCoefs = kwargs.pop('TagCatCoefs')
            if self._numTagCats[0] == 0 : self._tagCatCoefs = [ self._tagCatCoefs ]

        # cache integrals as a function of observables
        for dils in self._dilutions :
            for dil in dils :
                dil.setAttribute('CacheAndTrack')
                from ROOT import RooAbsReal, RooArgSet
                realObs = RooArgSet( [ obs._var for obs in dil.Observables() if isinstance( obs._var, RooAbsReal )  ]  )
                if len(realObs) > 0 :
                    print 'P2VV - INFO: TaggingParams.__init__(): invoking %s.parameterizeIntegral(%s)'\
                          % ( dil.GetName(), [ obs.GetName() for obs in realObs ] )
                    dil.setParameterizeIntegral(realObs)

        _util_conditionalObs_mixin.__init__( self, kwargs )
        _util_extConstraints_mixin.__init__( self, kwargs )
        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw ) :
        def raiseError(kw) : raise RuntimeError( 'TaggingParams.__getitem__(\'%s\'): need to specify tagging category' % kw )

        if kw == 'numTagCats' :
            return self._numTagCats[1] if self._numTagCats[0] == 0 else self._numTagCats
        if kw == 'tagCatCoefs' :
            return self._tagCatCoefs[0] if self._numTagCats[0] == 0 else self._tagCatCoefs
        if kw == 'dilution' :
            return self._dilutions[0][0] if self._numTagCats[0] == 0 and self._numTagCats[1] == 1 else raiseError(kw)
        if kw == 'dilutions' :
            return self._dilutions[1] if self._numTagCats[0] == 0 else self._dilutions
        if kw == 'ADilWTag' :
            return self._ADilWTags[0][0] if self._numTagCats[0] == 0 and self._numTagCats[1] == 1 else raiseError(kw)
        if kw == 'ADilWTags' :
            return self._ADilWTags[1] if self._numTagCats[0] == 0 else self._ADilWTags
        if kw == 'CEvenOdd' :
            return self._CEvenOdds[0][0] if self._numTagCats[0] == 0 and self._numTagCats[1] == 1 else raiseError(kw)
        if kw == 'CEvenOdds' :
            return self._CEvenOdds[0] if self._numTagCats[0] == 0 else self._CEvenOdds
        if kw in [ 'avgCEven', 'avgCOdd' ] :
            return self._CEvenOdds[0][0][kw] if self._numTagCats[0] == 0 and self._numTagCats[1] == 1 else raiseError(kw)
        if kw in [ 'avgCEvens', 'avgCOdds' ] :
            return [ CEvenOdd[ kw[:-1] ] for CEvenOdd in self._CEvenOdds[0] ] if self._numTagCats[0] == 0 \
                   else [ [ CEvenOdd[ kw[:-1] ] for CEvenOdd in CEvenOdds ] for CEvenOdds in self._CEvenOdds ]

        return getattr( self, '_' + kw )


class Trivial_TaggingParams( TaggingParams ) :
    def __init__( self ) :
        from RooFitWrappers import ConstVar
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        TaggingParams.__init__( self
                              , Dilutions = [ ConstVar( Name = 'one',  Value = 1 ) ]
                              , ADilWTags = [ ConstVar( Name = 'zero', Value = 0 ) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class WTag_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'wTag', kwargs, Title = 'B wrong-tag probability', Value = 0.25, MinMax = ( 0., 0.5 ) )

        self._check_extraneous_kw( kwargs )
        from RooFitWrappers import FormulaVar, ConstVar
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        TaggingParams.__init__( self
                              , Dilutions = [ FormulaVar(  'tagDilution', '1. - 2*@0 '
                                                         , [ self._wTag ], Title = 'Tagging dilution' ) ]
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class LinearEstWTag_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'estWTag',    kwargs, Title = 'Estimated wrong-tag probability',         Value = 0.25,  MinMax = ( 0.,  0.5 ) )
        self._parseArg( 'p0',         kwargs, Title = 'p0  tagging parameter',                   Value = 0.392, MinMax = ( 0.,  0.5 ) )
        self._parseArg( 'p1',         kwargs, Title = 'p1  tagging parameter',                   Value = 1.035, MinMax = ( 0.8, 1.2 ) )
        self._parseArg( 'avgEstWTag', kwargs, Title = 'Average estimated wrong-tag probability', Value = 0.391, MinMax = ( 0.,  0.5 )
                       , Constant = True
                      )

        constraints = [ ]
        if kwargs.pop( 'p0Constraint', None ) :
            from RooFitWrappers import Pdf, ConstVar
            from ROOT import RooGaussian as Gaussian
            constraints.append( Pdf(  Name = self._p0.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._p0
                                                    , ConstVar( Name = 'p0_mean',  Value = 0.392 )
                                                    , ConstVar( Name = 'p0_sigma', Value = 0.009 )
                                                   ]
                                   )
                              )
            self._p0['Error'] = 0.009

        if kwargs.pop( 'p1Constraint', None ) :
            from RooFitWrappers import Pdf, ConstVar
            from ROOT import RooGaussian as Gaussian
            constraints.append( Pdf(  Name = self._p1.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._p1
                                                    , ConstVar( Name = 'p1_mean',  Value = 1.035 )
                                                    , ConstVar( Name = 'p1_sigma', Value = 0.024 )
                                                   ]
                                   )
                              )
            self._p1['Error'] = 0.024

        self._check_extraneous_kw( kwargs )
        from RooFitWrappers import CalibratedDilution, ConstVar
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        TaggingParams.__init__(  self
                                , Dilutions    = [ CalibratedDilution(  Name       = 'tagDilution'
                                                                      , EstWTag    = self._estWTag
                                                                      , AvgEstWTag = self._avgEstWTag
                                                                      , P0         = self._p0
                                                                      , P1         = self._p1
                                                                     )
                                                 ]
                                , ADilWTags    = [ ConstVar( Name = 'zero', Value = 0) ]
                                , CEvenOdds    = [ Trivial_CEvenOdd() ]
                                , Conditionals = [ self._estWTag ]
                                , Constraints  = constraints
                              )

class Dilution_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'dilution', kwargs, Title = 'Tagging dilution', Value = 0.5, MinMax = ( 0., 1. ) )

        self._check_extraneous_kw( kwargs )
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        from RooFitWrappers import ConstVar
        TaggingParams.__init__( self
                              , Dilutions = [ self._dilution ]
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class WTagsCoefAsyms_TaggingParams( TaggingParams ) :
    """flavour tagging parameters with wrong-tag parameters and asymmetry coefficients
    """

    def __init__( self, **kwargs ) :
        # get wrong-tag parameters
        if 'wTag' in kwargs and 'wTagBar' in kwargs :
            # wrong tag + wrong tag bar
            self._parseArg( 'wTag',    kwargs, Title = 'B wrong-tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
            self._parseArg( 'wTagBar', kwargs, Title = 'Bbar wrong-tag probability', Value = 0.25, MinMax = ( 0., 0.5 ) )
        else :
            # average wrong tag + wrong tag asymmetry
            self._parseArg( 'WTag',  kwargs, Title = 'Average wrong-tag probability',   Value = 0.25, MinMax = (  0., 0.5 ) )
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
                                                   , [ self._wTag, self._wTagBar ], Title = 'Tagging dilution' ) ]
                                   , ADilWTags = [ FormulaVar(  'ADilWTag',    '(@0 - @1) / (1. - @0 - @1)'
                                                   , [ self._wTag, self._wTagBar ], Title = 'Dilution/wrong-tag asymmetry' ) ]
                                   , CEvenOdds = [ CEvenOdd ]
                                  )
        else :
            TaggingParams.__init__( self
                                   , Dilutions = [ FormulaVar(  'tagDilution', '1. - 2. * @0'
                                                   , [ self._WTag ],              Title = 'Tagging dilution' ) ]
                                   , ADilWTags = [ FormulaVar(  'ADilWTag',    '2. * @0 * @1 / (1. - 2. * @0)'
                                                   , [ self._WTag, self._AWTag ], Title = 'Dilution/wrong-tag asymmetry' ) ]
                                   , CEvenOdds = [ CEvenOdd ]
                                  )


class CatDilutionsCoefAsyms_TaggingParams( TaggingParams ) :
    """flavour tagging parameters with two tags, tagging categories, dilution parameters and asymmetry coefficients
    """

    def __init__( self, **kwargs ) :
        # get number of tagging categories
        if 'NumTagCats0' in kwargs and 'NumTagCats1' in kwargs : numTagCats = ( kwargs.pop('NumTagCats0'), kwargs.pop('NumTagCats1') )
        else                                                   : numTagCats = ( 0,                         kwargs.pop('NumTagCats')  )

        if numTagCats[0] < 0 or numTagCats[1] < 1 :
            raise KeyError('CatDilutionsCoefAsyms_TaggingParams: number of tagging categories must be greater than or equal to one')

        if numTagCats[0] == 0 : namePF = ''
        else                  : namePF = '1_'

        # initialize parameter lists
        tagCatCoefs             = [ [ ] for cat in range( numTagCats[0] if numTagCats[0] > 0 else 1 ) ]
        dilutions               = ( [ ], [ ] )
        ADilWTags               = ( [ ], [ ] )
        CEvenOdds               = [ [ ] for cat in range( numTagCats[0] if numTagCats[0] > 0 else 1 ) ]
        self._singleTagCatCoefs = ( [ ], [ ] )
        self._ATagEffVals       = ( [ ], [ ] )

        # get dilution parameters category 0
        startTagCat0 = 1
        if numTagCats[0] > 0 :
            if 'tagDilution0_0' in kwargs and 'ADilWTag0_0' in kwargs :
                startTagCat0 = 0
            else :
                from RooFitWrappers import ConstVar
                dilutions[0].append( ConstVar( Name = 'tagDilution0_0', Value = 0. ) )
                ADilWTags[0].append( ConstVar( Name = 'ADilWTag0_0',    Value = 0. ) )

        startTagCat1 = 1
        if 'tagDilution%s0' % namePF in kwargs and 'ADilWTag%s0' % namePF in kwargs :
            startTagCat1 = 0
        else :
            from RooFitWrappers import ConstVar
            dilutions[1].append( ConstVar( Name = 'tagDilution%s0' % namePF, Value = 0. ) )
            ADilWTags[1].append( ConstVar( Name = 'ADilWTag%s0'    % namePF, Value = 0. ) )

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
                                                             else { 'Name' : 'avgCEvenSum', 'Value' : avgCEvenSum }
                                                , avgCOdd  = avgCOddSum if isinstance( avgCOddSum, RooObject )   \
                                                             else { 'Name' : 'avgCOddSum', 'Value' : avgCOddSum }
                                               )

        CEvenOdds[0].append(CEvenOddSum)

        # get tagging category coefficients and asymmetries (assume factorization)
        from RooFitWrappers import ConstVar, FormulaVar, Product
        if numTagCats[0] > 1 :
            for index0 in range( 1, numTagCats[0] ) :
                self._parseArg(  'tagCatCoef0_%d' % index0, kwargs, ContainerList = self._singleTagCatCoefs[0]
                               , Title = 'Tagging category coefficient 0 %d' % index0
                               , Value = 1. / float( numTagCats[0] )
                               , MinMax = ( 0., 1. )
                              )

                if hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                    ATagEffVal = kwargs.pop( 'ATagEff0_%d' % index0, 0. )
                    if isinstance( ATagEffVal, RooObject ) : ATagEffVal = ATagEffVal.getVal()
                    self._ATagEffVals[0].append( ATagEffVal )

            self._singleTagCatCoefs[0].insert( 0, FormulaVar(  'tagCatCoef0_0'
                                                             , '1-@' + '-@'.join( str(i) for i in range( numTagCats[0] - 1 ) )
                                                             , self._singleTagCatCoefs[0]
                                                            )
                                             )

            if hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                tagCatProd0 = 0.
                for tagCatCoefVal, ATagEffVal in zip( self._singleTagCatCoefs[0][ 1 : ], self._ATagEffVals[0] ) :
                    tagCatProd0 -= tagCatCoefVal * ATagEffVal
                self._ATagEffVals[0].insert( 0, tagCatProd0 / self._singleTagCatCoefs[0][0] )

        else :
            self._ATagEffVals[0].append(0)

        if numTagCats[1] > 1 :
            for index1 in range( 1, numTagCats[1] ) :
                self._parseArg(  'tagCatCoef%s%d' % ( namePF, index1 ), kwargs, ContainerList = self._singleTagCatCoefs[1]
                               , Title = 'Tagging category coefficient 1 %d' % index1
                               , Value = 1. / float( numTagCats[1] )
                               , MinMax = ( 0., 1. )
                              )

                if hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                    ATagEffVal = kwargs.pop( 'ATagEff%s%d' % ( namePF, index1 ), 0. )
                    if isinstance( ATagEffVal, RooObject ) : ATagEffVal = ATagEffVal.getVal()
                    self._ATagEffVals[1].append( ATagEffVal )

            if numTagCats[0] > 0 :
                self._singleTagCatCoefs[1].insert( 0, FormulaVar(  'tagCatCoef%s0' % namePF
                                                                 , '1-@' + '-@'.join( str(i) for i in range( numTagCats[1] - 1 ) )
                                                                 , self._singleTagCatCoefs[1]
                                                                )
                                                 )

                if hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                    tagCatProd0 = 0.
                    for tagCatCoefVal, ATagEffVal in zip( self._singleTagCatCoefs[1][ 1 : ], self._ATagEffVals[1] ) :
                        tagCatProd0 -= tagCatCoefVal * ATagEffVal
                    self._ATagEffVals[1].insert( 0, tagCatProd0 / self._singleTagCatCoefs[1][0] )

            else :
                self._singleTagCatCoefs[1].insert( 0, None )
                self._ATagEffVals[1].insert( 0, None )

        else :
            self._ATagEffVals[1].append(0)
            self._singleTagCatCoefs[1][0] = ConstVar( Name = 'tagCatCoef%s0' % namePF, Value = 1. )

        # loop over tagging categories
        for index0 in range( numTagCats[0] if numTagCats[0] > 0 else 1 ) :
            if index0 >= startTagCat0 :
                # get dilution parameters categories 0
                self._parseArg(  'tagDilution0_%d' % index0, kwargs, ContainerList = dilutions[0]
                               , Title = 'Tagging dilution 0 %d' % index0
                               , Value = float(index0) / float( numTagCats[0] )
                               , MinMax = ( 0., 1. )
                              )
                self._parseArg(  'ADilWTag0_%d' % index0, kwargs, ContainerList = ADilWTags[0]
                               , Title = 'Dilution/wrong-tag asymmetry 0 %d' % index0
                               , Value = 0.
                               , MinMax = ( -1., 1. )
                              )

            for index1 in range( numTagCats[1] ) :
                if index0 == 0 and index1 >= startTagCat1 :
                    # get dilution parameters categories 1
                    self._parseArg(  'tagDilution%s%d' % ( namePF, index1 ), kwargs, ContainerList = dilutions[1]
                                   , Title = 'Tagging dilution 1 %d' % index1
                                   , Value = float(index1) / float( numTagCats[1] )
                                   , MinMax = ( 0., 1. )
                                  )
                    self._parseArg(  'ADilWTag%s%d' % ( namePF, index1 ), kwargs, ContainerList = ADilWTags[1]
                                   , Title = 'Dilution/wrong-tag asymmetry 1 %d' % index1
                                   , Value = 0.
                                   , MinMax = ( -1., 1. )
                                  )

                if index0 > 0 or index1 > 0 :
                    if numTagCats[0] > 0 or numTagCats[1] > 1 :
                        # build tagging category coefficient
                        if self._singleTagCatCoefs[0] :
                            if self._singleTagCatCoefs[1] :
                                tagCatCoef = Product(  'tagCatCoef%d-%d' % ( index0, index1 )
                                                     , [ self._singleTagCatCoefs[0][index0], self._singleTagCatCoefs[1][index1] ]
                                                    )
                            else :
                                tagCatCoef = self._singleTagCatCoefs[0][index0]
                        else :
                             tagCatCoef = self._singleTagCatCoefs[1][index1]
                        tagCatCoefs[index0].append(tagCatCoef)

                    # get average even and average odd coefficients
                    if ( 'CEvenOdd%d-%d' % ( index0, index1 ) if namePF else 'CEvenOdd%d' % index1 ) in kwargs :
                        # a coefficients object is specified
                        CEvenOdd = kwargs.pop( 'CEvenOdd%d-%d' % ( index0, index1 ) if namePF else 'CEvenOdd%d' % index1 )

                    else :
                        if ( 'AvgCEven%d-%d' % ( index0, index1 ) if namePF else 'AvgCEven%d' % index1 ) in kwargs\
                                and ( 'AvgCOdd%d-%d' % ( index0, index1 ) if namePF else 'AvgCOdd%d' % index1 ) in kwargs :
                            # (values for the) coefficients are specified
                            avgCEven = kwargs.pop( 'AvgCEven%d-%d' % ( index0, index1 ) if namePF else 'AvgCEven%d' % index1 )
                            avgCOdd  = kwargs.pop( 'AvgCOdd%d-%d'  % ( index0, index1 ) if namePF else 'AvgCOdd%d'  % index1 )

                            if not isinstance( avgCEven, RooObject ) and not isinstance( avgCOdd, RooObject )\
                                    and hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                                avgCOdd  /= 1. + self._AProdVal * self._ANormVal
                                avgCEven /= 1. + self._AProdVal * self._ANormVal

                        elif hasattr(self, '_AProdVal') and hasattr(self, '_ANormVal') :
                            # values for the asymmetries are specified
                            avgCEven = 1. + self._AProdVal * self._ANormVal \
                                          + self._AProdVal * self._ATagEffVals[0][index0] + self._ANormVal * self._ATagEffVals[0][index0] \
                                          + self._AProdVal * self._ATagEffVals[1][index1] + self._ANormVal * self._ATagEffVals[1][index1] \
                                          + self._AProdVal * self._ANormVal * self._ATagEffVals[0][index0] * self._ATagEffVals[1][index1]
                            avgCOdd  = self._AProdVal + self._ANormVal + self._ATagEffVals[0][index0] + self._ATagEffVals[1][index1] \
                                       + self._AProdVal * self._ANormVal * self._ATagEffVals[0][index0] \
                                       + self._AProdVal * self._ANormVal * self._ATagEffVals[1][index1]

                            avgCEven /= 1. + self._AProdVal * self._ANormVal
                            avgCOdd  /= 1. + self._AProdVal * self._ANormVal

                        else :
                            # use values for tagging efficiency asymmetry = 0
                            avgCEven = CEvenOddSum['avgCEven'].getVal()
                            avgCOdd  = CEvenOddSum['avgCOdd'].getVal()

                        CEvenOdd = Coefficients_CEvenOdd(  avgCEven = avgCEven if isinstance( avgCEven, RooObject ) \
                                                               else { 'Name' : 'avgCEven%d-%d' % ( index0, index1 ), 'Value' : avgCEven }
                                                         , avgCOdd  = avgCOdd if isinstance( avgCOdd, RooObject )   \
                                                               else { 'Name' : 'avgCOdd%d-%d'  % ( index0, index1 ), 'Value' : avgCOdd  }
                                                        )

                    CEvenOdds[index0].append(CEvenOdd)

        # get conditional observables
        conditionals = kwargs.pop( 'Conditionals', [] )

        # get external constraints
        constraints = kwargs.pop( 'Constraints', [] )

        # check for remaining keyword arguments and initialize
        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__(  self
                               , NumTagCats   = numTagCats   if numTagCats[0] > 0 else numTagCats[1]
                               , TagCatCoefs  = tagCatCoefs  if numTagCats[0] > 0 else tagCatCoefs[0]
                               , Dilutions    = dilutions    if numTagCats[0] > 0 else dilutions[1]
                               , ADilWTags    = ADilWTags    if numTagCats[0] > 0 else ADilWTags[1]
                               , CEvenOdds    = CEvenOdds    if numTagCats[0] > 0 else CEvenOdds[0]
                               , Conditionals = conditionals
                               , Constraints  = constraints
                              )


###########################################################################################################################################
## Tagging categories classes ##
################################

def getTagCatParamsFromData( data, estWTagName, tagCats = [ ], numSigmas = 1., avgEstWTag = 0.39, P0    = 0.392, P1    = 1.04
                                                                                                , P0Err = 0.009, P1Err = 0.02
                                                                                                , AP0   = 0.,    AP1   = 0.
                           ) :
    assert data, 'getTagCatParamsFromData(): no data set found'
    assert estWTagName and data.get(0).find(estWTagName), 'getTagCatParamsFromData(): estimated wrong-tag probability not found in data'

    from RooFitWrappers import RooObject
    if isinstance( avgEstWTag, RooObject ) : avgEstWTag = avgEstWTag.getVal()
    if isinstance( P0,         RooObject ) : P0         = P0.getVal()
    if isinstance( P1,         RooObject ) : P1         = P1.getVal()
    if isinstance( AP0,        RooObject ) : AP0        = AP0.getVal()
    if isinstance( AP1,        RooObject ) : AP1        = AP1.getVal()

    tagCatsCalc = tagCats[ : ]

    etaMin = 0.
    if not tagCatsCalc :
        # get minimum estimated wrong-tag probability
        etaMin = 0.5
        for varSet in data : etaMin = min( etaMin, varSet.getRealValue(estWTagName) )

        # determine binning in estimated wrong-tag probability from data
        bin = 0
        binUpperEdges = [ 0.499999 ]
        while True :
            bin += 1

            # get high bin edge
            highEdge = binUpperEdges[ bin - 1 ]

            # calculate low bin edge
            if highEdge >= avgEstWTag : lowEdge = highEdge - 2. * ( P0Err + P1Err * (highEdge - avgEstWTag) ) / ( P1 / numSigmas + P1Err )
            else                      : lowEdge = highEdge - 2. * ( P0Err - P1Err * (highEdge - avgEstWTag) ) / ( P1 / numSigmas - P1Err )

            # set low bin edge
            binUpperEdges.append(lowEdge)

            # check if this is the last bin
            if lowEdge < etaMin : break

        # scale bin widths to match range of estimated wrong-tag probability
        if binUpperEdges[-2] - etaMin < etaMin - binUpperEdges[-1] : del binUpperEdges[-1]
        binScale = ( binUpperEdges[0] - etaMin ) / ( binUpperEdges[0] - binUpperEdges[-1] )
        tagCatsCalc = [ ( 'Untagged', 0, 0.500001 ) ]
        for bin in range( 1, len(binUpperEdges) ) :
            binUpperEdges[bin] = binUpperEdges[0] - ( binUpperEdges[0] - binUpperEdges[bin] ) * binScale
            tagCatsCalc.append( ( 'TagCat%d' % bin, bin, binUpperEdges[ bin - 1 ] ) )

    # determine tagging category parameters
    numTagCats = len(tagCatsCalc)
    numEvTot   = data.sumEntries()
    numEvCats  = [0]  * numTagCats
    sumEtaCats = [0.] * numTagCats
    for varSet in data :
        # get estimated wrong-tag probability for this event
        eta = varSet.getRealValue(estWTagName)

        # determine tagging category
        cat = -1
        for catIter in range(numTagCats) :
          if eta >= tagCatsCalc[catIter][2] : break
          cat += 1

        if cat < 0 : raise RuntimeError('getTagCatParamsFromData(): estimated wrong-tag probability out of range')

        # update number of events and sum of estimated wrong-tag probabilities
        numEvCats[cat]  += data.weight()
        sumEtaCats[cat] += eta * data.weight()

    # check number of events
    numEvTotCount = 0.
    for numEv in numEvCats : numEvTotCount += numEv
    assert abs( numEvTotCount - numEvTot ) < 1.e-10 * abs( numEvTotCount + numEvTot ),\
           'getTagCatParamsFromData(): counted number of events is not equal to number of events in data set'

    # update tagging category parameters
    for cat in range(numTagCats) :
        avgEtaCat = sumEtaCats[cat] / numEvCats[cat]
        tagCatsCalc[cat] = tagCatsCalc[cat][ : 3 ]\
                           + (  avgEtaCat
                              , P0 + P1 * ( avgEtaCat - avgEstWTag ), 0.
                              , numEvCats[cat] / numEvTot, 0.
                              , numEvCats[cat]
                             )

    # print tagging category binning to screen
    print 'P2VV - INFO: getTagCatParamsFromData(): tagging category binning:'
    print '    <eta> = %.3f   P0 = %.3f +- %.3f   P1 = %.3f +- %.3f    P0 asym. = %.3f    P1 asym. = %.3f'\
          % ( avgEstWTag, P0, P0Err, P1, P1Err, AP0, AP1 )
    print '    minimum eta = %.3f    average eta = %.3f' % ( etaMin, data.mean( data.get(0).find(estWTagName) ) )
    for bin, cat in enumerate(tagCatsCalc) :
        deltaEta = ( ( cat[2] - tagCatsCalc[bin + 1][2] ) / 2. ) if bin < len(tagCatsCalc) - 1 else ( ( cat[2] - etaMin ) / 2. )
        if cat[2] >= avgEstWTag : binRangeSig = P1 * deltaEta / ( P0Err + P1Err * ( cat[2] - avgEstWTag - deltaEta ) )
        else                    : binRangeSig = P1 * deltaEta / ( P0Err - P1Err * ( cat[2] - avgEstWTag - deltaEta ) )
        print '    {0:<10s}  :  {1:.3f} -- {2:.3f} ({3:.2f} sigma)  :  <eta> = {4:.3f}  <w> = {5:.3f}  efficiency = {6:.4f} ({7:.0f} events)'\
              .format(  cat[0], cat[2], tagCatsCalc[bin + 1][2] if bin < len(tagCatsCalc) - 1 else etaMin
                      , binRangeSig, cat[3], cat[4], cat[6], cat[8]
                     )

    return tagCatsCalc


class TaggingCategories( _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin ) :
    def __init__( self, **kwargs ) :
        self._numTagCats = kwargs.pop('NumTagCats')
        if hasattr( self, '_combined' ) and self._combined :
            self._parseArg(  'tagCat0', kwargs, ObjectType = 'Category', SingleArgKey = 'Name'
                           , Title = 'Tagging Category 0', Observable = True
                           , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, self._numTagCats[0] ) ]
                          )
            self._parseArg(  'tagCat1', kwargs, ObjectType = 'Category', SingleArgKey = 'Name'
                           , Title = 'Tagging Category 1', Observable = True
                           , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, self._numTagCats[1] ) ]
                          )

        else :
            self._parseArg(  'tagCat', kwargs, ObjectType = 'Category', SingleArgKey = 'Name'
                           , Title = 'Tagging Category', Observable = True
                           , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, self._numTagCats ) ]
                          )

        self._tagCatCoefs  = kwargs.pop( 'TagCatCoefs'  )
        self._ATagEffs     = kwargs.pop( 'ATagEffs'     )
        self._tagDilutions = kwargs.pop( 'TagDilutions' )
        self._ADilWTags    = kwargs.pop( 'ADilWTags'    )

        _util_conditionalObs_mixin.__init__( self, kwargs )
        _util_extConstraints_mixin.__init__( self, kwargs )
        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )

    def tagCatsDict( self ) :
        if hasattr( self, '_combined' ) and self._combined :
            dictList = [ ]
            for cat0 in range(self._numTagCats[0]) :
                if cat0 > 0 :
                    dictList += [  ( 'tagDilution0_%d' % cat0, self._tagDilutions[0][ cat0 - 1 ] )
                                 , ( 'ADilWTag0_%d'    % cat0, self._ADilWTags[0][ cat0 - 1 ]    )
                                ]

                for cat1 in range(self._numTagCats[1]) :
                    if cat0 > 0 or cat1 > 0 :
                        if cat0 == 0 :
                            dictList += [  ( 'tagDilution1_%d' % cat1, self._tagDilutions[1][ cat1 - 1 ] )
                                         , ( 'ADilWTag1_%d'    % cat1, self._ADilWTags[1][ cat1 - 1 ]    )
                                        ]

                        dictList += [  ( 'tagCatCoef%d-%d' % ( cat0, cat1 ), self._tagCatCoefs[cat0][cat1] )
                                     , ( 'ATagEff%d-%d'    % ( cat0, cat1 ), self._ATagEffs[cat0][cat1]    )
                                    ]

        else :
            dictList = [ ( 'tagCatCoef%d'  % ( cat + 1 ), coef     ) for cat, coef     in enumerate( self._tagCatCoefs[ 1 : ] ) ]\
                     + [ ( 'ATagEff%d'     % ( cat + 1 ), asym     ) for cat, asym     in enumerate( self._ATagEffs[ 1 : ]    ) ]\
                     + [ ( 'tagDilution%d' % ( cat + 1 ), dilution ) for cat, dilution in enumerate( self._tagDilutions       ) ]\
                     + [ ( 'ADilWTag%d'    % ( cat + 1 ), ADilWTag ) for cat, ADilWTag in enumerate( self._ADilWTags          ) ]

        return dict(  [ ( 'NumTagCats', self._numTagCats ) ]
                    + dictList
                    + [ ( 'Conditionals', self.conditionalObservables() ) ]
                    + [ ( 'Constraints',  self.externalConstraints()    ) ]
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

        catParams = [ getCatParam(params) for params in [ 'TagCatCoefs', 'ATagEffs', 'TagDilutions', 'ADilWTags' ] ]

        if numTagCats == 6 :
            # set default parameters for the six (standard) categories
            if not catParams[0] : catParams[0] = [ 0.15, 0.07, 0.03, 0.01, 0.003 ]
            if not catParams[1] : catParams[1] = 5 * [ 0. ]
            if not catParams[2] : catParams[2] = [ 0.20, 0.30, 0.46, 0.52, 0.76  ]
            if not catParams[3] : catParams[3] = 5 * [ 0. ]

        else :
            # loop over tagging categories and set default parameters
            for cat in range( numTagCats - 1 ) :
                if len(catParams[0]) == cat :
                    from math import pow
                    numCatsFrac = float(numTagCats) / 6.
                    tagCatCoef  = 0.15 / numCatsFrac * pow( 0.5, float(cat) / numCatsFrac )
                    catParams[0].append( tagCatCoef )

                if len(catParams[1]) == cat : catParams[1].append(0.)
                if len(catParams[2]) == cat : catParams[2].append( float(cat + 1) / float(numTagCats) )
                if len(catParams[3]) == cat : catParams[3].append(0.)

        # check for remaining arguments and initialize
        self._check_extraneous_kw( kwargs )
        TaggingCategories.__init__(  self, NumTagCats = numTagCats, tagCat = tagCat, TagCatCoefs = catParams[0], ATagEffs = catParams[1]
                                   , TagDilutions = catParams[2], ADilWTags = catParams[3]
                                  )


class Linear_TaggingCategories( TaggingCategories ) :
    def __init__( self, **kwargs ) :
        # get type of tagging
        tagType = 'SS' if kwargs.pop( 'SameSide', False ) else 'OS'

        # get tagging category variable (or its name)
        tagCat = kwargs.pop( 'tagCat', 'tagCat' + tagType )

        # estimated wrong-tag variable
        if 'estWTag' in kwargs :
            self._parseArg( 'estWTag', kwargs, Name = 'estWTag' + tagType, Title = 'Estimated wrong-tag probability'
                           , Value = 0.25, MinMax = ( 0., 0.5 ) )

        # get linear calibration parameters
        self._parseArg(  'avgEstWTag', kwargs, Name = 'avgEstWTag' + tagType
                       , Value = 0.391 if tagType == 'OS' else 0.324
                       , ObjectType = 'ConstVar' )
        self._parseArg(  'wTagP0',     kwargs, Name = 'wTagP0' + tagType, Title = 'Average wrong-tag parameter p_0'
                       , Value = 0.392 if tagType == 'OS' else 0.350
                       , Error = 0.009 if tagType == 'OS' else 0.015
                       , MinMax = (  0.,  0.5 )
                      )
        self._parseArg(  'wTagP1', kwargs, Name = 'wTagP1' + tagType, Title = 'Average wrong-tag parameter p_1'
                       , Value = 1.035 if tagType == 'OS' else 0.509
                       , Error = 0.024 if tagType == 'OS' else 0.155
                       , MinMax = (  0., 2. )
                      )
        self._parseArg(  'wTagAP0', kwargs, Name = 'wTagAP0' + tagType, Title = 'Wrong tag parameter p_0 asymmetry'
                       , Value = 0., MinMax = ( -1.,  1.  )
                      )
        self._parseArg(  'wTagAP1', kwargs, Name = 'wTagAP1' + tagType, Title = 'Wrong tag parameter p_1 asymmetry'
                       , Value = 0., MinMax = ( -1.,  1.  )
                      )

        # constrain calibration parameters
        constraints = [ ]
        wTagP0Constraint = kwargs.pop( 'wTagP0Constraint', None )
        if wTagP0Constraint :
            from RooFitWrappers import Pdf, ConstVar
            from ROOT import RooGaussian as Gaussian
            constraints.append( Pdf(  Name = self._wTagP0.GetName() + tagType + '_constraint', Type = Gaussian
                                    , Parameters = [  self._wTagP0
                                                    , ConstVar( Name = 'wTagP0' + tagType + '_mean',  Value = self._wTagP0.getVal()   )
                                                    , ConstVar( Name = 'wTagP0' + tagType + '_sigma', Value = self._wTagP0.getError() )
                                                   ]
                                   )
                              )

        wTagP1Constraint = kwargs.pop( 'wTagP1Constraint', None )
        if wTagP1Constraint :
            from RooFitWrappers import Pdf, ConstVar
            from ROOT import RooGaussian as Gaussian
            constraints.append( Pdf(  Name = self._wTagP1.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._wTagP1
                                                    , ConstVar( Name = 'wTagP1' + tagType + '_mean',  Value = self._wTagP1.getVal()   )
                                                    , ConstVar( Name = 'wTagP1' + tagType + '_sigma', Value = self._wTagP1.getError() )
                                                   ]
                                   )
                              )

        # get data set
        self._data = kwargs.pop( 'DataSet', None )
        self._etaName = kwargs.pop( 'estWTagName', self._estWTag.GetName() if hasattr( self, '_estWTag' ) else '' )
        if not self._etaName : self._data = None

        # get tagging category binning in estimated wrong-tag probability (eta)
        tagCats = kwargs.pop( 'TagCats', None )
        if tagCats == None :
            tagCats = [  ( 'Untagged', 0, 0.500001, 0.50, 0.50, 0., 0.65,  0. )
                       , ( 'Tagged',   1, 0.499999, 0.40, 0.40, 0., 0.35,  0. )
                      ] if hasattr( self, '_estWTag' ) else \
                      [  ( 'Untagged', 0, 0.500001, 0.50, 0.50, 0., 0.65,  0. )
                       , ( 'TagCat1',  1, 0.499999, 0.44, 0.44, 0., 0.24,  0. )
                       , ( 'TagCat2',  2, 0.38,     0.35, 0.35, 0., 0.062, 0. )
                       , ( 'TagCat3',  3, 0.31,     0.28, 0.28, 0., 0.032, 0. )
                       , ( 'TagCat4',  4, 0.24,     0.21, 0.21, 0., 0.012, 0. )
                       , ( 'TagCat5',  5, 0.17,     0.15, 0.14, 0., 0.004, 0. )
                      ]

        # get number of calibration parameter standard deviations for tagging category bins
        nSigmaTagBins = kwargs.pop( 'NumSigmaTagBins', 1. )

        # determine tagging category parameters from data
        if self._data and self._etaName :
            self._tagCats = getTagCatParamsFromData(  self._data, estWTagName = self._etaName, tagCats = tagCats, numSigmas = nSigmaTagBins
                                                    , avgEstWTag = self._avgEstWTag
                                                    , P0    = self._wTagP0.getVal(),   P1    = self._wTagP1.getVal()
                                                    , P0Err = self._wTagP0.getError(), P1Err = self._wTagP1.getError()
                                                    , AP0   = self._wTagAP0,           AP1   = self._wTagAP1
                                                   )

            from math import sqrt
            tagCatCoefs = [ dict(  Value = catPars[6]
                                 , Error = sqrt( catPars[6] / self._data.sumEntries() )
                                 , Constant = True
                                ) for catPars in self._tagCats
                          ]

        else :
            self._tagCats = tagCats
            tagCatCoefs = [ catPars[6] for catPars in self._tagCats ]

        ATagEffs    = [ catPars[7] for catPars in self._tagCats ]
        dilutions = [ ]
        ADilWTags = [ ]
        self._estWTags = [ ]
        from RooFitWrappers import CalibratedDilution
        for cat, catPars in enumerate( self._tagCats[ 1 : ] ) :
            if hasattr( self, '_estWTag' ) :
                estWTag = self._estWTag
            else :
                from RooFitWrappers import ConstVar
                estWTag = ConstVar( Name = 'estWTag%s%d' % ( tagType, cat + 1 ), Title = 'Estimated wrong-tag probability ' + tagType
                                   , Value = catPars[3] )
                self._estWTags.append(estWTag)

            dilutions.append( CalibratedDilution(  Name       = 'tagDilution%s%d' % ( tagType, cat + 1 )
                                                 , Title      = 'Tagging dilution %s %d' % ( tagType, cat + 1 )
                                                 , EstWTag    = estWTag
                                                 , AvgEstWTag = self._avgEstWTag
                                                 , P0         = self._wTagP0
                                                 , P1         = self._wTagP1
                                                )
                            )
            ADilWTags.append( CalibratedDilution(  Name       = 'ADilWTag%s%d' % ( tagType, cat + 1 )
                                                 , Title      = 'Dilution/wrong-tag asymmetry %s %d' % ( tagType, cat + 1 )
                                                 , EstWTag    = estWTag
                                                 , AvgEstWTag = self._avgEstWTag
                                                 , P0         = self._wTagP0
                                                 , P1         = self._wTagP1
                                                 , AP0        = self._wTagAP0
                                                 , AP1        = self._wTagAP1
                                                )
                            )

        # adjust errors on calibration parameters
        if not wTagP0Constraint : self._wTagP0.setError(  3. * self._wTagP0.getError() )
        if not wTagP1Constraint : self._wTagP1.setError( 17. * self._wTagP1.getError() )

        # check for remaining arguments and initialize
        self._check_extraneous_kw( kwargs )
        TaggingCategories.__init__( self, NumTagCats = len(self._tagCats), tagCat = tagCat
                                   , TagCatCoefs = tagCatCoefs, ATagEffs = ATagEffs, TagDilutions = dilutions, ADilWTags = ADilWTags
                                   , Conditionals = [ self._estWTag ] if hasattr( self, '_estWTag' ) else [ ]
                                   , Constraints = constraints
                                  )


    def traditionalCatRange( self, tradCat ) :
        tradCatEdges = [ ( 0.5, 0.5 ), ( 0.499, 0.38 ), ( 0.38, 0.31 ), ( 0.31, 0.24 ), ( 0.24, 0.17 ), ( 0.17, 0. ) ]
        catLowEdge  = 0
        catHighEdge = 0
        for cat in range(len(self._tagCats)) :
            highEdge = self._tagCats[cat][2]
            lowEdge  = self._tagCats[cat + 1][2] if bin < len(self._tagCats) - 1 else -0.0001
            if tradCatEdges[tradCat][0] < highEdge and tradCatEdges[tradCat][0] >= lowEdge : catHighEdge = cat
            if tradCatEdges[tradCat][1] < highEdge and tradCatEdges[tradCat][1] >= lowEdge : catLowEdge  = cat

        return ( catHighEdge, catLowEdge )


class Combined_TaggingCategories( TaggingCategories ) :
    def __init__( self, Categories0, Categories1, **kwargs ) :
        self._combined = True

        tagCatCoefs = [ [ 0. for iter1 in range( Categories1['numTagCats'] ) ] for iter0 in range( Categories0['numTagCats'] ) ]
        ATagEffs    = [ [ 0. for iter1 in range( Categories1['numTagCats'] ) ] for iter0 in range( Categories0['numTagCats'] ) ]
        if hasattr( Categories0, '_data' ) and Categories0['data'] and hasattr( Categories1, '_data' ) and Categories1['data']\
                and hasattr( Categories0, '_tagCats' ) and hasattr( Categories1, '_tagCats' ) :
            # get tagging category coefficients and asymmetries for combined categories from data
            assert Categories0['data'] == Categories1['data'],\
                    'Combined_TaggingCategories.__init__(): data sets for Categories0 and Categories1 are not the same'

            self._data     = Categories0['data']
            self._tagCats  = ( Categories0['tagCats'], Categories1['tagCats'] )
            self._etaName  = ( Categories0['etaName'], Categories1['etaName'] )

            # determine tagging category parameters
            etaMin = ( 0.5, 0.5 )
            numTagCats = ( Categories0['numTagCats'], Categories1['numTagCats'] )
            numEvTot   = self._data.sumEntries()
            numEvCats  = [ [ 0 for iter1 in range( numTagCats[1] ) ] for iter0 in range( numTagCats[0] ) ]
            sumEtaCats = (  [ [ 0. for iter1 in range( numTagCats[1] ) ] for iter0 in range( numTagCats[0] ) ]
                          , [ [ 0. for iter1 in range( numTagCats[1] ) ] for iter0 in range( numTagCats[0] ) ]
                         )
            for varSet in self._data :
                # get estimated wrong-tag probability for this event
                eta = ( varSet.getRealValue(self._etaName[0]), varSet.getRealValue(self._etaName[1]) )

                # update minimum estimated wrong-tag probability
                etaMin = ( min( etaMin[0], eta[0] ), min( etaMin[1], eta[1] ) )

                # determine tagging category
                bin0 = -1
                for catIter0 in range( numTagCats[0] ) :
                  if eta[0] >= self._tagCats[0][catIter0][2] : break
                  bin0 += 1
                if bin0 < 0 : raise RuntimeError('Combined_TaggingCategories.__init__(): estimated wrong-tag probability 0 out of range')

                bin1 = -1
                for catIter1 in range( numTagCats[1] ) :
                  if eta[1] >= self._tagCats[1][catIter1][2] : break
                  bin1 += 1
                if bin1 < 0 : raise RuntimeError('Combined_TaggingCategories.__init__(): estimated wrong-tag probability 1 out of range')

                # update number of events and sum of estimated wrong-tag probabilities
                numEvCats[bin0][bin1] += self._data.weight()
                sumEtaCats[0][bin0][bin1] += eta[0] * self._data.weight()
                sumEtaCats[1][bin0][bin1] += eta[1] * self._data.weight()

            # check number of events
            numEvTotCount = 0.
            for numEvs in numEvCats :
                for numEv in numEvs : numEvTotCount += numEv
            assert abs( numEvTotCount - numEvTot ) < 1.e-10 * abs( numEvTotCount + numEvTot ),\
              'Combined_TaggingCategories.__init__(): counted number of events (%.1f) is not equal to number of events in data set (%.1f)'\
              % ( numEvTotCount, numEvTot )

            # get tagging calibration parameters
            avgEtas = ( Categories0['avgEstWTag'].getVal(), Categories1['avgEstWTag'].getVal() )
            P0s     = ( Categories0['wTagP0'].getVal(),     Categories1['wTagP0'].getVal()     )
            P0Errs  = ( Categories0['wTagP0'].getError(),   Categories1['wTagP0'].getError()   )
            P1s     = ( Categories0['wTagP1'].getVal(),     Categories1['wTagP1'].getVal()     )
            P1Errs  = ( Categories0['wTagP1'].getError(),   Categories1['wTagP1'].getError()   )
            AP0s    = ( Categories0['wTagAP0'].getVal(),    Categories1['wTagAP0'].getVal()    )
            AP1s    = ( Categories0['wTagAP1'].getVal(),    Categories1['wTagAP1'].getVal()    )

            # set coefficients and print tagging category binning to screen
            print 'P2VV - INFO: Combined_TaggingCategories.__init__(): tagging category binning:'
            print '    <eta_0> = %.3f   P0_0 = %.3f +- %.3f   P1_0 = %.3f +- %.3f    P0_0 asym. = %.3f    P1_0 asym. = %.3f'\
                  % ( avgEtas[0], P0s[0], P0Errs[0], P1s[0], P1Errs[0], AP0s[0], AP1s[0] )
            print '    minimum eta_0 = %.3f    average eta_0 = %.3f'\
                  % ( etaMin[0], self._data.mean( self._data.get(0).find( self._etaName[0] ) ) )
            print '    <eta_1> = %.3f   P0_1 = %.3f +- %.3f   P1_1 = %.3f +- %.3f    P0_1 asym. = %.3f    P1_1 asym. = %.3f'\
                  % ( avgEtas[1], P0s[1], P0Errs[1], P1s[1], P1Errs[1], AP0s[1], AP1s[1] )
            print '    minimum eta_1 = %.3f    average eta_1 = %.3f'\
                  % ( etaMin[1], self._data.mean(self._data.get(0).find(self._etaName[1])) )

            from math import sqrt
            self._tagCatsComb = [ [ ( ) for iter1 in range( numTagCats[1] ) ] for iter0 in range( numTagCats[0] ) ]
            for bin0 in range( numTagCats[0] ) :
                for bin1 in range( numTagCats[1] ) :
                    if numEvCats[bin0][bin1] > 0. :
                        avgEtaCat = (  sumEtaCats[0][bin0][bin1] / numEvCats[bin0][bin1]
                                     , sumEtaCats[1][bin0][bin1] / numEvCats[bin0][bin1]
                                    )
                    else :
                        binWidth = (  self._tagCats[0][bin0][2] - (self._tagCats[0][bin0+1][2] if bin0 < numTagCats[0] - 1 else etaMin[0])
                                    , self._tagCats[1][bin1][2] - (self._tagCats[1][bin1+1][2] if bin1 < numTagCats[1] - 1 else etaMin[1])
                                   )
                        avgEtaCat = (  self._tagCats[0][bin0][2] - 0.5 * binWidth[0]
                                     , self._tagCats[1][bin1][2] - 0.5 * binWidth[1]
                                    )

                    cat = (  'Untagged' if bin0 == 0 and bin1 == 0 else 'TagCat%d_%d' % ( bin0, bin1 )
                           , bin0 * numTagCats[1] + bin1
                           , ( avgEtaCat[0], avgEtaCat[1] )
                           , (  P0s[0] + P1s[0] * ( avgEtaCat[0] - avgEtas[0] )
                              , P0s[1] + P1s[1] * ( avgEtaCat[1] - avgEtas[1] ) )
                           , ( 0., 0. )
                           , numEvCats[bin0][bin1] / numEvTot
                           , 0.
                           , numEvCats[bin0][bin1]
                          )

                    self._tagCatsComb[bin0][bin1] = cat
                    tagCatCoefs[bin0][bin1] = dict(  Value = cat[5]
                                                   , Error = sqrt( ( cat[5] / self._data.sumEntries() ) if cat[5] > 0. else 0. )
                                                   , Constant = True
                                                  )
                    ATagEffs[bin0][bin1] = cat[6]

                    print '    {0:<10s}  :  0  :  {1:.3f} -- {2:.3f}  :  <eta> = {3:.3f}  <w> = {4:.3f}  efficiency = {5:.4f} ({6:.0f} events)'\
                          .format( cat[0], self._tagCats[0][bin0][2], self._tagCats[0][bin0 + 1][2]\
                                   if bin0 < numTagCats[0] - 1 else etaMin[0], cat[2][0], cat[3][0], cat[5], cat[7] )
                    print '    {0:<10s}     1  :  {1:.3f} -- {2:.3f}  :  <eta> = {3:.3f}  <w> = {4:.3f}'\
                          .format( '',     self._tagCats[1][bin1][2], self._tagCats[1][bin1 + 1][2]\
                                   if bin1 < numTagCats[1] - 1 else etaMin[1], cat[2][1], cat[3][1]                 )

        else :
            # make tagging category coefficients and asymmetries for combined categories products of the original values
            from math import sqrt
            for bin0, ( coef0, asym0 ) in enumerate( zip( Categories0['tagCatCoefs'], Categories0['ATagEffs'] ) ) :
                for bin1, ( coef1, asym1 ) in enumerate( zip( Categories1['tagCatCoefs'], Categories1['ATagEffs'] ) ) :
                    tagCatCoefs[bin0][bin1] = coef0 * coef1
                    ATagEffs[bin0][bin1] = asym0 * asym1

        # check for remaining arguments and initialize
        self._check_extraneous_kw( kwargs )
        TaggingCategories.__init__( self, NumTagCats = ( numTagCats[0], numTagCats[1] )
                                   , tagCat0 = Categories0['tagCat'], tagCat1 = Categories1['tagCat']
                                   , TagCatCoefs = tagCatCoefs, ATagEffs = ATagEffs
                                   , TagDilutions = ( Categories0['tagDilutions'], Categories1['tagDilutions'] )
                                   , ADilWTags = ( Categories0['ADilWTags'], Categories1['ADilWTags'] )
                                   , Conditionals = Categories0.conditionalObservables() + Categories1.conditionalObservables()
                                   , Constraints  = Categories0.externalConstraints() + Categories1.externalConstraints()
                                  )


###########################################################################################################################################
## Tagging PDFs ##
##################

class Trivial_TagPdf( _util_parse_mixin ) :
    def pdf( self ) : return self._pdf

    def __init__( self, tagdecision, **kwargs ) :
        double =  set([ -1,    +1 ])==set([ i for i in tagdecision.states().itervalues() ] )
        triple =  set([ -1, 0, +1 ])==set([ i for i in tagdecision.states().itervalues() ] )
        assert triple or double

        namePF = kwargs.pop( 'NamePrefix', '' )
        if namePF : namePF += '_'
        if triple :
           self._parseArg( namePF + 'TagEff', kwargs, Title = namePF + 'Tagging efficiency', Value = 0.25, MinMax = ( 0., 1. ) )
        self._parseArg( namePF + 'ATagEff', kwargs, Title = namePF + 'Tagging asymmetry ', Value = 0., MinMax = ( -1., 1. ) )

        from RooFitWrappers import GenericPdf
        name = kwargs.pop( 'Name', namePF + 'Trivial_TagPdf' )
        if triple :
            self._pdf = GenericPdf( name, Formula = '(@0==0)*(1-@1)+(@0!=0)*@1*0.5*(1+@0*@2)'
                                   , Arguments = [  tagdecision
                                                  , getattr( self, '_%sTagEff'%namePF  )
                                                  , getattr( self, '_%sATagEff'%namePF )
                                                 ]
                                  )
        else :
            self._pdf = GenericPdf( name, Formula = '0.5*(1+@0*@1)', Arguments = [ tagdecision, getattr(self, '_' + namePF + 'ATagEff') ] )

        self._pdf.setAttribute("CacheAndTrack") ;
        for ( k, v ) in kwargs.iteritems() : setattr( self, '_' + k, v )


class BinnedTaggingPdf( _util_parse_mixin ) :
    def __init__( self, Name, tagCat, iTag, tagBinCoefs ) :
        from RooFitWrappers import BinnedPdf
        self._name        = Name
        self._tagCat      = tagCat
        self._iTag        = iTag
        self._tagBinCoefs = tagBinCoefs
        self._pdf         = BinnedPdf( Name, Categories = [ tagCat, iTag ], Coefficients = tagBinCoefs )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def pdf( self ) : return self._pdf

    def _init( self, Name, tagCat, iTag, kwargs ) :
        # set tagging category and tag variables
        self._tagCat = tagCat
        self._iTag   = iTag

        # get name prefix from arguments
        self._namePF = kwargs.pop( 'NamePrefix', '' )
        if self._namePF : self._namePF += '_'

        # get name
        self._name = self._namePF + Name

        # get number of tagging categories
        self._numTagCats = tagCat.numTypes()

        # get tagging bin names
        self._untCatName = kwargs.pop( 'UntaggedCatName', 'Untagged' )
        self._tagCatName = kwargs.pop( 'TaggedCatName',   'Tagged'   )
        self._BName      = kwargs.pop( 'BName',           'B'        )
        self._BbarName   = kwargs.pop( 'BbarName',        'Bbar'     )

        # get data
        self._data = kwargs.pop( 'Data', None )

        # get work space
        from RooFitWrappers import RooObject
        self._ws = RooObject().ws()

        if self._data :
            # get number of events in each tagging bin
            from RooFitWrappers import ArgSet
            self._tagTable = self._data.table( ArgSet( self._name + '_tagTableSet', [ self._tagCat, self._iTag ] ) )
            self._nUntB    = self._tagTable.get( '{%s;%s}' % ( self._untCatName, self._BName    ) )
            self._nUntBbar = self._tagTable.get( '{%s;%s}' % ( self._untCatName, self._BbarName ) )

            self._nTagB    = [ 0. for cat in range(self._numTagCats) ]
            self._nTagBbar = [ 0. for cat in range(self._numTagCats) ]
            for catIter in range( 1, tagCat.numTypes() ) :
                tagCatNameMult = ( '%s%d' % ( self._tagCatName, catIter ) ) if self._numTagCats > 2 else self._tagCatName
                self._nTagB[catIter]    += self._tagTable.get( '{%s;%s}' % ( tagCatNameMult, self._BName    ) )
                self._nTagBbar[catIter] += self._tagTable.get( '{%s;%s}' % ( tagCatNameMult, self._BbarName ) )
                self._nTagB[0]    += self._nTagB[catIter]
                self._nTagBbar[0] += self._nTagBbar[catIter]

        # check order of types in tag (default: B-Bbar, reversed: Bbar-B)
        self._tagRevOrder = False
        for iter, catType in enumerate(iTag) :
            if iter == 0 and catType.getVal() < 0 : self._tagRevOrder = True

        # get tagging category coefficients
        self._tagCatCoefs = kwargs.pop( 'TagCatCoefs', None )
        if self._tagCatCoefs :
            self._relativeCatCoefs = False if not kwargs.pop( 'RelativeCoefs', True ) else True

        else :
            # create tagging category coefficients
            if self._data : from math import sqrt
            self._relativeCatCoefs = False
            self._tagCatCoefs = [ ]
            for coefIter in range( self._numTagCats ) :
                if self._data :
                    coefVal = ( self._nUntB + self._nUntBbar ) if coefIter == 0\
                              else ( self._nTagB[coefIter] + self._nTagBbar[coefIter] )
                    coefErr = sqrt(coefVal)
                    coefVal /= self._data.sumEntries()
                    coefErr /= self._data.sumEntries()
                else :
                    coefVal = ( 1. - float(coefIter) / float(self._numTagCats) ) / float(self._numTagCats)
                    coefErr = 0.01

                self._parseArg( self._namePF + 'tagCatCoef%d' % coefIter, kwargs, ContainerList = self._tagCatCoefs
                               , Title    = 'Tagging category coefficient %d' % coefIter
                               , Value    = coefVal
                               , Error    = coefErr
                               , MinMax   = ( 0., 1. )
                               , Constant = True if self._data else False
                              )

            if not self._data :
                untCoefVal = 1.
                for coef in self._tagCatCoefs[ 1 : ] : untCoefVal -= coef.getVal()
                self._tagCatCoefs[0].setVal(untCoefVal)

        # get tagging category coefficient names
        self._tagCatCoefNames = [ coef if type(coef) == str else coef.GetName() for coef in self._tagCatCoefs ]


class TagUntag_BinnedTaggingPdf( BinnedTaggingPdf ) :
    def __init__( self, Name, tagCat, iTag, **kwargs ) :
        # initialize
        self._init( Name, tagCat, iTag, kwargs )

        if self._data : from math import sqrt

        if self._relativeCatCoefs :
            # asymmetry between this untagged coefficient and provided untagged coefficient
            if self._data :
                AUntVal = float( self._nUntB + self._nUntBbar )
                AUntErr = sqrt(AUntVal)
                AUntVal = AUntVal / self._data.sumEntries() / self._ws[ self._tagCatCoefNames[0] ].getVal() - 1.
                AUntErr /=  self._data.sumEntries() * self._ws[ self._tagCatCoefNames[0] ].getVal()
            else :
                AUntVal = 0.
                AUntErr = 0.01

            self._parseArg( self._namePF + 'AUntagged', kwargs, Title  = 'Untagged-tagged asymmetry in tagging category coefficients'
                           , Value = AUntVal, Error = AUntErr, MinMax = ( -1., 1. ) )
            AUntagged = getattr( self, '_' + self._namePF + 'AUntagged' )

        # B-Bbar asymmetry for untagged events
        AUntBBbarErr = 1. / sqrt( float( self._nUntB + self._nUntBbar ) ) if self._data else 0.01
        self._parseArg( self._namePF + 'AUntBBbar', kwargs, Title = 'Untagged B-Bbar asymmetry'
                       , Value = 0., Error = AUntBBbarErr, MinMax = ( -1., 1. ) , Constant = True )
        AUntBBbar = getattr( self, '_' + self._namePF + 'AUntBBbar' )

        # B-Bbar asymmetry for tagged events
        ATagBBbarVal = float( self._nTagB[0] - self._nTagBbar[0]) / float( self._nTagB[0] + self._nTagBbar[0] )\
                       if self._data else 0.
        ATagBBbarErr = 1. / sqrt( float( self._nTagB[0] + self._nTagBbar[0] ) ) if self._data else 0.01
        self._parseArg( self._namePF + 'ATagBBbar', kwargs, Title = 'Tagged B-Bbar asymmetry'
                       , Value = ATagBBbarVal, Error = ATagBBbarErr , MinMax = ( -1., 1. ) )
        ATagBBbar = getattr( self, '_' + self._namePF + 'ATagBBbar' )

        # tagging category asymmetry factors
        from RooFitWrappers import FormulaVar
        if self._relativeCatCoefs :
            self._untaggedBbarCoef = FormulaVar(  self._namePF + ( 'untaggedBCoef' if self._tagRevOrder else 'untaggedBbarCoef' )
                                                , '0.5*(1.+@0)*(1.%s@1)' % ( '+' if self._tagRevOrder else '-' )
                                                , [ AUntagged, AUntBBbar ]
                                               )
            self._taggedBCoef      = FormulaVar(  self._namePF + ( 'taggedBbarCoef' if self._tagRevOrder else 'taggedBCoef' )
                                                , '0.5*(1.-@2/(1.-@2)*@0)*(1.%s@1)' % ( '-' if self._tagRevOrder else '+' )
                                                , [ AUntagged, ATagBBbar, self._ws[ self._tagCatCoefNames[0] ] ]
                                               )
            self._taggedBbarCoef   = FormulaVar(  self._namePF + ( 'taggedBCoef' if self._tagRevOrder else 'taggedBbarCoef' )
                                                , '0.5*(1.-@2/(1.-@2)*@0)*(1.%s@1)' % ( '+' if self._tagRevOrder else '-' )
                                                , [ AUntagged, ATagBBbar, self._ws[ self._tagCatCoefNames[0] ] ]
                                               )
        else :
            self._untaggedBbarCoef = FormulaVar(  self._namePF + ( 'untaggedBCoef' if self._tagRevOrder else 'untaggedBbarCoef' )
                                                , '0.5*(1.%s@0)' % ( '+' if self._tagRevOrder else '-' )
                                                , [ AUntBBbar ]
                                               )
            self._taggedBCoef      = FormulaVar(  self._namePF + ( 'taggedBbarCoef' if self._tagRevOrder else 'taggedBCoef' )
                                                , '0.5*(1.%s@0)' % ( '-' if self._tagRevOrder else '+' )
                                                , [ ATagBBbar ]
                                               )
            self._taggedBbarCoef   = FormulaVar(  self._namePF + ( 'taggedBCoef' if self._tagRevOrder else 'taggedBbarCoef' )
                                                , '0.5*(1.%s@0)' % ( '+' if self._tagRevOrder else '-' )
                                                , [ ATagBBbar ]
                                               )

        # tagging bin coefficients
        tagBinCoefs = [ ]
        from RooFitWrappers import Product
        for binIter in range( 1, self._numTagCats * 2 ) :
            cat = binIter % self._numTagCats

            # product of tagging category coefficient and bin asymmetry factor
            tagBinCoefs.append( Product(  self._namePF + 'tagBinCoef%d' % binIter
                                        , [ self._ws[ self._tagCatCoefNames[cat] ], self._untaggedBbarCoef if cat == 0\
                                            else ( self._taggedBCoef if binIter < self._numTagCats else self._taggedBbarCoef ) ]
                                       )
                              )

        # initialize TaggingPdf
        self._check_extraneous_kw( kwargs )
        BinnedTaggingPdf.__init__( self, self._name, self._tagCat, self._iTag, tagBinCoefs )


class TagCats_BinnedTaggingPdf( BinnedTaggingPdf ) :
    def __init__( self, Name, tagCat, iTag, **kwargs ) :
        # initialize
        self._init( Name, tagCat, iTag, kwargs )

        if self._data : from math import sqrt

        from RooFitWrappers import FormulaVar
        if self._relativeCatCoefs :
            # tagging category asymmetries
            self._ATagCats = [ ]
            for catIter in range( 1, self._numTagCats ) :
                if self._data :
                    ACatVal = float( self._nTagB[catIter] + self._nTagBbar[catIter] )
                    ACatErr = sqrt(ACatVal)
                    ACatVal = ACatVal / self._data.sumEntries() / self._ws[ self._tagCatCoefNames[catIter] ].getVal() - 1.
                    ACatErr /= self._data.sumEntries() * self._ws[ self._tagCatCoefNames[catIter] ].getVal()
                else :
                    ACatVal = 0.
                    ACatErr = 0.01

                self._parseArg( self._namePF + 'ATagCat%s' % catIter, kwargs, ContainerList = self._ATagCats
                               , Title  = 'Tagging category coefficient asymmetry %d' % catIter
                               , Value = ACatVal
                               , Error = ACatErr
                               , MinMax = ( -1., 1. )
                               , Constant = True if self._data else False
                              )

            self._ATagCats = [ FormulaVar(  self._namePF + 'ATagCat0'
                                          , '-1./@0*(%s)' % '+'.join( '@%d*@%d' % ( cat, cat + self._numTagCats - 1 )\
                                                                      for cat in range( 1, self._numTagCats ) )
                                          , [ self._ws[coefName] for coefName in self._tagCatCoefNames ] + self._ATagCats[ : ]
                                          , Title = 'Tagging category coefficient asymmetry 0'
                                         )
                             ] + self._ATagCats

        # B-Bbar asymmetries
        self._ABBbars = [ ]
        for catIter in range(self._numTagCats) :
            if catIter != 0 and self._data :
                ABBbarVal = float(self._nTagB[catIter] - self._nTagBbar[catIter]) / float(self._nTagB[catIter] + self._nTagBbar[catIter])
                ABBbarErr = 1. / sqrt( float( self._nTagB[catIter] + self._nTagBbar[catIter] ) )
            else :
                ABBbarVal = 0.
                ABBbarErr = 0.01

            self._parseArg(  self._namePF + 'ABBbar%d' % catIter, kwargs, ContainerList = self._ABBbars
                           , Title    = 'B-Bbar asymmetry tagging category %d' % catIter
                           , Value    = -ABBbarVal if self._tagRevOrder else ABBbarVal
                           , Error    = ABBbarErr
                           , MinMax   = ( -1., 1. )
                           , Constant = True if self._data else False
                          )

        # tagging bin coefficients
        tagBinCoefs = [ ]
        for binIter in range( 1, self._numTagCats * 2 ) :
            cat = binIter % self._numTagCats
            if self._relativeCatCoefs :
                tagBinCoefs.append( FormulaVar(  self._namePF + 'tagBinCoef%d' % binIter
                                               , '0.5*@0*(1+@1)*(1%s@2)' % ( '+' if binIter < self._numTagCats else '-' )
                                               , [ self._ws[ self._tagCatCoefNames[cat] ], self._ATagCats[cat], self._ABBbars[cat] ]
                                              )
                                  )
            else :
                tagBinCoefs.append( FormulaVar(  self._namePF + 'tagBinCoef%d' % binIter
                                               , '0.5*@0*(1%s@1)' % ( '+' if binIter < self._numTagCats else '-' )
                                               , [ self._ws[ self._tagCatCoefNames[cat] ], self._ABBbars[cat] ]
                                              )
                                  )

        # initialize TaggingPdf
        self._check_extraneous_kw( kwargs )
        BinnedTaggingPdf.__init__( self, self._name, self._tagCat, self._iTag, tagBinCoefs )
