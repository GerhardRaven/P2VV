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
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar,ConstVar
        self._parseArg( 'wTag', kwargs, Title = 'B wrong tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__( self
                              , Dilutions = [ FormulaVar(  'tagDilution', '1. - 2*@0 '
                                                         , [ self._wTag ], Title = 'Average tagging dilution' ) ]
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class Trivial_PerEventTaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar,ConstVar
        self._parseArg( 'wTag', kwargs, Title = 'B wrong tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
        self._parseArg( 'p0', kwargs, Title = 'p0  tagging parameter',    Value = 0.384, Constant = True, MinMax = ( 0., 0.5 ) )
        self._parseArg( 'p1', kwargs, Title = 'p1  tagging parameter',    Value = 1.037, Constant = True, MinMax = ( 0.8, 1.2 ) )
        self._parseArg( 'etaAverage', kwargs, Title = 'etaAverage  tagging parameter',    Value = 0.379, Constant = True, MinMax = ( 0., 0.5 ) )
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__( self
                              , Dilutions = [ FormulaVar(  'tagDilution', '1. - 2. * ( @2 + @3 * ( @0 - @1 ) ) '
                                                         , [ self._wTag, self._etaAverage, self._p0, self._p1 ], Title = 'Calibrated Per-Event Dilution' ) ]
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class Trivial_Dilution( TaggingParams ) :
    def __init__( self, **kwargs ) :
        from P2VVParameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        from RooFitWrappers import ConstVar
        TaggingParams.__init__( self
                              , Dilutions = kwargs.pop('Dilution')
                              , ADilWTags = [ ConstVar( Name = 'zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )
        self._check_extraneous_kw( kwargs )

class WTagsCoefAsyms_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar

        if 'NumTagCats' not in kwargs :
            # only one tagging category
            self._parseArg( 'wTag',    kwargs, Title = 'B wrong tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
            self._parseArg( 'wTagBar', kwargs, Title = 'Bbar wrong tag probability', Value = 0.25, MinMax = ( 0., 0.5 ) )

            if 'CEvenOdd' in kwargs :
                CEvenOdd = kwargs.pop('CEvenOdd')
            else :
                from P2VVParameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
                if 'AProd' in kwargs and 'ANorm' in kwargs :
                    AProdVal = kwargs.pop('AProd')
                    ANormVal = kwargs.pop('ANorm')
                    avgCOdd  = ( AProdVal + ANormVal ) / ( 1. + AProdVal * ANormVal )
                elif 'AvgCEven' in kwargs and 'AvgCOdd' in kwargs :
                    avgCEven = kwargs.pop('AvgCEven')
                    avgCOdd  = kwargs.pop('AvgCOdd') / avgCEven
                else :
                    avgCOdd = 0.
                CEvenOdd = Coefficients_CEvenOdd( avgCEven = 1., avgCOdd = avgCOdd )

            self._check_extraneous_kw( kwargs )
            TaggingParams.__init__( self
                                   , Dilutions = [ FormulaVar(  'tagDilution', '1. - @0 - @1'
                                                              , [ self._wTag, self._wTagBar ], Title = 'Average tagging dilution' ) ]
                                   , ADilWTags = [ FormulaVar(  'ADilWTag',    '(@0 - @1) / (1. - @0 - @1)'
                                                              , [ self._wTag, self._wTagBar ], Title = 'Dilution/wrong tag asymmetry' ) ]
                                   , CEvenOdds = [ CEvenOdd ]
                                  )

        else :
            # multiple tagging categories
            numTagCats = kwargs.pop('NumTagCats')
            if numTagCats < 1 :
                raise KeyError('WTagsCoefAsyms_TaggingParams: number of tagging categories must be greater than or equal to one')

            from RooFitWrappers import ConstVar, RealVar, FormulaVar
            tagCatCoefs    = [ ]
            dilutions      = [ ]
            ADilWTags      = [ ]
            CEvenOdds      = [ ]
            self._wTags    = [ ]
            self._wTagBars = [ ]

            # get wrong tag parameters category 0
            if 'wTag0' in kwargs and 'wTagBar0' in kwargs :
                startTagCat = 0
            else :
                startTagCat = 1
                self._wTags.append(None)
                self._wTagBars.append(None)
                dilutions.append( ConstVar( Name = 'tagDilution0', Value = 0. ) )
                ADilWTags.append( ConstVar( Name = 'ADilWTag0',    Value = 0. ) )

            # get average even and average odd coefficients of category 0
            from P2VVParameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
            if 'CEvenOddSum' in kwargs :
                CEvenOddSum = kwargs.pop('CEvenOddSum')
            elif 'AvgCEvenSum' in kwargs and 'AvgCOddSum' in kwargs :
                avgCEvenSum = kwargs.pop('AvgCEvenSum')
                avgCOddSum  = kwargs.pop('AvgCOddSum')
                CEvenOddSum = Coefficients_CEvenOdd( avgCEven = avgCEvenSum, avgCOdd = avgCOddSum )
            else :
                if 'AProd' in kwargs and 'ANorm' in kwargs :
                    self._AProdVal = kwargs.pop('AProd')
                    self._ANormVal = kwargs.pop('ANorm')
                    avgCOddSum  = ( self._AProdVal + self._ANormVal ) / ( 1. + self._AProdVal * self._ANormVal )
                else :
                    avgCOddSum = 0.
                CEvenOddSum = Coefficients_CEvenOdd(  avgCEven = ConstVar( Name = 'avgCEvenSum', Value = 1. )
                                                    , avgCOdd  = RealVar(  'avgCOddSum',  Value = avgCOddSum, MinMax = ( -2., 2. ) )
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
                        CEvenOdd = kwargs.pop('CEvenOdd%d' % index)
                    elif 'AvgCEven%d' % index in kwargs and 'AvgCOdd%d' % index in kwargs :
                        avgCEven = kwargs.pop('AvgCEven%d' % index)
                        avgCOdd  = kwargs.pop('AvgCOdd%d'  % index)
                        CEvenOdd = Coefficients_CEvenOdd( avgCEven = avgCEven, avgCOdd = avgCOdd )
                    else :
                        if 'ATagEff%d' % index in kwargs and hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                            ATagEffVal = kwargs.pop('ATagEff%d' % index)
                            avgCEven = 1. + self._AProdVal * self._ANormVal + self._AProdVal * ATagEffVal + self._ANormVal * ATagEffVal
                            avgCOdd  = self._AProdVal + self._ANormVal + ATagEffVal + self._AProdVal * self._ANormVal * ATagEffVal
                            avgCEven /= 1. + self._AProdVal * self._ANormVal
                            avgCOdd  /= 1. + self._AProdVal * self._ANormVal
                        else :
                            avgCEven = 1. / ( 1. + self._AProdVal * self._ANormVal )
                            avgCOdd  = 0.
                        CEvenOdd = Coefficients_CEvenOdd(  avgCEven = RealVar( 'avgCEven%d' % index, Value = avgCEven, MinMax = ( 0., 2.) )
                                                         , avgCOdd  = RealVar( 'avgCOdd%d'  % index, Value = avgCOdd,  MinMax = (-2., 2.) )
                                                        )
                    CEvenOdds.append(CEvenOdd)

                # get wrong tag parameters
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

            # check for remaining keyword arguments and initialize
            self._check_extraneous_kw( kwargs )
            TaggingParams.__init__(  self, NumTagCats = numTagCats
                                   , TagCatCoefs = tagCatCoefs
                                   , Dilutions   = dilutions
                                   , ADilWTags   = ADilWTags
                                   , CEvenOdds   = CEvenOdds
                                  )



from P2VVParameterizations.GeneralUtils import _util_parse_mixin
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
