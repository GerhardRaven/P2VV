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
                              , ADilWTags = [ ConstVar('zero',Value = 0) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

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

            from RooFitWrappers import RealVar, FormulaVar
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
                dilutions.append( RealVar( 'tagDilution0', Title = 'Average tagging dilution 0',     Value = 0. ) )
                ADilWTags.append( RealVar( 'ADilWTag0',    Title = 'Dilution/wrong tag asymmetry 0', Value = 0. ) )

            # get average even and average odd coefficients of category 0
            from P2VVParameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
            if 'CEvenOddSum' in kwargs :
                CEvenOddSum = kwargs.pop('CEvenOddSum')
            else :
                if 'AProd' in kwargs and 'ANorm' in kwargs :
                    self._AProdVal = kwargs.pop('AProd')
                    self._ANormVal = kwargs.pop('ANorm')
                    avgCOddSum  = ( self._AProdVal + self._ANormVal ) / ( 1. + self._AProdVal * self._ANormVal )
                elif 'AvgCEvenSum' in kwargs and 'AvgCOddSum' in kwargs :
                    avgCEvenSum = kwargs.pop('AvgCEvenSum')
                    avgCOddSum  = kwargs.pop('AvgCOddSum') / avgCEvenSum
                else :
                    avgCOddSum = 0.
                CEvenOddSum = Coefficients_CEvenOdd( avgCEven = 1., avgCOdd = avgCOddSum )

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
                    else :
                        if 'ATagEffVal%d' % index in kwargs and hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                            ATagEffVal = kwargs.pop('ATagEffVal%d' % index)
                            avgCOdd = ( self._AProdVal + self._ANormVal + ATagEffVal + self._AProdVal*self._ANormVal*ATagEffVal )\
                                      / ( 1. + self._AProdVal*self._ANormVal + self._AProdVal*ATagEffVal + self._ANormVal*ATagEffVal )
                        elif 'AvgCEven%d' % index in kwargs and 'AvgCOdd%d' % index in kwargs :
                            avgCEven = kwargs.pop('AvgCEven%d' % index)
                            avgCOdd  = kwargs.pop('AvgCOdd%d'  % index) / avgCEven
                        else :
                            avgCOdd = 0.
                        CEvenOdd = Coefficients_CEvenOdd( avgCEven = 1., avgCOdd = avgCOdd )

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

