###########################################################################################################################################
## P2VVParameterizations.FlavourTagging: Flavour tagging parameters                                                                      ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin

# initial values for tagging calibration parameters in fit
P0OSVal     =  0.38
P0OSErr     =  0.01
DelP0OSVal  =  0.014
DelP0OSErr  =  0.001
P1OSVal     =  1.00
P1OSErr     =  0.06
DelP1OSVal  =  0.07
DelP1OSErr  =  0.01

P0SSVal     =  0.44
P0SSErr     =  0.01
DelP0SSVal  = -0.016
DelP0SSErr  =  0.002
P1SSVal     =  1.0
P1SSErr     =  0.1
DelP1SSVal  =  0.02
DelP1SSErr  =  0.02

# values for tagging calibration parameters from tagging calibration: correspond to specific data sample!
from P2VV.Imports import extConstraintValues

( P0OSConstrVal, P0OSConstrErr, avgEtaOSVal ) = extConstraintValues.getSetVal( 'P0OS',    (  0.379,  0.011, 0.379 ) )
( DelP0OSConstrVal, DelP0OSConstrErr )        = extConstraintValues.getSetVal( 'DelP0OS', (  0.0137, 0.0012 ) )
( P1OSConstrVal,    P1OSConstrErr    )        = extConstraintValues.getSetVal( 'P1OS',    (  1.00,   0.06   ) )
( DelP1OSConstrVal, DelP1OSConstrErr )        = extConstraintValues.getSetVal( 'DelP1OS', (  0.070,  0.012  ) )

( P0SSConstrVal, P0SSConstrErr, avgEtaSSVal ) = extConstraintValues.getSetVal( 'P0SS',    (  0.437, 0.008, 0.437 ) )
( DelP0SSConstrVal, DelP0SSConstrErr )        = extConstraintValues.getSetVal( 'DelP0SS', ( -0.016, 0.002 ) )
( P1SSConstrVal,    P1SSConstrErr    )        = extConstraintValues.getSetVal( 'P1SS',    (  1.00,  0.12  ) )
( DelP1SSConstrVal, DelP1SSConstrErr )        = extConstraintValues.getSetVal( 'DelP1SS', (  0.015, 0.019 ) )

from ROOT import RooNumber
RooInf = RooNumber.infinity()

###########################################################################################################################################
## Tagging tools ##
###################

def getTagCatParamsFromData( data, estWTagName, tagCats = [ ], numSigmas = 1., SameSide = False, WeightVarName = ''
                            , avgEstWTag = None, P0 = None, P1 = None, P0Err = None, P1Err = None, AP0 = None, AP1 = None ) :
    assert data, 'getTagCatParamsFromData(): no data set found'
    assert estWTagName and data.get().find(estWTagName), 'getTagCatParamsFromData(): estimated wrong-tag probability not found in data'
    if WeightVarName : assert data.get().find(WeightVarName), 'getTagCatParamsFromData(): weight variable not found in data'

    if avgEstWTag == None : avgEstWTag = avgEtaOSVal                         if not SameSide else avgEtaSSVal
    if P0         == None : P0         = P0OSConstrVal                       if not SameSide else P0SSConstrVal
    if P1         == None : P1         = P1OSConstrVal                       if not SameSide else P1SSConstrVal
    if P0Err      == None : P0Err      = P0OSConstrErr                       if not SameSide else P0SSConstrErr
    if P1Err      == None : P1Err      = P1OSConstrErr                       if not SameSide else P1SSConstrErr
    if AP0        == None : AP0        = (DelP0OSConstrVal/2./P0OSConstrVal) if not SameSide else (DelP0SSConstrVal/2./P0SSConstrVal)
    if AP1        == None : AP1        = (DelP1OSConstrVal/2./P1OSConstrVal) if not SameSide else (DelP1SSConstrVal/2./P1SSConstrVal)

    from P2VV.RooFitWrappers import RooObject
    if isinstance( avgEstWTag, RooObject ) : avgEstWTag = avgEstWTag.getVal()
    if isinstance( P0,         RooObject ) : P0         = P0.getVal()
    if isinstance( P1,         RooObject ) : P1         = P1.getVal()
    if isinstance( AP0,        RooObject ) : AP0        = AP0.getVal()
    if isinstance( AP1,        RooObject ) : AP1        = AP1.getVal()

    tagCatsCalc = tagCats[ : ]

    etaMin = 0.5
    if not tagCatsCalc :
        # get minimum estimated wrong-tag probability
        for varSet in data : etaMin = min( etaMin, varSet.getRealValue(estWTagName) )

        # determine binning in estimated wrong-tag probability from data
        bin = 0
        binUpperEdges = [ 0.499999999 ]
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
        tagCatsCalc = [ ( 'Untagged', 0, 0.5000001 ) ]
        for bin in range( 1, len(binUpperEdges) ) :
            binUpperEdges[bin] = binUpperEdges[0] - ( binUpperEdges[0] - binUpperEdges[bin] ) * binScale
            tagCatsCalc.append( ( 'TagCat%d' % bin, bin, binUpperEdges[ bin - 1 ] ) )

    # determine tagging category parameters
    numTagCats = len(tagCatsCalc)
    numEvCats  = [0]  * numTagCats
    sumEtaCats = [0.] * numTagCats
    for varSet in data :
        # get estimated wrong-tag probability for this event
        eta = varSet.getRealValue(estWTagName)
        if eta < etaMin : etaMin = eta

        # determine tagging category
        cat = -1
        for catIter in range(numTagCats) :
          if eta >= tagCatsCalc[catIter][2] : break
          cat += 1

        if cat < 0 : raise RuntimeError('getTagCatParamsFromData(): estimated wrong-tag probability out of range')

        # update number of events and sum of estimated wrong-tag probabilities
        weight = varSet.getRealValue(WeightVarName) if WeightVarName else data.weight()
        numEvCats[cat]  += weight
        sumEtaCats[cat] += eta * weight

    # get total number of events
    numEvTot = sum(numEvCats)
    etaMean  = sum(sumEtaCats) / ( numEvTot if numEvTot else 1. )

    # update tagging category parameters
    for cat in range(numTagCats) :
        avgEtaCat = ( sumEtaCats[cat] / numEvCats[cat] ) if numEvCats[cat] != 0. else 0.
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
    print '    minimum eta = %.3f    average eta = %.3f' % ( etaMin, etaMean )
    for bin, cat in enumerate(tagCatsCalc) :
        deltaEta = ( ( cat[2] - tagCatsCalc[bin + 1][2] ) / 2. ) if bin < len(tagCatsCalc) - 1 else ( ( cat[2] - etaMin ) / 2. )
        if cat[2] >= avgEstWTag : binRangeSig = P1 * deltaEta / ( P0Err + P1Err * ( cat[2] - avgEstWTag - deltaEta ) )
        else                    : binRangeSig = P1 * deltaEta / ( P0Err - P1Err * ( cat[2] - avgEstWTag - deltaEta ) )
        print '    {0:<10s}  :  {1:.3f} -- {2:.3f} ({3:.2f} sigma)  :  <eta> = {4:.3f}  <w> = {5:.3f}  efficiency = {6:.4f} ({7:.0f} events)'\
              .format(  cat[0], cat[2], tagCatsCalc[bin + 1][2] if bin < len(tagCatsCalc) - 1 else etaMin
                      , binRangeSig, cat[3], cat[4], cat[6], cat[8]
                     )

    return tagCatsCalc


class getTaggingPerformance( ) :
    pass


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

        if self._numTagCats[0] > 0 or self._numTagCats[1] > 1 :
            self._tagCatCoefs = kwargs.pop('TagCatCoefs')
            if self._numTagCats[0] == 0 : self._tagCatCoefs = [ self._tagCatCoefs ]

        # cache integrals as a function of observables
        for dils in self._dilutions + self._ADilWTags :
            for dil in dils : dil.setAttribute('CacheAndTrack')

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
            return self._dilutions[1][0] if self._numTagCats[0] == 0 and self._numTagCats[1] == 1 else raiseError(kw)
        if kw == 'dilutions' :
            return self._dilutions[1] if self._numTagCats[0] == 0 else self._dilutions
        if kw == 'ADilWTag' :
            return self._ADilWTags[1][0] if self._numTagCats[0] == 0 and self._numTagCats[1] == 1 else raiseError(kw)
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
        from P2VV.Parameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        TaggingParams.__init__( self
                              , Dilutions = [ self._parseArg( 'one',  { }, Value = 1, ObjectType = 'ConstVar' ) ]
                              , ADilWTags = [ self._parseArg( 'zero', { }, Value = 0, ObjectType = 'ConstVar' ) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class WTag_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'wTag', kwargs, Title = 'B wrong-tag probability', Value = 0.25, MinMax = ( 0., 0.5 ) )

        self._check_extraneous_kw( kwargs )
        from P2VV.Parameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        TaggingParams.__init__( self
                              , Dilutions = [ self._parseArg( 'tagDilution', kwargs, Formula = '1. - 2*@0 ', Arguments = [ self._wTag ]
                                             , ObjectType = 'FormulaVar', Title = 'Tagging dilution' ) ]
                              , ADilWTags = [ self._parseArg( 'zero', kwargs, Value = 0, ObjectType = 'ConstVar' ) ]
                              , CEvenOdds = [ Trivial_CEvenOdd() ]
                              )

class LinearEstWTag_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'estWTag',    kwargs, Title = 'Estimated wrong-tag probability',         Value = avgEtaOSVal, MinMax = (0., 0.5) )
        self._parseArg( 'avgEstWTag', kwargs, Title = 'Average estimated wrong-tag probability', Value = avgEtaOSVal, MinMax = (0., 0.5)
                       , Constant = True )
        self._parseArg( 'p0', kwargs, Title = 'p0  tagging parameter', Value = P0OSVal, Error = P0OSErr, MinMax = ( 0.,  1. ) )
        self._parseArg( 'p1', kwargs, Title = 'p1  tagging parameter', Value = P1OSVal, Error = P1OSErr, MinMax = ( 0.5, 2. ) )

        constraints = set()
        p0Constr = kwargs.pop( 'p0Constraint', None )
        if type(p0Constr) == str and p0Constr == 'fixed' :
            self._p0.setConstant(True)
            self._p0.setVal(P0OSConstrVal)
            self._p0.setError(P0OSConstrErr)

        elif p0Constr :
            from P2VV.RooFitWrappers import Pdf
            from ROOT import RooGaussian as Gaussian
            constraints.add( Pdf(  Name = self._p0.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._p0
                                                    , self._parseArg( 'p0_mean',  kwargs, Value = P0OSConstrVal, ObjectType = 'ConstVar' )
                                                    , self._parseArg( 'p0_sigma', kwargs, Value = P0OSConstrErr, ObjectType = 'ConstVar' )
                                                   ]
                                   )
                              )

        p1Constr = kwargs.pop( 'p1Constraint', None )
        if type(p1Constr) == str and p1Constr == 'fixed' :
            self._p1.setConstant(True)
            self._p1.setVal(P1OSConstrVal)
            self._p1.setError(P1OSConstrErr)

        elif p1Constr :
            from P2VV.RooFitWrappers import Pdf
            from ROOT import RooGaussian as Gaussian
            constraints.add( Pdf(  Name = self._p1.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._p1
                                                    , self._parseArg( 'p1_mean',  kwargs, Value = P1OSConstrVal, ObjectType = 'ConstVar' )
                                                    , self._parseArg( 'p1_sigma', kwargs, Value = P1OSConstrErr, ObjectType = 'ConstVar' )
                                                   ]
                                   )
                              )

        self._check_extraneous_kw( kwargs )
        from P2VV.Parameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        TaggingParams.__init__(  self
                                , Dilutions    = [ self._parseArg( 'tagDilution', kwargs, 'CalibratedDilution'
                                                                  , EstWTag    = self._estWTag
                                                                  , AvgEstWTag = self._avgEstWTag
                                                                  , P0         = self._p0
                                                                  , P1         = self._p1
                                                                 )
                                                 ]
                                , ADilWTags    = [ self._parseArg( 'zero', kwargs, Value = 0, ObjectType = 'ConstVar' ) ]
                                , CEvenOdds    = [ Trivial_CEvenOdd() ]
                                , Conditionals = [ self._estWTag ]
                                , Constraints  = constraints
                              )

class Dilution_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        self._parseArg( 'dilution', kwargs, Title = 'Tagging dilution', Value = 0.5, MinMax = ( 0., 1. ) )

        self._check_extraneous_kw( kwargs )
        from P2VV.Parameterizations.BBbarAsymmetries import Trivial_CEvenOdd
        TaggingParams.__init__( self
                              , Dilutions = [ self._dilution ]
                              , ADilWTags = [ self._parseArg( 'zero', kwargs, Value = 0, ObjectType = 'ConstVar' ) ]
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
            from P2VV.RooFitWrappers import RooObject
            if 'AvgCEven' in kwargs and 'AvgCOdd' in kwargs :
                # (values for the) coefficients are specified
                avgCEven = kwargs.pop('AvgCEven')
                avgCOdd  = kwargs.pop('AvgCOdd')

            elif 'AProd' in kwargs and 'ANorm' in kwargs :
                # values for the production and normalization asymmetries are specified
                self._AProdVal = kwargs.pop('AProd')
                self._ANormVal = kwargs.pop('ANorm')
                if isinstance( self._AProdVal, RooObject ) : self._AProdVal = self._AProdVal.getVal()
                if isinstance( self._ANormVal, RooObject ) : self._ANormVal = self._ANormVal.getVal()

                avgCEven = 1. + self._AProdVal * self._ANormVal
                avgCOdd  = self._AProdVal + self._ANormVal

            else :
                # use default values
                avgCEven = 1.
                avgCOdd  = 0.

            if not isinstance( avgCEven, RooObject ) and not isinstance( avgCOdd, RooObject ) :
                avgCOdd  = dict( Name = 'avgCOdd',  Value = avgCOdd / avgCEven, Constant = True )
                avgCEven = dict( Name = 'avgCEven', Value = 1., ObjectType = 'ConstVar' )

            # create coefficients object
            from P2VV.Parameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
            CEvenOdd = Coefficients_CEvenOdd( avgCEven = avgCEven, avgCOdd = avgCOdd )

        # initialize and check for remaining arguments
        if hasattr( self, 'wTag' ) and hasattr( self, 'wTagBar' ) :
            TaggingParams.__init__( self
                                   , Dilutions = [ self._parseArg( 'tagDilution', kwargs, Title = 'Tagging dilution'
                                                                  , Coefficients = [ self._wTag, self._wTagBar ]
                                                                  , ObjectType = 'ComplementCoef' ) ]
                                   , ADilWTags = [ self._parseArg( 'ADilWTag', kwargs, Formula = '(@0 - @1) / (1. - @0 - @1)'
                                                                  , Arguments = [ self._wTag, self._wTagBar ], ObjectType = 'FormulaVar'
                                                                  , Title = 'Dilution/wrong-tag asymmetry' ) ]
                                   , CEvenOdds = [ CEvenOdd ]
                                  )
        else :
            TaggingParams.__init__( self
                                   , Dilutions = [ self._parseArg( 'tagDilution', kwargs, Formula = '1. - 2. * @0'
                                                   , Arguments = [ self._WTag ], ObjectType = 'FormulaVar', Title = 'Tagging dilution' ) ]
                                   , ADilWTags = [ self._parseArg( 'ADilWTag', kwargs, Formula = '2. * @0 * @1 / (1. - 2. * @0)'
                                                                  , Arguments = [ self._WTag, self._AWTag ], ObjectType = 'FormulaVar'
                                                                  , Title = 'Dilution/wrong-tag asymmetry' ) ]
                                   , CEvenOdds = [ CEvenOdd ]
                                  )
        self._check_extraneous_kw( kwargs )


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
                dilutions[0].append( self._parseArg( 'tagDilution0_0', kwargs, Value = 0., ObjectType = 'ConstVar' ) )
                ADilWTags[0].append( self._parseArg( 'ADilWTag0_0',    kwargs, Value = 0., ObjectType = 'ConstVar' ) )

        startTagCat1 = 1
        if 'tagDilution%s0' % namePF in kwargs and 'ADilWTag%s0' % namePF in kwargs :
            startTagCat1 = 0
        else :
            dilutions[1].append( self._parseArg( 'tagDilution%s0' % namePF, kwargs, Value = 0., ObjectType = 'ConstVar' ) )
            ADilWTags[1].append( self._parseArg( 'ADilWTag%s0'    % namePF, kwargs, Value = 0., ObjectType = 'ConstVar' ) )

        # get average even and average odd coefficients of category 0
        if 'CEvenOddSum' in kwargs :
            # a coefficients object is specified
            CEvenOddSum = kwargs.pop('CEvenOddSum')

        else :
            from P2VV.RooFitWrappers import RooObject
            if 'AvgCEvenSum' in kwargs and 'AvgCOddSum' in kwargs :
                # (values for the) coefficients are specified
                avgCEvenSum = kwargs.pop('AvgCEvenSum')
                avgCOddSum  = kwargs.pop('AvgCOddSum')

            elif 'AProd' in kwargs and 'ANorm' in kwargs :
                # values for the production and normalization asymmetries are specified
                self._AProdVal = kwargs.pop('AProd')
                self._ANormVal = kwargs.pop('ANorm')
                if isinstance( self._AProdVal, RooObject ) : self._AProdVal = self._AProdVal.getVal()
                if isinstance( self._ANormVal, RooObject ) : self._ANormVal = self._ANormVal.getVal()

                avgCEvenSum = 1. + self._AProdVal * self._ANormVal
                avgCOddSum  = self._AProdVal + self._ANormVal

            else :
                # use default values
                avgCEvenSum = 1.
                avgCOddSum  = 0.

            if not isinstance( avgCEvenSum, RooObject ) and not isinstance( avgCOddSum, RooObject ) :
                avgCOddSum  = dict( Name = 'avgCOddSum',  Value = avgCOddSum / avgCEvenSum, Constant = True )
                avgCEvenSum = dict( Name = 'avgCEvenSum', Value = 1., ObjectType = 'ConstVar' )

            from P2VV.Parameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
            CEvenOddSum = Coefficients_CEvenOdd( avgCEven = avgCEvenSum, avgCOdd = avgCOddSum )

        CEvenOdds[0].append(CEvenOddSum)

        # get tagging category coefficients and asymmetries (assume factorization)
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

            self._singleTagCatCoefs[0].insert( 0, self._parseArg( 'tagCatCoef0_0', kwargs, Title = 'Tagging category coefficient 0 0'
                                                                 , Coefficients = self._singleTagCatCoefs[0]
                                                                 , ObjectType = 'ComplementCoef' ) )

            if hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                tagCatProd0 = 0.
                for tagCatCoef, ATagEffVal in zip( self._singleTagCatCoefs[0][ 1 : ], self._ATagEffVals[0] ) :
                    tagCatProd0 -= tagCatCoef.getVal() * ATagEffVal
                self._ATagEffVals[0].insert( 0, tagCatProd0 / self._singleTagCatCoefs[0][0].getVal() )

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
                self._singleTagCatCoefs[1].insert( 0, self._parseArg( 'tagCatCoef%s0' % namePF, kwargs, ObjectType = 'ComplementCoef'
                                                                     , Title = 'Tagging category coefficient 1 0'
                                                                     , Coefficients = self._singleTagCatCoefs[1]
                                                                    )
                                                 )

                if hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                    tagCatProd0 = 0.
                    for tagCatCoef, ATagEffVal in zip( self._singleTagCatCoefs[1][ 1 : ], self._ATagEffVals[1] ) :
                        tagCatProd0 -= tagCatCoef.getVal() * ATagEffVal
                    self._ATagEffVals[1].insert( 0, tagCatProd0 / self._singleTagCatCoefs[1][0].getVal() )

            else :
                self._singleTagCatCoefs[1].insert( 0, None )
                self._ATagEffVals[1].insert( 0, None )

        else :
            self._ATagEffVals[1].append(0)
            self._singleTagCatCoefs[1].insert( 0, self._parseArg( 'tagCatCoef%s0' % namePF, kwargs, Value = 1., ObjectType = 'ConstVar' ) )

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
                                tagCatCoef = self._parseArg( 'tagCatCoef%d-%d' % ( index0, index1 ), kwargs, ObjectType = 'Product'
                                                            , Arguments = [ self._singleTagCatCoefs[0][index0]
                                                                           , self._singleTagCatCoefs[1][index1] ]
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

                        elif hasattr(self, '_AProdVal') and hasattr(self, '_ANormVal') :
                            # values for the asymmetries are specified
                            avgCEven = 1. + self._AProdVal * self._ANormVal \
                                          + self._AProdVal * self._ATagEffVals[0][index0] + self._ANormVal * self._ATagEffVals[0][index0] \
                                          + self._AProdVal * self._ATagEffVals[1][index1] + self._ANormVal * self._ATagEffVals[1][index1] \
                                          + self._AProdVal * self._ANormVal * self._ATagEffVals[0][index0] * self._ATagEffVals[1][index1]
                            avgCOdd  = self._AProdVal + self._ANormVal + self._ATagEffVals[0][index0] + self._ATagEffVals[1][index1] \
                                       + self._AProdVal * self._ANormVal * self._ATagEffVals[0][index0] \
                                       + self._AProdVal * self._ANormVal * self._ATagEffVals[1][index1]

                        else :
                            # use values for tagging efficiency asymmetry = 0
                            avgCEven = CEvenOddSum['avgCEven'].getVal()
                            avgCOdd  = CEvenOddSum['avgCOdd'].getVal()

                        if not isinstance( avgCEven, RooObject ) and not isinstance( avgCOdd, RooObject )\
                                and hasattr( self, '_AProdVal' ) and hasattr( self, '_ANormVal' ) :
                            avgCEven = dict( Name = 'avgCEven%d-%d' % ( index0, index1 )
                                            , Value = avgCOdd  / ( 1. + self._AProdVal * self._ANormVal ), Constant = True )
                            avgCOdd  = dict( Name = 'avgCOdd%d-%d'  % ( index0, index1 )
                                            , Value = avgCEven / ( 1. + self._AProdVal * self._ANormVal ), Constant = True )

                        CEvenOdd = Coefficients_CEvenOdd( avgCEven = avgCEven, avgCOdd = avgCOdd )

                    CEvenOdds[index0].append(CEvenOdd)

        # get conditional observables
        conditionals = kwargs.pop( 'Conditionals', [] )

        # get external constraints
        constraints = kwargs.pop( 'Constraints', set() )

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

class TaggingCategories( _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin ) :
    def __init__( self, **kwargs ) :
        self._numTagCats = kwargs.pop('NumTagCats')
        if hasattr( self, '_combined' ) and self._combined :
            self._parseArg(  'tagCat0', kwargs, ObjectType = 'Category', SingleArgKey = 'Name'
                           , Title = 'Tagging Category 0', Observable = True
                           , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, self._numTagCats[0] ) ]
                          )
            assert self._tagCat0.numTypes() == self._numTagCats[0]\
                   , 'TaggingCategories: tagging category "%s" has %d states, while %d were expected'\
                     % ( self._tagCat0.GetName(), self._tagCat0.numTypes(), self._numTagCats[0] )

            self._parseArg(  'tagCat1', kwargs, ObjectType = 'Category', SingleArgKey = 'Name'
                           , Title = 'Tagging Category 1', Observable = True
                           , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, self._numTagCats[1] ) ]
                          )
            assert self._tagCat1.numTypes() == self._numTagCats[1]\
                   , 'TaggingCategories: tagging category "%s" has %d states, while %d were expected'\
                     % ( self._tagCat1.GetName(), self._tagCat1.numTypes(), self._numTagCats[1] )

        else :
            self._parseArg(  'tagCat', kwargs, ObjectType = 'Category', SingleArgKey = 'Name'
                           , Title = 'Tagging Category', Observable = True
                           , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, self._numTagCats ) ]
                          )
            assert self._tagCat.numTypes() == self._numTagCats\
                   , 'TaggingCategories: tagging category "%s" has %d states, while %d were expected'\
                     % ( self._tagCat.GetName(), self._tagCat.numTypes(), self._numTagCats )

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
                    + [ ( 'Conditionals', self.ConditionalObservables() ) ]
                    + [ ( 'Constraints',  self.ExternalConstraints()    ) ]
                   )


class Independent_TaggingCategories( TaggingCategories ) :
    def __init__( self, **kwargs ) :
        # get tagging category parameters
        tagCats = kwargs.pop( 'TagCats', [ ] )

        # get number of tagging categories
        if not tagCats and 'NumTagCats' not in kwargs :
            raise KeyError('Independent_TaggingCategories: did not find "NumTagCats" argument')
        else :
            numTagCats = kwargs.pop( 'NumTagCats', len(tagCats) )

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

        # set default category parameters
        catParams = [ getCatParam(params) for params in [ 'TagCatCoefs', 'ATagEffs', 'TagDilutions', 'ADilWTags' ] ]
        for cat in range( numTagCats - 1 ) :
            catDil = 1. - 2. * tagCats[ cat + 1 ][4] if tagCats else float(cat + 1) / float(numTagCats)
            catCoef = tagCats[ cat + 1 ][6] if tagCats else 1.
            if len(catParams[0]) == cat :
                if not tagCats :
                    from math import pow
                    numCatsFrac = float(numTagCats) / 6.
                    catCoef  = 0.15 / numCatsFrac * pow( 0.5, float(cat) / numCatsFrac )
                catParams[0].append(catCoef)

            if len(catParams[1]) == cat : catParams[1].append(0.)
            if len(catParams[2]) == cat : catParams[2].append(catDil)
            if len(catParams[3]) == cat : catParams[3].append(0.)
        if len(catParams[0]) < numTagCats : catParams[0].insert( 0, 1. - sum(catParams[0]) )
        if len(catParams[1]) < numTagCats : catParams[1].insert( 0, 0. )

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
        if 'estWTag' in kwargs and kwargs['estWTag'] :
            self._parseArg( 'estWTag', kwargs, Name = 'estWTag' + tagType, Title = 'Estimated wrong-tag probability'
                           , Value = 0.25, MinMax = ( 0., 0.5 ) )
        else :
            self._estWTag = kwargs.pop( 'estWTag', None )

        # get name of estimated wrong-tag probability variable
        self._etaName = kwargs.pop( 'estWTagName', self._estWTag.GetName() if self._estWTag else '' )

        # set values for calibration parameters
        self._calVals = {    'avgEta'   : avgEtaOSVal      if tagType == 'OS' else avgEtaSSVal
                           , 'P0'       : P0OSConstrVal    if tagType == 'OS' else P0SSConstrVal
                           , 'P0Err'    : P0OSConstrErr    if tagType == 'OS' else P0SSConstrErr
                           , 'P1'       : P1OSConstrVal    if tagType == 'OS' else P1SSConstrVal
                           , 'P1Err'    : P1OSConstrErr    if tagType == 'OS' else P1SSConstrErr
                           , 'DelP0'    : DelP0OSConstrVal if tagType == 'OS' else DelP0SSConstrVal
                           , 'DelP0Err' : DelP0OSConstrErr if tagType == 'OS' else DelP0SSConstrErr
                           , 'DelP1'    : DelP1OSConstrVal if tagType == 'OS' else DelP1SSConstrVal
                           , 'DelP1Err' : DelP1OSConstrErr if tagType == 'OS' else DelP1SSConstrErr
                           , 'AP0'      : DelP0OSConstrVal / 2. / P0OSConstrVal if tagType == 'OS'\
                                          else DelP0SSConstrVal / 2. / P0SSConstrVal
                           , 'AP1'      : DelP1OSConstrVal / 2. / P1OSConstrVal if tagType == 'OS'\
                                          else DelP1SSConstrVal / 2. / P1SSConstrVal
                          }

        # get linear calibration parameters
        self._parseArg(  'avgEstWTag', kwargs, Name = 'avgEstWTag' + tagType
                       , Value = avgEtaOSVal if tagType == 'OS' else avgEtaSSVal
                       , ObjectType = 'ConstVar' )
        self._parseArg(  'wTagP0',     kwargs, Name = 'wTagP0' + tagType, Title = 'Average wrong-tag parameter p_0'
                       , Value = P0OSVal if tagType == 'OS' else P0SSVal
                       , Error = P0OSErr if tagType == 'OS' else P0SSErr
                       , MinMax = ( 0., 1. )
                      )
        self._parseArg(  'wTagP1', kwargs, Name = 'wTagP1' + tagType, Title = 'Average wrong-tag parameter p_1'
                       , Value = P1OSVal if tagType == 'OS' else P1SSVal
                       , Error = P1OSErr if tagType == 'OS' else P1SSErr
                       , MinMax = ( 0.5, 2. )
                      )
        self._parseArg(  'wTagDelP0', kwargs, Name = 'wTagDelP0' + tagType, Title = 'Wrong tag parameter p_0 difference'
                       , Value = DelP0OSVal if tagType == 'OS' else DelP0SSVal
                       , Error = DelP0OSErr if tagType == 'OS' else DelP0SSErr
                       , MinMax = ( -0.5, 0.5 )
                      )
        self._parseArg(  'wTagDelP1', kwargs, Name = 'wTagDelP1' + tagType, Title = 'Wrong tag parameter p_1 difference'
                       , Value = DelP1OSVal if tagType == 'OS' else DelP1SSVal
                       , Error = DelP1OSErr if tagType == 'OS' else DelP1SSErr
                       , MinMax = ( -0.5, 0.5 )
                      )

        self._wTagAP0 = self._parseArg( 'wTagAP0' + tagType, kwargs, Formula = '@0/2./@1', Arguments = [ self._wTagDelP0, self._wTagP0 ]
                                       , ObjectType = 'FormulaVar' )
        self._wTagAP1 = self._parseArg( 'wTagAP1' + tagType, kwargs, Formula = '@0/2./@1', Arguments = [ self._wTagDelP1, self._wTagP1 ]
                                       , ObjectType = 'FormulaVar' )

        # constrain calibration parameters
        constraints = set()
        wTagP0Constraint = kwargs.pop( 'wTagP0Constraint', None )
        if type(wTagP0Constraint) == str and wTagP0Constraint == 'fixed' :
            self._wTagP0.setConstant(True)
            self._wTagP0.setVal(self._calVals['P0'])
            self._wTagP0.setError(self._calVals['P0Err'])
            self._wTagDelP0.setConstant(True)
            self._wTagDelP0.setVal(self._calVals['DelP0'])
            self._wTagDelP0.setError(self._calVals['DelP0Err'])

        elif wTagP0Constraint :
            from P2VV.RooFitWrappers import Pdf
            from ROOT import RooGaussian as Gaussian
            constraints.add( Pdf(  Name = self._wTagP0.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._wTagP0
                                                    , self._parseArg( 'wTagP0' + tagType + '_mean', kwargs
                                                                      , Value = self._calVals['P0'], ObjectType = 'ConstVar' )
                                                    , self._parseArg( 'wTagP0' + tagType + '_sigma', kwargs
                                                                      , Value = self._calVals['P0Err'], ObjectType = 'ConstVar' )
                                                   ]
                                   )
                              )
            constraints.add( Pdf(  Name = self._wTagDelP0.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._wTagDelP0
                                                    , self._parseArg( 'wTagDelP0' + tagType + '_mean', kwargs
                                                                      , Value = self._calVals['DelP0'], ObjectType = 'ConstVar' )
                                                    , self._parseArg( 'wTagDelP0' + tagType + '_sigma', kwargs
                                                                      , Value = self._calVals['DelP0Err'], ObjectType = 'ConstVar' )
                                                   ]
                                   )
                              )

        wTagP1Constraint = kwargs.pop( 'wTagP1Constraint', None )
        if type(wTagP1Constraint) == str and wTagP1Constraint == 'fixed' :
            self._wTagP1.setConstant(True)
            self._wTagP1.setVal(self._calVals['P1'])
            self._wTagP1.setError(self._calVals['P1Err'])
            self._wTagDelP1.setConstant(True)
            self._wTagDelP1.setVal(self._calVals['DelP1'])
            self._wTagDelP1.setError(self._calVals['DelP1Err'])

        elif wTagP1Constraint :
            from P2VV.RooFitWrappers import Pdf
            from ROOT import RooGaussian as Gaussian
            constraints.add( Pdf(  Name = self._wTagP1.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._wTagP1
                                                    , self._parseArg( 'wTagP1' + tagType + '_mean', kwargs
                                                                     , Value = self._calVals['P1'], ObjectType = 'ConstVar' )
                                                    , self._parseArg( 'wTagP1' + tagType + '_sigma', kwargs
                                                                     , Value = self._calVals['P1Err'], ObjectType = 'ConstVar' )
                                                   ]
                                   )
                              )
            constraints.add( Pdf(  Name = self._wTagDelP1.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._wTagDelP1
                                                   , self._parseArg( 'wTagDelP1' + tagType + '_mean', kwargs
                                                                    , Value = self._calVals['DelP1'], ObjectType = 'ConstVar' )
                                                   , self._parseArg( 'wTagDelP1' + tagType + '_sigma', kwargs
                                                                    , Value = self._calVals['DelP1Err'], ObjectType = 'ConstVar' )
                                                   ]
                                   )
                              )

        # get tagging category binning in estimated wrong-tag probability (eta)
        self._tagCats = kwargs.pop( 'TagCats', None )
        if not self._tagCats :
            self._tagCats = [  ( 'Untagged', 0, 0.500000001, 0.50, 0.50, 0., 0.65,  0. )
                             , ( 'Tagged',   1, 0.499999999, 0.40, 0.40, 0., 0.35,  0. )
                            ] if self._estWTag else \
                            [  ( 'Untagged', 0, 0.500000001, 0.50, 0.50, 0., 0.65,  0. )
                             , ( 'TagCat1',  1, 0.499999999, 0.44, 0.44, 0., 0.24,  0. )
                             , ( 'TagCat2',  2, 0.38,      0.35, 0.35, 0., 0.062, 0. )
                             , ( 'TagCat3',  3, 0.31,      0.28, 0.28, 0., 0.032, 0. )
                             , ( 'TagCat4',  4, 0.24,      0.21, 0.21, 0., 0.012, 0. )
                             , ( 'TagCat5',  5, 0.17,      0.15, 0.14, 0., 0.004, 0. )
                            ]

        # get number of calibration parameter standard deviations for tagging category bins
        nSigmaTagBins = kwargs.pop( 'NumSigmaTagBins', 1. )

        # initialize category parameters
        tagCatCoefs    = [ catPars[6] for catPars in self._tagCats ]
        ATagEffs       = [ catPars[7] for catPars in self._tagCats ]
        dilutions      = [ ]
        ADilWTags      = [ ]
        self._estWTags = [ ]
        for cat, catPars in enumerate( self._tagCats[ 1 : ] ) :
            if self._estWTag :
                estWTag = self._estWTag
            else :
                estWTag = self._parseArg( 'estWTag%s%d' % ( tagType, cat + 1 ), kwargs, ObjectType = 'ConstVar'
                                         , Title = 'Estimated wrong-tag probability ' + tagType, Value = catPars[3] )
                self._estWTags.append(estWTag)

            dilutions.append( self._parseArg( 'tagDilution%s%d' % ( tagType, cat + 1 ), kwargs, ObjectType = 'CalibratedDilution'
                                             , Title      = 'Tagging dilution %s %d' % ( tagType, cat + 1 )
                                             , EstWTag    = estWTag
                                             , AvgEstWTag = self._avgEstWTag
                                             , P0         = self._wTagP0
                                             , P1         = self._wTagP1
                                            )
                            )
            ADilWTags.append( self._parseArg( 'ADilWTag%s%d' % ( tagType, cat + 1 ), kwargs, ObjectType = 'CalibratedDilution'
                                             , Title      = 'Dilution/wrong-tag asymmetry %s %d' % ( tagType, cat + 1 )
                                             , EstWTag    = estWTag
                                             , AvgEstWTag = self._avgEstWTag
                                             , P0         = self._wTagP0
                                             , P1         = self._wTagP1
                                             , AP0        = self._wTagAP0
                                             , AP1        = self._wTagAP1
                                            )
                            )

        # check for remaining arguments and initialize
        self._check_extraneous_kw( kwargs )
        TaggingCategories.__init__( self, NumTagCats = len(self._tagCats), tagCat = tagCat
                                   , TagCatCoefs = tagCatCoefs, ATagEffs = ATagEffs, TagDilutions = dilutions, ADilWTags = ADilWTags
                                   , Conditionals = [ self._estWTag ] if self._estWTag else [ ]
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
            avgEtas = ( Categories0['calVals']['avgEta'], Categories1['calVals']['avgEta'] )
            P0s     = ( Categories0['calVals']['P0'],     Categories1['calVals']['P0']     )
            P0Errs  = ( Categories0['calVals']['P0Err'],  Categories1['calVals']['P0Err']  )
            P1s     = ( Categories0['calVals']['P1'],     Categories1['calVals']['P1']     )
            P1Errs  = ( Categories0['calVals']['P1Err'],  Categories1['calVals']['P1Err']  )
            AP0s    = ( Categories0['calVals']['AP0'],    Categories1['calVals']['AP0']    )
            AP1s    = ( Categories0['calVals']['AP1'],    Categories1['calVals']['AP1']    )

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
                                   , Conditionals = Categories0.ConditionalObservables() | Categories1.ConditionalObservables()
                                   , Constraints  = Categories0.ExternalConstraints() | Categories1.ExternalConstraints()
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

        namePF = self.getNamePrefix(kwargs)
        if triple :
           self._parseArg( 'TagEff', kwargs, Title = 'Tagging efficiency', Value = 0.25, MinMax = ( 0., 1. ) )
        self._parseArg( 'ATagEff', kwargs, Title = 'Tagging asymmetry ', Value = 0., MinMax = ( -1., 1. ) )

        from P2VV.RooFitWrappers import GenericPdf
        name = kwargs.pop( 'Name', namePF + 'Trivial_TagPdf' )
        if triple :
            self._pdf = GenericPdf( name, Formula = '(@0==0)*(1-@1)+(@0!=0)*@1*0.5*(1+@0*@2)'
                                   , Arguments = [  tagdecision
                                                  , getattr( self, '_TagEff'  )
                                                  , getattr( self, '_ATagEff' )
                                                 ]
                                  )
        else :
            self._pdf = GenericPdf( name, Formula = '0.5*(1+@0*@1)', Arguments = [ tagdecision, getattr( self, '_' + 'ATagEff' ) ] )

        self._pdf.setAttribute("CacheAndTrack") ;
        for ( k, v ) in kwargs.iteritems() : setattr( self, '_' + k, v )


class BinnedTaggingPdf( _util_parse_mixin ) :
    def __init__( self, Name, tagBinCoefs, tagCat0, tagCat1, iTag0, iTag1 ) :
        from P2VV.RooFitWrappers import BinnedPdf
        self._name        = Name
        self._tagBinCoefs = tagBinCoefs
        self._tagCats     = [ tagCat0, tagCat1 ]
        self._iTags       = [ iTag0, iTag1 ]
        self._pdf         = BinnedPdf( Name, Categories = [ tagCat0, tagCat1, iTag0, iTag1 ], Coefficients = tagBinCoefs )\
                            if tagCat0 and iTag0 else BinnedPdf( Name, Categories = [ tagCat1, iTag1 ], Coefficients = tagBinCoefs )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def pdf( self ) : return self._pdf

    def _init( self, Name, tagCat0, tagCat1, iTag0, iTag1, kwargs ) :
        # set tagging category and tag variables
        self._tagCats = [ tagCat0, tagCat1 ]
        self._iTags   = [ iTag0, iTag1 ]

        # get name prefix from arguments
        namePF = self.getNamePrefix(kwargs)

        # get name
        self._name = namePF + Name

        # get number of tagging categories
        self._numTagCats = [ self._tagCats[0].numTypes() if self._tagCats[0] else 0, self._tagCats[1].numTypes() ]

        # get tagging bin names
        self._untCatName = kwargs.pop( 'UntaggedCatName', 'Untagged' )
        self._tagCatName = kwargs.pop( 'TaggedCatName',   'Tagged'   )
        self._BName      = kwargs.pop( 'BName',           'B'        )
        self._BbarName   = kwargs.pop( 'BbarName',        'Bbar'     )

        # get data
        self._data = kwargs.pop( 'Data', None )

        # get work space
        from P2VV.RooFitWrappers import RooObject
        self._ws = RooObject().ws()

        if self._data :
            # get number of events in each tagging bin
            self._nBB          = [ [ 0. for cat1 in range( self._numTagCats[1] ) ] for cat0 in range( max( 1, self._numTagCats[0] ) ) ]
            self._nBBbar       = [ [ 0. for cat1 in range( self._numTagCats[1] ) ] for cat0 in range( max( 1, self._numTagCats[0] ) ) ]
            self._nBbarB       = [ [ 0. for cat1 in range( self._numTagCats[1] ) ] for cat0 in range( max( 1, self._numTagCats[0] ) ) ]
            self._nBbarBbar    = [ [ 0. for cat1 in range( self._numTagCats[1] ) ] for cat0 in range( max( 1, self._numTagCats[0] ) ) ]
            self._nBBTag       = [ 0., 0., 0., 0. ]
            self._nBBbarTag    = [ 0., 0., 0., 0. ]
            self._nBbarBTag    = [ 0., 0., 0., 0. ]
            self._nBbarBbarTag = [ 0., 0., 0., 0. ]

            from P2VV.RooFitWrappers import ArgSet
            tagVars = [ self._tagCats[0], self._tagCats[1], self._iTags[0], self._iTags[1] ] if self._tagCats[0]\
                      else [ self._tagCats[1], self._iTags[1] ]
            self._tagTable = self._data.table( ArgSet( self._name + '_tagTableSet', tagVars ) )
            for cat0Iter in range( max( 1, self._numTagCats[0] ) ) :
                tagCat0Name = self._untCatName if cat0Iter == 0 else self._tagCatName if self._numTagCats[0] <= 2\
                              else '%s%d' % ( self._tagCatName, cat0Iter )

                for cat1Iter in range( self._numTagCats[1] ) :
                    tagCat1Name = self._untCatName if cat1Iter == 0 else self._tagCatName if self._numTagCats[1] <= 2\
                                  else '%s%d' % ( self._tagCatName, cat1Iter )

                    if not self._tagCats[0] :
                        self._nBB[0][cat1Iter]        = self._tagTable.get( '{%s;%s}' % ( tagCat1Name, self._BName    ) )
                        self._nBbarBbar[0][cat1Iter]  = self._tagTable.get( '{%s;%s}' % ( tagCat1Name, self._BbarName ) )

                        self._nBBTag[ 0 if cat1Iter == 0 else 1 ]       += self._nBB[0][cat1Iter]
                        self._nBbarBbarTag[ 0 if cat1Iter == 0 else 1 ] += self._nBbarBbar[0][cat1Iter]

                    else :
                        self._nBB[cat0Iter][cat1Iter]       = self._tagTable.get( '{%s;%s;%s;%s}'\
                                                              % ( tagCat0Name, tagCat1Name, self._BName, self._BName       ) )
                        self._nBbarB[cat0Iter][cat1Iter]    = self._tagTable.get( '{%s;%s;%s;%s}'\
                                                              % ( tagCat0Name, tagCat1Name, self._BbarName, self._BName    ) )
                        self._nBBbar[cat0Iter][cat1Iter]    = self._tagTable.get( '{%s;%s;%s;%s}'\
                                                              % ( tagCat0Name, tagCat1Name, self._BName, self._BbarName    ) )
                        self._nBbarBbar[cat0Iter][cat1Iter] = self._tagTable.get( '{%s;%s;%s;%s}'\
                                                              % ( tagCat0Name, tagCat1Name, self._BbarName, self._BbarName ) )

                        if   cat0Iter == 0 and cat1Iter == 0 : tagIndex = 0
                        elif cat0Iter == 0 and cat1Iter  > 0 : tagIndex = 1
                        elif cat0Iter  > 0 and cat1Iter == 0 : tagIndex = 2
                        else                                 : tagIndex = 3
                        self._nBBTag[tagIndex]       += self._nBB[cat0Iter][cat1Iter]
                        self._nBBbarTag[tagIndex]    += self._nBBbar[cat0Iter][cat1Iter]
                        self._nBbarBTag[tagIndex]    += self._nBbarB[cat0Iter][cat1Iter]
                        self._nBbarBbarTag[tagIndex] += self._nBbarBbar[cat0Iter][cat1Iter]

        # check order of types in tag (default: B-Bbar, reversed: Bbar-B)
        self._tagRevOrder = False
        for iter, catType in enumerate( self._iTags[1] ) :
            if iter == 0 and catType.getVal() < 0 : self._tagRevOrder = True
        if self._tagCats[0] :
            tagRevOrder0 = False
            for iter, catType in enumerate( self._iTags[0] ) :
                if iter == 0 and catType.getVal() < 0 : tagRevOrder0 = True
            if tagRevOrder0 != self._tagRevOrder :
                raise RuntimeError('BinnedTaggingPdf._init(): order of tags 0 is different from order of tags 1')

        # get tagging category coefficients
        self._relativeCatCoefs = False if not kwargs.pop( 'RelativeCoefs', False ) else True
        self._tagCatCoefs = kwargs.pop( 'TagCatCoefs', None )
        if not self._tagCatCoefs :
            # create tagging category coefficients
            if self._data : from math import sqrt
            self._tagCatCoefs0 = [ ]
            self._tagCatCoefs1 = [ ]
            self._tagCatCoefs = [ [] for cat in range( max( 1, self._numTagCats[0] ) ) ]

            if self._tagCats[0] :
                # coefficients for flavour tag 0
                for coef0Iter in range( 1, self._numTagCats[0] ) :
                    if self._data :
                        coefVal = 0.
                        for coef1Iter in range( self._numTagCats[1] ) :
                            coefVal += self._nBB[coef0Iter][coef1Iter] + self._nBbarB[coef0Iter][coef1Iter]\
                                       + self._nBBbar[coef0Iter][coef1Iter] + self._nBbarBbar[coef0Iter][coef1Iter]
                        if coefVal < 0. : coefVal = 0.
                        coefErr = sqrt(coefVal) if coefVal > 100. else 10.
                        coefVal /= self._data.sumEntries()
                        coefErr /= self._data.sumEntries()

                    else :
                        coefVal = 1. / float( self._numTagCats[0] )
                        coefErr = 0.01

                    self._parseArg(  'tagCatCoef0_%d' % coef0Iter, kwargs, ContainerList = self._tagCatCoefs0
                                   , Title    = 'Tagging categories 0 coefficient %d' % coef0Iter
                                   , Value    = coefVal
                                   , Error    = coefErr
                                   , MinMax   = ( 0., 1. )
                                   , Constant = True if self._data else False
                                  )

                # coefficient for category 0: "one minus the sum of other categories"
                self._tagCatCoefs0 = [ self._parseArg( 'tagCatCoef0_0', kwargs, Title = 'Tagging categories 0 coefficient 0'
                                                      , Coefficients = self._tagCatCoefs0, ObjectType = 'ComplementCoef' )
                                     ] + self._tagCatCoefs0

            # coefficients for flavour tag 1
            for coef1Iter in range( 1, self._numTagCats[1] ) :
                if self._data :
                    coefVal = 0.
                    for coef0Iter in range( max( 1, self._numTagCats[0] ) ) :
                        coefVal += self._nBB[coef0Iter][coef1Iter] + self._nBbarB[coef0Iter][coef1Iter]\
                                   + self._nBBbar[coef0Iter][coef1Iter] + self._nBbarBbar[coef0Iter][coef1Iter]
                    if coefVal < 0. : coefVal = 0.
                    coefErr = sqrt(coefVal) if coefVal > 100. else 10.
                    coefVal /= self._data.sumEntries()
                    coefErr /= self._data.sumEntries()

                else :
                    coefVal = 1. / float( self._numTagCats[0] )
                    coefErr = 0.01

                self._parseArg(  'tagCatCoef1_%d' % coef1Iter, kwargs, ContainerList = self._tagCatCoefs1
                               , Title    = 'Tagging categories 0 coefficient %d' % coef1Iter
                               , Value    = coefVal
                               , Error    = coefErr
                               , MinMax   = ( 0., 1. )
                               , Constant = True if self._data else False
                              )

            # coefficient for category 0: "one minus the sum of other categories"
            self._tagCatCoefs1 = [ self._parseArg( 'tagCatCoef1_0', kwargs, Title = 'Tagging categories 1 coefficient 0'
                                                  , Coefficients = self._tagCatCoefs1, ObjectType = 'ComplementCoef' )
                                 ] + self._tagCatCoefs1

            # tagging category coefficients
            for coef0Iter in range( max( 1, self._numTagCats[0] ) ) :
                for coef1Iter in range( self._numTagCats[1] ) :
                    if not self._tagCats[0] :
                        # one flavour tag
                        self._tagCatCoefs[coef0Iter].append( self._tagCatCoefs1[coef1Iter] )
                    else :
                        # two flavour tags
                        self._tagCatCoefs[coef0Iter].append( self._parseArg( 'tagCatCoef%d-%d' % ( coef0Iter, coef1Iter ), kwargs
                                                            , Arguments = [ self._tagCatCoefs0[coef0Iter], self._tagCatCoefs1[coef1Iter] ]
                                                            , Title = 'Tagging category coefficient %d-%d' % ( coef0Iter, coef1Iter )
                                                            , ObjectType = 'Product' )
                                                           )

        # get tagging category coefficient names
        self._tagCatCoefNames = [ [ coef if type(coef) == str else coef.GetName() for coef in coefs ] for coefs in self._tagCatCoefs ]


class TagUntag_BinnedTaggingPdf( BinnedTaggingPdf ) :
    def __init__( self, Name, tagCatOS, tagCatSS, iTagOS, iTagSS, **kwargs ) :
        # initialize BinnedTaggingPdf
        if   tagCatOS and iTagOS and     tagCatSS and     iTagSS : self._init( Name, tagCatOS, tagCatSS, iTagOS, iTagSS, kwargs )
        elif tagCatOS and iTagOS and not tagCatSS and not iTagSS : self._init( Name, None,     tagCatOS, None,   iTagOS, kwargs )
        elif tagCatSS and iTagSS and not tagCatOS and not iTagOS : self._init( Name, None,     tagCatSS, None,   iTagSS, kwargs )
        else : raise RuntimeError('TagUntag_BinnedTaggingPdf(): no tagging categories/tags specified')

        if self._data : from math import sqrt

        namePF = self.getNamePrefix(kwargs)

        if self._relativeCatCoefs :
            ###############################################################################################################################
            ## tagging category coefficient asymmetries relative to specified coefficient ##
            ################################################################################
            if self._data :
                AUntVal = max( 0., self._nBBTag[0] + self._nBbarBbarTag[0] + self._nBBbarTag[0] + self._nBbarBTag[0] )
                AUntErr = sqrt(AUntVal) if AUntVal > 100. else 10.
                AUntVal = AUntVal / self._data.sumEntries() / self._ws[ self._tagCatCoefNames[0][0] ].getVal() - 1.
                AUntErr /=  self._data.sumEntries() * self._ws[ self._tagCatCoefNames[0][0] ].getVal()
            else :
                AUntVal = 0.
                AUntErr = 0.01

            self._parseArg(  'AUntag', kwargs
                           , Title = 'Untagged asymmetry in tagging category coefficients'
                           , Value = AUntVal, Error = AUntErr, MinMax = ( -1., 1. )
                          )

            if self._tagCats[0] :
                # two flavour tags
                if namePF + 'tagCatCoefSSTag' not in self._ws :
                    self._tagCatCoefSSTag = self._parseArg(  namePF + 'tagCatCoefSSTag', kwargs
                                                           , Formula = '+'.join( '@%d' % cat for cat in range( self._numTagCats[1] - 1 ) )
                                                           , Arguments = [ self._ws[ self._tagCatCoefNames[0][cat] ]\
                                                                          for cat in range( 1, self._numTagCats[1] ) ]
                                                           , ObjectType = 'FormulaVar'
                                                          )
                if namePF + 'tagCatCoefOSTag' not in self._ws :
                    self._tagCatCoefOSTag = self._parseArg(  namePF + 'tagCatCoefOSTag', kwargs
                                                           , Formula = '+'.join( '@%d' % cat for cat in range( self._numTagCats[0] - 1 ) )
                                                           , Arguments = [ self._ws[ self._tagCatCoefNames[cat][0] ]\
                                                                          for cat in range( 1, self._numTagCats[0] ) ]
                                                           , ObjectType = 'FormulaVar'
                                                          )

                if self._data :
                    ASSTagVal = max( 0., self._nBBTag[1] + self._nBbarBbarTag[1] + self._nBBbarTag[1] + self._nBbarBTag[1] )
                    ASSTagErr = sqrt(ASSTagVal) if ASSTagVal > 100. else 10.
                    ASSTagVal = ASSTagVal / self._data.sumEntries() / self._tagCatCoefSSTag.getVal() - 1.
                    ASSTagErr /=  self._data.sumEntries() * self._tagCatCoefSSTag.getVal()

                    AOSTagVal = max( 0., self._nBBTag[2] + self._nBbarBbarTag[2] + self._nBBbarTag[2] + self._nBbarBTag[2] )
                    AOSTagErr = sqrt(AOSTagVal) if AOSTagVal > 100. else 10.
                    AOSTagVal = AOSTagVal / self._data.sumEntries() / self._tagCatCoefOSTag.getVal() - 1.
                    AOSTagErr /=  self._data.sumEntries() * self._tagCatCoefOSTag.getVal()


                else :
                    ASSTagVal = 0.
                    ASSTagErr = 0.01
                    AOSTagVal = 0.
                    AOSTagErr = 0.01

                self._parseArg(  'ASSTag', kwargs
                               , Title = 'Same-side tagged asymmetry in tagging category coefficients'
                               , Value = ASSTagVal, Error = ASSTagErr, MinMax = ( -1., 1. )
                              )
                self._parseArg(  'AOSTag', kwargs
                               , Title = 'Opposite-side tagged asymmetry in tagging category coefficients'
                               , Value = AOSTagVal, Error = AOSTagErr, MinMax = ( -1., 1. )
                              )

        ###################################################################################################################################
        ## B-Bbar asymmetries ##
        ########################

        self._dilutions = kwargs.pop( 'Dilutions', None )
        if not self._tagCats[0] :
            # one flavour tag
            if self._dilutions : raise RuntimeError('TagUntag_BinnedTaggingPdf(): dilutions cannot be used in a PDF with one flavour tag')
            numB    = max( 0., self._nBBTag[1] + self._nBbarBTag[1] )
            numBbar = max( 0., self._nBBbarTag[1] + self._nBbarBbarTag[1] )

            ABBbarVal = ( numB - numBbar ) / ( numB + numBbar ) if numB > 0. or numBbar > 0. else 0.
            ABBbarErr = 1. / sqrt( numB + numBbar )             if numB + numBbar > 100. else 0.1

            self._parseArg(  'ABBbarTag', kwargs
                           , Title    = 'B-Bbar asymmetry tagged category'
                           , Value    = ABBbarVal
                           , Error    = ABBbarErr
                           , MinMax   = ( -1., 1. )
                           , Constant = False
                          )

        else :
            # two flavour tags
            if self._dilutions :
                from P2VV.RooFitWrappers import RooObject
                objVal = lambda obj : obj.getVal() if hasattr( obj, 'getVal' ) or isinstance( obj, RooObject ) else obj
                self._ATagsVals = [ [ objVal(dilOS) * objVal(dilSS) for dilSS in self._dilutions[1] ] for dilOS in self._dilutions[0] ]

            # calculate asymmetry values for tagged categories
            if self._data :
                # SS tagged
                numBSS    = max( 0., self._nBBTag[1] + self._nBbarBTag[1] )
                numBbarSS = max( 0., self._nBBbarTag[1] + self._nBbarBbarTag[1] )

                ABBbarSSVal = ( numBSS - numBbarSS ) / ( numBSS + numBbarSS ) if numBSS > 0. or numBbarSS > 0. else 0.
                ABBbarSSErr = 1. / sqrt( numBSS + numBbarSS )                 if numBSS + numBbarSS > 100. else 0.1

                # OS tagged
                numBOS    = max( 0., self._nBBTag[2] + self._nBBbarTag[2] )
                numBbarOS = max( 0., self._nBbarBTag[2] + self._nBbarBbarTag[2] )

                ABBbarOSVal = ( numBOS - numBbarOS ) / ( numBOS + numBbarOS ) if numBOS > 0. or numBbarOS > 0. else 0.
                ABBbarOSErr = 1. / sqrt( numBOS + numBbarOS )                 if numBOS + numBbarOS > 100. else 0.1

                # OS and SS tagged
                numSame     = max( 0., self._nBBTag[3] + self._nBbarBbarTag[3] )
                numOpp      = max( 0., self._nBBbarTag[3] + self._nBbarBTag[3] )
                numBSame    = max( 0., self._nBBTag[3] )
                numBbarSame = max( 0., self._nBbarBbarTag[3] )
                numBOpp     = max( 0., self._nBBbarTag[3] )
                numBbarOpp  = max( 0., self._nBbarBTag[3] )

                ATagsVal      = ( numSame - numOpp ) / ( numSame + numOpp )         if numSame > 0. or numOpp > 0. else 0.
                ATagsErr      = 1. / sqrt( numSame + numOpp )                       if numSame + numOpp > 100. else 0.1
                ABBbarSameVal = (numBSame - numBbarSame) / (numBSame + numBbarSame) if numBSame > 0. or numBbarSame > 0. else 0.
                ABBbarSameErr = 1. / sqrt( numBSame + numBbarSame )                 if numBSame + numBbarSame > 100. else 0.1
                ABBbarOppVal  = ( numBOpp - numBbarOpp ) / ( numBOpp + numBbarOpp ) if numBOpp > 0. or numBbarOpp > 0. else 0.
                ABBbarOppErr  = 1. / sqrt( numBOpp + numBbarOpp )                   if numBOpp + numBbarOpp > 100. else 0.1

            else :
                # no data provided
                ABBbarSSVal   = 0.
                ABBbarSSErr   = 0.1
                ABBbarOSVal   = 0.
                ABBbarOSErr   = 0.1
                ATagsVal      = 0.
                ATagsErr      = 0.1
                ABBbarSameVal = 0.
                ABBbarSameErr = 0.1
                ABBbarOppVal  = 0.
                ABBbarOppErr  = 0.1

            # create asymmetry variables
            self._parseArg(  'ABBbarSSTag', kwargs
                           , Title    = 'B-Bbar asymmetry same-side tagged'
                           , Value    = ABBbarSSVal
                           , Error    = ABBbarSSErr
                           , MinMax   = ( -1., 1. )
                           , Constant = False
                          )
            self._parseArg(  'ABBbarOSTag', kwargs
                           , Title    = 'B-Bbar asymmetry opposite-side tagged'
                           , Value    = ABBbarOSVal
                           , Error    = ABBbarOSErr
                           , MinMax   = ( -1., 1. )
                           , Constant = False
                          )
            if self._dilutions :
                self._ATags = [ [ ] for i in range( len(self._ATagsVals) ) ]
                for OSIter, ATagsVals in enumerate(self._ATagsVals) :
                    for SSIter, ATagsVal in enumerate(ATagsVals) :
                        if OSIter == 0 or SSIter == 0 :
                            self._ATags[OSIter].append(None)
                        else :
                            self._parseArg(  'ATags%d-%d' % ( OSIter, SSIter ), kwargs, ContainerList = self._ATags[OSIter]
                                           , Title    = 'Tags asymmetry category %d-%d' % ( OSIter, SSIter )
                                           , Value    = ATagsVal
                                           , Error    = ATagsErr
                                           , MinMax   = ( -1., 1. )
                                           , Constant = False
                                          )
            else :
                self._parseArg(  'ATags', kwargs
                               , Title    = 'Tags asymmetry'
                               , Value    = ATagsVal
                               , Error    = ATagsErr
                               , MinMax   = ( -1., 1. )
                               , Constant = False
                              )
            self._parseArg(  'ABBbarSameTag', kwargs
                           , Title    = 'B-Bbar asymmetry same tags'
                           , Value    = ABBbarSameVal
                           , Error    = ABBbarSameErr
                           , MinMax   = ( -1., 1. )
                           , Constant = False
                          )
            self._parseArg(  'ABBbarOppTag', kwargs
                           , Title    = 'B-Bbar asymmetry opposite tags'
                           , Value    = ABBbarOppVal
                           , Error    = ABBbarOppErr
                           , MinMax   = ( -1., 1. )
                           , Constant = False
                          )

        ###################################################################################################################################
        ## tagging category asymmetry factors ##
        ########################################
        if not self._tagCats[0] :
            # one flavour tag
            if self._relativeCatCoefs :
                self._taggedBCoef    = self._parseArg( 'taggedBCoef', kwargs, Formula = '0.5*(1.-@2/(1.-@2)*@0)*(1.+@1)'
                                                      , Arguments = [self._AUntag, self._ABBbarTag, self._ws[self._tagCatCoefNames[0][0]]]
                                                      , ObjectType = 'FormulaVar'
                                                     )
                self._taggedBbarCoef = self._parseArg( 'taggedBbarCoef', kwargs, Formula = '0.5*(1.-@2/(1.-@2)*@0)*(1.-@1)'
                                                      , Arguments = [self._AUntag, self._ABBbarTag, self._ws[self._tagCatCoefNames[0][0]]]
                                                      , ObjectType = 'FormulaVar'
                                                     )
            else :
                self._taggedBCoef    = self._parseArg( 'taggedBCoef', kwargs, Formula = '0.5*(1.+@0)', Arguments = [ self._ABBbarTag ]
                                                      , ObjectType = 'FormulaVar' )
                self._taggedBbarCoef = self._parseArg( 'taggedBbarCoef', kwargs, Formula = '0.5*(1.-@0)', Arguments = [ self._ABBbarTag ]
                                                      , ObjectType = 'FormulaVar' )

        else :
            # two flavour tags
            if self._relativeCatCoefs :
                self._SSTaggedBCoef    = self._parseArg( 'SSTaggedBCoef', kwargs, Formula = '0.25*(1.+@0)*(1.+@1)'
                                                        , Arguments = [ self._ASSTag, self._ABBbarSSTag ], ObjectType = 'FormulaVar' )
                self._SSTaggedBbarCoef = self._parseArg( 'SSTaggedBbarCoef', kwargs, Formula = '0.25*(1.+@0)*(1.-@1)'
                                                         , Arguments = [ self._ASSTag, self._ABBbarSSTag ], ObjectType = 'FormulaVar' )
                self._OSTaggedBCoef    = self._parseArg( 'OSTaggedBCoef', kwargs, Formula = '0.25*(1.+@0)*(1.+@1)'
                                                         , Arguments = [ self._AOSTag, self._ABBbarOSTag ], ObjectType = 'FormulaVar' )
                self._OSTaggedBbarCoef = self._parseArg( 'OSTaggedBbarCoef', kwargs, Formula = '0.25*(1.+@0)*(1.-@1)'
                                                        , Arguments = [ self._AOSTag, self._ABBbarOSTag ], ObjectType = 'FormulaVar' )
                if not self._dilutions :
                    self._taggedBBCoef       = self._parseArg(  'taggedBBCoef', kwargs
                                                              , Formula = '0.25*(1.-(@5*@0+@6*@1+@7*@2)/(1.-@5-@6-@7))*(1.+@3)*(1.+@4)'
                                                              , Arguments = [ self._AUntag, self._AOSTag, self._ASSTag, self._ATags
                                                                             , self._ABBbarSameTag, self._ws[ self._tagCatCoefNames[0][0] ]
                                                                             , self._tagCatCoefOSTag, self._tagCatCoefSSTag ]
                                                              , ObjectType = 'FormulaVar'
                                                             )
                    self._taggedBbarBbarCoef = self._parseArg(  'taggedBbarBbarCoef', kwargs
                                                              , Formula = '0.25*(1.-(@5*@0+@6*@1+@7*@2)/(1.-@5-@6-@7))*(1.+@3)*(1.-@4)'
                                                              , Arguments = [ self._AUntag, self._AOSTag, self._ASSTag, self._ATags
                                                                             , self._ABBbarSameTag, self._ws[ self._tagCatCoefNames[0][0] ]
                                                                             , self._tagCatCoefOSTag, self._tagCatCoefSSTag ]
                                                              , ObjectType = 'FormulaVar'
                                                             )
                    self._taggedBBbarCoef    = self._parseArg(  'taggedBBbarCoef', kwargs
                                                              , Formula = '0.25*(1.-(@5*@0+@6*@1+@7*@2)/(1.-@5-@6-@7))*(1.-@3)*(1.+@4)'
                                                              , Arguments = [ self._AUntag, self._AOSTag, self._ASSTag, self._ATags
                                                                             , self._ABBbarOppTag, self._ws[ self._tagCatCoefNames[0][0] ]
                                                                             , self._tagCatCoefOSTag, self._tagCatCoefSSTag ]
                                                              , ObjectType = 'FormulaVar'
                                                             )
                    self._taggedBbarBCoef    = self._parseArg(  'taggedBbarBCoef', kwargs
                                                              , Formula = '0.25*(1.-(@5*@0+@6*@1+@7*@2)/(1.-@5-@6-@7))*(1.-@3)*(1.-@4)'
                                                              , Arguments = [ self._AUntag, self._AOSTag, self._ASSTag, self._ATags
                                                                             , self._ABBbarOppTag, self._ws[ self._tagCatCoefNames[0][0] ]
                                                                             , self._tagCatCoefOSTag, self._tagCatCoefSSTag ]
                                                              , ObjectType = 'FormulaVar'
                                                             )
            else :
                self._SSTaggedBCoef    = self._parseArg( 'SSTaggedBCoef', kwargs, Formula = '0.25*(1.+@0)'
                                                        , Arguments = [ self._ABBbarSSTag ], ObjectType = 'FormulaVar' )
                self._SSTaggedBbarCoef = self._parseArg( 'SSTaggedBbarCoef', kwargs, Formula = '0.25*(1.-@0)'
                                                        , Arguments = [ self._ABBbarSSTag ], ObjectType = 'FormulaVar' )
                self._OSTaggedBCoef    = self._parseArg( 'OSTaggedBCoef', kwargs, Formula = '0.25*(1.+@0)'
                                                        , Arguments = [ self._ABBbarOSTag ], ObjectType = 'FormulaVar' )
                self._OSTaggedBbarCoef = self._parseArg( 'OSTaggedBbarCoef', kwargs, Formula = '0.25*(1.-@0)'
                                                        , Arguments = [ self._ABBbarOSTag ], ObjectType = 'FormulaVar' )
                if not self._dilutions :
                    self._taggedBBCoef       = self._parseArg( 'taggedBBCoef', kwargs, Formula = '0.25*(1.+@0)*(1.+@1)'
                                                              , Arguments = [self._ATags, self._ABBbarSameTag], ObjectType = 'FormulaVar' )
                    self._taggedBbarBbarCoef = self._parseArg( 'taggedBbarBbarCoef', kwargs, Formula = '0.25*(1.+@0)*(1.-@1)'
                                                              , Arguments = [self._ATags, self._ABBbarSameTag], ObjectType = 'FormulaVar' )
                    self._taggedBBbarCoef    = self._parseArg( 'taggedBBbarCoef', kwargs, Formula = '0.25*(1.-@0)*(1.+@1)'
                                                              , Arguments = [self._ATags, self._ABBbarOppTag], ObjectType = 'FormulaVar' )
                    self._taggedBbarBCoef    = self._parseArg( 'taggedBbarBCoef', kwargs, Formula = '0.25*(1.-@0)*(1.-@1)'
                                                              , Arguments = [self._ATags, self._ABBbarOppTag], ObjectType = 'FormulaVar' )

        ###################################################################################################################################
        ## tagging bin coefficients ##
        ##############################
        tagBinCoefs = [ ]
        if not self._tagCats[0] :
            # one flavour tag: loop over #cats * 2 bins
            for binIter in range( self._numTagCats[1] * 2 ) :
                cat = binIter % self._numTagCats[1]
                Bbar = bool( ( binIter / self._numTagCats[1] ) % 2 )
                if self._tagRevOrder : Bbar = not Bbar

                # create product of tagging category coefficient and bin asymmetry factor
                if cat > 0 :
                    # tagged coefficient
                    tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, ObjectType = 'Product'
                                                       , Arguments = [ self._ws[ self._tagCatCoefNames[0][cat] ]
                                                                      , self._taggedBCoef if not Bbar else self._taggedBbarCoef ]
                                                      )
                                      )
                elif self._relativeCatCoefs :
                    # untagged coefficient (relative)
                    tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, Formula = '0.5*@0*(1.+@1)'
                                                       , Arguments = [ self._ws[ self._tagCatCoefNames[0][0] ], self._AUntag ]
                                                       , ObjectType = 'FormulaVar'
                                                      )
                                      )
                else :
                    # untagged coefficient (absolute)
                    tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, Formula = '0.5*@0'
                                                       , Arguments = [ self._ws[ self._tagCatCoefNames[0][0] ] ]
                                                       , ObjectType = 'FormulaVar'
                                                      )
                                      )

        else :
            # two flavour tags: loop over #cats0 * #cats1 * 2 * 2 bins
            for binIter in range( self._numTagCats[0] * self._numTagCats[1] * 4 ) :
                cat   = binIter % ( self._numTagCats[0] * self._numTagCats[1] )
                cat0  = cat % self._numTagCats[0]
                cat1  = cat / self._numTagCats[0]
                Bbar0 = bool( ( binIter / self._numTagCats[0] / self._numTagCats[1] )     % 2 )
                Bbar1 = bool( ( binIter / self._numTagCats[0] / self._numTagCats[1] / 2 ) % 2 )
                if self._tagRevOrder :
                    Bbar0 = not Bbar0
                    Bbar1 = not Bbar1

                # create product of tagging category coefficient and bin asymmetry factor
                if cat0 > 0 and cat1 > 0 :
                    # tagged coefficient
                    if self._dilutions and self._relativeCatCoefs :
                        tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs
                                                           , Formula = '0.25*@0*(1.-(@6*@1+@7*@2+@8*@3)/(1.-@6-@7-@8))*(1.%s@4)*(1.%s@5)'\
                                                                       % ( '-' if Bbar0 != Bbar1 else '+', '-' if Bbar0 else '+' )
                                                           , Arguments = [ self._ws[ self._tagCatCoefNames[cat0][cat1] ], self._AUntag
                                                                          , self._AOSTag, self._ASSTag, self._ATags[cat0][cat1]
                                                                          , self._ABBbarOppTag if Bbar0 != Bbar1 else self._ABBbarSameTag
                                                                          , self._ws[ self._tagCatCoefNames[0][0] ], self._tagCatCoefOSTag
                                                                          , self._tagCatCoefSSTag ]
                                                           , ObjectType = 'FormulaVar'
                                                          )
                                          )
                    elif self._dilutions :
                        tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs
                                                           , Formula = '0.25*@0*(1.%s@1)*(1.%s@2)'\
                                                                       % ( '-' if Bbar0 != Bbar1 else '+', '-' if Bbar0 else '+' )
                                                           , Arguments = [ self._ws[ self._tagCatCoefNames[cat0][cat1] ]
                                                                          , self._ATags[cat0][cat1], self._ABBbarOppTag if Bbar0 != Bbar1\
                                                                                                     else self._ABBbarSameTag ]
                                                           , ObjectType = 'FormulaVar'
                                                          )
                                          )
                    else :
                        tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, ObjectType = 'Product'
                                                           , Arguments = [  self._ws[ self._tagCatCoefNames[cat0][cat1] ]
                                                                          , self._taggedBBCoef if not Bbar0 and not Bbar1\
                                                                            else self._taggedBbarBbarCoef if Bbar0 and Bbar1\
                                                                            else self._taggedBBbarCoef if not Bbar0\
                                                                            else self._taggedBbarBCoef
                                                                         ]
                                                          )
                                          )
                elif cat0 > 0 or cat1 > 0 :
                    # OS tagged/SS tagged coefficient
                    tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, ObjectType = 'Product'
                                                       , Arguments = [  self._ws[ self._tagCatCoefNames[cat0][cat1] ]
                                                                      , self._SSTaggedBCoef if cat1 > 0 and not Bbar1\
                                                                        else self._SSTaggedBbarCoef if cat1 > 0\
                                                                        else self._OSTaggedBCoef if not Bbar0\
                                                                        else self._OSTaggedBbarCoef
                                                                     ]
                                                      )
                                      )
                elif self._relativeCatCoefs :
                    # untagged coefficient (relative)
                    tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, Formula = '0.25*@0*(1.+@1)'
                                                       , Arguments = [ self._ws[ self._tagCatCoefNames[0][0] ], self._AUntag ]
                                                       , ObjectType = 'FormulaVar'
                                                      )
                                      )
                else :
                    # untagged coefficient (absolute)
                    tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, Formula = '0.25*@0'
                                                       , Arguments = [ self._ws[ self._tagCatCoefNames[0][0] ] ]
                                                       , ObjectType = 'FormulaVar'
                                                      )
                                      )

        # initialize TaggingPdf
        self._check_extraneous_kw( kwargs )
        BinnedTaggingPdf.__init__( self, self._name, tagBinCoefs, self._tagCats[0], self._tagCats[1], self._iTags[0], self._iTags[1] )

        # force integral = 1
        self._pdf.setForceUnitIntegral(True)


class TagCats_BinnedTaggingPdf( BinnedTaggingPdf ) :
    def __init__( self, Name, tagCatOS, tagCatSS, iTagOS, iTagSS, **kwargs ) :
        # initialize BinnedTaggingPdf
        if   tagCatOS and iTagOS and     tagCatSS and     iTagSS : self._init( Name, tagCatOS, tagCatSS, iTagOS, iTagSS, kwargs )
        elif tagCatOS and iTagOS and not tagCatSS and not iTagSS : self._init( Name, None,     tagCatOS, None,   iTagOS, kwargs )
        elif tagCatSS and iTagSS and not tagCatOS and not iTagOS : self._init( Name, None,     tagCatSS, None,   iTagSS, kwargs )
        else : raise RuntimeError('TagUntag_BinnedTaggingPdf(): no tagging categories/tags specified')

        if self._data : from math import sqrt

        namePF = self.getNamePrefix(kwargs)

        if self._relativeCatCoefs :
            ###############################################################################################################################
            ## tagging category coefficient asymmetries relative to specified coefficient ##
            ################################################################################
            self._ATagCats = [ [ ] for i in range( max( 1, self._numTagCats[0] ) ) ]

            # loop over tagging categories to calculate asymmetry values
            for cat0Iter in range( max( 1, self._numTagCats[0] ) ) :
                for cat1Iter in range( self._numTagCats[1] ) :
                    # don't calculate for category (0, 0): this will be the sum of the other categories
                    if cat0Iter == 0 and cat1Iter == 0 : continue

                    # get asymmetry value and error
                    if self._data :
                        ACatVal = max( 0., self._nBB[cat0Iter][cat1Iter] + self._nBBbar[cat0Iter][cat1Iter]
                                  + self._nBbarB[cat0Iter][cat1Iter] + self._nBbarBbar[cat0Iter][cat1Iter] )
                        ACatErr = sqrt(ACatVal) if ACatVal > 100. else 10.
                        ACatVal = ACatVal / self._data.sumEntries() / self._ws[self._tagCatCoefNames[cat0Iter][cat1Iter]].getVal() - 1.\
                                  if self._ws[self._tagCatCoefNames[cat0Iter][cat1Iter]].getVal() > 0. else 0.
                        ACatErr /= self._data.sumEntries() * self._ws[self._tagCatCoefNames[cat0Iter][cat1Iter]].getVal()\
                                  if self._ws[self._tagCatCoefNames[cat0Iter][cat1Iter]].getVal() > 0. else 0.1

                    else :
                        ACatVal = 0.
                        ACatErr = 0.1

                    # create asymmetry variable
                    self._parseArg( 'ATagCat%s-%s' % ( cat0Iter, cat1Iter ), kwargs, ContainerList = self._ATagCats[cat0Iter]
                                   , Title    = 'Tagging category coefficient asymmetry %d-%d' % ( cat0Iter, cat1Iter )
                                   , Value    = ACatVal
                                   , Error    = ACatErr
                                   , MinMax   = ( -RooInf, RooInf )
                                   , Constant = False
                                  )

            # create variable for category (0, 0): "minus the weighted sum of other categories"
            numCats = max( 1, self._numTagCats[0] ) * self._numTagCats[1]
            ATagArgs = [ self._ws[coefName] for coefNames in self._tagCatCoefNames for coefName in coefNames ]
            for cat0 in range( max( 1, self._numTagCats[0] ) ) : ATagArgs += self._ATagCats[cat0][ : ]
            self._ATagCats[0] = [ self._parseArg( 'ATagCat0', kwargs
                                                 , Formula = '-1./@0*(%s)' % '+'.join( '@%d*@%d' % ( cat, cat + numCats - 1 )\
                                                                                       for cat in range( 1, numCats ) )
                                                 , Arguments = ATagArgs
                                                 , ObjectType = 'FormulaVar'
                                                 , Title = 'Tagging category coefficient asymmetry 0'
                                                )
                                ] + self._ATagCats[0]

        ###################################################################################################################################
        ## B-Bbar asymmetries ##
        ########################
        if self._tagCats[0] :
            self._ATags       = [ [ ] for i in range( max( 1, self._numTagCats[0] ) ) ]
            self._ABBbarsSame = [ [ ] for i in range( max( 1, self._numTagCats[0] ) ) ]
            self._ABBbarsOpp  = [ [ ] for i in range( max( 1, self._numTagCats[0] ) ) ]
        else :
            self._ABBbars = [ ]

        # loop over tagging categories to calculate asymmetry values
        for cat0Iter in range( max( 1, self._numTagCats[0] ) ) :
            for cat1Iter in range( self._numTagCats[1] ) :
                if not ( cat0Iter == 0 and cat1Iter == 0 ) and self._data :
                    # calculate asymmetry
                    if cat0Iter == 0 or cat1Iter == 0 :
                        # one of the two flavours is "untagged"
                        ATagsVal = 0.
                        ATagsErr = 0.01
                        if cat0Iter == 0 :
                            numB    = max( 0., self._nBB[cat0Iter][cat1Iter] + self._nBbarB[cat0Iter][cat1Iter] )
                            numBbar = max( 0., self._nBBbar[cat0Iter][cat1Iter] + self._nBbarBbar[cat0Iter][cat1Iter] )
                        else :
                            numB    = max( 0., self._nBB[cat0Iter][cat1Iter] + self._nBBbar[cat0Iter][cat1Iter] )
                            numBbar = max( 0., self._nBbarB[cat0Iter][cat1Iter] + self._nBbarBbar[cat0Iter][cat1Iter] )

                        ABBbarSameVal = ( numB - numBbar ) / ( numB + numBbar ) if numB > 0. or numBbar > 0. else 0.
                        ABBbarSameErr = 1. / sqrt( numB + numBbar )             if numB + numBbar > 100. else 0.1
                        ABBbarOppVal  = -ABBbarSameVal if cat0Iter == 0 else ABBbarSameVal
                        ABBbarOppErr  = ABBbarSameErr

                    else :
                        # both flavours are tagged
                        numSame     = max( 0., self._nBB[cat0Iter][cat1Iter] + self._nBbarBbar[cat0Iter][cat1Iter] )
                        numOpp      = max( 0., self._nBBbar[cat0Iter][cat1Iter] + self._nBbarB[cat0Iter][cat1Iter] )
                        numBSame    = max( 0., self._nBB[cat0Iter][cat1Iter] )
                        numBbarSame = max( 0., self._nBbarBbar[cat0Iter][cat1Iter] )
                        numBOpp     = max( 0., self._nBBbar[cat0Iter][cat1Iter] )
                        numBbarOpp  = max( 0., self._nBbarB[cat0Iter][cat1Iter] )

                        ATagsVal      = ( numSame - numOpp ) / ( numSame + numOpp )         if numSame > 0. or numOpp > 0. else 0.
                        ATagsErr      = 1. / sqrt( numSame + numOpp )                       if numSame + numOpp > 100. else 0.1
                        ABBbarSameVal = (numBSame - numBbarSame) / (numBSame + numBbarSame) if numBSame > 0. or numBbarSame > 0. else 0.
                        ABBbarSameErr = 1. / sqrt( numBSame + numBbarSame )                 if numBSame + numBbarSame > 100. else 0.1
                        ABBbarOppVal  = ( numBOpp - numBbarOpp ) / ( numBOpp + numBbarOpp ) if numBOpp > 0. or numBbarOpp > 0. else 0.
                        ABBbarOppErr  = 1. / sqrt( numBOpp + numBbarOpp )                   if numBOpp + numBbarOpp > 100. else 0.1

                else :
                    # no data provided or the untagged-untagged category (0, 0)
                    ATagsVal      = 0.
                    ATagsErr      = 0.1
                    ABBbarSameVal = 0.
                    ABBbarSameErr = 0.1
                    ABBbarOppVal  = 0.
                    ABBbarOppErr  = 0.1

                # create asymmetry variables
                if self._tagCats[0] :
                    # two flavour tags
                    if cat0Iter == 0 and cat1Iter == 0 :
                        # both flavours are untagged
                        self._ATags[cat0Iter].append(None)
                        self._ABBbarsSame[cat0Iter].append(None)
                        self._ABBbarsOpp[cat0Iter].append(None)

                    elif cat0Iter == 0 :
                        # only flavour 1 is tagged
                        self._ATags[cat0Iter].append(None)
                        self._parseArg(  'ABBbarSS%d' % cat1Iter, kwargs, ContainerList = self._ABBbarsSame[cat0Iter]
                                       , Title    = 'B-Bbar asymmetry same-side tagging category %d' % cat1Iter
                                       , Value    = ABBbarSameVal
                                       , Error    = ABBbarSameErr
                                       , MinMax   = ( -1., 1. )
                                       , Constant = False
                                      )
                        self._ABBbarsOpp[cat0Iter].append(None)

                    elif cat1Iter == 0 :
                        # only flavour 0 is tagged
                        self._ATags[cat0Iter].append(None)
                        self._parseArg(  'ABBbarOS%d' % cat0Iter, kwargs, ContainerList = self._ABBbarsSame[cat0Iter]
                                       , Title    = 'B-Bbar asymmetry opposite-side tagging category %d' % cat0Iter
                                       , Value    = ABBbarSameVal
                                       , Error    = ABBbarSameErr
                                       , MinMax   = ( -1., 1. )
                                       , Constant = False
                                      )
                        self._ABBbarsOpp[cat0Iter].append(None)

                    else :
                        # both flavours are tagged
                        self._parseArg(  'ATags%d-%d' % ( cat0Iter, cat1Iter ), kwargs, ContainerList = self._ATags[cat0Iter]
                                       , Title    = 'Tags asymmetry tagging category %d-%d' % ( cat0Iter, cat1Iter )
                                       , Value    = ATagsVal
                                       , Error    = ATagsErr
                                       , MinMax   = ( -1., 1. )
                                       , Constant = False
                                      )
                        self._parseArg(  'ABBbarSameTag%d-%d' % (cat0Iter, cat1Iter), kwargs, ContainerList = self._ABBbarsSame[cat0Iter]
                                       , Title    = 'B-Bbar asymmetry same tags tagging category %d-%d' % ( cat0Iter, cat1Iter )
                                       , Value    = ABBbarSameVal
                                       , Error    = ABBbarSameErr
                                       , MinMax   = ( -1., 1. )
                                       , Constant = False
                                      )
                        self._parseArg(  'ABBbarOppTag%d-%d' % ( cat0Iter, cat1Iter ), kwargs, ContainerList = self._ABBbarsOpp[cat0Iter]
                                       , Title    = 'B-Bbar asymmetry opposite tags tagging category %d-%d' % ( cat0Iter, cat1Iter )
                                       , Value    = ABBbarOppVal
                                       , Error    = ABBbarOppErr
                                       , MinMax   = ( -1., 1. )
                                       , Constant = False
                                      )
                else :
                    # one flavour tag
                    self._parseArg(  'ABBbar%d' % cat1Iter, kwargs, ContainerList = self._ABBbars
                                   , Title    = 'B-Bbar asymmetry tagging category %d' % cat1Iter
                                   , Value    = ABBbarSameVal
                                   , Error    = ABBbarSameErr
                                   , MinMax   = ( -1., 1. )
                                   , Constant = False
                                  )

        ###################################################################################################################################
        ## tagging bin coefficients ##
        ##############################
        tagBinCoefs = [ ]
        if not self._tagCats[0] :
            # one flavour tag: loop over #cats * 2 bins
            for binIter in range( self._numTagCats[1] * 2 ) :
                cat = binIter % self._numTagCats[1]
                Bbar = bool( ( binIter / self._numTagCats[1] ) % 2 )
                if self._tagRevOrder : Bbar = not Bbar
                if self._relativeCatCoefs :
                    # create bin coefficient relative to specified category coefficient
                    tagBinForm = '0.5*@0*(1+@1)*(1%s@2)' % '-' if Bbar  else '+'
                    tagBinArgs = [ self._ws[ self._tagCatCoefNames[0][cat] ], self._ATagCats[0][cat], self._ABBbars[cat] ]

                else :
                    # create bin coefficient with specified category coefficient
                    tagBinForm = '0.5*@0*(1%s@1)' % '-' if Bbar else '+'
                    tagBinArgs = [ self._ws[ self._tagCatCoefNames[0][cat] ], self._ABBbars[cat] ]

                tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, Formula = tagBinForm, Arguments = tagBinArgs
                                                   , ObjectType = 'FormulaVar' ) )

        else :
            # two flavour tags: loop over #cats0 * #cats1 * 2 * 2 bins
            for binIter in range( self._numTagCats[0] * self._numTagCats[1] * 4 ) :
                cat   = binIter % ( self._numTagCats[0] * self._numTagCats[1] )
                cat0  = cat % self._numTagCats[0]
                cat1  = cat / self._numTagCats[0]
                Bbar0 = bool( ( binIter / self._numTagCats[0] / self._numTagCats[1] )     % 2 )
                Bbar1 = bool( ( binIter / self._numTagCats[0] / self._numTagCats[1] / 2 ) % 2 )
                if self._tagRevOrder :
                    Bbar0 = not Bbar0
                    Bbar1 = not Bbar1

                if self._relativeCatCoefs :
                    # create bin coefficient relative to specified category coefficient
                    if cat0 == 0 and cat1 == 0 :
                        # both flavours are untagged
                        tagBinForm = '0.25*@0*(1+@1)'
                        tagBinArgs = [ self._ws[ self._tagCatCoefNames[cat0][cat1] ], self._ATagCats[cat0][cat1] ]

                    elif cat0 == 0 or cat1 == 0 :
                        # only one flavour is tagged
                        tagBinForm = '0.25*@0*(1+@1)*(1%s@2)'\
                                     % ( '+' if ( ( cat0 == 0 and not Bbar1 ) or ( cat1 == 0 and not Bbar0 ) ) else '-' )
                        tagBinArgs = [  self._ws[ self._tagCatCoefNames[cat0][cat1] ]
                                      , self._ATagCats[cat0][cat1]
                                      , self._ABBbarsSame[cat0][cat1]
                                     ]

                    else :
                        # both flavours are tagged
                        tagBinForm = '0.25*@0*(1+@1)*(1%s@2)*(1%s@3)' % ( '+' if Bbar0 == Bbar1 else '-', '+' if not Bbar0 else '-' )
                        tagBinArgs = [  self._ws[ self._tagCatCoefNames[cat0][cat1] ]
                                      , self._ATagCats[cat0][cat1]
                                      , self._ATags[cat0][cat1]
                                      , self._ABBbarsSame[cat0][cat1] if Bbar0 == Bbar1 else self._ABBbarsOpp[cat0][cat1]
                                     ]

                else :
                    if cat0 == 0 and cat1 == 0 :
                        # both flavours are untagged
                        tagBinForm = '0.25*@0'
                        tagBinArgs = [ self._ws[ self._tagCatCoefNames[cat0][cat1] ] ]

                    elif cat0 == 0 or cat1 == 0 :
                        # only one flavour is tagged
                        tagBinForm = '0.25*@0*(1%s@1)' % ( '+' if ( ( cat0 == 0 and not Bbar1 ) or ( cat1 == 0 and not Bbar0 ) ) else '-' )
                        tagBinArgs = [ self._ws[ self._tagCatCoefNames[cat0][cat1] ], self._ABBbarsSame[cat0][cat1] ]

                    else :
                        # both flavours are tagged
                        tagBinForm = '0.25*@0*(1%s@1)*(1%s@2)' % ( '+' if Bbar0 == Bbar1 else '-', '+' if not Bbar0 else '-' )
                        tagBinArgs = [  self._ws[ self._tagCatCoefNames[cat0][cat1] ]
                                      , self._ATags[cat0][cat1]
                                      , self._ABBbarsSame[cat0][cat1] if Bbar0 == Bbar1 else self._ABBbarsOpp[cat0][cat1]
                                     ]

                tagBinCoefs.append( self._parseArg( 'tagBinCoef%d' % binIter, kwargs, Formula = tagBinForm, Arguments = tagBinArgs
                                   , ObjectType = 'FormulaVar' ) )

        # initialize TaggingPdf
        self._check_extraneous_kw( kwargs )
        BinnedTaggingPdf.__init__( self, self._name, tagBinCoefs, self._tagCats[0], self._tagCats[1], self._iTags[0], self._iTags[1] )

        # force integral = 1
        self._pdf.setForceUnitIntegral(True)
