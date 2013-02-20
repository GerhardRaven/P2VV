###########################################################################################################################################
## P2VVParameterizations.TimePDFs: Parameterizations of PDFs that depend on decay time                                                   ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

class BDecayBasisCoefficients :
    def __init__(self, **kwargs ) :
        for i in ['sin','cos','sinh','cosh' ] : setattr(self,i,kwargs.pop(i))
        if kwargs : raise KeyError('unknown keyword arguments: %s' % kwargs )
    def __getitem__(self,kw) :
        return getattr(self,kw)

class JpsiphiBTagDecayBasisCoefficients( BDecayBasisCoefficients ) :
    def __init__( self, AngFuncs, Amplitudes, CPParams, Order ) :
        def combine( tCoefType, angFuncs, amplitudes, CPParams, iIndex, jIndex ) :
            from P2VV.RooFitWrappers import ConstVar, FormulaVar, Product, Addition
            one   = ConstVar( Name = 'one',   Value =  1 )
            minus = ConstVar( Name = 'minus', Value = -1 )

            # define functions which return Re(Ai* Aj), Im( Ai* Aj)
            # TODO: replace by Addition & Product
            Re = lambda Ai, Aj : FormulaVar( 'Re_amps_%s_%s' % ( Ai, Aj ), '@0*@2+@1*@3', [ Ai.Re, Ai.Im, Aj.Re, Aj.Im ] )
            Im = lambda Ai, Aj : FormulaVar( 'Im_amps_%s_%s' % ( Ai, Aj ), '@0*@3-@1*@2', [ Ai.Re, Ai.Im, Aj.Re, Aj.Im ] )

            # define functions which return the coefficients of the time dependence
            def CPVDec( termInd, iInd, jInd, amps, CPPar ) :
                # get parameters for CP violation in decay
                assert termInd in ( 'plus', 'min', 'Re', 'Im' )
                if CPPar.CPVInDecay() : return CPPar.R( termInd, iInd, jInd )

                # only CP violation in mixing (or lambdas for all polarizations identical):
                # R^+  =     ( 1 + eta_i * eta_j ) / 2
                # R^-  =     ( 1 - eta_i * eta_j ) / 2
                # R^Re =     ( eta_i + eta_j )     / 2
                # R^Im = i * ( eta_i - eta_j )     / 2
                if termInd == 'plus' : return ( None if amps[iInd].CP != amps[jInd].CP else one,                                 None )
                if termInd == 'min'  : return ( None if amps[iInd].CP == amps[jInd].CP else one,                                 None )
                if termInd == 'Re'   : return ( None if amps[iInd].CP != amps[jInd].CP else one if amps[iInd].CP > 0 else minus, None )
                if termInd == 'Im'   : return ( None, None if amps[iInd].CP == amps[jInd].CP else one if amps[iInd].CP > 0 else minus )

            def tCoefTerm( name, dec, mix, sign ) :
                # build one of the two terms within the time function coefficient:
                # + R^+ * 1,  + R^- * C,  + R^- * 1,  + R^+ * C,  + R^Re * D,  + R^Im * S,  + R^Im * D,  - R^Re * S
                if not dec or not mix : return None

                if dec == minus : sign = -sign
                if mix == minus : sign = -sign
                sign = minus if sign < 0 else None
                if dec in [ one, minus ] : dec = None
                if mix in [ one, minus ] : mix = None

                facs = [ fac for fac in ( sign, dec, mix ) if fac ]
                return Product( name, facs ) if len(facs) > 1 else facs[0] if len(facs) == 1 else one

            def tCoef( tCType, iInd, jInd, amps, CPPar ) :
                # build the coefficient of the time function:
                # cosh: + R^+  * 1 + R^-  * C
                # cos:  + R^-  * 1 + R^+  * C
                # sinh: + R^Re * D + R^Im * S
                # sin:  + R^Im * D - R^Re * S
                assert tCType in ( 'cosh', 'sinh', 'cos', 'sin' )

                CPVDecInds = [ 'plus', 'min' ] if tCType in ( 'cosh', 'cos' ) else [ 'Re', 'Im' ]
                if tCType in ( 'cos', 'sin' ) : CPVDecInds = [ CPVDecInds[1], CPVDecInds[0] ]
                CPVDecs = [ CPVDec( CPVDecInd, iInd, jInd, amps, CPPar ) for CPVDecInd in CPVDecInds ]

                CPVMixs = [ one, CPPar.C() ] if tCType in ( 'cosh', 'cos' ) else [ CPPar.D(), CPPar.S() ]
                signs = [ +1, +1 ] if tCType in ( 'cosh', 'cos', 'sinh' ) else [ +1, -1 ]

                name = '%s_%s_%s' % ( iInd , jInd, tCType )
                terms = [ (  tCoefTerm( 'Re_%s%d_%s_%s' % ( tCType, ind, iInd, jInd ), dec[0], mix, sign )
                           , tCoefTerm( 'Im_%s%d_%s_%s' % ( tCType, ind, iInd, jInd ), dec[1], mix, sign )
                          ) for ind, ( dec, mix, sign ) in enumerate( zip( CPVDecs, CPVMixs, signs ) )
                        ]
                assert any( comp for term in terms for comp in term )

                add = lambda name, args : Addition( name, args ) if len(args) > 1 else args[0] if len(args) == 1 else None
                return (  add( 'Re_%s_%s_%s' % ( tCType, iInd, jInd ), [ term[0] for term in terms if term[0] ] )
                        , add( 'Im_%s_%s_%s' % ( tCType, iInd, jInd ), [ term[1] for term in terms if term[1] ] )
                       )

            # create complex amplitude, time and angular factors
            ( reAmps, imAmps ) = ( Re( amplitudes[iIndex], amplitudes[jIndex] ), Im( amplitudes[iIndex], amplitudes[jIndex] ) )\
                             if iIndex != jIndex else ( amplitudes[iIndex].Mag2, None )
            ( reTime, imTime ) = tCoef( tCoefType, iIndex, jIndex, amplitudes, CPParams )
            ( reAng, imAng ) = angFuncs[ ( iIndex, jIndex ) ]

            # (real part of) product of amplitude, time and angular (complex) factors:
            # Re(abc) = Re(ab)Re(c) - Im(ab)Im(c) = Re(a)Re(b)Re(c) - Im(a)Im(b)Re(c) - Re(a)Im(b)Im(c) - Im(a)Re(b)Im(c)
            # (there is a minus in case there are 2 imaginary contributions)
            prod = lambda name, args : [ Product( name, args ) ] if all(args) else [ ]
            return prod( 'ReReRe_%s_%s_%s' % ( tCoefType, amplitudes[iIndex], amplitudes[jIndex] ), [        reAmps , reTime, reAng ] ) \
                 + prod( 'ImImRe_%s_%s_%s' % ( tCoefType, amplitudes[iIndex], amplitudes[jIndex] ), [ minus, imAmps , imTime, reAng ] ) \
                 + prod( 'ReImIm_%s_%s_%s' % ( tCoefType, amplitudes[iIndex], amplitudes[jIndex] ), [ minus, reAmps , imTime, imAng ] ) \
                 + prod( 'ImReIm_%s_%s_%s' % ( tCoefType, amplitudes[iIndex], amplitudes[jIndex] ), [ minus, imAmps , reTime, imAng ] )

        try : # this requires python 2.7 or later...
            from itertools import combinations_with_replacement as cwr
        except:
            from P2VV.Compatibility import cwr

        args = dict()
        from P2VV.RooFitWrappers import Addition
        for tCoefType in [ 'cosh', 'sinh', 'cos', 'sin' ] :
            # NOTE: 'Amplitudes' must be traversed 'in order' : A0, Apar, Aperp, AS, so we cannot use Amplitudes.keys() out of the box
            args[ tCoefType ] = Addition( '%sCoef' % tCoefType
                                         , [ term for ( i, j ) in cwr( Order, 2 )\
                                                  for term in combine( tCoefType, AngFuncs, Amplitudes, CPParams, i, j ) ]
                                        )

        BDecayBasisCoefficients.__init__( self, **args )

class JpsiphiBDecayBasisCoefficients( BDecayBasisCoefficients ) :
    def __init__(self,  angFuncs, Amplitudes,CP, itag, dilution, order ) :
        def combine( name, afun, A, CPparams, tag, i, j) :
            # TODO: deal with tag = None: create the untagged PDF in that case!
            from P2VV.RooFitWrappers import ConstVar, FormulaVar, Product
            plus  = ConstVar( Name = 'plus',  Value =  1 )
            minus = ConstVar( Name = 'minus', Value = -1 )
            if type(CP['C'])==ConstVar and CP['C'].getVal() == 0 : 
                Norm = [ ]
            else :
                Norm = [ FormulaVar('Norm','1.0/(1.0+sign(@0)*@1)',[tag,CP['C']] ) ]
                Norm[0].setAttribute("CacheAndTrack")
            # define functions which return Re(Conj(Ai) Aj), Im( Conj(Ai) Aj)
            # TODO: replace by Addition & Product... why? (only parameters)
            Re        = lambda ai, aj  : FormulaVar('Re_c_%s_%s'%(ai,aj),'@0*@2+@1*@3',[ai.Re,ai.Im,aj.Re,aj.Im])
            Im        = lambda ai, aj  : FormulaVar('Im_c_%s_%s'%(ai,aj),'@0*@3-@1*@2',[ai.Re,ai.Im,aj.Re,aj.Im])
            # define functions which return the coefficients that define the time-dependence...
            _minus_if = lambda b, x : [ minus ] + x if b else  x 
            coef = { 'cosh' : lambda ai,aj,CP : ( plus if ai.CP == aj.CP else  CP['C']
                                                , None )
                   , 'cos'  : lambda ai,aj,CP : ( tag  if ai.CP != aj.CP else Product('Re_%s_%s_cos' %(ai,aj),                              [ tag, CP['C'] ] ) 
                                                , None )
                   , 'sinh' : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP < 0     ,  [      CP['D'] ] ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP < aj.CP ,  [      CP['S'] ] )) )
                   , 'sin'  : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP > 0     ,  [ tag, CP['S'] ] ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP < aj.CP ,  [ tag, CP['D'] ] )) )
                   }
            (c_re,c_im) = coef[name](A[i],A[j],CPparams)
            (f_re,f_im) = afun[(i,j)]
            (a_re,a_im) = ( Re(A[i],A[j]),Im(A[i],A[j]) )
            # NOTE: thes sign are just the obvious Re(a b c)  = Re(a)Re(b)Re(c) - Re(a)Im(b)Im(c) - Im(a)Re(b)Im(c) - Im(a)Im(b)Re(c),
            #       i.e. there is a minus in case there are 2 imaginary contributions
            prod = lambda name, args : [ Product(name, args) ] if all(args) else []
            s  = prod('ReReRe_%s_%s_%s'%(name,A[i],A[j]), [        a_re , c_re, f_re ]+Norm ) \
               + prod('ImImRe_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_im , c_im, f_re ]+Norm ) \
               + prod('ImReIm_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_im , c_re, f_im ]+Norm ) \
               + prod('ReImIm_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_re , c_im, f_im ]+Norm )
            assert len(s) == 1 # for now, coefficients are either real, or imaginary, but not both... (not true in general, but I'm lazy today ;-)
            return s[0]

        args = dict()
        from P2VV.RooFitWrappers import Addition,Product, RealCategory
        try : # this requires python 2.7 or later...
            from itertools import combinations_with_replacement as cwr
        except:
            from P2VV.Compatibility import cwr

        tag = Product('tag',( RealCategory('tag_real', itag ),dilution))
        for name in [ 'cosh', 'sinh', 'cos', 'sin' ] :
            # NOTE: 'Amplitudes'  must be traversed 'in order' -- so we cannot use Amplitudes.keys() out of the box, but use the
            #       specified order (which also selects which ones to include!)...
            args[ name ] = Addition( 'a_%s'% name, [ combine(name,angFuncs,Amplitudes,CP,tag,i,j) for (i,j) in cwr( order, 2 ) ] )

        BDecayBasisCoefficients.__init__( self, **args )


def JpsiphiBDecay( Name, time, tag, lifetimeParams, sigtres, tagging,  angles, amplitudes, CP, order ) :
    from P2VV.Parameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
    basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes, CP, tag, tagging['dilution'], order )
    from P2VV.RooFitWrappers import BDecay
    return  BDecay(  Name                   = Name
                   , time                   = time
                   , dm                     = lifetimeParams['dM']
                   , tau                    = lifetimeParams['MeanLifetime']
                   , dGamma                 = lifetimeParams['dGamma']
                   , resolutionModel        = sigtres['model']
                   , coshCoef               = basisCoefficients['cosh']
                   , cosCoef                = basisCoefficients['cos']
                   , sinhCoef               = basisCoefficients['sinh']
                   , sinCoef                = basisCoefficients['sin']
                   , ConditionalObservables = sigtres.conditionalObservables() + tagging.conditionalObservables()
                   , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                              + sigtres.externalConstraints()\
                                              + tagging.externalConstraints()
                  )

def JpsiphiBTagDecay( Name, time, tag, lifetimeParams, sigtres, tagging,  angles, amplitudes, CP, order ) :
    from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
    basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes, CP, order )
    from P2VV.RooFitWrappers import BTagDecay
    return  BTagDecay(  Name                   = Name
                      , time                   = time
                      , iTag                   = tag
                      , dm                     = lifetimeParams['dM']
                      , tau                    = lifetimeParams['MeanLifetime']
                      , dGamma                 = lifetimeParams['dGamma']
                      , resolutionModel        = sigtres['model']
                      , coshCoef               = basisCoefficients['cosh']
                      , cosCoef                = basisCoefficients['cos']
                      , sinhCoef               = basisCoefficients['sinh']
                      , sinCoef                = basisCoefficients['sin']
                      , dilution               = tagging['dilution']
                      , ADilWTag               = tagging['ADilWTag']
                      , avgCEven               = tagging['avgCEven']
                      , avgCOdd                = tagging['avgCOdd']
                      , ConditionalObservables = sigtres.conditionalObservables() + tagging.conditionalObservables()
                      , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                                 + sigtres.externalConstraints()\
                                                 + tagging.externalConstraints()
                     )




from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin
class TimePdf( _util_parse_mixin ) :
    def __init__( self, **kwargs ) :
        self._pdf = kwargs.pop('pdf')
        self._efficiency = kwargs.pop( 'Efficiency', None )
        if self._efficiency : self._pdf = self._efficiency * self._pdf
        for (k,v) in kwargs.iteritems() : setattr( self, '_' + k, v )

    def pdf(self) : return self._pdf

class LP2011_Background_Time( TimePdf ) :
    def __init__(self, time, resolutionModel, **kwargs) :
        Name = kwargs.pop('Name', self.__class__.__name__)
        self._ml_tau = self._parseArg( '%s_ml_tau' % Name, kwargs, Title = 'medium lifetime background ', Unit = 'ps'
                                      , Value = 0.152, Error = 0.003, MinMax = ( 0.1, 10. ) )
        self._ll_tau = self._parseArg( '%s_ll_tau' % Name, kwargs, Title = 'long lifetime background ',   Unit = 'ps'
                                      , Value = 1.06,  Error = 0.04,  MinMax = ( 0.1, 10. ) )
        self._fml = self._parseArg(    '%s_fml' % Name,    kwargs, Title = 'fraction medium lifetime background'
                                      , Value = 0.79,  Error = 0.01,  MinMax = ( 0.,   1. ) )

        from P2VV.RooFitWrappers import  SumPdf,Pdf
        from ROOT import RooDecay as Decay
        self._mlPdf = Pdf(  Name = Name + '_ml'
                          , Type = Decay
                          , Parameters = (time, self._ml_tau, resolutionModel, 'SingleSided')
                          , ConditionalObservables = resolutionModel.ConditionalObservables()
                          , ExternalConstraints = resolutionModel.ExternalConstraints()
                         )
        self._llPdf = Pdf(  Name = Name + '_ll'
                          , Type = Decay
                          , Parameters = (time, self._ll_tau, resolutionModel, 'SingleSided')
                          , ConditionalObservables = resolutionModel.ConditionalObservables()
                          , ExternalConstraints = resolutionModel.ExternalConstraints()
                         )

        efficiency = kwargs.pop( 'Efficiency', None )
        if efficiency :
            self._mlPdf = efficiency * self._mlPdf
            self._llPdf = efficiency * self._llPdf

        self._check_extraneous_kw( kwargs )
        TimePdf.__init__( self, pdf = SumPdf(  Name = Name
                                             , PDFs = ( self._mlPdf, self._llPdf)
                                             , Yields = { self._mlPdf.GetName() : self._fml }
                                            )
                        )
        self._efficiency = efficiency

class Single_Exponent_Time( TimePdf ) :
    def __init__(self,time,resolutionModel,**kwargs) :
        self._parseArg('tau', kwargs, Title = 'lifetime', Unit = 'ps', Value = 1.5, Error = 0.05, MinMax = ( 0.1, 10. ) )

        from P2VV.RooFitWrappers import Pdf
        from ROOT import RooDecay as Decay
        TimePdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name',self.__class__.__name__)
                                          , Type = Decay
                                          , Parameters = (time, self._tau,resolutionModel,'SingleSided')
                                          , ConditionalObservables = resolutionModel.ConditionalObservables()
                                          , ExternalConstraints = resolutionModel.ExternalConstraints())
                         )

class Prompt_Peak( TimePdf ) :
    def __init__(self, time, resolutionModel, **kwargs) :
        self._parseArg('prompt_tau', kwargs, Title = 'lifetime', Unit = 'ps', Value = 0, Constant = True)

        from P2VV.RooFitWrappers import Pdf
        from ROOT import RooDecay as Decay
        TimePdf.__init__(self, pdf = Pdf( Name = kwargs.pop('Name',self.__class__.__name__)
                                          , Type = Decay
                                          , Parameters = (time, self._prompt_tau, resolutionModel, 'SingleSided')
                                          , ConditionalObservables = resolutionModel.ConditionalObservables()
                                          , ExternalConstraints = resolutionModel.ExternalConstraints())
                         )
