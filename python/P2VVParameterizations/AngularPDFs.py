###########################################################################################################################################
## P2VVParameterizations.AngularPDFs: Parameterizations of PDFs that only depend on decay angles                                         ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################


class AngularPdfTerms ( list ) :
    def __init__( self, AngTerms ) :
        self += ( term for term in AngTerms )

    def __getitem__( self, keyWord ) :
        if type(keyWord) == str : return getattr( self, '_' + keyWord )
        else :                    return list.__getitem__(self, keyWord)

    def buildSumPdf( self, Name ) :
        from RooFitWrappers import RealSumPdf
        return RealSumPdf( Name, self )

class Coefficients_AngularPdfTerms ( AngularPdfTerms ) :
    def __init__( self, **kwargs ) :
        # get angular functions from kwargs
        from RooFitWrappers import __check_req_kw__
        __check_req_kw( 'AngFunctions', kwargs )
        self._angFuncs = kwargs.pop('AngFunctions')

        # get angular function coefficients
        if 'AngCoefficients' in kwargs :
            self._angCoefs = kwargs.pop('AngCoefficients')
            for key in self._angFuncs.keys() :
                if key not in self._angCoefs : raise KeyError('Coefficients_AngularPdfTerms: no coefficient %s found' % str(key))
        else :
            from RooFitWrappers import RealVar
            self._angCoefs = dict( ( k, (  RealVar( '%s_%s_ReCoef' %( k[0],k[1] ), Value = 1. ) if self._angFuncs[k][0] else None
                                        ,  RealVar( '%s_%s_ImCoef' %( k[0],k[1] ), Value = 1. ) if self._angFuncs[k][1] else None ) ) for k, func in self._angFuncs.iteritems() )

        # check if there are no arguments left
        if kwargs: raise KeyError('Coefficients_AngularPdfTerms: got unknown keyword%s: %s'\
                                        % ( '' if len(kwargs) == 1 else 's', kwargs ))

        # build angular terms
        from RooFitWrappers import ConstVar, Product
        minus = ConstVar('minus',  Value = -1  )
        newAngTerm = lambda func, coef, minSign :\
              [ Product( coef.GetName() + '_x_' + func.GetName(), [ minSign, coef, func ] if minSign else [ coef, func ] ) ]\
              if func and coef else [ ]

        angTerms = []
        for key in self._angFuncs.keys() :
            angTerms += newAngTerm( self._angFuncs[key][0], self._angCoefs[key][0], None  )
            angTerms += newAngTerm( self._angFuncs[key][1], self._angCoefs[key][1], minus )

        # initialize
        AngularPdfTerms.__init__( self, angTerms )

class Amplitudes_AngularPdfTerms ( Coefficients_AngularPdfTerms ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import __check_req_kw__
        __check_req_kw__('AmpNames',kwargs)
        __check_req_kw__('Amplitudes',kwargs)
        __check_req_kw__('AngFunctions',kwargs)

        try :   from itertools import combinations_with_replacement as cwr
        except: from compatibility import cwr

        # get amplitude names from arguments
        self._ampNames = kwargs.pop('AmpNames')

        # get amplitudes from arguments
        self._amplitudes = kwargs.pop('Amplitudes')
        for amp in self._ampNames : assert amp in self._amplitudes, 'Amplitudes_AngularPdfTerms: no amplitude \'%s\' found' % amp

        # get angular functions from arguments
        angFuncs = { }
        angFuncsArg = kwargs.pop('AngFunctions')
        for amp1, amp2 in cwr( self._ampNames, 2 ) :
            assert ( amp1, amp2 ) in angFuncsArg, 'Amplitudes_AngularPdfTerms: no angular function %s found' % str(( amp1, amp2 ))
            angFuncs[( amp1, amp2 )] = angFuncsArg[( amp1, amp2 )]

        # check if there are no arguments left
        if kwargs: raise KeyError('Amplitudes_AngularPdfTerms: got unknown keyword%s: %s'\
                                        % ( '' if len(kwargs) == 1 else 's', kwargs ))

        # build angular coefficients
        from RooFitWrappers import FormulaVar
        angCoefs = { }
        Re = lambda Ai, Aj : FormulaVar( 'Re_c_%s_%s' % ( Ai, Aj ), '@0*@2 + @1*@3', [ Ai.Re, Ai.Im, Aj.Re, Aj.Im ] )
        Im = lambda Ai, Aj : FormulaVar( 'Im_c_%s_%s' % ( Ai, Aj ), '@0*@3 - @1*@2', [ Ai.Re, Ai.Im, Aj.Re, Aj.Im ] )
        for amp1, amp2 in cwr( self._ampNames, 2 ) :
            angCoefs[( amp1, amp2 )] = (  Re( self._amplitudes[amp1], self._amplitudes[amp2] )
                                        , Im( self._amplitudes[amp1], self._amplitudes[amp2] ) )
        # initialize
        Coefficients_AngularPdfTerms.__init__( self, AngCoefficients = angCoefs, AngFunctions = angFuncs )
