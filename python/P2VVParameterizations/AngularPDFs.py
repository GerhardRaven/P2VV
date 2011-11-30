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
        for term in AngTerms : self.append(term)

    def __getitem__( self, keyWord ) :
        if type(keyWord) == str : return getattr( self, '_' + keyWord )
        else :                    return list.__getitem__(self, keyWord)

    def buildSumPdf( self, Name ) :
        from RooFitWrappers import RealSumPdf
        return RealSumPdf( Name, self )

class Coefficients_AngularPdfTerms ( AngularPdfTerms ) :
    def __init__( self, **kwargs ) :
        # get angular functions from kwargs
        if 'AngFunctions' in kwargs : self._angFuncs = kwargs.pop('AngFunctions')
        else : raise KeyError('Coefficients_AngularPdfTerms: no key word \'AngFunctions\' found')

        # get angular function coefficients
        if 'AngCoefficients' in kwargs :
            self._angCoefs = kwargs.pop('AngCoefficients')
            for key in self._angFuncs.keys() :
                if key not in self._angCoefs : raise KeyError('Coefficients_AngularPdfTerms: no coefficient %s found' % str(key))
        else :
            from RooFitWrappers import RealVar
            self._angCoefs = {}
            for key, func in self._angFuncs.iteritems() :
                self._angCoefs[key] = (  RealVar( key[0] + '_' + key[1] + '_ReCoef', Value = 1. ) if self._angFuncs[key][0] else None
                                       , RealVar( key[0] + '_' + key[1] + '_ImCoef', Value = 1. ) if self._angFuncs[key][1] else None )

        # check if there are no arguments left
        if len(kwargs): raise KeyError('Coefficients_AngularPdfTerms: got unknown keyword%s: %s'\
                                        % ( '' if len(kwargs) == 1 else 's', kwargs ))

        # build angular terms
        from RooFitWrappers import ConstVar, Product
        angTerms = []
        minus = ConstVar('minus',  Value = -1  )
        newAngTerm = lambda func, coef, minSign :\
              [ Product( coef.GetName() + '_x_' + func.GetName(), [ minSign, coef, func ] if minSign else [ coef, func ] ) ]\
              if func and coef else [ ]

        for key in self._angFuncs.keys() :
            angTerms += newAngTerm( self._angFuncs[key][0], self._angCoefs[key][0], None  )
            angTerms += newAngTerm( self._angFuncs[key][1], self._angCoefs[key][1], minus )

        # initialize
        AngularPdfTerms.__init__( self, angTerms )

class Amplitudes_AngularPdfTerms ( Coefficients_AngularPdfTerms ) :
    def __init__( self, **kwargs ) :
        try :   from itertools import combinations_with_replacement as cwr
        except: from compatibility import cwr

        # get amplitude names from arguments
        if 'AmpNames' in kwargs : self._ampNames = kwargs.pop('AmpNames')
        else : raise KeyError('Amplitudes_AngularPdfTerms: no key word \'AmpNames\' found')

        # get amplitudes from arguments
        if 'Amplitudes' in kwargs :
            self._amplitudes = kwargs.pop('Amplitudes')
            for amp in self._ampNames : assert amp in self._amplitudes, 'Amplitudes_AngularPdfTerms: no amplitude \'%s\' found' % amp
        else :
            raise KeyError('Amplitudes_AngularPdfTerms: no key word \'Amplitudes\' found')

        # get angular functions from arguments
        if 'AngFunctions' in kwargs :
            angFuncs = { }
            angFuncsArg = kwargs.pop('AngFunctions')
            for amp1, amp2 in cwr( self._ampNames, 2 ) :
                assert ( amp1, amp2 ) in angFuncsArg, 'Amplitudes_AngularPdfTerms: no angular function %s found' % str(( amp1, amp2 ))
                angFuncs[( amp1, amp2 )] = angFuncsArg[( amp1, amp2 )]
        else :
            raise KeyError('Amplitudes_AngularPdfTerms: no key word \'AngFunctions\' found')

        # check if there are no arguments left
        if len(kwargs): raise KeyError('Amplitudes_AngularPdfTerms: got unknown keyword%s: %s'\
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

