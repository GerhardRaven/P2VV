###########################################################################################################################################
## P2VVParameterizations.AngularPDFs: Parameterizations of PDFs that only depend on decay angles                                         ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin


class AngularPdfTerms ( list ) :
    def __init__( self, AngTerms ) :
        # set the angular terms
        self += ( term for term in AngTerms )

    def __getitem__( self, keyWord ) :
        # return an attribute if the key word is a string and an angular term if the key word is an integer
        if type(keyWord) == str : return getattr( self, '_' + keyWord )
        else :                    return list.__getitem__(self, keyWord)

    def buildSumPdf( self, Name ) :
        # build a RealSumPdf from the current angular terms
        from P2VV.RooFitWrappers import RealSumPdf
        return RealSumPdf( Name, self )

    def buildAddition( self, Name ) :
        # build an Addition function from the current angular terms
        from P2VV.RooFitWrappers import Addition
        return Addition( Name, self )

class Coefficients_AngularPdfTerms ( AngularPdfTerms ) :
    def __init__( self, **kwargs ) :
        # get angular functions from arguments (dictionary with complex ( Re, Im ) wrapped RooFit objects)
        from P2VV.RooFitWrappers import __check_req_kw__
        __check_req_kw__( 'AngFunctions', kwargs )
        self._angFuncs = kwargs.pop('AngFunctions')

        # get keys for angular terms
        self._keys = kwargs.pop( 'Keys', self._angFuncs.keys() )

        # get angular function coefficients (dictionary with complex ( Re, Im ) wrapped RooFit objects)
        if 'AngCoefficients' in kwargs :
            # get coefficients from arguments
            self._angCoefs = kwargs.pop('AngCoefficients')
            for key in self._keys :
                if key not in self._angCoefs : raise KeyError('Coefficients_AngularPdfTerms: no coefficient %s found' % str(key))
        else :
            # build a coefficient for each angular function
            from P2VV.RooFitWrappers import RealVar
            self._angCoefs = dict( ( key, (  RealVar('%s_%s_ReCoef' % (key[0], key[1]), Value = 1.) if self._angFuncs[key][0] else None
                                           , RealVar('%s_%s_ImCoef' % (key[0], key[1]), Value = 1.) if self._angFuncs[key][1] else None) )\
                                     for key, func in self._angFuncs.iteritems() )

        # check if there are arguments left
        if kwargs: raise KeyError('Coefficients_AngularPdfTerms: got unknown keyword%s: %s' % ( '' if len(kwargs) == 1 else 's', kwargs ) )

        # build angular terms ( list of products of coefficients and functions: Re[C*F] = Re(C)*Re(F) - Im(C)*Im(F) )
        from P2VV.RooFitWrappers import ConstVar, Product
        minus = ConstVar(Name = 'minus',  Value = -1  )
        newAngTerm = lambda func, coef, minSign :\
              [ Product( coef.GetName() + '_x_' + func.GetName(), [ minSign, coef, func ] if minSign else [ coef, func ] ) ]\
              if func and coef else [ ]

        angTerms = []
        for key in self._keys :
            angTerms += newAngTerm( self._angFuncs[key][0], self._angCoefs[key][0], None  )
            angTerms += newAngTerm( self._angFuncs[key][1], self._angCoefs[key][1], minus )

        # initialize
        AngularPdfTerms.__init__( self, angTerms )

class AngleBasis_AngularPdfTerms ( Coefficients_AngularPdfTerms ) :
    def __init__( self, **kwargs ) :
        # get angles from arguments
        from P2VV.RooFitWrappers import __check_req_kw__
        __check_req_kw__( 'Angles', kwargs )
        self._angles = kwargs.pop('Angles')

        # initialize coefficients and functions lists
        keys     = []
        angCoefs = {}
        angFuncs = {}

        # get coefficients and angular basis function indices from arguments
        hasC000 = False
        for coefName in sorted(kwargs.keys()) :
            # get coefficient's arguments
            coefArgs = kwargs[coefName]

            # check if this argument is a coefficient
            if coefName[0] != 'C' or type(coefArgs) != dict : continue

            # process arguments of coefficient
            indices = coefArgs.pop( 'Indices', None )
            if not indices : raise KeyError( 'AngleBasis_AngularPdfTerms: no indices found for coefficient %s' % coefName )
            if indices == (0,0,0) : hasC000 = True

            # get (RooFit) name of angular term
            cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
            RooFitCoefName = coefArgs.pop( 'Name', 'Cab%d%d%s' % ( indices[0], indices[1], cnvrtInd(indices[2]) ) )

            # create coefficient and angular function
            from P2VV.RooFitWrappers import RealVar, P2VVAngleBasis
            keys.append( ( coefName, None ) )
            angCoefs[ ( coefName, None ) ] = ( RealVar( RooFitCoefName, **coefArgs ), None )
            angFuncs[ ( coefName, None ) ] = (  P2VVAngleBasis( Name = 'AB', Angles = self._angles
                                                               , Indices = ( indices[0], 0, indices[1], indices[2] ) , FixedCoef = 1. )
                                              , None
                                             )

        # remove coefficients from arguments
        for key in keys : kwargs.pop(key[0])

        # add P_0 Y_0_0 term ( 0, 0, 0 )
        if not hasC000 :
            from P2VV.RooFitWrappers import ConstVar,P2VVAngleBasis
            keys.insert( 0, ( 'C000', None ) )
            angCoefs[ keys[0] ] = ( ConstVar( Name = 'Cab000', Value = 1. ), None )
            angFuncs[ keys[0] ] = ( P2VVAngleBasis( Name = 'AB', Angles = self._angles, Indices = (0, 0, 0, 0), FixedCoef = 1. ), None )

        # check if there are no arguments left
        if kwargs: raise KeyError('AngleBasis_AngularPdfTerms: got unknown keyword%s: %s'\
                                        % ( '' if len(kwargs) == 1 else 's', kwargs ))

        # initialize
        Coefficients_AngularPdfTerms.__init__( self, Keys = keys, AngCoefficients = angCoefs, AngFunctions = angFuncs )


class Amplitudes_AngularPdfTerms ( Coefficients_AngularPdfTerms ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import __check_req_kw__
        __check_req_kw__('AmpNames',kwargs)
        __check_req_kw__('Amplitudes',kwargs)
        __check_req_kw__('AngFunctions',kwargs)

        try :   from itertools import combinations_with_replacement as cwr
        except: from P2VV.Compatibility import cwr

        # get amplitude names from arguments
        self._ampNames = kwargs.pop('AmpNames')

        # get amplitudes from arguments
        self._amplitudes = kwargs.pop('Amplitudes')
        for amp in self._ampNames : assert amp in self._amplitudes, 'Amplitudes_AngularPdfTerms: no amplitude \'%s\' found' % amp

        # get keys for angular terms
        keys = [ key for key in cwr( self._ampNames, 2 ) ]

        # get angular functions from arguments
        angFuncs = { }
        angFuncsArg = kwargs.pop('AngFunctions')
        for key in keys :
            assert key in angFuncsArg, 'Amplitudes_AngularPdfTerms: no angular function %s found' % str(key)
            angFuncs[key] = angFuncsArg[key]

        # check if there are no arguments left
        if kwargs: raise KeyError('Amplitudes_AngularPdfTerms: got unknown keyword%s: %s'\
                                        % ( '' if len(kwargs) == 1 else 's', kwargs ))

        # build angular coefficients
        from P2VV.RooFitWrappers import FormulaVar
        angCoefs = { }
        Re = lambda Ai, Aj : FormulaVar( 'Re_c_%s_%s' % ( Ai, Aj ), '@0*@2 + @1*@3', [ Ai.Re, Ai.Im, Aj.Re, Aj.Im ] )
        Im = lambda Ai, Aj : FormulaVar( 'Im_c_%s_%s' % ( Ai, Aj ), '@0*@3 - @1*@2', [ Ai.Re, Ai.Im, Aj.Re, Aj.Im ] )
        for key in keys :
            angCoefs[key] = (  Re( self._amplitudes[ key[0] ], self._amplitudes[ key[1] ] )
                             , Im( self._amplitudes[ key[0] ], self._amplitudes[ key[1] ] ) )
        # initialize
        Coefficients_AngularPdfTerms.__init__( self, Keys = keys, AngCoefficients = angCoefs, AngFunctions = angFuncs )


class Uniform_Angles( _util_parse_mixin ) :
    def pdf(self) :
        return self._pdf        
    def __init__( self, angles, **kwargs ) :
        from P2VV.RooFitWrappers import UniformPdf
        self._pdf =  UniformPdf(Name = kwargs.pop('Name','Uniform_AnglesPdf'), Arguments = ( angles['phi'],angles['ctheta'],angles['cpsi'] ) )
        for (k,v) in kwargs.iteritems() :
            setattr(self,'_'+k,v)


class SPlot_Moment_Angles( object ) :
    def pdf(self, **kwargs ) :
        from P2VV.Utilities.DataMoments import RealMomentsBuilder
        mb = RealMomentsBuilder()
        mb.appendPYList( self._angles, kwargs.pop('Indices') )
        mb.compute( self._splot.data( kwargs.pop('Component') ) )
        if kwargs.pop('Verbose',False) : 
            mb.Print()
            print 'computed moments, creating PDF with name %s' % kwargs['Name']
        return mb.createPDF( **kwargs )

    def __init__(self, angles, splot ) :
        self._splot  = splot
        self._angles = angles

