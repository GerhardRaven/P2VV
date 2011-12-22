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
    def __init__( self, angFuncs, Amplitudes, CP, order ) :
        def combine( name, afun, A, CPparams, i, j ) :
            from RooFitWrappers import ConstVar, FormulaVar, Product
            plus  = ConstVar('plus', Value = 1)
            minus = ConstVar('minus',  Value = -1  )
            one   = plus
            # define functions which return Re(Conj(Ai) Aj), Im( Conj(Ai) Aj)
            # TODO: replace by Addition & Product... why? (only parameters)
            Re        = lambda ai, aj  : FormulaVar('Re_c_%s_%s'%(ai,aj),'@0*@2+@1*@3',[ai.Re,ai.Im,aj.Re,aj.Im])
            Im        = lambda ai, aj  : FormulaVar('Im_c_%s_%s'%(ai,aj),'@0*@3-@1*@2',[ai.Re,ai.Im,aj.Re,aj.Im])
            # define functions which return the coefficients that define the time-dependence...
            _minus_if = lambda b, x : [ minus, x ] if b else [ x ]
            coef = { 'cosh' : lambda ai,aj,CP : ( one  if ai.CP == aj.CP else CP['C']  # make a trivial product just to get the labels right???
                                                , None )
                   , 'cos'  : lambda ai,aj,CP : ( CP['C'] if ai.CP == aj.CP else one 
                                                , None )
                   , 'sinh' : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP > 0    , CP['D'] ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP < aj.CP, CP['S'] )) )
                   , 'sin'  : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP > 0    , CP['S'] ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP > aj.CP, CP['D'] )) )
                   }
            (c_re,c_im) = coef[name](A[i],A[j],CPparams)
            (f_re,f_im) = afun[(i,j)]
            (a_re,a_im) = ( Re(A[i],A[j]),Im(A[i],A[j]) )
            # NOTE: thes sign are just the obvious Re(a b c)  = Re(a)Re(b)Re(c) - Re(a)Im(b)Im(c) - Im(a)Re(b)Im(c) - Im(a)Im(b)Re(c),
            #       i.e. there is a minus in case there are 2 imaginary contributions
            prod = lambda name, args : [ Product(name, args) ] if all(args) else []
            s  = prod('ReReRe_%s_%s_%s'%(name,A[i],A[j]), [        a_re , c_re, f_re ] ) \
               + prod('ImImRe_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_im , c_im, f_re ] ) \
               + prod('ImReIm_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_im , c_re, f_im ] ) \
               + prod('ReImIm_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_re , c_im, f_im ] )
            assert len(s) == 1 # for now, coefficients are either real, or imaginary, but not both... (not true in general, but I'm lazy today ;-)
            return s[0]

        args = dict()
        from RooFitWrappers import Addition
        try : # this requires python 2.7 or later...
            from itertools import combinations_with_replacement as cwr
        except:
            from compatibility import cwr

        for name in [ 'cosh', 'sinh', 'cos', 'sin' ] :
            # NOTE: 'Amplitudes'  must be traversed 'in order' : A0, Apar, Aperp, AS -- so we cannot use Amplitudes.keys() out of the box...
            args[ name ] = Addition( 'a_%s'% name, [ combine(name,angFuncs,Amplitudes,CP,i,j) for (i,j) in cwr( order, 2 ) ] )

        BDecayBasisCoefficients.__init__( self, **args )

class JpsiphiBDecayBasisCoefficients( BDecayBasisCoefficients ) :
    def __init__(self,  angFuncs, Amplitudes,CP, tag, order ) :
        def combine( name, afun, A, CPparams, tag, i, j) :
            from RooFitWrappers import ConstVar, FormulaVar, Product
            plus  = ConstVar('plus', Value = 1)
            minus = ConstVar('minus',  Value = -1  )
            Norm = FormulaVar('Norm','1.0/(1.0+@0*@1)',[tag,CP['C']] )
            # define functions which return Re(Conj(Ai) Aj), Im( Conj(Ai) Aj)
            # TODO: replace by Addition & Product... why? (only parameters)
            Re        = lambda ai, aj  : FormulaVar('Re_c_%s_%s'%(ai,aj),'@0*@2+@1*@3',[ai.Re,ai.Im,aj.Re,aj.Im])
            Im        = lambda ai, aj  : FormulaVar('Im_c_%s_%s'%(ai,aj),'@0*@3-@1*@2',[ai.Re,ai.Im,aj.Re,aj.Im])
            # define functions which return the coefficients that define the time-dependence...
            _minus_if = lambda b, x : [ minus ] + x if b else  x 
            coef = { 'cosh' : lambda ai,aj,CP : ( Norm  if ai.CP == aj.CP else Product("Re_%s_%s_cosh", [ Norm, CP['C'] ] )
                                                , None )
                   , 'cos'  : lambda ai,aj,CP : ( Product('Re_%s_%s_cos'%(ai,aj), [tag, Norm ] + ( [ CP['C'] ] if ai.CP == aj.CP else [ ] )  )
                                                , None )
                   , 'sinh' : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP > 0     ,  [      Norm, CP['D'] ] ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sinh'%(ai,aj),  _minus_if( ai.CP < aj.CP ,  [      Norm, CP['S'] ] )) )
                   , 'sin'  : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('Re_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP > 0     ,  [ tag, Norm, CP['S'] ] ))
                                                , None if ai.CP == aj.CP else Product('Im_%s_%s_sin' %(ai,aj),  _minus_if( ai.CP > aj.CP ,  [ tag, Norm, CP['D'] ] )) )
                   }
            (c_re,c_im) = coef[name](A[i],A[j],CPparams)
            (f_re,f_im) = afun[(i,j)]
            (a_re,a_im) = ( Re(A[i],A[j]),Im(A[i],A[j]) )
            # NOTE: thes sign are just the obvious Re(a b c)  = Re(a)Re(b)Re(c) - Re(a)Im(b)Im(c) - Im(a)Re(b)Im(c) - Im(a)Im(b)Re(c),
            #       i.e. there is a minus in case there are 2 imaginary contributions
            prod = lambda name, args : [ Product(name, args) ] if all(args) else []
            s  = prod('ReReRe_%s_%s_%s'%(name,A[i],A[j]), [        a_re , c_re, f_re ] ) \
               + prod('ImImRe_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_im , c_im, f_re ] ) \
               + prod('ImReIm_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_im , c_re, f_im ] ) \
               + prod('ReImIm_%s_%s_%s'%(name,A[i],A[j]), [ minus, a_re , c_im, f_im ] )
            assert len(s) == 1 # for now, coefficients are either real, or imaginary, but not both... (not true in general, but I'm lazy today ;-)
            return s[0]

        args = dict()
        from RooFitWrappers import Addition
        try : # this requires python 2.7 or later...
            from itertools import combinations_with_replacement as cwr
        except:
            from compatibility import cwr

        for name in [ 'cosh', 'sinh', 'cos', 'sin' ] :
            # NOTE: 'Amplitudes'  must be traversed 'in order' -- so we cannot use Amplitudes.keys() out of the box, but use the
            #       specified order (which also selects which ones to include!)...
            args[ name ] = Addition( 'a_%s'% name, [ combine(name,angFuncs,Amplitudes,CP,tag,i,j) for (i,j) in cwr( order, 2 ) ] )

        BDecayBasisCoefficients.__init__( self, **args )


from P2VVParameterizations.GeneralUtils import _util_parse_mixin
class TimePdf( _util_parse_mixin ) :
    def __init__(self, **kwargs ) :
        for (k,v) in kwargs.iteritems() :
            setattr(self,'_'+k,v)
    def pdf(self) :
        return self._pdf

class LP2011_Background_Time( TimePdf ) :
    def __init__(self,time,resolutionModel,**kwargs) :
        self._parseArg('t_bkg_ml_tau', kwargs, Title = 'medium lifetime background ', Unit = 'ps', Value = 0.21, MinMax = (0.01,0.5) )
        self._parseArg('t_bkg_ll_tau', kwargs, Title = 'long lifetime background ', Unit = 'ps', Value = 1.92, MinMax = (0.5,2.5) )
        self._parseArg('t_bkg_fll',kwargs, Title = 'fraction long lifetime background', Value = 0.3, MinMax = (0., 1.) )
        from RooFitWrappers import  SumPdf,Pdf
        from ROOT import RooDecay as Decay
        Name = kwargs.pop('Name',self.__class__.__name__)
        ml = Pdf( kwargs.pop('t_bkg_ml',Name+'_t_bkg_ml'),Type =Decay, Parameters = (time,self._t_bkg_ml_tau,resolutionModel,'SingleSided'))
        ll = Pdf( kwargs.pop('t_bkg_ll',Name+'_t_bkg_ll'),Type =Decay, Parameters = (time,self._t_bkg_ll_tau,resolutionModel,'SingleSided'))
        TimePdf.__init__(self, pdf = SumPdf(Name = Name,   PDFs = (  ml, ll)  , Yields = { ml.GetName() : self._t_bkg_fll } ) )


class Single_Exponent_Time( TimePdf ) :
    def __init__(self,time,resolutionModel,**kwargs) :
        self._parseArg('t_sig_tau', kwargs, Title = 'lifetime', Unit = 'ps', Value = 1.5, MinMax = (0.5,2.5) )
        from RooFitWrappers import Pdf
        from ROOT import RooDecay as Decay
        TimePdf.__init__(self, pdf = Pdf(kwargs.pop('Name',self.__class__.__name__),Type =Decay, Parameters = (time, self._t_sig_tau,resolutionModel,'SingleSided')) )
