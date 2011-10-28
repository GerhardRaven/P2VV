from RooFitWrappers import *
from ROOT import gSystem
gSystem.Load('libP2VV.so')
ws = RooObject()
ws.setWorkspace( RooWorkspace('myworkspace') )

plus  = ConstVar('plus', Value = 1)
minus = ConstVar('minus',  Value = -1  )
one   = plus

args = dict()
from math import pi
cpsiAng   = RealVar(  'helcthetaK', Title = 'cosine of kaon polarization angle',   Observable = True,  MinMax=(-1., 1.))
cthetaAng = RealVar(  'helcthetaL', Title = 'cosine of lepton polarization angle', Observable = True,  MinMax=(-1., 1.))
phiAng    = RealVar(  'helphi',     Title = 'angle between decay planes',          Observable = True,  MinMax=(-pi, pi))
t         = RealVar(  't',          Title = 'decay time', Unit='ps',               Observable = True,  MinMax=(-5,14)  )

unbiased  = Category( 'unbiased',   Title = 'Unbiased trigger?',         Observable = True, States = { 'yes': 1 ,   'no' : 0 } )
decay     = Category( 'decaytype',  Title = 'J/psiX decay type',         Observable = True, States = { 'JpsiKplus' : 10, 'JpsiKmin' : 11, 'JpsiKstar0' : 20, 'JpsiKstar0bar': 21, 'Jpsiphi': 40 } )
iTag      = Category( 'iTag',       Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1 } )

from ROOT import RooTruthModel as TruthModel
args.update( { 'dm'     : RealVar( 'dm',        Title = 'delta m',       Unit = 'ps^{-1}', Observable = False, Value = 17.8 )  
             , 'tau'    : RealVar( 't_sig_tau', Title = 'mean lifetime', Unit = 'ps',      Observable = False, Value =  1.5, MinMax = (1.3, 1.8) )
             , 'dGamma' : RealVar( 'dGamma',    Title = 'dGamma',        Unit = 'ps^{-1}', Observable = False, Value = 0.05, MinMax = (-0.3,0.3) )
             , 'resolutionModel' : ResolutionModel( 'resModel', Type = TruthModel, Observables = [ t ] )
             , 'decayType' : 'SingleSided'
             } )

# TODO: package this in a seperate class, which gives access to CCP, DCP, and SCP, so we can hide the 'internal' parameterization
# TODO: explicitly export externally visible symbols .... For now, interal stuff is prefixed by an _
_phiCP      = RealVar( 'phiCP',    Title = 'CP violation param. phi_s',     Observable = False, Value = -0.2,  MinMax = (-2*pi,2*pi) )
_lambdaCPSq = RealVar( 'lambda^2', Title = 'CP violation param |lambda|^2', Observable = False, Value = 1)

class CPParam :
    def __init__(self,**kwargs) :
        for i in 'CDS' : setattr(self,i,kwargs.pop(i))

CP = CPParam( C = FormulaVar('C', '(1.-@0)/(1.+@0)',                   [ _lambdaCPSq ] )
            , D = FormulaVar('D', ' 2 * sqrt(@0) * cos(@1) / (1.+@0)', [ _lambdaCPSq, _phiCP ] )
            , S = FormulaVar('S', '-2 * sqrt(@0) * sin(@1) / (1.+@0)', [ _lambdaCPSq, _phiCP ] )
            )

tagDilution = RealVar( 'tagDilution', Title = 'Average Tagging Dilution',     Observable = False, Value = 1 )
ADilWTag    = RealVar( 'ADilWTag',    Title = 'dilution/wrong tag asymmetry', Observable = False, Value = 0 )


_AProd   = RealVar(    'AProd',    Title = 'production asymmetry',         Observable = False, Value = 0 )
_ATagEff = RealVar(    'ATagEff',  Title = 'tagging efficiency asymmetry', Observable = False, Value = 0 )
_ANorm   = Product(    'ANorm',   [minus,CP.C],  Title = 'normalization asymmetry' )

args.update( { 'avgCEven' : FormulaVar( 'avgCEven', '1. + @0*@1 + @0*@2 + @1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average even coefficients')
             , 'avgCOdd'  : FormulaVar( 'avgCOdd',     '@0 + @1 + @2 + @0*@1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average odd coefficients')
             } )


# construct amplitudes with polar parameters -- these are the 'externally visible' (expected) parameters -- 2*4=8 terms for 4 amplitudes
class Polar2_Amplitude :
    def __init__(self,name, r2, arg, CP ) :
        self.name = name
        self.Re = FormulaVar('Re%s'%name, 'sqrt(@0) * cos(@1)', [r2,arg], Title = 'Re(%s)'% name )
        self.Im = FormulaVar('Im%s'%name, 'sqrt(@0) * sin(@1)', [r2,arg], Title = 'Im(%s)'% name )
        self.CP = CP
    def __str__(self) : return self.name

# polar transversity amplitudes -- this is 'internal only'
_A0Mag2    = RealVar('A0Mag2',    Title = '|A0|^2',     Observable = False, Value = 0.556, MinMax = (0., 1.))
_A0Ph      = RealVar('delta0',    Title = 'delta_0',    Observable = False, Value = 0. )
_AperpMag2 = RealVar('AperpMag2', Title = '|A_perp|^2', Observable = False, Value = 0.233, MinMax = ( 0., 1.))
_AperpPh   = RealVar('deltaPerp', Title = 'delta_perp', Observable = False, Value = pi-2.91, MinMax=( -2. * pi, 2. * pi))
_AparMag2  = FormulaVar('AparMag2', '1. - @0 - @1', [_A0Mag2, _AperpMag2],  Title = '|A_par|^2' )
_AparPh    = RealVar('deltaPar',  Title = 'delta_par',  Observable = False, Value = 2.93, MinMax = ( -2. * pi, 2. * pi))
_ASMag2    = RealVar('ASMag2',    Title = '|A_S|^2',    Observable = False, Value = 0.05, MinMax=( 0., 1.))
_ASPh      = RealVar('deltaS',    Title = 'delta_S',    Observable = False, Value = 2.2, MinMax=( -2. * pi, 2. * pi))


    
# using transversity amplitudes and helicity angles
from math import sqrt
_ba = lambda  name,args : Addition(name, [ AngleBasis((cpsiAng,cthetaAng,phiAng) , *a) for a in args ] )

# the following two tables, define 'everything'....####################################################
angFuncs = { ('A0',   'A0')    :  ('Re', _ba('J_0020x0020_0', [(0, 0, 0,  0,  4.             )
                                                              ,(0, 0, 2,  0, -sqrt( 16. / 5.))
                                                              ,(2, 0, 0,  0,  8.             )
                                                              ,(2, 0, 2,  0, -sqrt( 64. / 5.))]))
           , ('Apar', 'Apar')  :  ('Re', _ba('J_22x002022_0', [(2, 2, 0,  0,  2.             )
                                                              ,(2, 2, 2,  0,  sqrt(  1. / 5.))
                                                              ,(2, 2, 2,  2, -sqrt(  3. / 5.))]))
           , ('Aperp','Aperp') :  ('Re', _ba('J_22x002022_1', [(2, 2, 0,  0,  2.             )
                                                              ,(2, 2, 2,  0,  sqrt(  1. / 5.))
                                                              ,(2, 2, 2,  2,  sqrt(  3. / 5.))]) )
           , ('A0',   'Apar')  :  ('Re', _ba('J_21x21_0',     [(2, 1, 2,  1,  sqrt( 24. / 5.))]))
           , ('A0',   'Aperp') :  ('Im', _ba('J_21x2m1_0',    [(2, 1, 2, -1, -sqrt( 24. / 5.))]))
           , ('Apar', 'Aperp') :  ('Im', _ba('J_22x2m2_0',    [(2, 2, 2, -2,  sqrt( 12. / 5.))])) 
           , ('AS',   'AS')    :  ('Re', _ba('J_00x0020_0',   [(0, 0, 0,  0,  4.             )
                                                              ,(0, 0, 2,  0, -sqrt( 16. / 5.))]))
           , ('A0',   'AS')    :  ('Re', _ba('J_10x0020_0',   [(1, 0, 0,  0,  sqrt(192.     ))
                                                              ,(1, 0, 2,  0, -sqrt(192. / 5.))]))
           , ('Apar', 'AS')    :  ('Re', _ba('J_11x21_0',     [(1, 1, 2,  1,  sqrt( 72. / 5.))]))
           , ('Aperp','AS')    :  ('Im', _ba('J_11x2m1_0',    [(1, 1, 2, -1,  sqrt( 72. / 5.))]))
           }

# in python 2.7 and later, could use collections.OrderedDict so we can traverse without having to think...
Amplitudes = { 'A0'    : Polar2_Amplitude( 'A0',    _A0Mag2,    _A0Ph,    +1 )
             , 'Apar'  : Polar2_Amplitude( 'Apar',  _AparMag2,  _AparPh,  +1 )
             , 'Aperp' : Polar2_Amplitude( 'Aperp', _AperpMag2, _AperpPh, -1 )
             , 'AS'    : Polar2_Amplitude( 'AS',    _ASMag2,    _ASPh,    -1 )
             }
#######################################################################################################

def combine( name, f, A, CPparams, i, j) :
    # define functions which return Re(Conj(Ai) Aj), Im( Conj(Ai) Aj)
    Real   = lambda ai, aj  : FormulaVar('Re_c_%s_%s'%(ai,aj),'@0*@2+@1*@3',[ai.Re,ai.Im,aj.Re,aj.Im])
    Imag   = lambda ai, aj  : FormulaVar('Im_c_%s_%s'%(ai,aj),'@0*@3-@1*@2',[ai.Re,ai.Im,aj.Re,aj.Im])
    # define functions which return the coefficients that define the time-dependence...
    coef = { 'cosh' : lambda ai,aj,CP : ( one  if ai.CP == aj.CP else CP.C, None )
           , 'cos'  : lambda ai,aj,CP : ( CP.C if ai.CP == aj.CP else one , None )
           , 'sinh' : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('PR_%s_%s_D'%(ai,aj), [ minus if ai.CP > 0     else plus,  CP.D ])
                                        , None if ai.CP == aj.CP else Product('PI_%s_%s_S'%(ai,aj), [ plus  if ai.CP > aj.CP else minus, CP.S ]) )
           , 'sin'  : lambda ai,aj,CP : ( None if ai.CP != aj.CP else Product('PR_%s_%s_S'%(ai,aj), [ minus if ai.CP > 0     else plus,  CP.S ])
                                        , None if ai.CP == aj.CP else Product('PI_%s_%s_D'%(ai,aj), [ minus if ai.CP > aj.CP else plus,  CP.D ]) )
        }
    (re,im) = coef[name](A[i],A[j],CPparams)
    # for now, coefficients are either real, or imaginary, but not both... (not true in general, but I'm lazy today ;-)
    assert not ( re and im )
    if f[(i,j)][0] == 'Re' : 
        if re : return Product('rr_%s_%s_%s'%(name,A[i],A[j]), [        Real(A[i],A[j]), re, f[(i,j)][1] ] ) # Re( z) = +Re(z)
        if im : return Product('ri_%s_%s_%s'%(name,A[i],A[j]), [ minus, Imag(A[i],A[j]), im, f[(i,j)][1] ] ) # Re(iz) = -Im(z)
    else  :                                            
        if re : return Product('ir_%s_%s_%s'%(name,A[i],A[j]), [        Imag(A[i],A[j]), re, f[(i,j)][1] ] ) # Im( z) = +Im(z)
        if im : return Product('ri_%s_%s_%s'%(name,A[i],A[j]), [        Real(A[i],A[j]), im, f[(i,j)][1] ] ) # Im(iz) = +Re(z)

for name in [ 'cosh', 'sinh', 'cos', 'sin' ] :
    from itertools import combinations_with_replacement as cwr
    # NOTE: 'Amplitudes'  must be traversed 'in order' : A0, Apar, Aperp, AS -- so we cannot use Amplitudes.keys() out of the box...
    args[ '%sCoef' % name ] = Addition( 'a_%s'% name, [ combine(name,angFuncs,Amplitudes,CP,i,j) for (i,j) in cwr( ['A0','Apar','Aperp','AS'], 2 ) ] )

# build PDF
args.update(  { 'time'     : t
              , 'iTag'     : iTag
              , 'dilution' : tagDilution
              , 'ADilWTag' : ADilWTag
              } )
pdf = BTagDecay( 'sig_pdf', args )

data = pdf.generate( [ cpsiAng,cthetaAng,phiAng,t,iTag, ] , 10000 )
pdf.fitTo(data)

#ws.ws().Print('V')
