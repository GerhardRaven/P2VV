from RooFitWrappers import *
from ROOT import gSystem
gSystem.Load('libP2VV.so')
ws = RooObject()
ws.setWorkspace( RooWorkspace('myworkspace') )

from math import pi
cpsiAng   = RealVar('helcthetaK', Title = 'cosine of kaon polarization angle',   Observable = True,  MinMax=(-1., 1.))
cthetaAng = RealVar('helcthetaL', Title = 'cosine of lepton polarization angle', Observable = True,  MinMax=(-1., 1.))
phiAng    = RealVar('helphi',     Title = 'angle between decay planes',          Observable = True,  MinMax=(-pi, pi))
t         = RealVar('t',          Title = 'decay time', Unit='ps',               Observable = True,  MinMax=(-5,14)  )

zero  = ConstVar('zero',   Value =  0  )
one   = ConstVar('one',    Value =  1  )
minus = ConstVar('minus',  Value = -1  )
half  = ConstVar('half',   Value =  0.5)

unbiased = Category( 'unbiased',  Title = 'Unbiased trigger?',         Observable = True, States = { 'yes': 1 ,   'no' : 0 } )
decay    = Category( 'decaytype', Title = 'J/psiX decay type',         Observable = True, States = { 'JpsiKplus' : 10, 'JpsiKmin' : 11, 'JpsiKstar0' : 20, 'JpsiKstar0bar': 21, 'Jpsiphi': 40 } )
iTag     = Category( 'iTag',      Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1 } )

from ROOT import RooTruthModel as TruthModel
resolution = ResolutionModel( 'resModel', Type = TruthModel, Observables = [ t ] )

dm        =  RealVar( 'dm',        Title = 'delta m',       Unit = 'ps^{-1}', Observable = False, Value = 17.8 )
BMeanLife =  RealVar( 't_sig_tau', Title = 'mean lifetime', Unit = 'ps',      Observable = False, Value =  1.5, MinMax = (1.3, 1.8) )
dGamma    =  RealVar( 'dGamma',    Title = 'dGamma',        Unit = 'ps^{-1}', Observable = False, Value = 0.05, MinMax = (-0.3,0.3) )

# TODO: package this in a seperate class, which gives access to CCP, DCP, and SCP, so we can hide the 'internal' parameterization
# TODO: explicitly export externally visible symbols .... For now, interal stuff is prefixed by an _
_phiCP      = RealVar( 'phiCP',    Title = 'CP violation param. phi_s',     Observable = False, Value = -0.2,  MinMax = (-2*pi,2*pi) )
_lambdaCPSq = RealVar( 'lambda^2', Title = 'CP violation param |lambda|^2', Observable = False, Value = 1)

class CPParam :
    def __init__(self,**kwargs) :
        self.C = kwargs.pop('C')
        self.D = kwargs.pop('D')
        self.S = kwargs.pop('S')

CP = CPParam( C = FormulaVar('C', '(1.-@0)/(1.+@0)',                   [ _lambdaCPSq ] )
            , D = FormulaVar('D', ' 2 * sqrt(@0) * cos(@1) / (1.+@0)', [ _lambdaCPSq, _phiCP ] )
            , S = FormulaVar('S', '-2 * sqrt(@0) * sin(@1) / (1.+@0)', [ _lambdaCPSq, _phiCP ] )
            )

tagDilution = RealVar( 'tagDilution', Title = 'Average Tagging Dilution', Observable = False, Value = 1 )
ADilWTag    = RealVar( 'ADilWTag',    Title = 'dilution/wrong tag asymmetry', Observable = False, Value = 0 )


_AProd   = RealVar(    'AProd',    Title = 'production asymmetry', Observable = False, Value = 0)
_ANorm   = FormulaVar( 'ANorm',    '-@0', [CP.C], Title = 'normalization asymmetry' )
_ATagEff = RealVar(    'ATagEff',  Title = 'tagging efficiency asymmetry', Observable = False, Value=0.)
avgCEven = FormulaVar( 'avgCEven', '1. + @0*@1 + @0*@2 + @1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average even coefficients')
avgCOdd  = FormulaVar( 'avgCOdd',     '@0 + @1 + @2 + @0*@1*@2', [_AProd, _ANorm, _ATagEff], Title = 'CP average odd coefficients')


# polar transversity amplitudes -- this is 'internal only'
_A0Mag2    = RealVar('A0Mag2',    Title = '|A0|^2',     Observable = False, Value = 0.556, MinMax = (0., 1.))
_A0Ph      = RealVar('delta0',    Title = 'delta_0',    Observable = False, Value = 0. )
_AperpMag2 = RealVar('AperpMag2', Title = '|A_perp|^2', Observable = False, Value = 0.233, MinMax = ( 0., 1.))
_AperpPh   = RealVar('deltaPerp', Title = 'delta_perp', Observable = False, Value = pi-2.91, MinMax=( -2. * pi, 2. * pi))
_AparMag2  = FormulaVar('AparMag2', '1. - @0 - @1', [_A0Mag2, _AperpMag2],  Title = '|A_par|^2' )
_AparPh    = RealVar('deltaPar',  Title = 'delta_par',  Observable = False, Value = 2.93, MinMax = ( -2. * pi, 2. * pi))
_ASMag2    = RealVar('ASMag2',    Title = '|A_S|^2',    Observable = False, Value = 0.05, MinMax=( 0., 1.))
_ASPh      = RealVar('deltaS',    Title = 'delta_S',    Observable = False, Value = 2.2, MinMax=( -2. * pi, 2. * pi))

# construct cartesian amplitudes with polar parameters -- these are the 'externally visible' (expected) parameters -- 2*4=8 terms for 4 amplitudes
class Carth_Amplitude :
    def __init__(self,x,y) :
        self.Re = x
        self.Im = y

Amplitudes = { 'A0'    : Carth_Amplitude( FormulaVar('ReA0',   'sqrt(@0) * cos(@1)', [_A0Mag2,    _A0Ph],    Title = 'Re(A_0)'     )
                                        , FormulaVar('ImA0',   'sqrt(@0) * sin(@1)', [_A0Mag2,    _A0Ph],    Title = 'Im(A_0)'     ))
             , 'Apar'  : Carth_Amplitude( FormulaVar('ReApar', 'sqrt(@0) * cos(@1)', [_AparMag2,  _AparPh],  Title = 'Re(A_par)'   )
                                        , FormulaVar('ImApar', 'sqrt(@0) * sin(@1)', [_AparMag2,  _AparPh],  Title = 'Im(A_par)'   ))
             , 'Aperp' : Carth_Amplitude( FormulaVar('ReAperp','sqrt(@0) * cos(@1)', [_AperpMag2, _AperpPh], Title = 'Re(A_perp)'  )
                                        , FormulaVar('ImAperp','sqrt(@0) * sin(@1)', [_AperpMag2, _AperpPh], Title = 'Im(A_perp)'  ))
             , 'AS'    : Carth_Amplitude( FormulaVar('ReAS',   'sqrt(@0) * cos(@1)', [_ASMag2,    _ASPh],    Title = 'Re(A_S)'     )
                                        , FormulaVar('ImAS',   'sqrt(@0) * sin(@1)', [_ASMag2,    _ASPh],    Title = 'Im(A_S)'     ))
             }
    
# these are the angular terms: 4x(4+1)/2 = 10
J_0020x0020_0 = { 'cosh' : FormulaVar('J_0020x0020_0_cosh',   '@0 * @0 + @1 * @1',       [Amplitudes['A0'].Re, Amplitudes['A0'].Im      ], Title = 'J_0020x0020_0 cosh coefficient') # +|A_0|^2 * 1
                , 'cos'  : FormulaVar('J_0020x0020_0_cos',   '(@0 * @0 + @1 * @1) * @2', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, CP.C], Title = 'J_0020x0020_0 cos  coefficient') # +|A_0|^2 * C
                , 'sinh' : FormulaVar('J_0020x0020_0_sinh', '-(@0 * @0 + @1 * @1) * @2', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, CP.D], Title = 'J_0020x0020_0 sinh coefficient') # -|A_0|^2 * D
                , 'sin'  : FormulaVar('J_0020x0020_0_sin',  '-(@0 * @0 + @1 * @1) * @2', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, CP.S], Title = 'J_0020x0020_0 sin  coefficient') # -|A_0|^2 * S
                }
J_22x002022_0 = { 'cosh' : FormulaVar('J_22x002022_0_cosh',   '@0 * @0 + @1 * @1',       [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im     ], Title = 'J_22x002022_0 cosh coefficient') # +|A_par|^2 * 1
                , 'cos'  : FormulaVar('J_22x002022_0_cos',   '(@0 * @0 + @1 * @1) * @2', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, CP.C], Title = 'J_22x002022_0 cos  coefficient') # +|A_par|^2 * C
                , 'sinh' : FormulaVar('J_22x002022_0_sinh', '-(@0 * @0 + @1 * @1) * @2', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, CP.D], Title = 'J_22x002022_0 sinh coefficient') # -|A_par|^2 * D
                , 'sin'  : FormulaVar('J_22x002022_0_sin',  '-(@0 * @0 + @1 * @1) * @2', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, CP.S], Title = 'J_22x002022_0 sin  coefficient') # -|A_par|^2 * S
                }
J_22x002022_1 = { 'cosh' : FormulaVar('J_22x002022_1_cosh',   '@0 * @0 + @1 * @1',       [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im     ], Title = 'J_22x002022_1 cosh coefficient') # +|A_perp|^2 * 1
                , 'cos'  : FormulaVar('J_22x002022_1_cos',   '(@0 * @0 + @1 * @1) * @2', [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.C], Title = 'J_22x002022_1 cos  coefficient') # +|A_perp|^2 * C
                , 'sinh' : FormulaVar('J_22x002022_1_sinh',  '(@0 * @0 + @1 * @1) * @2', [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.D], Title = 'J_22x002022_1 sinh coefficient') # +|A_perp|^2 * D
                , 'sin'  : FormulaVar('J_22x002022_1_sin',   '(@0 * @0 + @1 * @1) * @2', [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.S], Title = 'J_22x002022_1 sin  coefficient') # +|A_perp|^2 * S
                }
J_21x21_0     = { 'cosh' : FormulaVar('J_21x21_0_cosh',       '@0 * @2 + @1 * @3',       [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Apar'].Re, Amplitudes['Apar'].Im      ], Title = 'J_21x21_0 cosh coefficient') # +Re(A_0* A_par) * 1
                , 'cos'  : FormulaVar('J_21x21_0_cos',       '(@0 * @2 + @1 * @3) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, CP.C], Title = 'J_21x21_0 cos coefficient' ) # +Re(A_0* A_par) * C
                , 'sinh' : FormulaVar('J_21x21_0_sinh',     '-(@0 * @2 + @1 * @3) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, CP.D], Title = 'J_21x21_0 sinh coefficient') # -Re(A_0* A_par) * D
                , 'sin'  : FormulaVar('J_21x21_0_sin',      '-(@0 * @2 + @1 * @3) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, CP.S], Title = 'J_21x21_0 sin coefficient' ) # -Re(A_0* A_par) * S
                }
J_21x2m1      = { 'cosh' : FormulaVar('J_21x2m1_0_cosh',     '(@0 * @3 - @1 * @2) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.C], Title = 'J_21x2m1_0 cosh coefficient') # +Im(A_0* A_perp) * C
                , 'cos'  : FormulaVar('J_21x2m1_0_cos',       '@0 * @3 - @1 * @2',       [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im      ], Title = 'J_21x2m1_0 cos coefficient' ) # +Im(A_0* A_perp) * 1
                , 'sinh' : FormulaVar('J_21x2m1_0_sinh',     '(@0 * @2 + @1 * @3) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.S], Title = 'J_21x2m1_0 sinh coefficient') # +Re(A_0* A_perp) * S
                , 'sin'  : FormulaVar('J_21x2m1_0_sin',     '-(@0 * @2 + @1 * @3) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.D], Title = 'J_21x2m1_0 sin coefficient' ) # -Re(A_0* A_perp) * D
                }
J_22x2m2_0    = { 'cosh' : FormulaVar('J_22x2m2_0_cosh',     '(@0 * @3 - @1 * @2) * @4', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.C],   Title = 'J_22x2m2_0 cosh coefficient') # +Im(A_par* A_perp) * C
                , 'cos'  : FormulaVar('J_22x2m2_0_cos',       '@0 * @3 - @1 * @2',       [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im      ],   Title = 'J_22x2m2_0 cos coefficient' ) # +Im(A_par* A_perp) * 1
                , 'sinh' : FormulaVar('J_22x2m2_0_sinh',     '(@0 * @2 + @1 * @3) * @4', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.S],   Title = 'J_22x2m2_0 sinh coefficient') # +Re(A_par* A_perp) * S
                , 'sin'  : FormulaVar('J_22x2m2_0_sin',     '-(@0 * @2 + @1 * @3) * @4', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, CP.D],   Title = 'J_22x2m2_0 sin coefficient' ) # -Re(A_par* A_perp) * D
                }
J_00x0020_0   = { 'cosh' : FormulaVar('J_00x0020_0_cosh',     '@0 * @0 + @1 * @1',       [Amplitudes['AS'].Re, Amplitudes['AS'].Im     ], Title = 'J_00x0020_0 cosh coefficient') # +|A_S|^2 * 1
                , 'cos'  : FormulaVar('J_00x0020_0_cos',     '(@0 * @0 + @1 * @1) * @2', [Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.C], Title = 'J_00x0020_0 cos coefficient' ) # +|A_S|^2 * C
                , 'sinh' : FormulaVar('J_00x0020_0_sinh',    '(@0 * @0 + @1 * @1) * @2', [Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.D], Title = 'J_00x0020_0 sinh coefficient') # +|A_S|^2 * D
                , 'sin'  : FormulaVar('J_00x0020_0_sin',     '(@0 * @0 + @1 * @1) * @2', [Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.S], Title = 'J_00x0020_0 sin coefficient' ) # +|A_S|^2 * S
                }
J_10x0020_0   = { 'cosh' : FormulaVar('J_10x0020_0_cosh',    '(@0 * @2 + @1 * @3) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.C], Title = 'J_10x0020_0 cosh coefficient' ) # +Re(A_0* A_S) * C
                , 'cos'  : FormulaVar('J_10x0020_0_cos',      '@0 * @2 + @1 * @3',       [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im     ], Title = 'J_10x0020_0 cos coefficient'  ) # +Re(A_0* A_S) * 1
                , 'sinh' : FormulaVar('J_10x0020_0_sinh',   '-(@0 * @3 - @1 * @2) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.S], Title = 'J_10x0020_0 sinh coefficient' ) # -Im(A_0* A_S) * S
                , 'sin'  : FormulaVar('J_10x0020_0_sin',     '(@0 * @3 - @1 * @2) * @4', [Amplitudes['A0'].Re, Amplitudes['A0'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.D], Title = 'J_10x0020_0 sin coefficient'  ) # +Im(A_0* A_S) * D
                }
J_11x21_0     = { 'cosh' : FormulaVar('J_11x21_0_cosh',      '(@0 * @2 + @1 * @3) * @4', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.C],  Title = 'J_11x21_0 cosh coefficient')  # +Re(A_par* A_S) * C
                , 'cos'  : FormulaVar('J_11x21_0_cos',        '@0 * @2 + @1 * @3',       [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im     ],  Title = 'J_11x21_0 cos coefficient')   # +Re(A_par* A_S) * 1
                , 'sinh' : FormulaVar('J_11x21_0_sinh',     '-(@0 * @3 - @1 * @2) * @4', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.S],  Title = 'J_11x21_0 sinh coefficient')  # -Im(A_par* A_S) * S
                , 'sin'  : FormulaVar('J_11x21_0_sin',       '(@0 * @3 - @1 * @2) * @4', [Amplitudes['Apar'].Re, Amplitudes['Apar'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.D],  Title = 'J_11x21_0 sin coefficient')   # +Im(A_par* A_S) * D
                }
J_11x2m1_0    = { 'cosh' : FormulaVar('J_11x2m1_0_cosh',      '@0 * @3 - @1 * @2',       [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im     ],  Title = 'J_11x2m1_0 cosh coefficient' ) # +Im(A_perp* A_S) * 1
                , 'cos'  : FormulaVar('J_11x2m1_0_cos',      '(@0 * @3 - @1 * @2) * @4', [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.C],  Title = 'J_11x2m1_0 cos coefficient'  ) # +Im(A_perp* A_S) * C
                , 'sinh' : FormulaVar('J_11x2m1_0_sinh',     '(@0 * @3 - @1 * @2) * @4', [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.D],  Title = 'J_11x2m1_0 sinh coefficient' ) # +Im(A_perp* A_S) * D
                , 'sin'  : FormulaVar('J_11x2m1_0_sin',      '(@0 * @3 - @1 * @2) * @4', [Amplitudes['Aperp'].Re, Amplitudes['Aperp'].Im, Amplitudes['AS'].Re, Amplitudes['AS'].Im, CP.S],  Title = 'J_11x2m1_0 sin coefficient'  ) # +Im(A_perp* A_S) * S
                }

# TODO: put the above definitions 'inline' in dict 'coef' and push that into a function to hide the internals...
coef =  { ('A0',   'A0')    : J_0020x0020_0
        , ('A0',   'Apar')  : J_21x21_0 
        , ('A0',   'Aperp') : J_21x2m1
        , ('A0',   'AS')    : J_10x0020_0
        , ('Apar', 'Apar')  : J_22x002022_0
        , ('Apar', 'Aperp') : J_22x2m2_0
        , ('Apar', 'AS')    : J_11x21_0
        , ('Aperp','Aperp') : J_22x002022_1
        , ('Aperp','AS')    : J_11x2m1_0
        , ('AS',   'AS')    : J_00x0020_0
        }

from itertools import combinations_with_replacement
# TODO: 'Amplitudes'  must be traversed 'in order' : A0, Apar, Aperp, AS -- so we cannot use Amplitudes.keys() out of the box...
for (i,j) in combinations_with_replacement( ['A0','Apar','Aperp','AS'], 2 ) :
    print (i,j)
    zz = coef[ (i,j) ]

z = AngleBasis( (cpsiAng,cthetaAng,phiAng), 0,0,0,0,4.)


# using helicity angles
#angFuncs.append(('J_0020x0020_0', [(0, 0, 0,  0,  4.             ),
#                                   (0, 0, 2,  0, -sqrt( 16. / 5.)),
#                                   (2, 0, 0,  0,  8.             ),
#                                   (2, 0, 2,  0, -sqrt( 64. / 5.))]))
#angFuncs.append(('J_22x002022_0', [(2, 2, 0,  0,  2.             ),
#                                   (2, 2, 2,  0,  sqrt(  1. / 5.)),
#                                   (2, 2, 2,  2, -sqrt(  3. / 5.))]))
#angFuncs.append(('J_22x002022_1', [(2, 2, 0,  0,  2.             ),
#                                   (2, 2, 2,  0,  sqrt(  1. / 5.)),
#                                   (2, 2, 2,  2,  sqrt(  3. / 5.))]))
#angFuncs.append(('J_21x21_0',     [(2, 1, 2,  1,  sqrt( 24. / 5.))]))
#angFuncs.append(('J_21x2m1_0',    [(2, 1, 2, -1, -sqrt( 24. / 5.))]))
#angFuncs.append(('J_22x2m2_0',    [(2, 2, 2, -2,  sqrt( 12. / 5.))]))
#
#angFuncs.append(('J_00x0020_0', [(0, 0, 0,  0,  4.             ),
#                                 (0, 0, 2,  0, -sqrt( 16. / 5.))]))
#angFuncs.append(('J_10x0020_0', [(1, 0, 0,  0,  sqrt(192.     )),
#                                 (1, 0, 2,  0, -sqrt(192. / 5.))]))
#angFuncs.append(('J_11x21_0',   [(1, 1, 2,  1,  sqrt( 72. / 5.))]))
#angFuncs.append(('J_11x2m1_0',  [(1, 1, 2, -1,  sqrt( 72. / 5.))]))
#


ws.ws().Print('V')
