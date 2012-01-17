from math import sqrt, pi
from RooFitWrappers import *

indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')

# define observables
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5200, 5550), nBins =  48
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  16 )

t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.3,14.), nBins =  54 )
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.0, 0.15),  nBins =  50 )
eta_os  = RealVar('tagomega_os',      Title = 'estimated mistag OS',          Observable = True, MinMax = (0,0.50001),  nBins =  25)
#The peak at 0.5 seems to be shifted to -2 in the SS eta!
eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (-2.0001,0.50001),  nBins =  25)
iTag_os = Category( 'tagdecision_os', Title = 'initial state flavour tag OS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
#The peak at 0 seems to be shifted to -1000 in the SS tagdecision
iTag_ss = Category( 'tagdecision_ss', Title = 'initial state flavour tag SS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : -1000 } )
sel  = Category( 'sel',            Title = 'selection',                 Observable = True, States = { 'good': +1 } )
triggerdec = Category( 'triggerDecision',            Title = 'triggerdec',                 Observable = True, States = { 'triggered': +1 } )
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )
bkgcat = Category( 'bkgcat',            Title = 'bkgcat',                 Observable = True, States = { 'bkgcat0': 0, 'bkgcat10': 10 } )

####################
### Set time cut ###
####################

tcut = 0.3
t['MinMax'] = (tcut,14.)

from P2VVGeneralUtils import readData
#Read MC data
MCdata = readData('/data/bfys/dveijk/MC/2012/Bs2JpsiPhi_MC11a_ntupleB_for_fitting_20120109.root'
                  , dataSetName = 'DecayTree'
                  , NTuple = True
                  , observables = [ t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], iTag_os,eta_os, triggerdec,sel,bkgcat,mphi]
                  )

print 'Number of events', MCdata.numEntries()

#Time Resolution Model for MC
from P2VVParameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
tres = TimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tres.setConstant('.*')

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                       , deltaGamma = dict( Name = 'dGamma'
                                                            , Value = 0.060
                                                            )
                                       , deltaM = dict( Value = 17.8, MinMax = (16.5,18.5), Constant = False) 
                                       )

# define tagging parameter 
from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
tagging = TaggingParams( estWTag = eta_os ) # Constant = False, Constrain = True ) TODO!!!

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam(  phiCP      = dict( Name = 'phi_s'
                                              , Value = -0.04
                                              , MinMax = (-pi,pi))
                         , lambdaCPSq = dict( Value = 1., Constant = True )
                        )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
from P2VVParameterizations.DecayAmplitudes import JpsiVPolar_AmplitudeSet
amplitudes = JpsiVPolar_AmplitudeSet(  A0Mag2 = 0.60, A0Phase = 0
                                     , AperpMag2 = 0.16, AperpPhase = -0.17 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                     , AparPhase = 2.5
                                     , ASMag2 = dict( Value = 0.0, Constant = True )
                                     , ASPhase = dict( Value = 0.0, Constant = True )
                                     , PWaveNorm = False
                                    )

# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions
                                                  , amplitudes
                                                  , CP
                                                  , Product('tag',(iTag_os,tagging['dilution']))
                                                  , ['A0','Apar','Aperp','AS'] ) 


from RooFitWrappers import BDecay
MC_sig_t_angles = BDecay( Name      = 'MC_sig_t_angles'
                          , time      = t
                          , dm        = lifetimeParams['deltaM'] 
                          , tau       = lifetimeParams['MeanLifetime']
                          , dGamma    = lifetimeParams['deltaGamma'] 
                          , resolutionModel = tres.model()
                          , coshCoef  = basisCoefficients['cosh']
                          , cosCoef   = basisCoefficients['cos']
                          , sinhCoef  = basisCoefficients['sinh']
                          , sinCoef   = basisCoefficients['sin']
#                          , ConditionalObservables = ( eta_os, )
#                          , ConditionalObservables = ( eta_os, iTag_os, )
                          )

#####################################
### Angular acceptance correction ###
#####################################
MCpdf = MC_sig_t_angles

print 'Number of MC events', MCdata.numEntries()
allObs = MCpdf.getObservables( MCdata.get() )
print 'MCobservables:', [ i.GetName() for i in allObs ]
o = MCpdf.getObservables(MCdata.get() )

from P2VVGeneralUtils import RealMomentsBuilder
nset = angles.angles.values()

canomoms = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, MCpdf,nset) for v in angles.functions.itervalues() for i in v if i ) )
canomoms.compute(MCdata)
canomoms.Print(Scales = [1./(16.*sqrt(pi)),1./(16.*sqrt(pi)),1./(16.*sqrt(pi))])
canomoms.write('canomoms_tcut_%s.txt'%(str(tcut)))

from itertools import chain
#momindices = indices(3,3)
momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))

eff = RealMomentsBuilder()
eff.appendPYList( angles.angles, momindices, PDF = MCpdf, NormSet = nset)
eff.compute(MCdata)
eff.Print()
eff.write('effmoments_tcut_%s.txt'%(str(tcut)))

def calc_moments_from_series( c ) :
    return { 'Re_ang_A0_A0'       :   4*( c[(0,0,0)]+2*c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]+2*c[(2,2,0)]/5) - sqrt(3./20)*(c[(0,2,2)]+2*c[(2,2,2)]/5)  )
           , 'Re_ang_Apar_Apar'   :   4*( c[(0,0,0)]-  c[(2,0,0)]/5 + sqrt(1./20)*( c[(0,2,0)]-  c[(2,2,0)]/5) + sqrt(3./20)*(c[(0,2,2)]-  c[(2,2,2)]/5)  )
           , 'Re_ang_Aperp_Aperp' :   4*( c[(0,0,0)]-  c[(2,0,0)]/5 - sqrt(1./ 5)*( c[(0,2,0)]-  c[(2,2,0)]/5 ) )
           , 'Im_ang_Apar_Aperp'  :   4*sqrt(3./5.)*( c[(0,2,-1)] - c[(2,2,-1)]/5 )
           , 'Re_ang_A0_Apar'     :   4*sqrt(6./5.)* 3*pi/32 *( c[(1,2,-2)] - c[(3,2,-2)]/4 - 5*c[(5,2,-2)]/128  - 7*c[(7,2,-2)]/512 - 105*c[(9,2,-2)]/16384)
           , 'Im_ang_A0_Aperp'    :  -4*sqrt(6./5.)* 3*pi/32 *( c[(1,2, 1)] - c[(3,2, 1)]/4 - 5*c[(5,2, 1)]/128  - 7*c[(7,2, 1)]/512 - 105*c[(9,2, 1)]/16384)
           , 'Re_ang_AS_AS'       :   2*(2*c[(0,0,0)]+sqrt(1./5)*c[(0,2,0)]-sqrt(3./5)*c[(0,2,2)])
           , 'Re_ang_Apar_AS' :   12*sqrt(2./5.)*pi/8 *( c[(0,2,-2)] - c[(2,2,-2)]/8 - c[(4,2,-2)]/64 - 5*pi*c[(6,2,-2)]/1024 -35*pi*c[(8,2,-2)]/16384)
           , 'Im_ang_Aperp_AS'     :  -12*sqrt(2./5.)*pi/8 *( c[(0,2,1)] - c[(2,2,1)]/8 - c[(4,2,1)]/64 - 5*pi*c[(6,2,1)]/1024 -35*pi*c[(8,2,1)]/16384)
           , 'Re_ang_A0_AS'       :   (2./3)*(4*sqrt(3)*c[(1,0,0)]+2*sqrt(3./5)*c[(1,2,0)]-6*sqrt(1./5)*c[(1,2,2)])
           }

def compare_methods (canomoments,moments): #compute the 'canonical' moments given the Fourier series
    c = dict()
    for m in moments :
        c[ ( moments.basisFuncIndices()[m][0],moments.basisFuncIndices()[m][1],moments.basisFuncIndices()[m][2] ) ] = moments.coefficients()[m][0]
    xi_c = calc_moments_from_series( c )
    
    for name in xi_c.iterkeys() :
        print '%s : direct moment: %s ;  from Fourier series: %s ; ratio = %s ' % ( name, canomoments.coefficients()[name][0], xi_c[name], canomoments.coefficients()[name][0]/xi_c[name])    

compare_methods(canomoms,eff)
