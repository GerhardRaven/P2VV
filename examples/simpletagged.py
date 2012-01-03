from math import sqrt, pi
from RooFitWrappers import *

obj  = RooObject( workspace = 'workspace')
indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )

from P2VVGeneralUtils import numCPU
from ROOTDecorators import  ROOTversion as Rv

fitOpts = dict( NumCPU = numCPU() 
              , Timer=1
              , Save = True
              , Verbose = False
              , Optimize = True if Rv[1]<32 else 0 # do NOT optimize in 5.32 or later... ( Optimize = 1 only works on a single CPU, 2 doesn't work at all )
              , Minimizer = ('Minuit2','minimize')
              )

# define observables
t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.3, 14),    nBins =  54 )
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.0, 0.15),  nBins =  50 )
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5259, 5451), nBins =  48
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3025, 3169), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1012, 1028), nBins =  16 )
iTag = Category( 'tagdecision_os', Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
sel  = Category( 'sel', Title = 'selection', Observable = True, States = { 'good': +1 } )
eta  = RealVar('tagomega_os',      Title = 'estimated mistag',          Observable = True, MinMax = (0,0.50001),  nBins =  50)
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )

from P2VVGeneralUtils import readData
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
               , 'DecayTree'
               , True
               , [ iTag,eta, mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], st, sel ]
               )
tag  = RealVar('tagdilution', Title = 'estimated dilution', Observable = True, Formula = ( "@0*(1-2*@1)",[iTag,eta], data ), MinMax=(-1,1),nBins =  51)

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5368, MinMax = (5363,5372) ) )
cmb_m = Background_BMass( Name = 'cmb_m', mass = m, m_bkg_exp  = dict( Name = 'm_cmb_exp' ) )

# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tres.setConstant('.*')

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.65, deltaGamma = 0.10, deltaM = dict( Value = 17.8, MinMax = (16.5,18.5), Constant = True) )

from P2VVParameterizations.FlavourTagging import Trivial_Dilution
taggingParams = Trivial_Dilution( Dilution = tag ) # TODO: add calibration: Dilution = FormulaVar("calibrated_tag","...",{p0,p1,etaAverage,tag} )
# TODO: add external constraint terms for p0 and p1...

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP      = dict( Name = 'phi_s', Value = -0.04, MinMax = (-pi,pi), Constant = False )
                        , lambdaCPSq = ConstVar( Name = 'one', Value = 1 ) 
                        )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011( A0Mag2 = 0.50, A0Phase = 0
                                    , AperpMag2 = 0.25, AperpPhase = dict( Value = -0.17 ) # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                    , AparPhase = 2.5
                                    , ASMag2 = dict( Value = 0, Constant = True ) , ASPhase = dict( Value = 0, Constant = True ) 
                                    )
# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes,CP, tag, ['A0','Apar','Aperp','AS'] ) 

from RooFitWrappers import BDecay
sig_t_angles = BDecay( Name      = 'sig_t_angles'
                     , time      = t
                     , dm        = lifetimeParams['deltaM'] 
                     , tau       = lifetimeParams['MeanLifetime']
                     , dGamma    = lifetimeParams['deltaGamma'] 
                     , resolutionModel = tres.model()
                     , coshCoef  = basisCoefficients['cosh']
                     , cosCoef   = basisCoefficients['cos']
                     , sinhCoef  = basisCoefficients['sinh']
                     , sinCoef   = basisCoefficients['sin']
                     , ConditionalObservables = ( tag, )
                     )

### TODO: fit MC version of sig_t_angles on MC data and get efficiency,
###       then multiply sig_t_angles with the efficiency...
###       

#mcpdf = BDecay( ... )
#from P2VVGeneralUtils import RealMomentsBuilder
#eff = RealMomentsBuilder()
#eff.appendPYList( angles.angles, indices(4,2) , PDF = mcpdf, NormSet = mcpdf.getObservables( mcdata.get() ) )
#eff.compute(mcdata)
#
## replace signal pdf with efficiency corrected signal pdf
#sig_t_angles = eff * sig_t_angles

#####
from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
cmb_t = Background_Time( Name = 'cmb_t', time = t, resolutionModel = tres.model(), t_bkg_fll = dict( Name = 't_cmb_fll', Value = 0.61), t_bkg_ll_tau = dict( Name = 't_cmb_ll_tau', Value = 2.5, MinMax = (1.0,3.5)), t_bkg_ml_tau = dict( Name = 't_cmb_ml_tau', Value = 0.43) )

# Create components
ntot= data.numEntries()
(nsig,ncmb) = ( 0.9*ntot, 0.1*ntot) 
signal         = Component('signal',         ( sig_m.pdf(), sig_t_angles, ), Yield = ( nsig, 0, 1.1*nsig) )
cmb_background = Component('cmb_background', ( cmb_m.pdf(), cmb_t.pdf(),  ), Yield = ( ncmb, 0, 1.3*ncmb) )

# make sweighted dataset using J/psi phi mass
pdf_m = buildPdf((signal, cmb_background ), Observables = (m,), Name='pdf_m')
c_m = pdf_m.fitTo(data,**fitOpts)
for p in pdf_m.Parameters() : p.setConstant( not p.getAttribute('Yield') )
from P2VVGeneralUtils import SData
splot_m = SData(Pdf = pdf_m, Data = data, Name = 'MassSplot')

# create hist pdf for dilution
from RooFitWrappers import HistPdf
for c in [ signal, cmb_background ] :
    c += HistPdf(Name = "dilution_%s_pdf"%c.GetName(), Observables = [ tag ], Data = splot_m.data(c.GetName() ) )

# create PDF for angular background
from P2VVParameterizations.AngularPDFs import SPlot_Moment_Angles
mompdfBuilder = SPlot_Moment_Angles( angles.angles , splot_m )
cmb_background += mompdfBuilder.pdf( Component = cmb_background.GetName(), Indices = [ i for i in indices(3,4) ], Name = 'cmb_angles', MinSignificance = 0.5, Scale = sqrt(50.) )

# tagged fit
pdf   = buildPdf((signal, cmb_background, ), Observables = (m,t,tag)+tuple(angles.angles.itervalues()), Name='pdf')
c_pdf = pdf.fitTo(data,**fitOpts)
