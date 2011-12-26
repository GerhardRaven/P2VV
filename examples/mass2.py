from RooFitWrappers import *
obj  = RooObject( workspace = 'workspace')

# define observables
t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.5, 14),    nBins =  54 )
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.0, 0.15),  nBins =  50 )
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5259, 5451), nBins =  48 )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3025, 3169), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1012, 1028), nBins =  16 )
iTag = Category( 'tagdecision' , Title = 'initial state flavour tag',   Observable = True,  States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )

observables = [ iTag, mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], st ]

# TODO: add a wrapper that allows one to pass ranges as dict: m.setRanges( { 'leftsideband' : (m.getMin(),5330), 'signal' : (5330,5410), ... } )
#       even more fance: allow None in the upper or lower limit to indicate it should extend to getMax resp. getMin...
#       m.setRanges( { 'leftsideband'  : ( None, 5330 )
#                    , 'signal'        : ( 5330, 5410 )
#                    , 'rightsideband' : ( 5410, None ) } )
m.setRange('leftsideband', (m.getMin(),5330) )
m.setRange('signal',(5330,5410) )
m.setRange('rightsideband',(5410,m.getMax()) )

#unbiased data only: /data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_DTT_after_yuehongs_script_20111220.root
from P2VVGeneralUtils import readData
# data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
data = readData( '/tmp/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
               , 'DecayTree'
               , True
               , observables
               )
# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5368, MinMax = (5363,5372) ) )
psi_m = Background_BMass( Name = 'psi_m', mass = m, m_bkg_exp  = dict( Name = 'm_psi_exp' ) )
cmb_m = Background_BMass( Name = 'cmb_m', mass = m, m_bkg_exp  = dict( Name = 'm_cmb_exp' ) )

# J/psi mass pdf
from P2VVParameterizations.MassPDFs import Signal_PsiMass, Background_PsiMass
sig_mpsi = Signal_PsiMass(     Name = 'sig_mpsi', mass = mpsi )
bkg_mpsi = Background_PsiMass( Name = 'bkg_mpsi', mass = mpsi )

# sigma(t) pdf
from P2VVParameterizations.TimeResolution import Gamma_Sigmat
sig_st = Gamma_Sigmat( Name = 'sig_st', st = st, st_sig_gamma = dict( Name = 'st_sig_gamma', Value = 12. ), st_sig_beta = dict( Name = 'st_sig_beta', Value = 0.003 ) )
cmb_st = Gamma_Sigmat( Name = 'cmb_st', st = st, st_sig_gamma = dict( Name = 'st_cmb_gamma', Value = 9.6 ), st_sig_beta = dict( Name = 'st_cmb_beta', Value = 0.004 ) )
psi_st = Gamma_Sigmat( Name = 'psi_st', st = st, st_sig_gamma = dict( Name = 'st_psi_gamma', Value = 5.5 ), st_sig_beta = dict( Name = 'st_psi_beta', Value = 0.008 ) )

# phi mass pdf

# angle pdfs (background only)
from P2VVParameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
###TODO: this builds a RealSumPdf with coefficients one, and functions prod( C,fun ). Replace with coef C, functions fun
coefPDFTerms = AngleBasis_AngularPdfTerms(  Angles = angles.angles
                                          , **dict( (  ('C%d%d%d' % i ).replace('-','m')
                                                     , {  'Name'     : ( 'C_%d%d%d' % i ).replace('-','m')
                                                        , 'Value'    : 1. if i == (0,0,0) else 0.
                                                        , 'MinMax'   : ( 0.99, 1.01 )  if i == (0,0,0) else ( -0.4,+0.4 ) if i[1] == 0 else (-0.1,0.1)
                                                        , 'Indices'  : i
                                                        , 'Constant' : i == (0,0,0)
                                                       }
                                                    ) for i in indices(3,3)
                                                  )
                                         )
#TODO: adjust syntax to match the rest... move Name into __init__ as kw...
bkg_angles = coefPDFTerms.buildSumPdf('bkg_angles_PDF')
#from P2VVParameterizations.AngularPDFs import Uniform_Angles
#bkg_angles =  Uniform_Angles( angles = angles.angles )

# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tres.setConstant('.*')

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.68, deltaGamma = 0.05, deltaM = dict( Value = 17.8, MinMax = (16,19), Constant = True) )

from P2VVParameterizations.FlavourTagging import Trivial_TaggingParams
#taggingParams = Trivial_TaggingParams( wTag = eta ) # FormulaVar('wTag','@2 + @3*(@0-@1)',[eta,etaAverage,p0,p1] ) )
taggingParams = Trivial_TaggingParams( wTag = ConstVar( Name = 'half', Value = 0.5 ) ) # FormulaVar('wTag','@2 + @3*(@0-@1)',[eta,etaAverage,p0,p1] ) )

from math import pi
from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP      = dict( Name = 'phi_s', Value = -0.04, MinMax = (-pi,pi), Constant = True )
                        , lambdaCPSq = ConstVar( Name = 'one', Value = 1 ) 
                        )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011( A0Mag2 = 0.60, A0Phase = 0
                                    , AperpMag2 = 0.160, AperpPhase = -0.17
                                    , AparPhase = 2.5
                                    , ASMag2 = dict( Value = 0, Constant = True ) , ASPhase = dict( Value = 0, Constant = True ) 
                                    )
#amplitudes.setConstant('.*') # not all parameters appear in the untagged PDF...

# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

# TODO: should be able to write BTagDecay('mypdf', **lifetimeParams.BTagDecay() + **basisCoefficients.BTagDecay() + **taggingParams.BTagDecay() )
# TODO: unify keys left and right....
sig_t_angles_tag = BTagDecay( Name      = 'sig_t_angles_tag'
                            , time      = t
                            , iTag      = iTag
                            , dm        = lifetimeParams['deltaM'] 
                            , tau       = lifetimeParams['MeanLifetime']
                            , dGamma    = lifetimeParams['deltaGamma'] 
                            , resolutionModel = tres.model()
                            , coshCoef  = basisCoefficients['cosh']
                            , cosCoef   = basisCoefficients['cos']
                            , sinhCoef  = basisCoefficients['sinh']
                            , sinCoef   = basisCoefficients['sin']
                            , avgCEven  = taggingParams['avgCEven'] 
                            , avgCOdd   = taggingParams['avgCOdd']
                            , dilution  = taggingParams['dilution']
                            , ADilWTag  = taggingParams['ADilWTag']
                            )

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as Signal_Time, LP2011_Background_Time as Background_Time
sig_t = Signal_Time(     Name = 'sig_t', time = t, resolutionModel = tres.model(), t_sig_tau = dict( Value = 1.5, Name = 't_sig_tau', MinMax=(1.0,2.0) ) )
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = tres.model(), t_bkg_fll = dict( Name = 't_psi_fll'), t_bkg_ll_tau = dict( Name = 't_psi_ll_tau'), t_bkg_ml_tau = dict( Name = 't_psi_ml_tau') )
cmb_t = Background_Time( Name = 'cmb_t', time = t, resolutionModel = tres.model(), t_bkg_fll = dict( Name = 't_cmb_fll'), t_bkg_ll_tau = dict( Name = 't_cmb_ll_tau'), t_bkg_ml_tau = dict( Name = 't_cmb_ml_tau') )

# itag distribution (background only)
from P2VVParameterizations.FlavourTagging import Trivial_Background_Tag
bkg_tag = Trivial_Background_Tag( tagdecision = iTag, bkg_tag_delta = 0.0 )

# Create components
(ntot,nsig,fpsi) = (data.numEntries(), 23000, 0.43)
npsi = (ntot-nsig)*fpsi
ncmb = (ntot-nsig)*(1-fpsi)
signal         = Component('signal',         ( sig_m.pdf(),  sig_mpsi.pdf(), sig_st.pdf(), sig_t_angles_tag                       ), Yield = ( nsig, 0.9*nsig, 1.1*nsig) )
psi_background = Component('psi_background', ( psi_m.pdf(),  sig_mpsi.pdf(), sig_st.pdf(), psi_t.pdf(),  bkg_tag.pdf() ), Yield = ( npsi, 0.7*npsi, 1.3*npsi) )
cmb_background = Component('cmb_background', ( cmb_m.pdf(),  bkg_mpsi.pdf(), cmb_st.pdf(), cmb_t.pdf(),  bkg_tag.pdf() ), Yield = ( ncmb, 0.7*ncmb, 1.3*ncmb) )


def FitAndPlot( pdf, data, fitOpts = dict() ) :
    # Fit
    from P2VVGeneralUtils import numCPU
    from ROOTDecorators import  ROOTversion as Rv
    result = pdf.fitTo(data, NumCPU = numCPU() 
                           , Timer=1
                           , Save = True
                           , Verbose = False
                           , Optimize = True if Rv[1]<32 else 0 # do NOT optimize in 5.32 or later... ( Optimize = 1 only works on a single CPU, 2 doesn't work at all )
                           , Minimizer = ('Minuit2','minimize')
                           , **fitOpts
                           )

    # Plot: 
    from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
    canvas = dict()
    for rng in ( None, 'signal','leftsideband','rightsideband','leftsideband,rightsideband' ) :
        canvas[rng] = TCanvas('%s'%rng)
        obs = observables
        obs =  [ o for o in obs if o in pdf.Observables() ]
        obs =  [ o for o in obs if hasattr(o,'frame') ]
        for (p,o) in zip( canvas[rng].pads(len(obs)), obs ) :
            dataRng = dict( CutRange =        rng ) if rng else dict()
            pdfRng  = dict( ProjectionRange = rng ) if rng else dict()
            from P2VVGeneralUtils import plot
            plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                                , 'psi*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                                                , 'cmb*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                                }
                                 , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataRng )
                                 , pdfOpts  = dict( LineWidth = 2, **pdfRng )
                                 , logy = ( o == t )
                                 )
    return (result,canvas)

# Build, Fit and Plot PDFs
# TODO: add sPlots for various components
pdf_m = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi), Name='pdf_m')
c_m = FitAndPlot(pdf_m,data)
for p in pdf_m.Parameters() : p.setConstant( not p.getAttribute('Yield') )
#TODO: move more into Moment_Angles... and allow one to get RooRealVar (centered on the moment) as coefficients instead of ConstVar...
from ROOT import RooStats
splot_m = RooStats.SPlot("splotdata","splotdata",data,pdf_m._var, RooArgList( p._var for p in pdf_m.Parameters() if p.getAttribute('Yield') ) )
sdata   = splot_m.GetSDataSet()
from P2VVParameterizations.AngularPDFs import Moment_Angles
mompdfBuilder = Moment_Angles( angles.angles , splot_m.GetSDataSet() )
psi_background += mompdfBuilder.pdf( Component = 'psi_background', Indices = [ i for i in indices(3,3) ], Name = 'psi_angles' )
cmb_background += mompdfBuilder.pdf( Component = 'cmb_background', Indices = [ i for i in indices(3,3) ], Name = 'cmb_angles' )

pdf_mst = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,st), Name='pdf_mst')
c_mst = FitAndPlot(pdf_mst,data)
for p in pdf_mst.Parameters() : p.setConstant(True)

# first fit background angles in sidebands, then fix the background angle parameters
# TODO: compute moments -- that's a LOT faster...
#bkg_pdf  = buildPdf((cmb_background, psi_background), Observables = (m,mpsi)+tuple(angles.angles.itervalues()), Name='bkg_pdf')
#c_bkg = FitAndPlot( bkg_pdf, data, fitOpts = dict( Range = 'leftsideband,rightsideband' ) )
#for p in bkg_pdf.Parameters() : p.setConstant(True)

# full fit
pdf  = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,t,iTag)+tuple(angles.angles.itervalues()), Name='pdf')
c_pdf = FitAndPlot(pdf,data)

