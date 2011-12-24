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
coefPDFTerms = AngleBasis_AngularPdfTerms(  Angles = angles.angles
                                          , **dict( (  ('C%d%d%d' % i ).replace('-','m')
                                                     , {  'Name'     : ( 'C_%d%d%d' % i ).replace('-','m')
                                                        , 'Value'    : 1. if i == (0,0,0) else 0.
                                                        , 'MinMax'   : ( 0.99, 1.01 )  if i == (0,0,0) else ( -0.4,+0.4 ) if i[1]==0 else (-0.1,0.1)
                                                        , 'Indices'  : i
                                                        , 'Constant' : i == (0,0,0)
                                                       }
                                                    ) for i in indices(4,3)
                                                  )
                                         )
bkg_angles = coefPDFTerms.buildSumPdf('angCoefsPDF')
#from P2VVParameterizations.AngularPDFs import Uniform_Angles
#bkg_angles =  Uniform_Angles( angles = angles.angles )

# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tres.setConstant('.*')

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
signal         = Component('signal',         ( sig_m.pdf(),  sig_mpsi.pdf(),  sig_t.pdf(), sig_st.pdf(), bkg_angles ), Yield = ( nsig, 0.9*nsig, 1.1*nsig) )
psi_background = Component('psi_background', ( psi_m.pdf(),  sig_mpsi.pdf(),  psi_t.pdf(), sig_st.pdf(), bkg_angles ), Yield = ( npsi, 0.7*npsi, 1.3*npsi) )
cmb_background = Component('cmb_background', ( cmb_m.pdf(),  bkg_mpsi.pdf(),  cmb_t.pdf(), cmb_st.pdf(), bkg_angles ), Yield = ( ncmb, 0.7*ncmb, 1.3*ncmb) )

# Build PDF
pdf  = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,t,st), Name='pdf')
#pdf  = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,t)+tuple(angles.angles.itervalues()), Name='pdf')

# Fit
from P2VVGeneralUtils import numCPU, ROOTversion
(ROOTmajor,ROOTminor,ROOTpatch) = ROOTversion()
pdf.fitTo(data, NumCPU = numCPU() 
              , Timer=1
              , Verbose = False
              , Optimize = True if ROOTminor<32 else 0 # do NOT optimize in 5.32 or later... ( Optimize = 1 only works on a single CPU )
              , Minimizer = ('Minuit2','minimize'))

# Plot: TODO: define mass ranges for signal, sideband, and restrict plots to those... (see sig.py for an example)
from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
canvas = dict()
for rng in ( None, 'signal','leftsideband,rightsideband' ) :
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
