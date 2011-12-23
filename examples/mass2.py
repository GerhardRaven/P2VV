from RooFitWrappers import *
ws = RooObject( workspace = 'workspace')

t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.5, 14) )
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5259, 5451), nBins = 48 )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3025, 3169), nBins = 32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1020-8, 1020+8) )
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )

observables = [ mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'] ]

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

# phi mass pdf

# angle pdfs
#from P2VVParameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
#indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
#sig_angles = AngleBasis_AngularPdfTerms(  Angles = angles.angles
#                                          , **dict( (  'C%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
#                                                     , {  'Name'    : 'C_sig_ab%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
#                                                        , 'Value'   : 1.
#                                                        , 'MinMax'  : ( -1., +1. )
#                                                        , 'Indices' : inds
#                                                       }
#                                                    ) for inds in indices
#                                                  )
#                                         )
#mcpdf = coefPDFTerms.buildSumPdf('angCoefsPDF')
from P2VVParameterizations.AngularPDFs import Uniform_Angles
all_angles =  Uniform_Angles( angles = angles.angles )

# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t)
tres.setConstant('.*')

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as Signal_Time, LP2011_Background_Time as Background_Time
sig_t = Signal_Time(     Name = 'sig_t', time = t, resolutionModel = tres.model(), t_sig_tau = dict( Value = 1.5, Name = 't_sig_tau', MinMax=(1.0,2.0) ) )
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = tres.model() )
cmb_t = Background_Time( Name = 'cmb_t', time = t, resolutionModel = tres.model() )

# Create components
(ntot,nsig,fpsi) = (data.numEntries(), 20000, 0.4)
npsi = (ntot-nsig)*fpsi
ncmb = (ntot-nsig)*(1-fpsi)
signal         = Component('signal',         ( sig_m.pdf(),  sig_mpsi.pdf(),  sig_t.pdf(), all_angles.pdf() ), Yield = ( nsig, 0.8*nsig, 1.2*nsig) )
psi_background = Component('psi_background', ( psi_m.pdf(),  sig_mpsi.pdf(),  psi_t.pdf(), all_angles.pdf() ), Yield = ( npsi, 0.6*npsi, 1.4*npsi) )
cmb_background = Component('cmb_background', ( cmb_m.pdf(),  bkg_mpsi.pdf(),  cmb_t.pdf(), all_angles.pdf() ), Yield = ( ncmb, 0.6*ncmb, 1.4*ncmb) )

# Build PDF
pdf  = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,t)+tuple(angles.angles.itervalues()), Name='pdf')

# Fit
from P2VVGeneralUtils import numCPU
pdf.fitTo(data, NumCPU = numCPU() , Timer=1, Extended = True, Verbose = False,  Optimize=0, Minimizer = ('Minuit2','minimize'))

# Plot: TODO: define mass ranges for signal, sideband, and restrict plots to those... (see sig.py for an example)
from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
canvas = TCanvas('canvas', 'canvas', 1000, 500)
obs = [ o for o in observables if o in pdf.Observables() ]
for (p,o) in zip( canvas.pads(len(obs)), obs ) :
    from P2VVGeneralUtils import plot
    plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                        , 'psi*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                                        , 'cmb*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                        }
                         , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack )
                         , pdfOpts  = dict( LineWidth = 2 )
                         , logy = o == t
                         )
