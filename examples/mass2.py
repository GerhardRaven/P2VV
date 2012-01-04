from RooFitWrappers import *
#from ROOT import RooMsgService, RooFit
#msgService = RooMsgService.instance()
#msgService -= (RooFit.INFO,    RooFit.Topic(RooFit.Plotting))
#msgService += (RooFit.PROGRESS,RooFit.Topic(RooFit.Plotting))

obj  = RooObject( workspace = 'workspace')

# define observables
t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.5, 14),    nBins =  54 )
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.0, 0.15),  nBins =  50 )
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5259, 5451), nBins =  48
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3025, 3169), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1012, 1028), nBins =  16 )
iTag = Category( 'tagdecision_os', Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
eta  = RealVar('tagomega_os',      Title = 'estimated mistag',          Observable = True, MinMax = (0,0.50001),  nBins =  50)
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )

observables = [ iTag,eta, mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], st ]

#unbiased data only: /data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_DTT_after_yuehongs_script_20111220.root
from P2VVGeneralUtils import readData
#data = readData( '/tmp/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
               , dataSetName = 'DecayTree'
               , NTuple      = True
               , observables = observables
               )
#tag  = RealVar('tagdilution', Title = 'estimated dilution', Observable = True, Formula = ( "@0*(1-2*@1)",[iTag,eta], data ), MinMax=(-1,1),nBins =  51)
tag  = ConstVar(Name = 'tagdilution', Title = 'estimated dilution',  Value = 0)
#observables = [ tag, mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], st ]
observables = [ tag, mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], st ]

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
psi_st = Gamma_Sigmat( Name = 'psi_st', st = st, st_sig_gamma = dict( Name = 'st_psi_gamma', Value = 5.5 ), st_sig_beta = dict( Name = 'st_psi_beta', Value = 0.008 ) )
cmb_st = Gamma_Sigmat( Name = 'cmb_st', st = st, st_sig_gamma = dict( Name = 'st_cmb_gamma', Value = 9.6 ), st_sig_beta = dict( Name = 'st_cmb_beta', Value = 0.004 ) )

# phi mass pdf
# angle pdfs (background only)

# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tres.setConstant('.*')

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.65, deltaGamma = 0.10, deltaM = dict( Value = 17.8, MinMax = (16.5,18.5), Constant = True) )

from P2VVParameterizations.FlavourTagging import Trivial_Dilution
taggingParams = Trivial_Dilution( Dilution = tag ) # TODO: add calibration...

from math import pi
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
#amplitudes.setConstant('.*') # not all parameters appear in the untagged PDF...

# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes,CP, tag, ['A0','Apar','Aperp','AS'] ) 

# TODO: should be able to write BTagDecay('mypdf', **lifetimeParams.BTagDecay() + **basisCoefficients.BTagDecay() + **taggingParams.BTagDecay() )
# TODO: unify keys left and right....
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
                     #, ConditionalObservables = ( tag, )
                     )

from P2VVParameterizations.TimePDFs import Single_Exponent_Time as Signal_Time, LP2011_Background_Time as Background_Time
sig_t = Signal_Time(     Name = 'sig_t', time = t, resolutionModel = tres.model(), t_sig_tau = dict( Value = 1.5, Name = 't_sig_tau', MinMax=(1.0,2.0) ) )
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = tres.model(), t_bkg_fll = dict( Name = 't_psi_fll', Value = 0.47), t_bkg_ll_tau = dict( Name = 't_psi_ll_tau', Value = 1.3), t_bkg_ml_tau = dict( Name = 't_psi_ml_tau', Value = 0.21) )
cmb_t = Background_Time( Name = 'cmb_t', time = t, resolutionModel = tres.model(), t_bkg_fll = dict( Name = 't_cmb_fll', Value = 0.61), t_bkg_ll_tau = dict( Name = 't_cmb_ll_tau', Value = 2.9, MinMax = (2.5,3.5)), t_bkg_ml_tau = dict( Name = 't_cmb_ml_tau', Value = 0.43) )

# itag distribution (background only)
#from P2VVParameterizations.FlavourTagging import Trivial_Background_Tag
#bkg_tag = Trivial_Background_Tag( Name = 'bkg_tag', tagdecision = iTag, bkg_tag_eps = dict( Name = 'tag_bkg_eff', Value = 0.35 ), bkg_tag_delta = dict( Name = 'tag_bkg_delta', Value = 0 ) )
#sig_tag = Trivial_Background_Tag( Name = 'sig_tag', tagdecision = iTag, bkg_tag_eps = dict( Name = 'tag_sig_eff', Value = 0.35 ), bkg_tag_delta = dict( Name = 'tag_sig_delta', Value = 0 ) )


# Create components
(ntot,nsig,fpsi) = (data.numEntries(), 23000, 0.43)
npsi = (ntot-nsig)*fpsi
ncmb = (ntot-nsig)*(1-fpsi)
signal         = Component('signal',         ( sig_m.pdf(),  sig_mpsi.pdf(), sig_st.pdf(), sig_t_angles, ), Yield = ( nsig, 0.9*nsig, 1.1*nsig) )
psi_background = Component('psi_background', ( psi_m.pdf(),  sig_mpsi.pdf(), psi_st.pdf(), psi_t.pdf(),  ), Yield = ( npsi, 0.7*npsi, 1.3*npsi) )
cmb_background = Component('cmb_background', ( cmb_m.pdf(),  bkg_mpsi.pdf(), cmb_st.pdf(), cmb_t.pdf(),  ), Yield = ( ncmb, 0.7*ncmb, 1.3*ncmb) )


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
    for rng in ( None, 'signal','leftsideband,rightsideband' ) :
        canvas[rng] = TCanvas('%s'%rng)
        obs = observables
        obs =  [ o for o in obs if o in pdf.Observables() ]
        obs =  [ o for o in obs if hasattr(o,'frame') ]
        for (p,o) in zip( canvas[rng].pads(len(obs)), obs ) :
            dataOpts = dict( CutRange =        rng ) if rng else dict()
            pdfOpts  = dict( ProjectionRange = rng ) if rng else dict()
            if o == t  and tag in pdf.Observables() : 
                # why does prod( f(x,y)|y , g(y) ), when g(y) is a histogram, 
                # perform a numerical integral over y when plotting as a function
                # of x -- it could instead just sum over the histogram bins of g(y),
                # weighted with the binheights, i.e. sum_i g_i f(x,y_i) where 
                # g_i is g(y_i), with i = 1,Nbins....
                from ROOT import RooArgSet
                pdfOpts[ 'ProjWData' ] = ( RooArgSet(tag._var), data, True )
            from P2VVGeneralUtils import plot
            plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                                , 'psi*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                                                , 'cmb*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                                }
                                 , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
                                 , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
                                 , logy = ( o == t )
                                 )
    return (result,canvas)



# TODO: find a way (using sorted?) to specify the order of components & observables,
#       and too pass the value of 'c_opts' from outside this bit of code...
def splot( pdf, sdata ) :
    # switch off all yields, except current one
    from contextlib import contextmanager
    @contextmanager
    def __select_component( i, yields ):
        orig = dict( (j,j.getVal()) for j in yields )
        [ j.setVal(0) for j in orig.iterkeys() if j!=i ]
        try     : yield
        finally : [ j.setVal(v) for (j,v) in orig.iteritems() ]
    from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
    canvas = TCanvas(pdf.GetName() + '_splot')
    obs = [ o for o in pdf.Observables() if hasattr(o,'frame') and o not in sdata.usedObservables() ]
    for (p,o) in zip( canvas.pads(len(obs)), obs ) :
        # select yields
        _yields = [ p for p in pdf.Parameters() if p.getAttribute('Yield') ]
        # loop over components
        for (pp,i) in zip( p.pads(1,len(_yields)), _yields ) :
            # switch off all yields, except current one
            with __select_component( i, _yields ) :
                # plot both weighed data and PDF
                # TODO: add the same color coding as above...
                c_name = i.GetName()[2:]
                c_opts = { 'signal'             : dict( LineColor = kGreen )
                         , 'psi_background'     : dict( LineColor = kRed )
                         , 'cmb_background'     : dict( LineColor = kBlue )
                         }
                from P2VVGeneralUtils import plot
                plot( pp, o, sdata.data( c_name ), pdf, pdfOpts = c_opts[c_name] if c_name in c_opts else {})
    return canvas

# Build, Fit and Plot PDFs
pdf_m = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi), Name='pdf_m')
c_m = FitAndPlot(pdf_m,data)
for p in pdf_m.Parameters() : p.setConstant( not p.getAttribute('Yield') )
from P2VVGeneralUtils import SData
splot_m = SData(Pdf = pdf_m, Data = data, Name = 'MassSplot')

# get histograms for dilution...
#from RooFitWrappers import HistPdf
#for c in [ signal, psi_background,cmb_background ] :
#    c += HistPdf(Name = "dilution_%s_pdf"%c.GetName(), Observables = [ tag ], Data = splot_m.data(c.GetName() ) )
#
#pdf_di = buildPdf( (signal, cmb_background, psi_background), Observables = (tag,), Name = 'pdf_di')
#c_di = splot( pdf_di, splot_m)

from P2VVParameterizations.AngularPDFs import SPlot_Moment_Angles
mompdfBuilder = SPlot_Moment_Angles( angles.angles , splot_m )
#TODO: allow one to get RooRealVar (centered on the computed moment, with an +-N signma range) as coefficients instead of ConstVar...
indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
from math import sqrt
psi_background += mompdfBuilder.pdf( Component = psi_background.GetName(), Indices = [ i for i in indices(2,3) ], Name = 'psi_angles', MinSignificance = 0.5, Scale = sqrt(50.) )
cmb_background += mompdfBuilder.pdf( Component = cmb_background.GetName(), Indices = [ i for i in indices(3,4) ], Name = 'cmb_angles', MinSignificance = 0.5, Scale = sqrt(50.) )

#pdf_mst = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,st), Name='pdf_mst')
#c_mst = FitAndPlot(pdf_mst,data,sdata = splot_m)
#for p in pdf_mst.Parameters() : p.setConstant(True)

# untagged fit
pdf_untagged  = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,t)+tuple(angles.angles.itervalues()), Name='untagged_pdf')
CP.setConstant('.*')
amplitudes.setConstant('AperpPhase')
u_res = FitAndPlot(pdf_untagged, data )
u_spl = splot( pdf_untagged, splot_m )
# tagged fit
#pdf  = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,t,tag)+tuple(angles.angles.itervalues()), Name='pdf')
#c_pdf = FitAndPlot(pdf,data, sdata = splot_m)

