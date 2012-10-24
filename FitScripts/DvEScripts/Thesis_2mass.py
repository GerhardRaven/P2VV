from math import sqrt, pi, cos, sin
from RooFitWrappers import *

import RootStyle
from ROOT import (gROOT,gStyle,TStyle)
MyStyle = RootStyle.MyStyle()
gROOT.SetStyle(MyStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')
from ROOT import RooMsgService
#RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Integration))
RooMsgService.instance().getStream(1).removeTopic(RooFit.Caching)
RooMsgService.instance().getStream(1).removeTopic(RooFit.Eval)

from P2VVGeneralUtils import numCPU
fitOpts = dict( NumCPU = numCPU()
              , Timer=1
              , Save = True
#              , Verbose = True
#              , Minimizer = ('Minuit2','minimize')
              )

tmincut = 0.3

# define observables
m    = RealVar('mass',  Title = 'm(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5200, 5550), nBins =  50
               ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                              , 'signal'        : ( 5330, 5410 )
                              , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'm(#mu^{+}#mu^{-})',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32
               ,Ranges = { 'leftpsisideband' : ( None, 3070 )
                           , 'psisignal'     : ( 3070, 3130 )    
                           , 'rightpsisideband' : ( 3130, None)
                           }
               )
mphi = RealVar('mdau2', Title = 'm(K^{+}K^{-})',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  16 )
t    = RealVar('time',  Title = 'Decay time',    Unit = 'ps',  Observable = True, MinMax = (tmincut, 14),    nBins =  54 )
#Set the left boundary of sigmat to non-zero to prevent problems with integration when making plots. Checked the data: no events below 0.007.
st   = RealVar('sigmat',Title = '#sigma_{t}',     Unit = 'ps',  Observable = True, MinMax = (0.007, 0.12),  nBins =  20 )
eta_os  = RealVar('tagomega_os',      Title = 'estimated mistag OS',          Observable = True, MinMax = (0,0.50001),  nBins =  20)
#The peak at 0.5 seems to be shifted to -2 in the SS eta!
eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (0,0.50001),  nBins =  20)
iTag_os = Category( 'tagdecision_os', Title = 'initial state flavour tag OS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
#The peak at 0 seems to be shifted to -1000 in the SS tagdecision
iTag_ss = Category( 'tagdecision_ss', Title = 'initial state flavour tag SS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
sel  = Category( 'sel',            Title = 'selection',                 Observable = True, States = { 'good': +1 } )
triggerdec = Category( 'triggerDecision',            Title = 'triggerdec',                 Observable = True, States = { 'triggered': +1 } )
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )

from P2VVGeneralUtils import readData

#Read data
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
                 , dataSetName = 'DecayTree'
                 , NTuple = True
                 , observables = [ m, mpsi, mphi, t, st, eta_os, eta_ss, iTag_os, iTag_ss, sel, triggerdec, angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]
                 , Rename = 'Data_DecayTree'
                 )

#TODO: Fix bug: When Rename is on, data = nil, but the dataset is imported in the ws anyway, so this is a quick fix:
#data = obj.ws().data('Data_DecayTree')

print 'Number of events', data.numEntries()

data.table(iTag_os).Print('v')
data.table(iTag_ss).Print('v')

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5368, MinMax = (5363,5372) ) )
psi_m = Background_BMass( Name = 'psi_m', mass = m, m_bkg_exp  = dict( Name = 'm_psi_exp' ) )
cmb_m = Background_BMass( Name = 'cmb_m', mass = m, m_bkg_exp  = dict( Name = 'm_cmb_exp' ) )

# J/psi mass pdf
from P2VVParameterizations.MassPDFs import Signal_PsiMass, Background_PsiMass
sig_mpsi = Signal_PsiMass(     Name = 'sig_mpsi', mass = mpsi )
bkg_mpsi = Background_PsiMass( Name = 'bkg_mpsi', mass = mpsi )

#Signal: per-event error
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as DataTimeResolution
sigtres = DataTimeResolution( time = t, timeResSFConstraint = True, sigmat = st)
psitres = sigtres
cmbtres = sigtres

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                       , dGamma = dict( Name = 'dGamma'
                                                            , Value = 0.060
                                                            #, Blind = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
                                                            )
                                       , dM = dict( Value = 17.58, MinMax = (16.5,18.5), Constant = False)
                                       , dMConstraint = True
                                      )


# define tagging parameter 
# WARNING: we don't try to describe wtag, so when plotting you must use ProjWData for eta_os !!!
from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
tagging = TaggingParams( estWTag = eta_os, p0Constraint = True, p1Constraint = True )
tagging.addConditional(iTag_os)

from P2VVParameterizations.CPVParams import LambdaArg_CPParam
CP = LambdaArg_CPParam( phiCP      = dict( Name = 'phi_s'
                                              , Value = -0.04
                                              , MinMax = (-pi,pi)
                                              #, Blind =  ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
                                              )
                           )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet
amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.52, A0Phase = 0
                                              , AperpMag2 = 0.25, AperpPhase = 2.77 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                              , AparPhase = 3.2
                                              , f_S = dict( Value = 0.02, Constant = False )
                                              , ASPhase = dict( Value = 2.7, Constant = False )
                                             )

# need to specify order in which to traverse...
from RooFitWrappers import RealCategory
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes, CP, ['A0','Apar','Aperp','AS'] )

from RooFitWrappers import BTagDecay
sig_t_angles_iTag = BTagDecay(  Name                   = 'sig_t_angles_iTag'
                              , time                   = t
                              , iTag                   = iTag_os
                              , dm                     = lifetimeParams['dM']
                              , tau                    = lifetimeParams['MeanLifetime']
                              , dGamma                 = lifetimeParams['dGamma']
                              , resolutionModel        = sigtres['model']
                              , coshCoef               = basisCoefficients['cosh']
                              , cosCoef                = basisCoefficients['cos']
                              , sinhCoef               = basisCoefficients['sinh']
                              , sinCoef                = basisCoefficients['sin']
                              , dilution               = tagging['dilution']
                              , ADilWTag               = tagging['ADilWTag']
                              , avgCEven               = tagging['avgCEven']
                              , avgCOdd                = tagging['avgCOdd']
                              , ConditionalObservables = sigtres.conditionalObservables() + tagging.conditionalObservables()
                              , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                                         + sigtres.externalConstraints()\
                                                         + tagging.externalConstraints()
                             )


from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
psi_t = Background_Time( Name = 'psi_t', time = t, resolutionModel = psitres.model()
                         , t_bkg_fll    = dict( Name = 't_psi_fll',    Value = 0.21 )
                         , t_bkg_ll_tau = dict( Name = 't_psi_ll_tau', Value = 1.06, MinMax = (0.5,2.5) )
                         , t_bkg_ml_tau = dict( Name = 't_psi_ml_tau', Value = 0.15, MinMax = (0.01,0.5) )
                         )
cmb_t = Background_Time( Name = 'cmb_t', time = t, resolutionModel = cmbtres.model()
                         , t_bkg_fll = dict( Name = 't_cmb_fll', Value = 0.61)
                         , t_bkg_ll_tau = dict( Name = 't_cmb_ll_tau', Value = 1.5, MinMax = (0.1,2.5))
                         , t_bkg_ml_tau = dict( Name = 't_cmb_ml_tau', Value = 0.43) )

print "*****************************"
print [ i.GetName() for i in sig_t_angles_iTag.ConditionalObservables()  ]
print "*****************************"

print "*****************************"
print [ i.GetName() for i in psi_t._pdf.ConditionalObservables()  ]
print "*****************************"

print "*****************************"
print [ i.GetName() for i in cmb_t._pdf.ConditionalObservables()  ]
print "*****************************"

#####################################
### Angular acceptance correction ###
#####################################
from P2VVGeneralUtils import RealMomentsBuilder
#nset = angles.angles.values()

#canomoms = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, MCpdf,nset) for v in angles.functions.itervalues() for i in v if i ) )
#canomoms.compute(MCdata)
#canomoms.Print(Scales = [1./(16.*sqrt(pi)),1./(16.*sqrt(pi)),1./(16.*sqrt(pi))])

from itertools import chain
#momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))
#These are the relevant terms as found with MinSignificance>3
momindices = [(0,0,0),(0,2,0),(0,2,2),(2,0,0)]

eff = RealMomentsBuilder()
#Don't specify pdf and normset here, we're gonna read moments and not calculate any.
eff.appendPYList( angles.angles, momindices)
eff.read('/data/bfys/dveijk/DataJpsiPhi/2012/effmoments_tcut_%s.txt'%(str(tmincut)))

#Build Angular acceptance corrected PDF
sig_t_angles_iTag = eff * sig_t_angles_iTag

##############################
### Proper time acceptance ###
##############################
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/data/bfys/dveijk/DataJpsiPhi/2012/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root', Histogram = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')
sig_t_angles_iTag = acceptance * sig_t_angles_iTag
#sig_t_angles_iTag.setAttribute("NOCacheAndTrack")

# Create components
(ntot,nsig,fpsi) = (data.numEntries(), 23000, 0.7)
npsi = (ntot-nsig)*fpsi
ncmb = (ntot-nsig)*(1-fpsi)
signal         = Component('signal',         ( sig_m.pdf(),  sig_mpsi.pdf(), sig_t_angles_iTag), Yield = ( nsig, 0.5*nsig, 2.0*nsig) )
psi_background = Component('psi_background', ( psi_m.pdf(),  sig_mpsi.pdf(), psi_t.pdf(), { eta_os: None, iTag_os : None }), Yield = ( npsi, 0.5*npsi, 5.0*npsi) )
cmb_background = Component('cmb_background', ( cmb_m.pdf(),  bkg_mpsi.pdf(), cmb_t.pdf(), { eta_os: None, iTag_os : None }), Yield = ( ncmb, 0.5*ncmb, 5.0*ncmb) )

############################
### B MASS AND JPSI MASS ###
############################
pdf_m_mpsi = buildPdf((signal, psi_background, cmb_background), Observables = (m,mpsi), Name='pdf_m_mpsi')
result_m_mpsi = pdf_m_mpsi.fitTo(data, **fitOpts)

plot_m_mpsi = True

from ROOT import kDashed, kRed, kBlack, kGreen, kOrange, TCanvas, TLatex
from P2VVGeneralUtils import plot

if plot_m_mpsi:

    #plot on m frame
    canvas_m = TCanvas('bmass','bmass',972,600)
    canvas_m.SetFixedAspectRatio(True)
    plot( canvas_m
          , m
          , data
          , pdf_m_mpsi
          , xTitleOffset = 1.1
          , components = {'signal*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          ,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                          ,'cmb*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                          }
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack)
          , pdfOpts = dict( LineWidth = 3 )
          #, plotResidHist = True
          , logy = False
          , frameOpts = dict( Title = ''
                              , TitleOffset = (1.1,"xy")
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=70
                              ) 
          )

    #plot on mpsi frame
    canvas_mpsi = TCanvas('psimass','psimass',972,600)
    canvas_mpsi.SetFixedAspectRatio(True)
    plot( canvas_mpsi
          , mpsi
          , data
          , pdf_m_mpsi
          , xTitleOffset = 1.1
          , components = {'signal*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          ,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                          ,'cmb*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                          }
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack)
          , pdfOpts = dict( LineWidth = 3 )
          #, plotResidHist = True
          , logy = False
          , frameOpts = dict( Title = ''
                              , TitleOffset = (1.1,"xy")
                              , Object = ( TLatex(0.15,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=30
                              ) 
          )

###########################
### Calculate S-weights ###
###########################
# make sweighted dataset using Bs AND Jpsi mass
from P2VVGeneralUtils import SData

#for p in pdf_m_mpsi.Parameters() : p.setConstant( not p.getAttribute('Yield') )
splot_m = SData(Pdf = pdf_m_mpsi, Data = data, Name = 'MassSplot')

##################################
### Define components datasets ###
##################################
signaldata = splot_m.data('signal')
psibkgdata = splot_m.data('psi_background')
cmbbkgdata = splot_m.data('cmb_background')

##################
### BKG angles ###
##################
#Flat
if False:
    for i in angles.angles.itervalues():
        psi_background[i] = None
        cmb_background[i] = None
#Histograms
else:
    psi_background += HistPdf( Name = 'psibkg_angles'
                               , Observables = angles.angles.itervalues()
                               , Binning =  { angles.angles['cpsi']   : 5
                                              , angles.angles['ctheta'] : 7
                                              , angles.angles['phi' ]   : 9
                                              }
                               , Data  = psibkgdata
                               )
    cmb_background += HistPdf( Name = 'cmbbkg_angles'
                               , Observables = angles.angles.itervalues()
                               , Binning =  { angles.angles['cpsi']   : 5
                                              , angles.angles['ctheta'] : 7
                                              , angles.angles['phi' ]   : 18
                                              }
                               , Data  = cmbbkgdata
                               )

bkgangleplot = True

if bkgangleplot:
    psibkganglepdf = buildPdf((psi_background,), Observables = (angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']), Name='psibkganglepdf')
    cmbbkganglepdf = buildPdf((cmb_background,), Observables = (angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']), Name='cmbbkganglepdf')

    orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )
    obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']]

    psibkganglecanvas = TCanvas('psibkganglecanvas','psibkganglecanvas',972,600)
    psibkganglecanvas.SetFixedAspectRatio(True)        
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( psibkganglecanvas.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        from P2VVGeneralUtils import plot
        from ROOT import RooArgSet
        pdfOpts = dict()
        pdfOpts[ 'ProjWData' ] = ( RooArgSet(st._var, eta_os._var), data, True )
        plot( p, o, psibkgdata, psibkganglepdf
              , yTitleOffset = 1.2
              , xTitleOffset = 1.1
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack,)
              , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
              , plotResidHist = False
              , frameOpts = dict( Object = ( TLatex(0.30,.4,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                                  , Title = ""
                                  )
              , logy = ( o == t )
              )

    cmbbkganglecanvas = TCanvas('cmbbkganglecanvas','cmbbkganglecanvas',972,600)
    cmbbkganglecanvas.SetFixedAspectRatio(True)        
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( cmbbkganglecanvas.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        from P2VVGeneralUtils import plot
        from ROOT import RooArgSet
        pdfOpts = dict()
        pdfOpts[ 'ProjWData' ] = ( RooArgSet(st._var, eta_os._var), data, True )
        plot( p, o, cmbbkgdata, cmbbkganglepdf
              , yTitleOffset = 1.2
              , xTitleOffset = 1.1
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack,)
              , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
              , plotResidHist = False
              , frameOpts = dict( Object = ( TLatex(0.30,.18,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                                  , Title = ""
                                  )
              , logy = ( o == t )
              )

############
# BKG ONLY #
############

bkgtimeplot = True
if bkgtimeplot:
    psibkgtimepdf = buildPdf((psi_background,), Observables = (t,), Name = 'psibkgtimepdf')
    psibkgtimepdf.fitTo(psibkgdata,**fitOpts)

    canvas_psibkg_t = TCanvas('psibkg_t','psibkg_t',972,600)
    canvas_psibkg_t.SetFixedAspectRatio(True)
    plot( canvas_psibkg_t
          , t
          , psibkgdata
          , psibkgtimepdf
          , xTitleOffset = 1.1
          #, components = {'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange   )
          #                ,'cmb*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , pdfOpts = dict(  LineWidth = 3, ProjWData = ( psibkgdata.reduce(  ArgSet = psibkgtimepdf.ConditionalObservables() ), True ) )
          #, pdfOpts = dict( LineWidth = 3 )
          , plotResidHist = False
          , logy = False
          , frameOpts = dict( Title = ''
                              , TitleOffset = (1.2,'y')
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              , Bins=137 ) 
          )

    cmbbkgtimepdf = buildPdf((cmb_background,), Observables = (t,), Name = 'cmbbkgtimepdf')
    cmbbkgtimepdf.fitTo(cmbbkgdata,**fitOpts)

    canvas_cmbbkg_t = TCanvas('cmbbkg_t','cmbbkg_t',972,600)
    canvas_cmbbkg_t.SetFixedAspectRatio(True)
    plot( canvas_cmbbkg_t
          , t
          , cmbbkgdata
          , cmbbkgtimepdf
          , xTitleOffset = 1.1
          #, components = {'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange   )
          #                ,'cmb*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , pdfOpts = dict(  LineWidth = 3, ProjWData = ( cmbbkgdata.reduce(  ArgSet = cmbbkgtimepdf.ConditionalObservables() ), True ) )
          #, pdfOpts = dict( LineWidth = 3 )
          , plotResidHist = False
          , logy = False
          , frameOpts = dict( Title = ''
                              , TitleOffset = (1.2,'y')
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              , Bins=137 ) 
          )


pdf   = buildPdf((signal,psi_background,cmb_background), Observables = (m,mpsi,t,)+tuple(angles.angles.itervalues()), Name='fullpdf')

def search(fname,path) :
    import os
    for f in ( os.path.join(p,fname) for p in os.path.split(os.pathsep) ) :
        if os.path.exists(f) : return f
    return None

import os

fitname = 'cfit_2mass_acc'

paramfile = search(fitname+'params.txt',os.pathsep.join(['.','FitScripts']) )
if paramfile :
    print 'Reading fit result from %s' % paramfile
    fitset = pdf.getParameters(data)
    fitset.readFromFile(paramfile)

fit = True
if fit or not paramfile:
    cfitresult = pdf.fitTo(data, **fitOpts)
    cfitresult.writepars(fitname+'result',False)
    cfitresult.Print()
    fitset = pdf.getParameters(data)
    fitset.writeToFile(fitname+"params.txt")


assert False

print 120 * '='
print 'DvEFit: parameters in NLL:'
for par in pdf.getVariables() : par.Print()
print 120 * '='

##     #for rng in ( None, 'signal','leftsideband,rightsideband' ) :
##     for rng in ( None, ) :
##         canvas[rng] = TCanvas('%s'%rng,'%s'%rng)
##         obs =  [ o for o in pdf.Observables() if hasattr(o,'frame') ]
##         for (p,o) in zip( canvas[rng].pads(len(obs)), obs ) :
##             dataOpts = dict( CutRange =        rng ) if rng else dict()
##             pdfOpts  = dict( ProjectionRange = rng ) if rng else dict()
##             if o == t  and tag in pdf.Observables() : 
##                 # why does prod( f(x,y)|y , g(y) ), when g(y) is a histogram, 
##                 # perform a numerical integral over y when plotting as a function
##                 # of x -- it could instead just sum over the histogram bins of g(y),
##                 # weighted with the binheights, i.e. sum_i g_i f(x,y_i) where 
##                 # g_i is g(y_i), with i = 1,Nbins....
##                 from ROOT import RooArgSet
##                 pdfOpts[ 'ProjWData' ] = ( RooArgSet(tag._var), data, True )
##             from P2VVGeneralUtils import plot
##             plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
##                                                 , 'psi*'     : dict( LineColor = kRed,   LineStyle = kDashed )
##                                                 , 'cmb*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
##                                                 }
##                                  , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
##                                  , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
##                                  , logy = ( o == t )
##                                  )
##     return (result,canvas)



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
        _yields = [ y for y in pdf.Parameters() if y.getAttribute('Yield') ]
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


from P2VVGeneralUtils import SData
splot_m = SData(Pdf = pdf_m, Data = data, Name = 'MassSplot')

from P2VVParameterizations.AngularPDFs import SPlot_Moment_Angles
mompdfBuilder = SPlot_Moment_Angles( angles.angles , splot_m )
#TODO: allow one to get RooRealVar (centered on the computed moment, with an +-N signma range) as coefficients instead of ConstVar...
indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
from math import sqrt
psi_background += mompdfBuilder.pdf( Component = psi_background.GetName(), Indices = [ i for i in indices(2,2) ], Name = 'psi_angles', MinSignificance = 0.5, Scale = sqrt(50.) )
cmb_background += mompdfBuilder.pdf( Component = cmb_background.GetName(), Indices = [ i for i in indices(2,2) ], Name = 'cmb_angles', MinSignificance = 0.5, Scale = sqrt(50.) )


assert False
# untagged fit
pdf_untagged  = buildPdf((signal, cmb_background, psi_background), Observables = (m,mpsi,t)+tuple(angles.angles.itervalues()), Name='untagged_pdf')
CP.setConstant('.*')
amplitudes.setConstant('AperpPhase')
u_res = FitAndPlot(pdf_untagged, data )
u_spl = splot( pdf_untagged, splot_m )

