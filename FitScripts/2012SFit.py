from math import sqrt, pi, cos, sin
from RooFitWrappers import *
from ROOT import RooMsgService
from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
#RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Integration))
RooMsgService.instance().getStream(1).removeTopic(RooFit.Caching)
RooMsgService.instance().getStream(1).removeTopic(RooFit.Eval)

## import RootStyle
## from ROOT import (gROOT,gStyle,TStyle)
## MyStyle = RootStyle.MyStyle()
## gROOT.SetStyle(MyStyle.GetName())
## gROOT.ForceStyle()
## gStyle.UseCurrentStyle()   

indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')

from P2VVGeneralUtils import numCPU

fitOpts = dict( Timer = 1
              , NumCPU = numCPU()
              , Save = True
              , Minimizer = ('Minuit2','minimize')
#              , Verbose = True
              )

useHelicity = False
tmincut = 0.3
path =  '/data/bfys/dveijk/DataJpsiPhi/2012/'
path =  '/Users/raven/data/'
fpsi = 0.7
fitname = 'sfit' if fpsi==0 else 'sfit_mm'
fitname += '_hel' if useHelicity else '_tr'

# define observables
m    = RealVar('mass',  Title = 'M(J/#psi#varphi)', Unit = 'MeV/c^{2}', Observable = True, MinMax = (5200, 5550), nBins =  100
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  40 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  16 )
t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (tmincut, 14),    nBins =  54 )
#Set the left boundary of sigmat to non-zero to prevent problems with integrations (eg when making plots). Checked the data: no events below 0.007.
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.007, 0.12),  nBins =  15 )
eta_os  = RealVar('tagomega_os',      Title = 'estimated mistag OS',          Observable = True, MinMax = (0,0.50001),  nBins =  20)
#The peak at 0.5 seems to be shifted to -2 in the SS eta!
eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (0,0.50001),  nBins =  25)
iTag_os = Category( 'tagdecision_os', Title = 'initial state flavour tag OS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
#The peak at 0 seems to be shifted to -1000 in the SS tagdecision
iTag_ss = Category( 'tagdecision_ss', Title = 'initial state flavour tag SS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
sel  = Category( 'sel',            Title = 'selection',                 Observable = True, States = { 'good': +1 } )
triggerdec = Category( 'triggerDecision',            Title = 'triggerdec',                 Observable = True, States = { 'triggered': +1 } )
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles

if useHelicity :
    angles    = HelAngles( cpsi   = dict( Name = 'helcosthetaK', Title = 'cos(#theta_{K}',    nBins = 12)
                         , ctheta = dict( Name = 'helcosthetaL', Title = 'cos(#theta_{#ell}', nBins = 12)
                         , phi    = dict( Name = 'helphi',       Title = '#phi',              nBins = 12) )
else :
    angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 12 )
                        , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 12 )
                        , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 12 ) 
                        )

# add 20 bins for caching the normalization integral
for i in [ eta_os, st ] : i.setBins( 20 , 'cache' )

# define plot order
orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )

#Read data
from P2VVGeneralUtils import readData
data = readData(path + 'Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
                 , dataSetName = 'DecayTree'
                 , NTuple = True
                 , observables = [ m, mpsi, mphi, t, st, eta_os, eta_ss, iTag_os, iTag_ss, sel, triggerdec, angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]
                 , Rename = 'Data_DecayTree'
                 )

print 'Number of events', data.numEntries()

data.table(iTag_os).Print('v')
data.table(iTag_ss).Print('v')

####################
### Create components which will contain the various PDFs ###
####################

(nsig,nbkg) = (21000,10500)
signal         = Component('sgnl', Yield = ( nsig,          0, 1.4*nsig) )
background     = Component('cmb',  Yield = ( nbkg*(1-fpsi), 0, 1.4*nbkg*(1-fpsi)) )
psi            = Component('psi',  Yield = ( nbkg*fpsi,     0, 1.4*nbkg*fpsi ) )

# J/psi mass pdf
from P2VVParameterizations.MassPDFs import Signal_PsiMass, Background_PsiMass
sig_mpsi = Signal_PsiMass(     Name = 'sig_mpsi', mass = mpsi )
bkg_mpsi = Background_PsiMass( Name = 'cmb_mpsi', mass = mpsi )

signal     += sig_mpsi.pdf()
background += bkg_mpsi.pdf()
psi        += sig_mpsi.pdf()

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5365, MinMax = (5363,5372) ) )
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )
psi_m = Background_BMass( Name = 'psi_m', mass = m, m_bkg_exp  = dict( Name = 'm_psi_exp' ) )

signal     += sig_m.pdf()
background += bkg_m.pdf()
psi        += psi_m.pdf()

#Time Resolution Model
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as DataTimeResolution
tresdata = DataTimeResolution( time = t, timeResSFConstraint = True, sigmat = st)
signal     += { st: None }
background += { st: None }
psi        += { st: None }

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                       , dGamma = dict( Name = 'dGamma' , Value = 0.060)
                                       , dM = dict( Value = 17.58, MinMax = (16.5,18.5), Constant = False)
                                       , dMConstraint = True
                                     )

# define tagging parameter 
from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
tagging = TaggingParams( estWTag = eta_os, p0Constraint = True, p1Constraint = True )
tagging.addConditional(iTag_os)

from P2VVParameterizations.CPVParams import LambdaArg_CPParam, LambdaSqArg_CPParam
CP = LambdaArg_CPParam( phiCP      = dict( Name = 'phi_s' , Value = -0.04 , MinMax = (-pi,pi)))
#Fit for |lambda|^2
#CP = LambdaSqArg_CPParam( phiCP      = dict( Name = 'phi_s' , Value = -0.04 , MinMax = (-pi,pi))
#                           )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet
amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.52, A0Phase = 0
                                              , AperpMag2 = 0.25, AperpPhase = 2.77 
                                              , AparPhase = 3.2
                                              , f_S = dict( Value = 0.02, Constant = False )
                                              , ASPhase = dict( Value = 2.7, Constant = False )
                                             )
#from P2VVParameterizations.DecayAmplitudes import  JpsiVCarthesian_AmplitudeSet
#amplitudes = JpsiVCarthesian_AmplitudeSet( ReAS    = -0.21, ImAS    =  0.07
#                                         , ReApar  = -0.65, ImApar  = -0.14
#                                         , ReAperp = -0.65, ImAperp =  0.20 )
#

# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes, CP, ['A0','Apar','Aperp','AS'] )

from RooFitWrappers import BTagDecay
sig_t_angles_iTag = BTagDecay(  Name                   = 'sig_t_angles_iTag'
                              , time                   = t
                              , iTag                   = iTag_os
                              , dm                     = lifetimeParams['dM']
                              , tau                    = lifetimeParams['MeanLifetime']
                              , dGamma                 = lifetimeParams['dGamma']
                              , resolutionModel        = tresdata['model']
                              , coshCoef               = basisCoefficients['cosh']
                              , cosCoef                = basisCoefficients['cos']
                              , sinhCoef               = basisCoefficients['sinh']
                              , sinCoef                = basisCoefficients['sin']
                              , dilution               = tagging['dilution']
                              , ADilWTag               = tagging['ADilWTag']
                              , avgCEven               = tagging['avgCEven']
                              , avgCOdd                = tagging['avgCOdd']
                              , ConditionalObservables = tresdata.conditionalObservables() + tagging.conditionalObservables()
                              , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                                         + tresdata.externalConstraints()\
                                                         + tagging.externalConstraints()
                             )

signal += { eta_os : None }
if iTag_os in sig_t_angles_iTag.ConditionalObservables():
    from P2VVParameterizations.FlavourTagging import Trivial_TagPdf
    print 'Adding PDF for %s' % iTag_os.GetName()
    # TODO: count # of iTag_os = +1,0,-1 in sweighted data set in order to set the additonal
    #       parameters to their eventual fit values...
    signal += Trivial_TagPdf( iTag_os, NamePrefix = 'os_tag_' ).pdf()

#####################################
### Angular acceptance correction ###
#####################################
from P2VVGeneralUtils import RealMomentsBuilder
from itertools import chain
#momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))
#These are the relevant terms as found with MinSignificance>3
momindices = [(0,0,0),(0,2,0),(0,2,2),(2,0,0)]

eff = RealMomentsBuilder()
#Don't specify pdf and normset here, we're gonna read moments and not calculate any.
eff.appendPYList( angles.angles, momindices)
eff.read(path+'effmoments_tcut_%s.txt'%(tmincut))
eff.Print()

#Build Angular acceptance corrected PDF
sig_t_angles_iTag = eff * sig_t_angles_iTag

##############################
### Proper time acceptance ###
##############################
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
acceptance = Moriond2012_TimeAcceptance( time = t, Input = path +'BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root', Histogram = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_20bins')
#sig_t_angles_iTag = acceptance * sig_t_angles_iTag

# register with the signal component
signal += sig_t_angles_iTag

############
### SFIT ###
############
# make sweighted dataset. TODO: use mumu mass as well...
from P2VVGeneralUtils import SData, splot

pdf = buildPdf((signal,), Observables = angles.angles.values()+[t,iTag_os] , Name='pdf')
if fpsi==0 :
    masspdf = buildPdf((signal,background), Observables = (m,), Name = 'masspdf')
else :
    masspdf = buildPdf((signal,psi,background), Observables = (m,mpsi), Name = 'mass2pdf')

masspdf.fitTo(data,**fitOpts)


massplot = True

if massplot:
    from ROOT import kDashed, kRed, kGreen, kBlue, TCanvas, TLatex
    from P2VVGeneralUtils import plot
    obs = masspdf.Observables() 
    canvas = TCanvas()
    for pad,o in zip( canvas.pads( len(obs) ) , obs ) :
        plot( pad, o, data, masspdf, components = { 'sgnl*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen )
                                                  , 'bkg*'  : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed   ) 
                                                  , 'psi*'  : dict( LineStyle = kDashed, LineWidth=3, LineColor = kBlue  ) 
                                                  }
            , pdfOpts   = dict( LineWidth = 3 )
            , plotResidHist = True
            , dataOpts  = { 'MarkerSize' : 0.9,      'XErrorSize' : 0  }
            , frameOpts = dict( Object = ( TLatex(0.55,.8,"#splitline{LHCb preliminary}{#sqrt{s} = 7 TeV, L = 1.03 fb^{-1}}", NDC = True), )
                              #, Bins  = 50
                              , Title = ""
                              )

            )

for p in masspdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
splot_m = SData(Pdf = masspdf, Data = data, Name = 'MassSplot')


if True :
    # TODO: one canvas, overlay all components, normalized area, different color
    mcol = { 'sgnl' : kGreen, 'cmb' : kRed, 'psi' : kBlue }
    scanvas = TCanvas('splot_obs')
    obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']
            # ,iTag_os ,eta_os
            ,t  ,st
          ]
    from P2VVGeneralUtils import Sorter
    for p,o in zip( scanvas.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        print o.GetName()
        fr = o.frame( )
        for comp in splot_m.components() :
            if comp == 'sgnl' : continue
            hpdf = HistPdf( Name = 'hpdf_%s_%s' % (o.GetName(),comp), Observables = (o,), Data = splot_m.data(comp) )
            hpdf.plotOn(fr, LineColor = mcol[comp ] )
        fr.Draw( pad = p)
    scanvas.Print('%s_splot.pdf'% (fitname),'pdf')

assert False

def search(fname,path) :
    import os
    for f in ( os.path.join(p,fname) for p in os.path.split(os.pathsep) ) :
        if os.path.exists(f) : return f
    return None

import os
paramfile = search(fitname+'_params.txt',os.pathsep.join(['.','FitScripts']) )
if paramfile :
    print 'Reading fit result from %s' % paramfile
    fitset = pdf.getParameters(data)
    fitset.readFromFile(paramfile)

fit = True
if fit or not paramfile:
    sfitresult = pdf.fitTo( splot_m.data('sgnl'), SumW2Error = False, **fitOpts)
    sfitresult.Print()
    sfitresult.writepars(fitname+'result',False)
    fitset = pdf.getParameters(data)
    fitset.writeToFile(fitname)

########
# PLOT #
########

from P2VVGeneralUtils import plot

obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']
      # ,iTag_os ,eta_os
      ,t # ,st
      ]

canvas = dict()
for rng in ( None, ) :
    canvas[rng] = TCanvas('%s'%rng)
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( canvas[rng].pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        dataOpts = dict( CutRange =        rng ) if rng else dict()
        pdfOpts  = dict( ProjectionRange = rng ) if rng else dict()
        from P2VVGeneralUtils import plot
        from ROOT import RooArgSet
        pdfOpts[ 'ProjWData' ] = ( RooArgSet(st._var, eta_os._var), data, True )
        plot( p, o, splot_m.data('sgnl'), pdf
              , yTitleOffset = 1.5
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
              , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
              , plotResidHist = True
              , frameOpts = dict( Object = ( TLatex(0.15,.3,"#splitline{LHCb preliminary}{#sqrt{s} = 7 TeV, L = 1.03 fb^{-1}}", NDC = True), )
                                , Title = ""
                                )
              #Error for events with negative weight for negative log
              , logy = ( o == t )
              )
    canvas[rng].Print('%s_%s_fitproj.pdf'% (fitname,rng),'pdf')
assert False



#Turn this on when fit is fast with NumCPU = numCPU() working!!!
#######################
# Profile likelihoods #
#######################

pllvar = lifetimeParams['dM']

from ROOT import RooMinuit
#Need to implement conditionalobservables and externalconstraints here

nll = pdf.createNLL(splot_m.data('signal'), NumCPU = numCPU()
                    , Verbose = False
                    #, ConditionalObservables =
                    #, ExternalConstraints = 
                    )
minuit = RooMinuit(nll)
minuit.setProfile(1) #Timer =1
minimize = minuit.migrad() #Save = True
sfitresult = minuit.save()
pll = nll.createProfile(RooArgSet(pllvar)) #This creates a RooProfileLL object!
from ROOT import TCanvas
canvas = TCanvas()
frame = pllvar.frame()
pll.plotOn(frame)
frame.Draw()

