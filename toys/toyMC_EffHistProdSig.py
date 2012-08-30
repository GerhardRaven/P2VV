import sys
import os
from ToyMCUtils import Toy

toy = Toy()
parser = toy.parser()
parser.add_option("-a", dest = "acceptance", default = '',
                  type = 'string', action = 'store', help = 'use the acceptance')
parser.add_option("-p", dest = "pdf", default = 'sig_pdf',
                  action = 'store', help = 'which pdf to use')

(options, args) = toy.configure()

if options.pdf not in ['sig_pdf', 'comb_pdf']:
    print "PDF must be either sig_pdf or comb_pdf"
    sys.exit(-2)

    
from RooFitWrappers import *

w = RooObject( workspace = 'w' )
w = w.ws()

from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles = TrAngles( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
t      = RealVar('t', Title = 'decay time', Unit = 'ps', Observable = True, MinMax=(-3,14), nBins = 48)
iTag   = Category('tagdecision' , Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1 } ) # , 'untagged' : 0 } )
eta    = RealVar('eta', Title = 'estimated mis tag', Observable = False, Constant = True, Value = 0.3, MinMax=(0,0.5) )
mass   = RealVar('m', Observable = True, Unit = 'MeV/c^2', MinMax = (5200, 5550), nBins = 48)
observables = list( angles.angles.itervalues() ) + [ t,iTag, mass ]

for i in angles.angles.itervalues() : i.setBins(16)
mass.setRange('leftsideband', (mass.getMin(),5330) )
mass.setRange('signal',(5330,5410) )
mass.setRange('rightsideband',(5410,mass.getMax()) )

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP = { 'Name': 'HelloWorld', 'Value': -0.4, 'MinMax': (-3.2,3.2) }, lambdaCPSq = ConstVar(Name ='one',Value=1) )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiVPolar_AmplitudeSet
amplitudes = JpsiVPolar_AmplitudeSet( A0Mag2 = 0.60, A0Phase = 0
                                    , AperpMag2 = 0.160, AperpPhase = -0.17
                                    , AparPhase = 2.5
                                    , ASMag2 = { 'Value' : 0, 'Constant': True} , ASPhase = { 'Value': 0, 'Constant':True } )
#amplitudes.setConstant('.*AS.*',True)

#### Package from here until the "BTagDecay('name', args)" into a dedicated class/function...
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
# need to specify order in which to traverse...
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

#from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
#basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes,CP, iTag,  ['A0','Apar','Aperp','AS'] ) 

from P2VVParameterizations.FlavourTagging import WTag_TaggingParams
taggingParams = WTag_TaggingParams( wTag = eta ) # FormulaVar('wTag','@2 + @3*(@0-@1)',[eta,etaAverage,p0,p1] ) )

# now build the actual signal PDF...
from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.68, dGamma = 0.05, dM = dict( Value = 17.8, MinMax = (16,19), Constant = True) )

from P2VVParameterizations.TimeResolution import Truth_TimeResolution
args = { 'time'      : t
       , 'iTag'      : iTag
       , 'dm'        : lifetimeParams['dM']
       , 'tau'       : lifetimeParams['MeanLifetime']
       , 'dGamma'    : lifetimeParams['dGamma']
       , 'resolutionModel' : Truth_TimeResolution(time = t)['model']
       , 'coshCoef'  : basisCoefficients['cosh']
       , 'cosCoef'   : basisCoefficients['cos']
       , 'sinhCoef'  : basisCoefficients['sinh']
       , 'sinCoef'   : basisCoefficients['sin']
       , 'avgCEven'  : taggingParams['avgCEven'] 
       , 'avgCOdd'   : taggingParams['avgCOdd']
       , 'dilution'  : taggingParams['dilution']
       , 'ADilWTag'  : taggingParams['ADilWTag']
       }

# TODO: should be able to write BTagDecay('mypdf', **lifetimeParams.BTagDecay() + **basisCoefficients.BTagDecay() + **taggingParams.BTagDecay() )

# update resolution model, and build again...
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres =  TimeResolution(time = t)
tres.setConstant('.*')
args[ 'resolutionModel' ]  = tres.model()

sig_pdf = BTagDecay('sig_pdf', **args )

from ROOT import TFile
if options.acceptance:
    acc_file = TFile.Open(options.acceptance)
    _hist = acc_file.Get('eff_hist')
    eff_func = HistFunc('acceptance', Observables = [t], Histogram = _hist)
    sig_pdf = eff_func * sig_pdf

from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass
signal = Component('signal', (LP2011_Signal_Mass(mass = mass).pdf(), sig_pdf), Yield = (1000,0,15000))

from P2VVParameterizations.MassPDFs import LP2011_Background_Mass
from P2VVParameterizations.TimePDFs import LP2011_Background_Time
from P2VVParameterizations.FlavourTagging import Trivial_TagPdf
from P2VVParameterizations.AngularPDFs import Uniform_Angles
bkg  = Component('bkg',(  LP2011_Background_Mass( mass = mass ).pdf()
                       ,  LP2011_Background_Time( time = t , resolutionModel = tres.model()).pdf()
                       ,  Trivial_TagPdf( tagdecision = iTag, ATagEff = 0.3, NamePF = 'bkg' ).pdf()
                       ,  Uniform_Angles( angles = angles.angles ).pdf()
                       ), Yield = (4000,1000,15000) )

comb_pdf = buildPdf((signal, bkg), Observables = observables,  Name = 'comb_pdf')

w.addClassDeclImportDir("..")
w.addClassImplImportDir("..")
w.importClassCode()

# Get the PDF which we want to use
pdf = locals()[options.pdf]

obs = w.argSet(','.join([o['Name'] for o in observables]))

toy.run(Observables = obs, Pdf = pdf)

toy.write_output()
