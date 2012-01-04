import optparse
import sys
import os

parser = optparse.OptionParser(usage = '%prog')

parser.add_option("-o", "--output", dest = "output", default = 'toy.root',
                  type = 'string', help = "set output filename")
parser.add_option("-n", "--ntoys", dest = "ntoys", default = 100,
                  type = 'int', help = 'number of toys to run')
parser.add_option("--ncpu", dest = "ncpu", default = 4,
                  type = 'int', help = 'number of CPUs to use')
parser.add_option("-e", "--nevents", dest = "nevents", default = 10000,
                  type = 'int', help = 'number of events to generate')
parser.add_option("-a", dest = "acceptance", default = '',
                  type = 'string', action = 'store', help = 'use the acceptance')
parser.add_option("-p", dest = "pdf", default = 'comb_pdf',
                  action = 'store', help = 'which pdf to use')
parser.add_option("-s", "--snapshot", dest = "snapshot", default = '',
                  action = 'store', help = 'Extract a snapshot to current directory.')

(options, args) = parser.parse_args()

if options.pdf not in ['sig_pdf', 'comb_pdf']:
    print "PDF must be either sig_pdf or comb_pdf"
    sys.exit(-2)

if options.snapshot:
    import tarfile
    with tarfile.open(options.snapshot, 'r:bz2') as archive:
        python_dirs = []
        for member in archive.getmembers():
            if member.isfile() and os.path.exists(member.path):
                print "File %s already exists, skipping" % member.path
            else:
                archive.extract(member)
            if member.isdir() and member.path.endswith('python'):
                python_dirs.append(member.path)
        sys.path.extend(python_dirs)
    try:
        print 'Running with snapshot %(comment)s' % archive.pax_headers
    except KeyError:
        pass
    
from ROOT import RooDataHist, RooHistFunc
from RooFitWrappers import *
from ROOT import TH1F, TFile

w = RooObject( workspace = 'w' )
w = w.ws()

from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles = TrAngles( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
t      = RealVar('t', Title = 'decay time', Unit = 'ps', Observable = True, MinMax=(-3,14), nBins = 48)
iTag   = Category('tagdecision' , Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1 } ) # , 'untagged' : 0 } )
eta    = RealVar('eta', Title = 'estimated mis tag', Observable = False, Constant = True, Value = 0.3, MinMax=(0,0.5) )
mass   = RealVar('m', Observable = True, Unit = 'MeV/c^2', MinMax = (5200, 5550), nBins = 48)
observables = [ i for i in angles.angles.itervalues() ] + [ t,iTag, mass ]

for i in angles.angles.itervalues() : i.setBins(16)
mass.setRange('leftsideband', (mass.getMin(),5330) )
mass.setRange('signal',(5330,5410) )
mass.setRange('rightsideband',(5410,mass.getMax()) )

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP = { 'Name': 'HelloWorld', 'Value': -0.4, 'MinMax': (-3.2,3.2) }, lambdaCPSq = ConstVar(Name ='one',Value=1) )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011( A0Mag2 = 0.60, A0Phase = 0
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

from P2VVParameterizations.FlavourTagging import Trivial_TaggingParams
taggingParams = Trivial_TaggingParams( wTag = eta ) # FormulaVar('wTag','@2 + @3*(@0-@1)',[eta,etaAverage,p0,p1] ) )

# now build the actual signal PDF...
from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.68, deltaGamma = 0.05, deltaM = dict( Value = 17.8, MinMax = (16,19), Constant = True) )

from P2VVParameterizations.TimeResolution import Truth_TimeResolution
args = { 'time'      : t
       , 'iTag'      : iTag
       , 'dm'        : lifetimeParams['deltaM'] 
       , 'tau'       : lifetimeParams['MeanLifetime']
       , 'dGamma'    : lifetimeParams['deltaGamma'] 
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

from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass
signal = Component('signal', (LP2011_Signal_Mass(mass = mass).pdf(), sig_pdf), Yield = (1000,0,15000))

from P2VVParameterizations.MassPDFs import LP2011_Background_Mass
from P2VVParameterizations.TimePDFs import LP2011_Background_Time
from P2VVParameterizations.FlavourTagging import Trivial_Background_Tag
from P2VVParameterizations.AngularPDFs import Uniform_Angles
bkg  = Component('bkg',(  LP2011_Background_Mass( mass = mass ).pdf()
                       ,  LP2011_Background_Time( time = t , resolutionModel = tres.model()).pdf()
                       ,  Trivial_Background_Tag( tagdecision = iTag, bkg_tag_delta = 0.3 ).pdf()
                       ,  Uniform_Angles( angles = angles.angles ).pdf()
                       ), Yield = (4000,1000,15000) )

comb_pdf = buildPdf((signal, bkg), Observables = observables,  Name = 'comb_pdf')

##############################################
### Define acceptance function a la Wouter ###
##############################################
if options.acceptance:
    acc_file = TFile.Open(options.acceptance)
    acc_workspace = acc_file.Get('w')
    eff_hist = acc_workspace.data('eff_hist')
    eff_func = RooHistFunc('acceptance', 'time acceptance', RooArgSet(t._target_()), eff_hist)
    w.put(eff_func)
    w.factory("EffHistProd::acc_pdf(%s, acceptance)" % options.pdf)

w.addClassDeclImportDir("..")
w.addClassImplImportDir("..")
w.importClassCode()

# Get the PDF which included the acceptance
pdf = None
if options.acceptance:
    pdf= w.pdf('acc_pdf')
else:
    pdf = w.pdf(options.pdf)

# Get the bare observable objects and copy them
obs = w.argSet(','.join([o['Name'] for o in observables]))
pdf_params = pdf.getParameters(obs)
## for param in pdf_params:
##     if param.GetName() not in ['Gamma', 'deltaGamma']:
##         param.setConstant()
gen_params = pdf_params.snapshot(True)

# Make another ArgSet to put the fit results in
result_params = RooArgSet(pdf_params, "result_params")

# Some extra numbers of interest
NLL = RooRealVar('NLL', '-log(Likelihood', 1.)
ngen = RooRealVar('ngen', 'number of generated events', options.nevents)
seed = RooRealVar('seed', 'random seed', 0.)
result_params.add(NLL)
result_params.add(ngen)
result_params.add(seed)

# The dataset to store the results
result_data = RooDataSet('result_data', 'result_data', result_params)
data_params = result_data.get()

from ROOT import RooRandom
import struct, os

for i in range(options.ntoys):
    # Get a good random seed, set it and store it
    s = struct.unpack('I', os.urandom(4))[0]    
    RooRandom.randomGenerator().SetSeed(s)
    seed.setVal(s)

    # Reset pdf parameters to initial values. Note: this does not reset the estimated errors...
    pdf_params.assignValueOnly( gen_params ) 
    data = pdf.generate(obs, options.nevents)
    from ROOTDecorators import  ROOTversion as Rv
    fit_result = pdf.fitTo(data, RooFit.NumCPU(options.ncpu), RooFit.Save(True),
                           RooFit.Optimize( True if Rv[1]<32 else 0 ),
                           RooFit.Minos(False), RooFit.Minimizer('Minuit2'))
    if fit_result.status() != 0:
        print 'Fit result status = %s' % fit_result.status()
        continue
    NLL.setVal(fit_result.minNll())
    for result_param in result_params:
        data_param = data_params.find(result_param.GetName())
        data_param.setVal(result_param.getVal())
        # This sets a symmetric error, but since we don't run Minos, that's ok
        data_param.setError(result_param.getError())
    result_data.fill()

# Write the results to a file
output_file = TFile.Open(options.output, 'recreate')
output_file.WriteTObject(result_data, result_data.GetName())
output_file.WriteTObject(gen_params, 'gen_params')
output_file.Close()

## from ROOT import RooFit, RooGenFitStudy
## from ROOT import RooStudyManager

## gfs = RooGenFitStudy();
## obs = ','.join([o.GetName() for o in observables])
## gfs.setGenConfig("pdf", obs, RooFit.NumEvents(1000))
## gfs.setFitConfig("pdf", obs, RooFit.Minimizer("Minuit2"))

## mgr = RooStudyManager(w, gfs)


