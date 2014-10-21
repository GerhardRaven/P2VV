###########################################################################################################################################
## script parameters ##
#######################

numEvents     = 10000  # -1 to read dataset from file
dataFilePath  = 'testBTagDecay.root'
plotsFilePath = 'testBTagDecay.pdf'

phis         = 0.8
lambdas      = 1.
wrongTagProb = 0.1


###########################################################################################################################################
## build simplified B_s^0 -> J/psi phi signal PDF ##
####################################################

# dictionary for observables
observables = { }

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# B lifetime
from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
lifetimeParams = LifetimeParams()

# time resolution
from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
resModel = TimeResolution( time = dict( Value = 0., MinMax = ( 0., 5. ) ) )
observables['time'] = resModel.parameter('time')

# CP violation parameters
from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
lambdaCP = CPParam()
lambdaCP.parameter('phiCP').setVal(phis)
lambdaCP.parameter('lambdaCP').setVal(lambdas)

# angular functions
from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
angleFuncs = AngleFuncs()
for obsName in [ 'cpsi', 'ctheta', 'phi' ] : observables[obsName] = angleFuncs.angles[obsName]

# decay amplitudes
ampNames = [ 'A0', 'Apar', 'Aperp', 'AS' ]
from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
amplitudes = Amplitudes()
amplitudes.parameter('C_SP').setConstant(True)

# coefficients for time functions
from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, ampNames )

# tagging observables and parameters
from P2VV.RooFitWrappers import Category
observables['iTag'] = Category( Name = 'iTag', Title = 'flavour tag', States = { 'B' : +1, 'Bbar' : -1 } )
from P2VV.Parameterizations.FlavourTagging import WTag_TaggingParams as TaggingParams
taggingParams = TaggingParams()
taggingParams.parameter('wTag').setVal(wrongTagProb)

# build signal PDF
pdfArgs = dict(  time            = observables['time']
               , iTag            = observables['iTag']
               , tau             = lifetimeParams['MeanLifetime']
               , dGamma          = lifetimeParams['dGamma']
               , dm              = lifetimeParams['dM']
               , dilution        = taggingParams['dilution']
               , ADilWTag        = taggingParams['ADilWTag']
               , avgCEven        = taggingParams['avgCEven']
               , avgCOdd         = taggingParams['avgCOdd']
               , coshCoef        = timeBasisCoefs['cosh']
               , sinhCoef        = timeBasisCoefs['sinh']
               , cosCoef         = timeBasisCoefs['cos']
               , sinCoef         = timeBasisCoefs['sin']
               , resolutionModel = resModel['model']
              )
from P2VV.RooFitWrappers import BTagDecay
pdf = BTagDecay( 'sig_t_angles', **pdfArgs )

# collect python garbage
import gc
gc.collect()

# print PDF
pdf.Print()
pdf.getVariables().Print('v')


###########################################################################################################################################
## generate and fit ##
######################

# set maximum value of PDF for generating angles with accept/reject (is 0.1 a good value?!!)
pdf.setMaxVal(0.1)

from ROOT import TFile
if numEvents > 0 :
    # generate data
    print 'generating %d events' % numEvents
    data = pdf.generate( [ observables[name] for name in [ 'time', 'iTag', 'cpsi', 'ctheta', 'phi' ] ], numEvents, Name = 'JpsiphiData' )
    dataFile = TFile.Open( dataFilePath, 'RECREATE' )
    dataFile.Append(data)
    from ROOT import TObject
    dataFile.Write( dataFilePath, TObject.kOverwrite )
    dataFile.Close()

else :
    # read data from file
    dataFile = TFile.Open(dataFilePath)
    data = dataFile.Get('JpsiphiData')
    print 'read %d events from file' % data.numEntries()
    dataFile.Close()

# fit data
print 'fitting data'
fitResult = pdf.fitTo( data, Minimizer = 'Minuit2', NumCPU = 8, Timer = True, Save = True, Offset = True )
print 'fit results:'
fitResult.PrintSpecial( text = True, LaTeX = False, normal = False )


###########################################################################################################################################
## plot data and PDFs ##
########################

# set LHCb plot style
from P2VV.Load import LHCbStyle
from ROOT import kFullDotLarge

# create plots
frames = [ observables['time'].frame(30), observables['time'].frame(60) ]
data.plotOn( frames[0], MarkerStyle = kFullDotLarge, MarkerSize = 0.7, LineWidth = 2 )
data.plotOn( frames[1], Asymmetry = observables['iTag'], MarkerStyle = kFullDotLarge, MarkerSize = 0.6, LineWidth = 2 )
pdf.plotOn( frames[0], LineWidth = 3 )
pdf.plotOn( frames[1], Asymmetry = observables['iTag'], LineWidth = 3 )

frames[0].GetYaxis().SetTitle( 'Candidates / (%.2g ps)' % frames[0].GetXaxis().GetBinWidth(1) )
frames[1].GetYaxis().SetTitle('B/#bar{B}-tag asymmetry')
frames[0].SetTitleOffset( 1.2, 'Y' )
frames[1].SetTitleOffset( 1.2, 'Y' )

# draw plots
from ROOT import TCanvas
canvs = [ ]
for it, frame in enumerate(frames) :
    canvs.append( TCanvas( 'canv%d' % it ) )
    canvs[-1].SetLeftMargin(0.18)
    canvs[-1].SetRightMargin(0.05)
    canvs[-1].SetBottomMargin(0.18)
    canvs[-1].SetTopMargin(0.05)
    frame.Draw()
    canvs[-1].Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len(frames) - 1 else '' ) )
