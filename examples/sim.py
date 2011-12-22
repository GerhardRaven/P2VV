from RooFitWrappers import *
from ROOT import RooDecay as Decay,             \
                 RooExponential as Exponential, \
                 RooGaussian as Gaussian,       \
                 RooTruthModel as TruthModel,   \
                 RooGaussModel as GaussModel

from ROOT import TFile

# setup (singleton) workspace
from ROOT import RooWorkspace
ws = RooObject( workspace = 'myworkspace' )

# now (consistently!) create/declare observables
m = RealVar('m',Observable=True,Unit='MeV/c^2',MinMax=(5000,6000))
t = RealVar('t',Observable=True,MinMax=(-1,14),Unit='ps',Value=0)
#c = Category('tagdecision',{+1: 'B', -1 : 'Bbar'}, Value='B',  Observable=True )
#tau = RealVar('sig_tau',Observable=False,Blinded=('UnblindUniform','blindingString', 1.),MinMax=(1,2) )

signal = Component('signal', Yield = (3000,100,6000) )

res_mean  = RealVar( 'res_mean',    Unit = 'ps',   Value = 0,  MinMax = ( -10, 10 ) )
res_sigma = RealVar( 'res_sigma',  Unit = '1/ps', Value = 50, MinMax = (  20, 60 ) )
sig_res   = ResolutionModel( 'sig_res', Type = GaussModel,  Parameters = [t, res_mean, res_sigma ] )
sig_tau   = RealVar( 'sig_tau',  Unit = 'ps', Value = 1.5, MinMax = ( 1., 2. ) )
signal   += Pdf( 'time', Type = Decay,  Parameters = (t, sig_tau, sig_res,  'SingleSided' ) )

mass_mean  = RealVar( 'mass_mean',   Unit = 'MeV', Value = 5300, MinMax = ( 5200, 5800 ) )
mass_sigma = RealVar( 'mass_sigma',  Unit = 'MeV', Value = 15, MinMax = ( 10, 20 ) )
signal    += Pdf( 'mass', Type = Gaussian,  Parameters = (m, mass_mean, mass_sigma ) )

# create signal and background



background = Component('background',  Yield= (3000,1000,6000))

background_c = RealVar( 'background_c',  Unit = '1/MeV', Value = -0.0004)
# TODO: auto mangle name??
background += Pdf( 'background_m',  Type = Exponential, Parameters = (m, background_c ) )

background_tau = RealVar( 'background_tau',  Unit = 'ps', Value = 0.4, MinMax = ( 0.1, 0.9 ) )
background_res = ResolutionModel( 'background_res', Type = TruthModel, Parameters = [ t ] )

# TODO: auto mangle name??
background += Pdf( 'background_t', Type = Decay,  Parameters = (t, background_tau, background_res, 'SingleSided' ) )

# background = Component('background', ( background_t, background_m ), Yield= (3000,1000,6000))

pdf = buildPdf( (background,signal) , Observables = (m,t), Name='pdf' )

##########################################

data = pdf.generate((m,t), 2000)

#rootFile = TFile.Open("data.root", "recreate")
#rootFile.WriteTObject(data, "data")
#rootFile.Close()

pdf.fitTo(data)

from P2VVGeneralUtils import plot
from ROOT import TCanvas, RooFit, kDashed
sigcolor = RooFit.kGreen
bkgcolor = RooFit.kRed
lw = 2
ms = 0.4
xe = 0
dashed = kDashed
c = TCanvas()
for ( cc, obs, logy ) in zip( c.pads( 1, 2 ), ( m, t ), ( False, True ) ) :
    plot(  cc.cd(), obs, data, pdf
         , components = { 'signal*'     : { 'LineColor' : sigcolor, 'LineStyle' : dashed }
                        , 'background*' : { 'LineColor' : bkgcolor, 'LineStyle' : dashed }
                        }
         , plotResidHist = True, logy = logy
         , dataOpts = { 'XErrorSize' : xe, 'MarkerSize' : ms }
        )

# create a continuous tagging variable...

# split in categories

