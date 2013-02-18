from P2VV.RooFitWrappers import *
from ROOT import RooDecay as Decay,             \
                 RooExponential as Exponential, \
                 RooGaussian as Gaussian,       \
                 RooTruthModel as TruthModel,   \
                 RooGaussModel as GaussModel


# setup (singleton) workspace
from ROOT import RooWorkspace
ws = RooObject( workspace = 'myworkspace' )

# now (consistently!) create/declare observables
m = RealVar(Name = 'm',Observable=True,Unit='MeV/c^2',MinMax=(5000,6000))
t = RealVar(Name = 't',Observable=True,MinMax=(-1,14),Unit='ps')

signal     = Component( Name = 'signal',      Yield = (3000, 100,6000) )
background = Component( Name = 'background',  Yield = (3000,1000,6000) )

res_mean  = RealVar( Name = 'res_mean',   Unit = 'ps',   Value = 0,  MinMax = ( -10, 10 ) )
res_sigma = RealVar( Name = 'res_sigma',  Unit = '1/ps', Value = 50, MinMax = (  20, 60 ) )
sig_res   = ResolutionModel( Name = 'sig_res', Type = GaussModel,  Parameters = [t, res_mean, res_sigma ] )
sig_tau   = RealVar( Name = 'sig_tau',    Unit = 'ps', Value = 1.5, MinMax = ( 1., 2. ) )
signal   += Pdf( Name = 'time', Type = Decay,  Parameters = (t, sig_tau, sig_res,  'SingleSided' ) )

mass_mean  = RealVar( Name = 'mass_mean',  Unit = 'MeV', Value = 5300, MinMax = ( 5200, 5800 ) )
mass_sigma = RealVar( Name = 'mass_sigma', Unit = 'MeV', Value = 15,   MinMax = (   10,   20 ) )
signal    += Pdf( Name = 'mass', Type = Gaussian,  Parameters = (m, mass_mean, mass_sigma ) )


background_c = RealVar( Name = 'background_c',  Unit = '1/MeV', Value = -0.0004)
background += Pdf( Name = 'background_m',  Type = Exponential, Parameters = (m, background_c ) )

background_tau = RealVar( Name =  'background_tau',  Unit = 'ps', Value = 0.4, MinMax = ( 0.1, 0.9 ) )
background_res = ResolutionModel( Name = 'background_res', Type = TruthModel, Parameters = [ t ] )
background += Pdf( Name = 'background_t', Type = Decay,  Parameters = (t, background_tau, background_res, 'SingleSided' ) )


pdf = buildPdf( (background,signal) , Observables = (m,t), Name='pdf' )

##########################################

data = pdf.generate((m,t), 2000)

#from ROOT import TFile
#rootFile = TFile.Open("data.root", "recreate")
#rootFile.WriteTObject(data, "data")
#rootFile.Close()

pdf.fitTo(data)

from P2VV.GeneralUtils import plot
from ROOT import TCanvas, RooFit, kDashed, kGreen, kRed
( sigcolor, bkgcolor) = (kGreen,kRed )
(lw,ms,xe) = (2,0.4,0)
c = TCanvas()
for ( cc, obs, logy ) in zip( c.pads( 1, 2 ), ( m, t ), ( False, True ) ) :
    plot(  cc.cd(), obs, data, pdf
         , components = { 'signal*'     : { 'LineColor' : sigcolor, 'LineStyle' : kDashed }
                        , 'background*' : { 'LineColor' : bkgcolor, 'LineStyle' : kDashed }
                        }
         , plotResidHist = True, logy = logy
         , dataOpts = { 'XErrorSize' : xe, 'MarkerSize' : ms }
        )
