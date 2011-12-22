from RooFitWrappers import *
from ROOT import TFile

# setup (singleton) workspace
ws = RooObject( workspace = 'myworkspace' )

# now (consistently!) create/declare observables
m = RealVar('data_m',Observable=True,Unit='MeV/c^2',MinMax=(5000,6000))
t = RealVar('data_t',Observable=True,MinMax=(-1,14),Unit='ps',Value=0)
#c = Category('tagdecision',{+1: 'B', -1 : 'Bbar'}, Value='B',  Observable=True )
#tau = RealVar('sig_tau',Observable=False,Blinded=('UnblindUniform','blindingString', 1.),MinMax=(1,2) )

res_mean = RealVar( 'res_mean',    Unit = 'ps',   Value = 0,  MinMax = ( -10, 10 ) )
res_sigma = RealVar( 'res_sigma',  Unit = '1/ps', Value = 50, MinMax = (  20, 60 ) )

sig_res = ResolutionModel( 'sig_res', Type = 'RooGaussModel', Parameters = [ t, res_mean, res_sigma ] )
sig_tau = RealVar( 'sig_tau',  Unit = 'ps', Value = 1.5, MinMax = ( 1., 2. ) )
sig_t = Pdf( 'time', Type = 'Decay', Parameters = ( t, sig_tau, sig_res, 'SingleSided') )

mass_mean = RealVar( 'mass_mean',  Unit = 'MeV', Value = 5300, MinMax = ( 5200, 5800 ) )
mass_sigma = RealVar( 'mass_sigma',  Unit = 'MeV', Value = 15, MinMax = ( 10, 20 ) )
sig_m = Pdf( 'mass', Type = 'Gaussian',  Parameters = ( m, mass_mean, mass_sigma ) )

# create signal and background
signal = Component('signal',  Yield=(3000,100,6000) )
signal += sig_m
signal += sig_t

background = Component('background', Yield = (3000,1000,6000) )

background_c = RealVar( 'background_c',  Unit = '1/MeV', Value = -0.0004)
background += Pdf( 'background_m',  Type = 'Exponential', Parameters = ( m, background_c, ) )

background_tau = RealVar( 'background_tau',  Unit = 'ps', Value = 0.4, MinMax = ( 0.1, 0.9 ) )
background_res = ResolutionModel( 'background_res', Type = 'RooTruthModel', Parameters = [ t ] )

background += Pdf( 'background_t', Type = 'Decay', Parameters = ( t, background_tau, background_res , 'SingleSided', ) )

pdf = buildPdf( (background,signal) , Observables = (m,t), Name='pdf' )

##########################################

data = pdf.generate((m,t), 2000)

rootFile = TFile.Open("data.root", "recreate")
rootFile.WriteTObject(data, "data")
rootFile.Close()
