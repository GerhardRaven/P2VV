from RooFitWrappers import *

from ROOT import TFile

# setup (singleton) workspace
from ROOT import RooWorkspace
ws = RooObject()
ws.setWorkspace( RooWorkspace("myworkspace") )

# now (consistently!) create/declare observables
m = RealVar('data_m',Observable=True,Unit='MeV/c^2',MinMax=(5000,6000))
t = RealVar('data_t',Observable=True,MinMax=(-1,14),Unit='ps',Value=0)
#c = Category('tagdecision',{+1: 'B', -1 : 'Bbar'}, Value='B',  Observable=True )
#tau = RealVar('sig_tau',Observable=False,Blinded=('UnblindUniform','blindingString', 1.),MinMax=(1,2) )

res_mean = RealVar( 'res_mean',   Observable = False, Unit = 'ps',   Value = 0,  MinMax = ( -10, 10 ) )
res_sigma = RealVar( 'res_sigma', Observable = False, Unit = '1/ps', Value = 50, MinMax = (  20, 60 ) )

sig_res = ResolutionModel( 'sig_res', Type = 'RooGaussModel', Observables = [ t ],
                           Parameters = [ res_mean, res_sigma ] )

sig_tau = RealVar( 'sig_tau', Observable = False, Unit = 'ps', Value = 1.5, MinMax = ( 1., 2. ) )

sig_t = Pdf( 'time', Type = 'Decay', Observables = ( t, ), Parameters = ( sig_tau, ),
             Options = ( 'SingleSided', ), ResolutionModel = sig_res )

mass_mean = RealVar( 'mass_mean', Observable = False, Unit = 'MeV', Value = 5300, MinMax = ( 5200, 5800 ) )
mass_sigma = RealVar( 'mass_sigma', Observable = False, Unit = 'MeV', Value = 15, MinMax = ( 10, 20 ) )

sig_m = Pdf( 'mass', Type = 'Gaussian', Observables = ( m, ), Parameters = ( mass_mean, mass_sigma ) )

# create signal and background
signal = Component('signal')
signal.setYield(3000,100,6000)
## signal[m] = 'Gaussian(m,5300,15)'
signal[m] = sig_m
## signal[t] = 'Decay(t,sig_tau[1.5,1.0,2.0],TruthModel(t),SingleSided)'
signal[t] = sig_t

background = Component('background')
background.setYield(3000,1000,6000)

background_c = RealVar( 'background_c', Observable = False, Unit = '1/MeV', Value = -0.0004)
background[ m ] = Pdf( 'background', Observables = ( m, ), Type = 'Exponential',
                       Parameters = ( background_c, ) )

background_tau = RealVar( 'background_tau', Observable = False, Unit = 'ps', Value = 0.4, MinMax = ( 0.1, 0.9 ) )
background_res = ResolutionModel( 'background_res', Type = 'RooTruthModel', Observables = [ t ] )

#background[t] = 'Decay(t,bkg_tau[0.4,0.1,0.9],TruthModel(t),SingleSided)'
background[ t ] = Pdf( 'background', Type = 'Decay', Observables = ( t, ),
                       Parameters = ( background_tau, ), Options = ( 'SingleSided', ),
                       ResolutionModel = background_res )
# Exponential(m,-0.004)
#background[m,t] = 'PROD(Exponential(m,-0.004),Decay(t,bkg_tau[0.4,0.1,0.9],TruthModel(t),SingleSided))'

pdf = buildPdf( (background,signal) , Observables = (m,t), Name='pdf' )

##########################################

data = pdf.generate((m,t), 2000)

rootFile = TFile.Open("data.root", "recreate")
rootFile.WriteTObject(data, "data")
rootFile.Close()
