from RooFitWrappers import *
from ROOT import gSystem
gSystem.Load('libP2VV.so')

from ROOT import RooMsgService

RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

ws = RooObject()
ws.setWorkspace( RooWorkspace('myworkspace') )
ws.ws().importClassCode("RooBTagDecay",True);
ws.ws().importClassCode("RooP2VVAngleBasis",True);

from math import pi
t = RealVar('t', Title = 'decay time', Unit='ps', Observable = True, MinMax=(-5, 14))
m = RealVar('m', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5155, 5450))
mpsi = RealVar('mpsi', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3170))
observables = [t, m, mpsi]

# now build the actual signal PDF...
args = { 'tau'       : RealVar( 't_sig_tau', Title = 'mean lifetime', Unit = 'ps',       Value =  1.5,  MinMax = (1.3, 1.8) )
       , 'resolutionModel' : ResolutionModel( 'resModel', Type = TruthModel, Observables = [ t ] )
       , 'decayType' : 'SingleSided' 
       , 'time'      : t } 

mcpdf = Decay('mc_pdf',  args)

# update resolution model, and build again...
from parameterizations import ResolutionModelLP2011
tres = ResolutionModelLP2011(t).Model
args[ 'resolutionModel' ]  = tres

sig_tau = Decay('B_tau', args)

m_mean  = RealVar('m_mean',  Observable = False, Unit = 'MeV', Value = 5300, MinMax = (5200, 5800))
m_sigma = RealVar('m_sigma', Observable = False, Unit = 'MeV', Value = 15, MinMax = (10, 20))
sig_m = Pdf('B_mass', Type = Gaussian, Observables = (m,), Parameters = (m_mean, m_sigma ))

mpsi_mean  = RealVar('mpsi_mean',  Observable = False, Unit = 'MeV', Value = 3097, MinMax = (3070, 3110))
mpsi_sigma = RealVar('mpsi_sigma', Observable = False, Unit = 'MeV', Value = 10, MinMax = (5, 20))
sig_mpsi = Pdf('Jpsi_mass', Type = Gaussian, Observables = (mpsi,), Parameters = (msi_mean, mpsi_sigma))

# Create signal component
signal = Component('signal')
signal.setYield(10000,100,15000)
signal[m] = sig_m
signal[mpsi] = sig_mpsi
signal[t] = sig_tau

# Create background component
background = Component('background')
background.setYield(5000,1000,15000)
background_m_c = RealVar( 'background_m_c', Observable = False, Unit = '1/MeV', Value = -0.0004, MinMax(-0.1, -0.00001))
background_mpsi_c = RealVar( 'background_mpsi_c', Observable = False, Unit = '1/MeV', Value = -0.0004, MinMax(-0.1, -0.00001))
background[m] = Pdf('background_m', Observables = (m,), Type = Exponential, Parameters = (background_m_c,))
background[mpsi] = Pdf('background_mpsi', Observables = (m,), Type = Exponential, Parameters = (background_mpsi_c,))

background_tau = RealVar('background_tau', Observable = False, Unit = 'ps', Value = 0.4, MinMax = ( 0.1, 0.9 ) )

#background[t] = 'Decay(t,bkg_tau[0.4,0.1,0.9],TruthModel(t),SingleSided)'
background[t] = Pdf('background', Type = Decay, Observables = (t,), Parameters = (background_tau, tres, 'SingleSided'))

pdf = buildPdf((background,signal), observables = (m,mpsi,t), name='pdf')

#l = RooArgSet()
#pdf.branchNodeServerList( l )
#for i in l : i.setAttribute( "CacheAndTrack" )

pdf.Print("t")

print 'generating data'
data = pdf.generate(observables , 10000)

print 'fitting data'
from ROOT import RooCmdArg
NumCPU = RooCmdArg(RooFit.NumCPU(4))
pdf.fitTo(data, NumCPU, RooFit.Timer(1))
