###############################################################################
## compAngEfficiencyReal:                                                    ##
##   P2VV example script for computing efficiency moments in real data       ##
##                                                                           ##
## * decay channels: B0->J/psiK* or B_s0->J/psiphi                           ##
## * writes moments to ascii file                                            ##
## * assumes that $P2VVROOT/python is in $PYTHONPATH                         ##
## * assumes that $P2VVROOT/lib is in $LD_LIBRARYPATH                        ##
##                                                                           ##
## authors:                                                                  ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                ##
##                                                                           ##
###############################################################################

# specify decay mode ('Bd2JpsiKstar' or 'Bs2Jpsiphi')
#mode = 'Bd2JpsiKstar'
mode = 'Bs2Jpsiphi'

# efficiency moments file path
momentsFilePath = 'effMoments'

# data set name and file
dataSetName  = 'MyTree'
dataFilePath = '/data/bfys/dveijk/MC/ReducedMCNTuple.root'

###############################################################################
import P2VV, P2VVConfiguration, P2VVModelBuilders
from ROOT import RooFit, TCanvas
from math import sqrt, sin, cos

# load the P2VV library
P2VV.loadP2VVLib()

# create P2VV configuration object
config = P2VVConfiguration.getP2VVConfig(mode, ['onlySignal',
    'effType=angular'])

# custom RooFit variable settings
if config.value('anglesType')[0] == 'trans' :
  config['cpsiAng'].set(name = 'trcospsi')
  config['cthetaAng'].set(name = 'trcostheta')
  config['phiAng'].set(name = 'trphi')
else :
  config['cpsiAng'].set(name = 'helcosthetaK')
  config['cthetaAng'].set(name = 'helcosthetaL')
  config['phiAng'].set(name = 'helphi')

config['BLifetime'].set(name = 't', min = -2., max = 20.)
config['iTag'].set(name = 'tagdecision')
if mode == 'Bd2JpsiKstar' :
  config['fTag'].set(name = 'qrec')

# allow only positive values for the true lifetime: truth matched events
config.addSetting('trueBLifetime', P2VVConfiguration.RooRealSetting('TRUEt',
    'true B lifetime (ps)', 'par', 0., 0., 20.))

if mode == 'Bd2JpsiKstar' :
  ReHp = 0.159 * cos(1.563) / 0.775
  ImHp = 0.159 * sin(1.563) / 0.775
  ReHm = 0.612 * cos(2.712) / 0.775
  ImHm = 0.612 * sin(2.712) / 0.775
  config['ReApar'].set(val  = (ReHp + ReHm) / sqrt(2.))
  config['ImApar'].set(val  = (ImHp + ImHm) / sqrt(2.))
  config['ReAperp'].set(val = (ReHp - ReHm) / sqrt(2.))
  config['ImAperp'].set(val = (ImHp - ImHm) / sqrt(2.))
  config['Gamma'].set(val = 0.655737)
  config['dGamma'].set(val = 0.)
  config['dm'].set(val = 0.507)
  config['lambdaCPSq'].set(val = 1.)
  config['phiCP'].set(val = 0.)
elif mode == 'Bs2Jpsiphi' :
  config['ReApar'].set(val  = 0.49 * cos( 2.5)  / 0.775)
  config['ImApar'].set(val  = 0.49 * sin( 2.5)  / 0.775)
  config['ReAperp'].set(val = 0.40 * cos(-0.17) / 0.775)
  config['ImAperp'].set(val = 0.40 * sin(-0.17) / 0.775)
  config['Gamma'].set(val = 0.679348)
  config['dGamma'].set(val = 0.0599979)
  config['dm'].set(val = 17.8)
  config['lambdaCPSq'].set(val = 1.)
  config['phiCP'].set(val = -0.04)

# declare RooFit variables and store them in RooWorkspace
config.declareRooVars()

# get workspace
ws = config.workspace()

# build the PDF
pdf = P2VVModelBuilders.getP2VVPDF(config)

# print contents of RooWorkspace to screen
#config.workspace().Print()

# create data set from NTuple file(s)
dataObs = ws.set('angles').clone('dataObs')
dataObs.add(ws.arg(config['trueBLifetime'].name()))
data = P2VV.readData(dataFilePath, dataSetName, True, dataObs)

print 'compAngEfficiencyReal: %d events in data set' % data.numEntries()

# get sets of observables
angles = ws.set('angles')
marginalObs = ws.set('observables').clone('margObs')
marginalObs.remove(angles)
dataObs = data.get()

print 'compAngEfficiencyReal: angles:'
angles.Print()
print 'compAngEfficiencyReal: marginal observables:'
marginalObs.Print()
print 'compAngEfficiencyReal: observables in data set:'
dataObs.Print()

# integrate over all variables but the angles
angPDF = pdf.createProjection(marginalObs)

# compute efficiency moments
effBuilder = config.modelBuilder('efficiency')
effBuilder.buildEffBasis()
effBuilder.computeEffMoments(data, angPDF, dataObs)

# print efficiency moments
effBuilder.printEffMoments()
effBuilder.writeEffMoments(momentsFilePath)

