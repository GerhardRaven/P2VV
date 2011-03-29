###############################################################################
## compAngEfficiencyToy:                                                     ##
##   P2VV example script for computing efficiency moments in toy or EvtGen   ##
##   data (efficiency = 1)                                                   ##
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

import P2VV, P2VVConfiguration, P2VVModelBuilders, P2VVPlots
from ROOT import RooDataSet, RooFit, TCanvas, TChain, TFile
from math import sqrt, sin, cos

###############################################################################
# specify decay mode ('Bd2JpsiKstar' or 'Bs2Jpsiphi')
mode = 'Bd2JpsiKstar'
#mode = 'Bs2Jpsiphi'

# efficiency moments file path
momentsFilePath = 'effMoments'

# data set name and file
dataSetName  = mode[3:] + 'Data'
dataFilePath = dataSetName + '.root'
#dataFilePath = '/data/bfys/jleerdam/Bd2JpsiKst/EvtGen/Gauss-11144000-*.root'
#dataFilePath = '/data/bfys/jleerdam/Bs2Jpsiphi/EvtGen/Gauss-13144001-*.root'

# generate events?
generate = True
nEvents = 50000

# read events from NTuple or RooDataset
NTuple = False

###############################################################################
# load the P2VV library
P2VV.loadP2VVLib()

# create P2VV configuration object
config = P2VVConfiguration.getP2VVConfig(mode, ['onlySignal', 'noKSWave'])

# custom settings
config['cpsiAng'].set(name = 'hel_cthetak')
config['cthetaAng'].set(name = 'hel_cthetal')
config['phiAng'].set(name = 'hel_phi')
config['BLifetime'].set(name = 't', min = 0., max = 4.)
config['iTag'].set(name = 'tagInitial')
if mode == 'Bd2JpsiKstar' :
  config['fTag'].set(name = 'tagFinal')
config['misTag'].set(realType = 'par', val = 0., min = '', max = '')

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
elif mode == 'Bs2Jpsiphi' :
  config['ReApar'].set(val  = 0.49 * cos( 2.5)  / 0.775)
  config['ImApar'].set(val  = 0.49 * sin( 2.5)  / 0.775)
  config['ReAperp'].set(val = 0.40 * cos(-0.17) / 0.775)
  config['ImAperp'].set(val = 0.40 * sin(-0.17) / 0.775)
  config['Gamma'].set(val = 0.679348)
  config['dGamma'].set(val = 0.0599979)
  config['dm'].set(val = 17.8)
  config['phiCP'].set(val = -0.04)

# declare RooFit variables and store them in RooWorkspace
config.declareRooVars()

# get workspace
ws = config.workspace()

# build the B0->J/psiK* PDF
pdf = P2VVModelBuilders.getP2VVPDF(config)

# print contents of RooWorkspace to screen
#config.workspace().Print()

if generate :
  # generate events
  print 'compAngEfficiency: generating %d events' % nEvents
  if config.value('BDecayClass') == 'RooBDecay' : P2VV.registerMultiCatGen()
  data = pdf.generate(ws.set('observables'), nEvents)

  # write events to file
  print "compAngEfficiency: writing RooDataSet '%s' to file '%s'"\
      % (dataSetName, dataFilePath)
  file = TFile.Open(dataFilePath, 'RECREATE')
  data.Write(dataSetName)
  file.Close()

elif NTuple :
  # create data set from NTuple file(s)
  print "compAngEfficiency: reading NTuple(s) '%s' from file(s) '%s'"\
      % (dataSetName, dataFilePath)
  files = TChain(dataSetName)
  files.Add(dataFilePath)
  data = RooDataSet(dataSetName, dataSetName, files, ws.set('observables'))

else :
  # get data set from file
  print "compAngEfficiency: reading RooDataset '%s' from file '%s'"\
      % (dataSetName, dataFilePath)
  file = TFile.Open(dataFilePath, 'READ')
  data = file.Get(dataSetName)
  file.Close()

print 'compAngEfficiency: %d events in data set' % data.numEntries()

# adjust efficiency settings
config['effBasisType'].setValue('angular')
config['angEffBasisFuncs'].setValue((4, 4))

# get sets of observables
angles = ws.set('angles')
marginalObs = ws.set('observables').clone('margObs')
marginalObs.remove(angles)
dataObs = data.get()

print 'compAngEfficiency: angles:'
angles.Print()
print 'compAngEfficiency: marginal observables:'
marginalObs.Print()
print 'compAngEfficiency: observables in data set:'
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

