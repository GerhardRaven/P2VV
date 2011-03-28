###############################################################################
## compAngEfficiency: P2VV example script for computing efficiency moments   ##
##                                                                           ##
## decay channels: B0->J/psiK* or B_s0->J/psiphi                             ##
##                                                                           ##
## * assumes that $P2VVROOT/python is in $PYTHONPATH                         ##
## * assumes that $P2VVROOT/lib is in $LD_LIBRARYPATH                        ##
##                                                                           ##
## authors:                                                                  ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                ##
##                                                                           ##
###############################################################################

import P2VV, P2VVConfiguration, P2VVModelBuilders, P2VVPlots
from ROOT import RooDataSet, RooFit, TCanvas, TChain, TFile

###############################################################################
# specify decay mode ('Bd2JpsiKstar' or 'Bs2Jpsiphi')
mode = 'Bd2JpsiKstar'
#mode = 'Bs2Jpsiphi'

# data set name and file
#dataFilePath = '/data/bfys/jleerdam/Bs2Jpsiphi/EvtGen/Gauss-13144008-*.root'
dataSetName  = mode[3:] + 'Data'
dataFilePath = dataSetName + '.root'

# generate events?
generate = True
nEvents = 1000000

# read events from NTuple or RooDataset
NTuple = False

###############################################################################
# load the P2VV library
P2VV.loadP2VVLib()

# create P2VV configuration object
config = P2VVConfiguration.getP2VVConfig(mode, ['onlySignal'])

# custom settings
config['cpsiAng'].set(name = 'hel_cthetak')
config['cthetaAng'].set(name = 'hel_cthetal')
config['phiAng'].set(name = 'hel_phi')
config['BLifetime'].set(name = 't', min = 0., max = 4.)
config['iTag'].set(name = 'tagInitial')
if mode == 'Bd2JpsiKstar' :
  config['fTag'].set(name = 'tagFinal')
config['misTag'].set(realType = 'par', val = 0., min = '', max = '')

config['ReApar'].set(val = -0.6)
config['ImApar'].set(val = -0.1)
config['ReAperp'].set(val = -0.6)
config['ImAperp'].set(val = 0.1)
config['ReAS'].set(val = -0.2)
config['ImAS'].set(val = 0.3)

if mode == 'Bd2JpsiKstar' :
  config['dm'].set(min = -1., max = 2.)
elif mode == 'Bs2Jpsiphi' :
  config['dm'].set(min = 13., max = 23)
  config['phiCP'].set(val = -0.7)

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
  print 'computEfficiency: generating %d events' % nEvents
  if config.value('BDecayClass') == 'RooBDecay' : P2VV.registerMultiCatGen()
  data = pdf.generate(ws.set('observables'), nEvents)

  # write events to file
  print "computEfficiency: writing RooDataSet '%s' to file '%s'"\
      % (dataSetName, dataFilePath)
  file = TFile.Open(dataFilePath, 'RECREATE')
  data.Write(dataSetName)
  file.Close()

elif NTuple :
  # create data set from NTuple file(s)
  print "computEfficiency: reading NTuple(s) '%s' from file(s) '%s'"\
      % (dataSetName, dataFilePath)
  files = TChain(dataSetName)
  files.Add(dataFilePath)
  data = RooDataSet(dataSetName, dataSetName, files, ws.set('observables'))

else :
  # get data set from file
  print "computEfficiency: reading RooDataset '%s' from file '%s'"\
      % (dataSetName, dataFilePath)
  file = TFile.Open(dataFilePath, 'READ')
  data = file.Get(dataSetName)
  file.Close()

print 'compAngEfficiency: %d events in data set' % data.numEntries()

# adjust efficiency settings
config['effBasisType'].setValue('angular')
config['angEffBasisFuncs'].setValue((3, 3))

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

