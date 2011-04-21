###############################################################################
## fitJpsiVSignal:                                                           ##
##   P2VV example script for a fit to signal data (signal only, no           ##
##   experimental effects)                                                   ##
##                                                                           ##
## * decay channels: B0->J/psiK* or B_s0->J/psiphi                           ##
## * assumes that $P2VVROOT/python is in $PYTHONPATH                         ##
## * assumes that $P2VVROOT/lib is in $LD_LIBRARYPATH                        ##
##                                                                           ##
## authors:                                                                  ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                ##
##                                                                           ##
###############################################################################

# specify decay mode ('Bd2JpsiKstar' or 'Bs2Jpsiphi')
mode = 'Bd2JpsiKstar'
#mode = 'Bs2Jpsiphi'

# plots file
plotsFile = mode[3:] + 'Plots.ps'

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
import P2VV, P2VVConfiguration, P2VVModelBuilders, P2VVPlots
from ROOT import RooDataSet, RooFit, TCanvas, TChain, TFile

# load the P2VV library
P2VV.loadP2VVLib()

# set RooFit output
P2VV.setRooFitOutput()

# create P2VV configuration object
config = P2VVConfiguration.getP2VVConfig(mode, ['onlySignal'])
    # , 'transAngles', 'ampsType=polar', 'lambdaCPType=polar', 'noKSWave'
    # , 'RooBDecay', 'tResModel=3Gauss'])

# custom settings
if config.value('anglesType')[0] == 'trans' :
  config['cpsiAng'].set(name = 'tr_cpsi')
  config['cthetaAng'].set(name = 'tr_ctheta')
  config['phiAng'].set(name = 'tr_phi')
else :
  config['cpsiAng'].set(name = 'hel_cthetak')
  config['cthetaAng'].set(name = 'hel_cthetal')
  config['phiAng'].set(name = 'hel_phi')

config['BLifetime'].set(name = 't', min = 0., max = 4.)
config['iTag'].set(name = 'tagInitial')
if mode == 'Bd2JpsiKstar' :
  config['fTag'].set(name = 'tagFinal')
config['misTag'].set(realType = 'par', val = 0., min = '', max = '')

if config.value('ampsType') == 'polar' :
  # A_par^2 = 1 - A_0^2 - A_perp^2 :: Im(A_0) = 0
  config['A0Mag2'].set(val = 0.60)
  config['AperpMag2'].set(val = 0.16)
  config['ASMag2'].set(val = 0.05)
  config['AparPh'].set(val = 2.5 )
  config['AperpPh'].set(val = -0.2)
  config['ASPh'].set(val = 2.2)
else :
  # Re(A_0) = 1 :: Im(A_0) = 0
  config['ReApar'].set(val = -0.51)
  config['ImApar'].set(val = 0.38)
  config['ReAperp'].set(val = 0.51)
  config['ImAperp'].set(val = -0.10)
  config['ReAS'].set(val = -0.17)
  config['ImAS'].set(val = 0.23)

if mode == 'Bd2JpsiKstar' :
  config['dm'].set(min = -1., max = 2.)
elif mode == 'Bs2Jpsiphi' :
  config['dm'].set(min = 13., max = 23)
  if config.value('lambdaCPType') == 'polar' :
    config['phiCP'].set(val = -0.2)
    config['lambdaCPSq'].set(val = 1.)
  else :
    config['ReLambdaCP'].set(val = 0.980)
    config['ImLambdaCP'].set(val = 0.199)

# declare RooFit variables and store them in RooWorkspace
config.declareRooVars()

# get workspace
ws = config.workspace()

# build the PDF
pdf = P2VVModelBuilders.getP2VVPDF(config)

# print contents of RooWorkspace to screen
#config.workspace().Print()

if generate :
  # generate events
  print 'fitJpsiVSignal: generating %d events' % nEvents
  if config.value('BDecayClass') == 'RooBDecay' : P2VV.registerMultiCatGen()
  data = pdf.generate(ws.set('observables'), nEvents)

  # write events to file
  print "fitJpsiVSignal: writing RooDataSet '%s' to file '%s'"\
      % (dataSetName, dataFilePath)
  file = TFile.Open(dataFilePath, 'RECREATE')
  data.Write(dataSetName)
  file.Close()

elif NTuple :
  # create data set from NTuple file(s)
  print "fitJpsiVSignal: reading NTuple(s) '%s' from file(s) '%s'"\
      % (dataSetName, dataFilePath)
  files = TChain(dataSetName)
  files.Add(dataFilePath)
  data = RooDataSet(dataSetName, dataSetName, files, ws.set('observables'))

else :
  # get data set from file
  print "fitJpsiVSignal: reading RooDataset '%s' from file '%s'"\
      % (dataSetName, dataFilePath)
  file = TFile.Open(dataFilePath, 'READ')
  data = file.Get(dataSetName)
  file.Close()

print 'fitJpsiVSignal: %d events in data set' % data.numEntries()

# fit data
fitResult = pdf.fitTo(data, RooFit.Minos(False), RooFit.Hesse(False),
    RooFit.NumCPU(8), RooFit.Save())

# print polar (cartesian) amplitudes if 'ampsType' is cartesian (polar)
P2VV.convertAmplitudes(config, fitResult, True)

# get tags
itName = config['iTag'].name()
itPlus = config['iTag'].catTypesDict()[1]
itMin  = config['iTag'].catTypesDict()[-1]
if mode == 'Bd2JpsiKstar' :
  ftName = config['fTag'].name()
  ftPlus = config['fTag'].catTypesDict()[1]
  ftMin  = config['fTag'].catTypesDict()[-1]

# set plot style
P2VVPlots.setP2VVPlotStyle()

# lifetime plots
tCanv = TCanvas('tCanv', 'Lifetime')
tCanv.Divide(2, 2)

if mode == 'Bd2JpsiKstar' :
  P2VVPlots.plot(config, 'BLifetime', tCanv.cd(1), data, pdf,
      xTitle = 't (ps)',
      frameOpts = [RooFit.Title('Lifetime - B / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'BLifetime', tCanv.cd(2), data, pdf,
      xTitle = 't (ps)',
      frameOpts = [RooFit.Title('Lifetime - B / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftMin)])
  P2VVPlots.plot(config, 'BLifetime', tCanv.cd(3), data, pdf,
      xTitle = 't (ps)',
      frameOpts = [RooFit.Title('Lifetime - #bar{B} / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'BLifetime', tCanv.cd(4), data, pdf,
      xTitle = 't (ps)',
      frameOpts = [RooFit.Title('Lifetime - #bar{B} / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftMin)])

elif mode == 'Bs2Jpsiphi' :
  P2VVPlots.plot(config, 'BLifetime', tCanv.cd(1), data, pdf,
      xTitle = 't (ps)',
      frameOpts = [RooFit.Title('Lifetime - CP Average'),
                   RooFit.Bins(70)],
      dataOpts  = [], pdfOpts   = [])
  P2VVPlots.plot(config, 'BLifetime', tCanv.cd(3), data, pdf,
      xTitle = 't (ps)',
      frameOpts = [RooFit.Title('Lifetime - B'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus)])
  P2VVPlots.plot(config, 'BLifetime', tCanv.cd(4), data, pdf,
      xTitle = 't (ps)',
      frameOpts = [RooFit.Title('Lifetime - #bar{B}'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin)])

tCanv.Print(plotsFile + '(', 'ps')

# cos(theta_K) plots
cpsiCanv = TCanvas('cpsiCanv', 'cos(theta_K)')
cpsiCanv.Divide(2, 2)

if mode == 'Bd2JpsiKstar' :
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(1), data, pdf,
      xTitle = 'cos(#theta_{K})',
      frameOpts = [RooFit.Title('cos(#theta_{K}) - B / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(2), data, pdf,
      xTitle = 'cos(#theta_{K})',
      frameOpts = [RooFit.Title('cos(#theta_{K}) - B / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftMin)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(3), data, pdf,
      xTitle = 'cos(#theta_{K})',
      frameOpts = [RooFit.Title('cos(#theta_{K}) - #bar{B} / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(4), data, pdf,
      xTitle = 'cos(#theta_{K})',
      frameOpts = [RooFit.Title('cos(#theta_{K}) - #bar{B} / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftMin)])

elif mode == 'Bs2Jpsiphi' :
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(1), data, pdf,
      xTitle = 'cos(#theta_{K})',
      frameOpts = [RooFit.Title('cos(#theta_{K}) - CP Average'),
                   RooFit.Bins(70)],
      dataOpts  = [], pdfOpts   = [])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(3), data, pdf,
      xTitle = 'cos(#theta_{K})',
      frameOpts = [RooFit.Title('cos(#theta_{K}) - B'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(4), data, pdf,
      xTitle = 'cos(#theta_{K})',
      frameOpts = [RooFit.Title('cos(#theta_{K}) - #bar{B}'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin)])

cpsiCanv.Print(plotsFile, 'ps')

# cos(theta_l) plots
cthetaCanv = TCanvas('cthetaCanv', 'cos(theta_l)')
cthetaCanv.Divide(2, 2)

if mode == 'Bd2JpsiKstar' :
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(1), data, pdf,
      xTitle = 'cos(#theta_{l})',
      frameOpts = [RooFit.Title('cos(#theta_{l}) - B / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(2), data, pdf,
      xTitle = 'cos(#theta_{l})',
      frameOpts = [RooFit.Title('cos(#theta_{l}) - B / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftMin)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(3), data, pdf,
      xTitle = 'cos(#theta_{l})',
      frameOpts = [RooFit.Title('cos(#theta_{l}) - #bar{B} / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(4), data, pdf,
      xTitle = 'cos(#theta_{l})',
      frameOpts = [RooFit.Title('cos(#theta_{l}) - #bar{B} / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftMin)])

elif mode == 'Bs2Jpsiphi' :
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(1), data, pdf,
      xTitle = 'cos(#theta_{l})',
      frameOpts = [RooFit.Title('cos(#theta_{l}) - CP Average'),
                   RooFit.Bins(70)],
      dataOpts  = [], pdfOpts   = [])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(3), data, pdf,
      xTitle = 'cos(#theta_{l})',
      frameOpts = [RooFit.Title('cos(#theta_{l}) - B'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(4), data, pdf,
      xTitle = 'cos(#theta_{l})',
      frameOpts = [RooFit.Title('cos(#theta_{l}) - #bar{B}'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin)])

cthetaCanv.Print(plotsFile, 'ps')

# phi plots
phiCanv = TCanvas('phiCanv', 'phi')
phiCanv.Divide(2, 2)

if mode == 'Bd2JpsiKstar' :
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(1), data, pdf,
      xTitle = '#phi',
      frameOpts = [RooFit.Title('#phi - B / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(2), data, pdf,
      xTitle = '#phi',
      frameOpts = [RooFit.Title('#phi - B / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftMin)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(3), data, pdf,
      xTitle = '#phi',
      frameOpts = [RooFit.Title('#phi - #bar{B} / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(4), data, pdf,
      xTitle = '#phi',
      frameOpts = [RooFit.Title('#phi - #bar{B} / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftMin)])

elif mode == 'Bs2Jpsiphi' :
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(1), data, pdf,
      xTitle = '#phi',
      frameOpts = [RooFit.Title('#phi - CP Average'),
                   RooFit.Bins(70)],
      dataOpts  = [], pdfOpts   = [])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(3), data, pdf,
      xTitle = '#phi',
      frameOpts = [RooFit.Title('#phi - B'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(4), data, pdf,
      xTitle = '#phi',
      frameOpts = [RooFit.Title('#phi - #bar{B}'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin)])

phiCanv.Print(plotsFile + ')', 'ps')

