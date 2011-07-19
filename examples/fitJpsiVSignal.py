###############################################################################
## fitJpsiVSignal:                                                           ##
##   P2VV example script for a fit to B->J/psiV signal data                  ##
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

# amplitude values
A0Mag2    = 0.60
A0Ph      = 0.
AparMag2  = 0.24
AparPh    = 2.5
AperpMag2 = 0.16
AperpPh   = -0.17
ASMag2    = 0.05
ASPh      = 2.2

# lambda
phiCP      = -0.2
lambdaCPSq = 1.

# B lifetime and mixing
BLifetimeError = 0.05
if mode == 'Bd2JpsiKstar' :
  Gamma  = 0.65
  dGamma = 0.
  dm     = 0.51
elif mode == 'Bs2Jpsiphi' :
  Gamma  = 0.68
  dGamma = 0.05
  dm     = 17.8

# flavour tags
AProd   = 0.
ANorm   = -(1. - lambdaCPSq) / (1. + lambdaCPSq)
wTag    = 0.35
wTagBar = 0.35


###############################################################################
from math import sqrt, sin, cos
import P2VV, P2VVConfiguration, P2VVModelBuilders, P2VVPlots
from ROOT import RooFit, TCanvas

# load the P2VV library
P2VV.loadP2VVLib()

# set RooFit output
P2VV.setRooFitOutput()

# create P2VV configuration object
config = P2VVConfiguration.getP2VVConfig(mode, ['onlySignal','KSWave=include'])
  # additional options:
  # 'anglesType=trans'      : use transversity angles
  # 'KSWave=includeEven'    : a CP even S-wave in B_s -> J/psi phi
  # 'ampsType=transPolar'   : use polar transversity amplitudes
  # 'lambdaCPType=polar'    : use a polar parametrisation of lambda
  # 'RooBDecay'             : use RooBDecay instead of RooBTagDecay
  # 'tResModel=3Gauss'      : use a time resolution model with three Gaussians
  # 'allowITagZero'         : allow 0 value for init. state tag (untagged)
  # 'asymType=coefficients' : fit for the "sum" of the CP asymmetries
  # 'tagType=categories5'   : split into 5 tagging categories
  # 'effType=angular'       : split into 5 tagging categories

# custom settings
if config.value('anglesType')[0] == 'trans' :
  config['cpsiAng'].set(name = 'tr_cpsi')
  config['cthetaAng'].set(name = 'tr_ctheta')
  config['phiAng'].set(name = 'tr_phi')
else :
  config['cpsiAng'].set(name = 'hel_cthetak')
  config['cthetaAng'].set(name = 'hel_cthetal')
  config['phiAng'].set(name = 'hel_phi')

config['BLifetime'].set(name = 't', min = -0.5, max = 5.)
config['iTag'].set(name = 'tagInitial')
if mode == 'Bd2JpsiKstar' :
  config['fTag'].set(name = 'tagFinal')

if config.value('ampsType') == 'transPolar' :
  # A_par^2 = 1 - A_0^2 - A_perp^2 :: Im(A_0) = 0
  config['A0Mag2'].set(val = A0Mag2)
  config['AperpMag2'].set(val = AperpMag2)
  config['AparPh'].set(val = AparPh)
  config['AperpPh'].set(val = AperpPh)
  if config.value('KSWave')[:7] == 'include' :
    config['ASMag2'].set(val = ASMag2)
    config['ASPh'].set(val = ASPh)
elif config.value('ampsType') == 'transCartesian' :
  # Re(A_0) = 1 :: Im(A_0) = 0
  config['ReApar'].set(val = sqrt(AparMag2 / A0Mag2) * cos(AparPh))
  config['ImApar'].set(val = sqrt(AparMag2 / A0Mag2) * sin(AparPh))
  config['ReAperp'].set(val = sqrt(AperpMag2 / A0Mag2) * cos(AperpPh))
  config['ImAperp'].set(val = sqrt(AperpMag2 / A0Mag2) * sin(AperpPh))
  if config.value('KSWave')[:7] == 'include' :
    config['ReAS'].set(val = sqrt(ASMag2 / A0Mag2) * cos(ASPh))
    config['ImAS'].set(val = sqrt(ASMag2 / A0Mag2) * sin(ASPh))

config['Gamma'].set(val = Gamma)
config['dGamma'].set(val = dGamma)
if mode == 'Bd2JpsiKstar' :
  config['dm'].set(min = -1., max = 2.)
  if config.value('lambdaCPType') == 'polar' :
    config['lambdaCPSq'].set(val = lambdaCPSq)
  else :
    config['ReLambdaCP'].set(val = sqrt(lambdaCPSq))
    config['ImLambdaCP'].set(val = 0.)
elif mode == 'Bs2Jpsiphi' :
  config['dm'].set(min = 13., max = 23)
  if config.value('lambdaCPType') == 'polar' :
    config['phiCP'].set(val = phiCP)
    config['lambdaCPSq'].set(val = lambdaCPSq)
  else :
    config['ReLambdaCP'].set(val = sqrt(lambdaCPSq) * cos(-phiCP))
    config['ImLambdaCP'].set(val = sqrt(lambdaCPSq) * sin(-phiCP))

config['wTag'].set(val = wTag, min=0., max=0.5)
config['wTagBar'].set(val = wTagBar, min=0., max=0.5)
if config.value('asymType') == 'coefficients' :
  config['avgCOdd'].set(val = (AProd + ANorm) / (1. + AProd * ANorm))
else :
  config['AProd'].set(val = AProd)

if 'Gauss' in config.value('tResModel') :
  config['BLifetimeError'].set(val = BLifetimeError)

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
  P2VV.writeData(dataFilePath, dataSetName, data, NTuple)

else :
  # get data from file
  data = P2VV.readData(dataFilePath, dataSetName, NTuple,
      ws.set('observables'))

print 'fitJpsiVSignal: %d events in data set' % data.numEntries()

# fit data
fitResult = pdf.fitTo(data, RooFit.Minos(False), RooFit.Hesse(False),
    RooFit.NumCPU(4), RooFit.Save())

# print polar (cartesian) amplitudes if 'ampsType' is cartesian (polar)
if config.value('ampsType') == 'transCartesian'\
    or config.value('ampsType') == 'transPolar' :
  P2VV.convertAmplitudes(config, fitResult, True)

# print polar (cartesian) lambda if 'lambdaCPType' is cartesian (polar)
if mode == 'Bs2Jpsiphi' and (config.value('lambdaCPType') == 'cartesian'\
    or (config.value('lambdaCPType') == 'polar'\
    and config['lambdaCPSq'].minValue()\
    != config['lambdaCPSq'].maxValue())) :
  P2VV.convertLambda(config, fitResult, True)

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

# cos(theta_K)/cos(psi_tr) plots
if config.value('anglesType')[0] == 'trans' :
  angleName = 'cos(#psi_{tr})'
else :
  angleName = 'cos(#theta_{K})'

cpsiCanv = TCanvas('cpsiCanv', 'cos(psi)')
cpsiCanv.Divide(2, 2)

if mode == 'Bd2JpsiKstar' :
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(1), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(2), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftMin)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(3), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B} / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(4), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B} / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftMin)])

elif mode == 'Bs2Jpsiphi' :
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(1), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - CP Average'),
                   RooFit.Bins(70)],
      dataOpts  = [], pdfOpts   = [])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(3), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus)])
  P2VVPlots.plot(config, 'cpsiAng', cpsiCanv.cd(4), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B}'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin)])

cpsiCanv.Print(plotsFile, 'ps')

# cos(theta_l)/cos(theta_tr) plots
if config.value('anglesType')[0] == 'trans' :
  angleName = 'cos(#theta_{tr})'
else :
  angleName = 'cos(#theta_{l})'

cthetaCanv = TCanvas('cthetaCanv', 'cos(theta)')
cthetaCanv.Divide(2, 2)

if mode == 'Bd2JpsiKstar' :
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(1), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(2), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftMin)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(3), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B} / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(4), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B} / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftMin)])

elif mode == 'Bs2Jpsiphi' :
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(1), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - CP Average'),
                   RooFit.Bins(70)],
      dataOpts  = [], pdfOpts   = [])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(3), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus)])
  P2VVPlots.plot(config, 'cthetaAng', cthetaCanv.cd(4), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B}'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin)])

cthetaCanv.Print(plotsFile, 'ps')

# phi_hel/phi_tr plots
if config.value('anglesType')[0] == 'trans' :
  angleName = '#phi_{tr}'
else :
  angleName = '#phi_{hel}'

phiCanv = TCanvas('phiCanv', 'phi')
phiCanv.Divide(2, 2)

if mode == 'Bd2JpsiKstar' :
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(1), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(2), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == 1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus),
                   RooFit.Slice(ws.cat(ftName), ftMin)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(3), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B} / Oscillated'),
                   RooFit.Bins(30)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == 1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftPlus)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(4), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B} / Unoscillated'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1 && %s == -1' % (itName, ftName))],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin),
                   RooFit.Slice(ws.cat(ftName), ftMin)])

elif mode == 'Bs2Jpsiphi' :
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(1), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - CP Average'),
                   RooFit.Bins(70)],
      dataOpts  = [], pdfOpts   = [])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(3), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - B'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == 1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itPlus)])
  P2VVPlots.plot(config, 'phiAng', phiCanv.cd(4), data, pdf,
      xTitle = angleName,
      frameOpts = [RooFit.Title(angleName + ' - #bar{B}'),
                   RooFit.Bins(50)],
      dataOpts  = [RooFit.Cut('%s == -1' % itName)],
      pdfOpts   = [RooFit.Slice(ws.cat(itName), itMin)])

phiCanv.Print(plotsFile + ')', 'ps')

