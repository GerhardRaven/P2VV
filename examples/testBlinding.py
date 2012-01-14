from ROOT import *
ws = RooWorkspace('ws')

ws.factory('{xObs[0., -1., 2.], mPar_blind[0.5, -1., 2.], sPar_blind[2., 1., 3.]}')
ws.factory('RooGaussian::gPDF(xObs, mPar_blind, sPar_blind)')

print '---------------'
print 'Generating data'
print '---------------'
data = ws.pdf('gPDF').generate(ws.argSet('xObs'), 1000000)

print '-----------------'
print 'fitting with gPDF'
print '-----------------'
ws.arg('mPar_blind').setVal(0.5)
ws.arg('mPar_blind').setError(0.)
ws.arg('sPar_blind').setVal(2.)
ws.arg('sPar_blind').setError(0.)

ws.pdf('gPDF').fitTo(data, RooFit.Minos(False), RooFit.Hesse(False))

xFrame = ws.arg('xObs').frame()
data.plotOn(xFrame)
ws.pdf('gPDF').plotOn(xFrame)

print '-----------------------'
print 'fitting with gPDF_blind'
print '-----------------------'
ws.arg('mPar_blind').setVal(0.5)
ws.arg('mPar_blind').setError(0.)
ws.arg('sPar_blind').setVal(2.)
ws.arg('sPar_blind').setError(0.)

ws.factory("RooUnblindUniform::mPar('Blaat', 1., mPar_blind)")
ws.factory("RooUnblindUniform::sPar('Boehhhh', 1., sPar_blind)")
ws.factory('RooGaussian::gPDF_blind(xObs, mPar, sPar)')
ws.pdf('gPDF_blind').fitTo(data, RooFit.Minos(False), RooFit.Hesse(False))

xFrame_blind = ws.arg('xObs').frame()
data.plotOn(xFrame_blind)
ws.pdf('gPDF_blind').plotOn(xFrame_blind)

print '-------------'
print 'drawing plots'
print '-------------'
from P2VVLoad import ROOTStyle
canv = TCanvas('canv', '')
canv.Divide(2, 2)

canv.cd(1)
xFrame.Draw()
canv.cd(2)
xFrame_blind.Draw()

