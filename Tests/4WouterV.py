from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi
from array import array

from RooFitDecorators import *

name = 'TestSimFit'
wsfile = TFile('ToySimWS.root')

ws = wsfile.Get('ws')

pdf =  ws['sig_pdf']

ws.defineSet("observables","t,trcospsi,trcostheta,trphi,m,tagdecision,fitcat")
ras = ws.set('observables')
data = pdf.generate(ras,2)

assert False
