from ROOT import *
gSystem.Load('libp2vv')
from ModelBuilders import buildMomentPDF
from itertools import count,product
from math import pi

fname = 'duitsedata.root'
dataName = 'Bs2JpsiPhi'

f = TFile(fname)
tree = f.Get(dataName)

w = RooWorkspace("w")
w.factory("{trcospsi[-1,1],trcostheta[-1,1],trphi[%s,%s]}"%(-pi,pi))

data = RooDataSet('data','data',tree,w.argSet('trcospsi,trcostheta,trphi'))

ab = abasis(w,'trcospsi','trcostheta','trphi')
moments = []
#  Warning: the Y_lm are orthonormal, but the P_i are orthogonal, with dot product 2/(2*i+1)
for (i,l) in product(range(4),range(12)) :
      moments += [ Moment( ab('mom',i,0,l,m,1.), float(2*i+1)/2 ) for m in range(-l,l+1) ]
pdf = buildMomentPDF( w, "bkg_angles_pdf", data, moments )

c = TCanvas()
c.Divide(3,1)
for (v,i) in zip( w.argSet('trcospsi,trcostheta,trphi'), count(1) )  :
    c.cd(i)
    frame = v.frame()
    data.plotOn(frame)
    pdf.plotOn(frame)
    frame.Draw()
c.Flush()
