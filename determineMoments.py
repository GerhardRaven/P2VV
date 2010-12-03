from ROOT import *
gSystem.Load('libp2vv')
from ModelBuilders import buildMomentPDF
from itertools import count,product

fname = 'p2vv_7.root'
dataName = 'pdfData'
wsname = 'w'

f = TFile(fname)
w = f.Get(wsname)
data = w.data(dataName)

ab = abasis(w,'trcospsi','trcostheta','trphi')
moments = []
#  Warning: the Y_lm are orthonormal, but the P_i are orthogonal, with dot product 2/(2*i+1)
for (i,l) in product(range(5),range(5)) :
      moments += [ Moment( ab('mom',i,0,l,m,1.), float(2*i+1)/2 ) for m in range(-l,l+1) ]
pdf = buildMomentPDF( w, "bkg_angles_pdf", data, moments )

c = TCanvas()
c.Divide(3,1)
for (v,i) in zip( ['trcospsi','trcostheta','trphi'], count(1) )  :
    c.cd(i)
    frame = w.var(v).frame()
    data.plotOn(frame)
    pdf.plotOn(frame)
    frame.Draw()
c.Flush()
