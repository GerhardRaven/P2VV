from ROOT import *
gSystem.Load('libp2vv')
from ModelBuilders import buildMomentPDF,apybasis
from itertools import count,product
from math import pi

fname = 'Bs2JpsiPhiTuple.root'
dataName = 'dataset'

f = TFile(fname)
tree = f.Get(dataName)

w = RooWorkspace("w")
w.factory("{helcosthetaK[-1,1],helcosthetaL[-1,1],helphi[%s,%s]}"%(-pi,pi))
w.factory("{trcospsi[-1,1],trcostheta[-1,1],trphi[%s,%s]}"%(-pi,pi))
helangles =  w.argSet('helcosthetaK,helcosthetaL,helphi')
trangles  =  w.argSet('trcospsi,trcostheta,trphi')

#angles = trangles
angles = helangles

data = RooDataSet('data','data',tree,angles)

ab = apybasis(w,angles)
moments = []
#  Warning: the Y_lm are orthonormal, but the P_i are orthogonal, with dot product 2/(2*i+1)
for (i,l) in product(range(3),range(21)) :
      # transversity needs many more trphi moments
      # moments += [ Moment( ab.build('mom',i,0,l,m,1.), float(2*i+1)/2 ) for m in range(max(-l,-8),min(l+1,9)) ]
      # helicity has phi almost flat...
      moments += [ Moment( ab.build('mom',i,0,l,m,1.), float(2*i+1)/2 ) for m in range(max(-l,-2),min(l+1,3)) ]
pdf = buildMomentPDF( w, "bkg_angles_pdf", data, moments )

c = TCanvas()
c.Divide(3,2)
for (v,i) in zip( angles, count(1) ) :
    c.cd(i)
    frame = v.frame()
    data.plotOn(frame)
    pdf.plotOn(frame)
    frame.Draw()

    others = RooArgList( angles )
    others.remove( v )
    hist = pdf.createHistogram( others.names() )
    pdf.fillHistogram( hist,others,1., RooArgSet(v))
    c.cd(3+i)
    hist.Draw('COLZ')

c.Flush()

