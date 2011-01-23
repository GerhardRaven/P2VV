import RooFitDecorators
from itertools import product,count
from ModelBuilders import *
from ROOT import *
gSystem.Load("libp2vv.so")

class efficiency :
    def __init__(self,*args) : 
        (self.cpsi,self.ctheta,self.phi) = args[0] if len(args)==1 else args
    def accept(self) :
        x = self.cpsi.getVal()
        y = self.ctheta.getVal()
        from math import cos
        z = -1+2*cos( self.phi.getVal() )
        from random import random
        return random() > ( x*x*y/( 1 if z<0 else 1-z ) )

fname="p2vv_9.root"
pdfName = "jpsiphipdf"
dataName = "jpsiphipdfData"
workspaceName = "w"

f = TFile(fname)
w = f.Get(workspaceName) 
pdf = w[pdfName]
data = w[dataName] 
allObs = pdf.getObservables( data.get() )

#angles = w.argSet("transversityangles")
angles = w.set("helicityangles")

#replace input by inefficient data
eps = efficiency( angles )
inEffData = RooDataSet( "inEffData","inEffData", allObs )
for event in  data :
    allObs.assignValueOnly( event )
    if eps.accept() : inEffData.add( allObs )
data = inEffData;
          
# define the moments used to describe the efficiency
# for this, we need the PDF used to generate the data
marginalObs = pdf.getObservables( data.get() )
marginalObs.remove( angles )
# marginalize pdf over 'the rest' so we get the normalization of the moments right...
pdf_marginal = pdf.createProjection(marginalObs)
moments = []

ab = abasis(w,angles)
for (i,l) in product(range(4),range(4)) :
    # if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
    # Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
    moments += [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, pdf_marginal, allObs ) for m in range(-l,l+1) ]

# loop over all data, determine moments
pdf_eff = buildEffMomentsPDF(w, "_eff", pdf, data, moments)


c = TCanvas()
c.Divide(3,1);
for (i,var) in enumerate(angles) :
     c.cd(1+i) # start argument of enumerate is python >= 2.6...
     plot = var.frame()
     data.plotOn(plot); 
     pdf.plotOn(plot,RooFit.LineColor(kBlue))
     pdf_eff.plotOn(plot,RooFit.LineColor(kRed))
     plot.Draw()
c.Flush()

