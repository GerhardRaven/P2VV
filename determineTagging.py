from ROOT import *
from ModelBuilders import declareObservables,buildTagging

fname = 'Bs2JpsiPhiTuple.root'
dataName = 'dataset'

f = TFile(fname)
tree = f.Get(dataName)

w = RooWorkspace("w")
declareObservables(w)

tags = w.argSet('tagdecision,tagomega') 
data = RooDataSet('data','data',tree,tags)

(tagcat,pdf) = buildTagging(w,'sigtag',[0.25,0.35,0.45])

pdf.fitTo(data)

plot1 = w['tagomega'].frame()
plot2 = w['tagomega'].frame()
data.plotOn( plot1 )
binning = pdf.getBinning() # need to make sure it is not garbage collected too early...
data.plotOn( plot2, RooFit.Binning( binning ),RooFit.MarkerColor(kRed)  )
pdf.plotOn( plot1 )
pdf.plotOn( plot2 )

c = TCanvas()
c.Divide(1,2)
c.cd(1)
plot1.Draw()
c.cd(2)
plot2.Draw()
