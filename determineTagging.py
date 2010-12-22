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

### TODO: split the tagomega distribution (&parameterization!) according to tagdecision
###       write dilution used in the fit as Dilution_i(tagomega) ( = 1 - 2 * (alpha_ij P_j(w) )  for i = b,bbar. 
###       make overall PDF conditional on tagomega distribution...(something like RooEfficiecy( tagdecision 
### TODO: add a 'parametricstep'-like PDF which uses the (relative) efficiency in each
###       bin as parameterization -- the will simplify the interpretation of the parameters
###       esp. in the case of variable binwidth.  Also allow a user to decide which is the
###       'remainder' bin (eg. would like to use the 'untagged' i.e. tagomega....
### ALTERNATIVE:
###       split in categories (using RooThresholdCategory), assign one b, one bbar mistag per catagorie
###       and add a Multinomial (with parameters split for b,bbar) to account for the efficiency...
###       I suspect that Multionial | RooThresholdCategory == RooParametricStep....

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
