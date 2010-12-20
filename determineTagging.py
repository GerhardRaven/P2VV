from ROOT import *
from ModelBuilders import declareObservables

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

tagcat = RooThresholdCategory("tagcat","tagcat",w['tagomega'],"untagged",0)
pdf = RooThresholdPdf("tageff","tageff",w['tagomega'])
effs = []
for (name,upper) in [ ("tag1", 0.2), ("tag2",0.3), ("tag3", 0.45) ] :
    eff = RooRealVar( name + "_eff", name + "_eff", 0.2, 0., 1.)
    effs += [ eff ]
    tagcat.addThreshold(upper,name)
    pdf.addThreshold(upper, eff )


# TODO: replace by multinomial....
#x = RooEfficiency("tageff","tageff",eps,tagcat,"untagged")
#x.fitTo(data)
pdf.fitTo(data)
plot = w['tagomega'].frame()
data.plotOn( plot )
pdf.plotOn( plot )
plot.Draw()
