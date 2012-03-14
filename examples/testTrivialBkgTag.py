from RooFitWrappers import *
RooObject(workspace='mine')
iTag      = Category( 'tagdecision' , Title = 'initial state flavour tag',   Observable = True,  States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
from P2VVParameterizations.FlavourTagging import Trivial_TagPdf
p = Trivial_TagPdf(tagdecision=iTag, NamePrefix = 'bkg')
pdf = p.pdf()
data = pdf.generate(iTag,100000)
t = data.table(iTag)
t.Print("V")
