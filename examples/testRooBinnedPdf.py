###########################################################################################################################################
## testRooBinnedPdf:                                                                                                              ##
##   test/demonstrate the behaviour of the class RooBinnedPdf                                                                     ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                                                                            ##
##                                                                                                                                       ##
###########################################################################################################################################

# plots file
plotsFile = 'BinnedPdfPlots.ps'

###########################################################################################################################################

# load the P2VV library
from P2VVLoad import P2VVLibrary

# create a base variable and its binning
from array import array
from ROOT import RooRealVar, RooBinning
baseVar = RooRealVar( 'baseVar', 'baseVar', 0., -1., 1. )
binBoundaries = array( 'd', [ -1., -0.5, -0.25, -0.1, 0.3, 1. ] )
bins = RooBinning( 5, binBoundaries, 'baseVarBins' )
baseVar.setBinning( bins, 'baseVarBins' )

# create N - 1 bin coefficients
from ROOT import RooArgList
c1 = RooRealVar( 'c1', 'c1', 0.40, 0., 2. )
c2 = RooRealVar( 'c2', 'c2', 0.15, 0., 2. )
c3 = RooRealVar( 'c3', 'c3', 0.10, 0., 2. )
c4 = RooRealVar( 'c4', 'c4', 0.20, 0., 2. )
cList = RooArgList( c1, c2, c3, c4 )

# create a step function
from ROOT import RooBinnedPdf
mmn = RooBinnedPdf( 'mmn', 'mmn', baseVar, 'baseVarBins', cList, 0, 0 )

# create a Gaussian PDF
from ROOT import RooGaussian
mean  = RooRealVar(  'mean', 'mean',   0.  )
sigma = RooRealVar(  'sigma', 'sigma', 0.4 )
gauss = RooGaussian( 'gauss', 'gauss', baseVar, mean, sigma )

# create a product of the step function and the PDF
from ROOT import RooProdPdf
pdf = RooProdPdf( 'pdf', 'pdf', mmn, gauss )

# generate some output
print 'bin centres:'
for i in range( bins.numBins()  ) : print i, bins.binCenter(i)
print '\nbin widths:'
for i in range( bins.numBins()  ) : print i, bins.binWidth(i)
print '\ncoefficients:'
for i in range( cList.getSize() ) : print i, cList.at(i).getVal()

# plot the step function, the PDF and PDF data
from P2VVLoad import RooFitOutput
from P2VVLoad import LHCbStyle

from ROOT import TCanvas
canv = TCanvas()
canv.Divide(2, 2)

print 'plotting the step PDF vs the base variable'
canv.cd(1)
frame1 = baseVar.frame()
mmn.plotOn(frame1)
frame1.Draw()

print 'plotting the step PDF vs the third coefficient'
canv.cd(2)
frame2 = c3.frame()
mmn.plotOn(frame2)
frame2.Draw()

print 'generating baseVar data'
from ROOT import RooArgSet
data3 = pdf.generate(RooArgSet(baseVar), 100000)

print 'plotting the PDF vs the base variable'
canv.cd(3)
frame3 = baseVar.frame()
data3.plotOn(frame3)
pdf.plotOn(frame3)
frame3.Draw()

print 'generating c3 data'
data4 = pdf.generate(RooArgSet(c3), 100000)

print 'plotting the PDF vs the third coefficient'
canv.cd(4)
frame4 = c3.frame()
data4.plotOn(frame4)
pdf.plotOn(frame4)
frame4.Draw()

canv.Print(plotsFile)

