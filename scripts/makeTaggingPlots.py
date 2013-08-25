# script paramters
plotsFilePath   = 'plots/taggingPlots.ps'
dataSetFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets_4KKMassBins_noTagCats.root'
fullDataSetName = 'JpsiKK_splotdata'
sigDataSetName  = 'JpsiKK_splotdata_weighted_sigMass'
cbkgDataSetName = 'JpsiKK_splotdata_weighted_cbkgMass'
numEstWTagBins  = 50

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

# read data sets
from P2VV.Utilities.DataHandling import readData
fullData = readData( filePath = dataSetFilePath, dataSetName = fullDataSetName, NTuple = False )
sigData  = readData( filePath = dataSetFilePath, dataSetName = sigDataSetName,  NTuple = False )
cbkgData = readData( filePath = dataSetFilePath, dataSetName = cbkgDataSetName, NTuple = False )

# get observables
from P2VV.RooFitWrappers import RealVar, Category
estWTagOS = RealVar('tagomega_os')
estWTagSS = RealVar('tagomega_ss')
tagCatOS  = Category('tagCatP2VVOS')
tagCatSS  = Category('tagCatP2VVSS')

# build PDFs for estimated wrong-tag probabilities
taggedDataOS = sigData.reduce( '%s > 0' % tagCatOS.GetName() )
taggedDataSS = sigData.reduce( '%s > 0' % tagCatSS.GetName() )
                                                                                                                     
from P2VV.RooFitWrappers import HistPdf
tagPdfOS = HistPdf(  Name = 'sig_bkg_estWTagOS'
                   , Observables = [ estWTagOS ]
                   , Binning = { estWTagOS : numEstWTagBins }
                   , Data = taggedDataOS
                  )
                                                                                                                     
tagPdfSS = HistPdf(  Name = 'sig_bkg_estWTagSS'
                   , Observables = [ estWTagSS ]
                   , Binning = { estWTagSS : numEstWTagBins }
                   , Data = taggedDataSS
                  )
                                                                                                                     
# get normalization correction for tagged events
untagFracOS    = fullData.table(tagCatOS).getFrac('Untagged')
untagFracSigOS = sigData.table(tagCatOS).getFrac('Untagged')
untagFracBkgOS = cbkgData.table(tagCatOS).getFrac('Untagged')
untagFracSS    = fullData.table(tagCatSS).getFrac('Untagged')
untagFracSigSS = sigData.table(tagCatSS).getFrac('Untagged')
untagFracBkgSS = cbkgData.table(tagCatSS).getFrac('Untagged')
                                                                                                                     
# plot estimated wrong-tag probabilities for signal and for background
from P2VV.Load import LHCbStyle
from P2VV.Utilities.Plotting import plot
from ROOT import TCanvas, kBlue, kFullDotLarge
canvsOS = [ TCanvas( 'estWTagCanvOS%d' % it ) for it in range(3) ]
for ( canv, data, nBins, norm ) in zip(  canvsOS
                                      , [ fullData, sigData, cbkgData ]
                                      , 3 * [ numEstWTagBins ]
                                      , [ 1. - untagFracOS, 1. - untagFracSigOS, 1. - untagFracBkgOS ]
                                     ) :
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)
                                                                                                                     
    #plot(  canv, estWTagOS, data, tagPdfOS
    plot(  canv, estWTagOS,  data, yScale = ( 0., None ), xTitleOffset = 1.10, yTitleOffset = 1.15
         , xTitle     = '#eta^{OS}'
         , yTitle     = 'Candidates / %.2f' % ( 0.499999 / float(nBins) )
         , frameOpts  = dict( Bins = nBins, Range = ( 0., 0.499999 ), Name = estWTagOS.GetName() )
         , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.7, LineWidth = 3 )
         #, pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm )
        )

canvsSS = [ TCanvas( 'estWTagCanvSS%d' % it ) for it in range(3) ]
for ( canv, data, nBins, norm ) in zip(  canvsSS
                                      , [ fullData, sigData, cbkgData ]
                                      , 3 * [ numEstWTagBins ]
                                      , [ 1. - untagFracSS, 1. - untagFracSigSS, 1. - untagFracBkgSS ]
                                     ) :
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)
                                                                                                                     
    #plot(  canv, estWTagSS, data, tagPdfSS
    plot(  canv, estWTagSS,  data, yScale = ( 0., None ), xTitleOffset = 1.10, yTitleOffset = 1.15
         , xTitle     = '#eta^{SSK}'
         , yTitle     = 'Candidates / %.2f' % ( 0.499999 / float(nBins) )
         , frameOpts  = dict( Bins = nBins, Range = ( 0., 0.499999 ), Name = estWTagSS.GetName() )
         , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.7, LineWidth = 3 )
         #, pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm )
        )

for it, canv in enumerate( canvsOS + canvsSS ) :
    canv.Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len( canvsOS + canvsSS ) - 1 else '' ) )
