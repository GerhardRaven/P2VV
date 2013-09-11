# plot parameters
dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_4KKMassBins_2TagCats.root'
dataSetName = 'JpsiKK_cbkgSWeight'
plotsFilePath = 'bkgAnglePlots.ps'
nBins = dict( cpsi = 20, ctheta = 40, phi = 20 )
angleTitles = dict( cpsi = 'cos#kern[0.1]{#theta_{K}}', ctheta = 'cos#kern[0.1]{#theta_{#mu}}', phi = '#varphi_{h} [rad]' )
fitOpts = dict( NumCPU = 8, Optimize = 2, Timer = True, Minimizer = 'Minuit2', Strategy = 1, Offset = True )

# workspace
from P2VV.RooFitWrappers import RooObject, RealVar
ws = RooObject( workspace = 'angleWorkspace' ).ws()

# read data and get observables
from P2VV.Utilities.DataHandling import readData
cbkgData = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
observables = { }
for key, name in [ ( 'cpsi', 'helcosthetaK' ), ( 'ctheta', 'helcosthetaL' ), ( 'phi', 'helphi' ) ] : observables[key] = RealVar(name)

# background angles PDF
angPdfContainer = { }
from P2VV.Parameterizations.FullPDFs import buildBs2JpsiphiCombBkgAnglesPdf
cbkgAnglePDF = buildBs2JpsiphiCombBkgAnglesPdf( angPdfContainer, observables = observables, bkgAnglePdfType = 'hybrid'
                                               , cbkgData = cbkgData, transAngles = False, numAngleBins = ( 10, 24, 5 )
                                               , fitOptions = fitOpts )

# set LHCb plot style
from P2VV.Load import LHCbStyle

# create canvases
from ROOT import TCanvas
canvs = [  TCanvas( 'bkgAngles1D_ctk' )
         , TCanvas( 'bkgAngles1D_ctl' )
         , TCanvas( 'bkgAngles1D_phi' )
         , TCanvas( 'bkgAngles2D_ctk_ctl' )
         , TCanvas( 'bkgAngles2D_ctk_phi' )
         , TCanvas( 'bkgAngles2D_phi_ctl' )
        ]

# plot background angles
from P2VV.Utilities.Plotting import plot
from ROOT import kBlue
for ( name, canv ) in zip( [ 'cpsi', 'ctheta', 'phi' ], canvs[ : 3 ] ) :
    plot(  canv, observables[name], cbkgData, cbkgAnglePDF, xTitle = angleTitles[name], yTitleOffset = 1.5
         , frameOpts  = dict( Bins = nBins[name], Title = observables[name].GetTitle() )
         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3  )
        )

# plot 2-D angular distributions
from P2VV.Utilities.Plotting import _P2VVPlotStash
for angle0, angle1, canv in [ ( 'ctheta', 'cpsi', canvs[3] ), ( 'phi', 'cpsi', canvs[4] ), ( 'ctheta', 'phi', canvs[5] ) ] :
    bkgAngHist = cbkgData.createHistogram( observables[angle0]._var, observables[angle1]._var, nBins[angle0], nBins[angle1] )
    _P2VVPlotStash.append(bkgAngHist)
    bkgAngHist.SetStats(False)
    bkgAngHist.SetTitle( '%s vs. %s' % ( angleTitles[angle0], angleTitles[angle1] ) )
    bkgAngHist.SetMinimum(0.)
    bkgAngHist.GetXaxis().SetTitle( angleTitles[angle0] )
    bkgAngHist.GetYaxis().SetTitle( angleTitles[angle1] )
    bkgAngHist.GetXaxis().SetLabelOffset(0.01)
    bkgAngHist.GetYaxis().SetLabelOffset(0.008)
    bkgAngHist.GetXaxis().SetTitleOffset(1.8)
    bkgAngHist.GetYaxis().SetTitleOffset(1.8)
    bkgAngPad = canv.cd()
    bkgAngPad.SetLeftMargin(0.08)
    bkgAngPad.SetRightMargin(0.05)
    bkgAngPad.SetBottomMargin(0.05)
    bkgAngPad.SetTopMargin(0.05)
    bkgAngHist.Draw('lego2')

# print canvases to file
for it, canv in enumerate(canvs) : canv.Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len(canvs) - 1 else '' ) )
