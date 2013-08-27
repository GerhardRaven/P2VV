# script parameters
#dataSetFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2011Reco12_4KKMassBins_2TagCats_tagUntag.root'
dataSetFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets2012Reco14_4KKMassBins_2TagCats_tagUntag.root'
#dataSetFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSetsMC11a_4KKMassBins_2TagCats_tagUntag.root'
dataSetName     = 'JpsiKK_sigSWeight'
dataSetCats     = [  ( 'UU', 'tagCatP2VVOS0_tagCatP2VVSS0' )
                   , ( 'UT', 'tagCatP2VVOS0_tagCatP2VVSS1' )
                   , ( 'TU', 'tagCatP2VVOS1_tagCatP2VVSS0' )
                   , ( 'TT', 'tagCatP2VVOS1_tagCatP2VVSS1' )
                  ]
plotsFilePath   = 'plots/meanAcceptance.ps'

varNames = dict( time = 'time', ctk = 'helcosthetaK', ctl = 'helcosthetaL', phi = 'helphi' )
xTitles = dict(  time = 'Decay time [ps]'
               , ctk  = 'cos#kern[0.1]{#theta_{K}}'
               , ctl  = 'cos#kern[0.1]{#theta_{#mu}}'
               , phi  = '#varphi_{h} [rad]'
              )
yTitles = dict(  time = 'Candidates / (%.2g ps)'
               , ctk  = 'Candidates / %.2g'
               , ctl  = 'Candidates / %.2g'
               , phi  = 'Candidates / (%.2g#kern[0.3]{#pi} rad)'
              )
asymRange = ( -0.3, +0.3 )
numBins = 15

# workspace
from P2VV.RooFitWrappers import RooObject
worksp = RooObject( workspace = 'JpsiKK' ).ws()

# read data sets from file
dataSets = { }
from itertools import product
from P2VV.Utilities.DataHandling import readData
for cat in dataSetCats :
    dataSets[cat[0]] = readData( filePath = dataSetFilePath, dataSetName = dataSetName + '_' + cat[1],  NTuple = False )

# create histograms
from math import exp, log
from array import array
timeRange = ( 0.3, 12., numBins )
timeBinLog = log( timeRange[1] / timeRange[0] ) / float( timeRange[2] )
timeBins = array( 'd', [ timeRange[0] * exp( float(it) * timeBinLog ) for it in range( timeRange[2] + 1 ) ] )

from P2VV.Load import LHCbStyle
from ROOT import TH1D
from math import pi
TH1D.SetDefaultSumw2(True)
timeHists = dict( [ ( cat[0], TH1D( '%sHistTime' % cat[0], '%sHistTime' % cat[0], timeRange[2], timeBins ) ) for cat in dataSetCats ] )
ctkHists  = dict( [ ( cat[0], TH1D( '%sHistCtk'  % cat[0], '%sHistCtk'  % cat[0], numBins, -1., +1. ) ) for cat in dataSetCats ] )
ctlHists  = dict( [ ( cat[0], TH1D( '%sHistCtl'  % cat[0], '%sHistCtl'  % cat[0], numBins, -1., +1. ) ) for cat in dataSetCats ] )
phiHists  = dict( [ ( cat[0], TH1D( '%sHistPhi'  % cat[0], '%sHistPhi'  % cat[0], numBins, -pi, +pi ) ) for cat in dataSetCats ] )

# fill histograms
for cat in dataSetCats :
    numEv = dataSets[cat[0]].sumEntries()
    #numEv = dataSets[cat[0]].sumEntries('hlt1_excl_biased_dec==1')
    for obsSet in dataSets[cat[0]] :
        #if obsSet.getCatIndex('hlt1_excl_biased_dec') != 1 : continue

        timeVal = obsSet.getRealValue(varNames['time'])
        timeWidthRatio = timeHists[cat[0]].GetBinWidth(1) / timeHists[cat[0]].GetBinWidth( timeHists[cat[0]].FindBin(timeVal) )

        timeHists[cat[0]].Fill( timeVal,                              dataSets[cat[0]].weight() / numEv * timeWidthRatio )
        ctkHists[cat[0]].Fill(  obsSet.getRealValue(varNames['ctk']), dataSets[cat[0]].weight() / numEv )
        ctlHists[cat[0]].Fill(  obsSet.getRealValue(varNames['ctl']), dataSets[cat[0]].weight() / numEv )
        phiHists[cat[0]].Fill(  obsSet.getRealValue(varNames['phi']), dataSets[cat[0]].weight() / numEv )

# draw histograms and asymmetries
from math import sqrt
from ROOT import TCanvas, kFullDotLarge, TGraphAsymmErrors as Graph
canvs = [ ]
UUObsSet = dataSets['UU'].get()
for var, hists in zip( [ 'time', 'ctk', 'ctl', 'phi' ], [ timeHists, ctkHists, ctlHists, phiHists ] ) :
    # get histogram minimum and maximum
    histMin = min( [ hist.GetMinimum() for hist in hists.itervalues() if hist.GetMinimum() > 0. ] )
    histMax = max( [ hist.GetMaximum() for hist in hists.itervalues() ] )

    # get horizontal values and errors for asymmetries
    binMeans    = array( 'd', [ ] )
    binErrsLow  = array( 'd', [ ] )
    binErrsHigh = array( 'd', [ ] )
    for bin in range( 1, hists['UU'].GetNbinsX() + 1 ) :
        binLow  = hists['UU'].GetBinLowEdge(bin)
        binHigh = hists['UU'].GetBinLowEdge(bin) + hists['UU'].GetBinWidth(bin)
        binMeans.append( dataSets['UU'].mean( UUObsSet.find(varNames[var]), '{0:s}>={1:.6g} && {0:s}<{2:.6g}'\
                                                                            .format( varNames[var], binLow, binHigh ) ) )
        binErrsLow.append( binMeans[-1] - binLow )
        binErrsHigh.append( binHigh - binMeans[-1] )

    for cat in dataSetCats :
        # canvas for distribution
        canvs.append( TCanvas( hists[cat[0]].GetName() ) )
        canvs[-1].SetLeftMargin(0.18)
        canvs[-1].SetRightMargin(0.05)
        canvs[-1].SetBottomMargin(0.18)
        canvs[-1].SetTopMargin(0.05)
        if var == 'time' : canvs[-1].SetLogy(True)

        # distribution histogram settings
        binWidth = hists[cat[0]].GetBinWidth(1) / ( pi if var == 'phi' else 1. )
        hists[cat[0]].GetXaxis().SetTitle( xTitles[var] )
        hists[cat[0]].GetYaxis().SetTitle( yTitles[var] % binWidth )
        hists[cat[0]].GetXaxis().SetTitleOffset(1.15)
        hists[cat[0]].GetYaxis().SetTitleOffset( 1.00 if var == 'time' else 1.15 )
        hists[cat[0]].SetMinimum( 0.1 * histMin if var == 'time' else 0. )
        hists[cat[0]].SetMaximum( ( 2. if var == 'time' else 1.1 ) * histMax )
        hists[cat[0]].SetMarkerStyle(kFullDotLarge)
        hists[cat[0]].SetMarkerSize(0.7)
        hists[cat[0]].SetLineWidth(3)

        # draw distribution
        hists[cat[0]].Draw('E1')

        if var == 'time' :
            # also draw distribution on double log-scale for decay time
            canvs.append( TCanvas( hists[cat[0]].GetName() + '_log' ) )
            canvs[-1].SetLeftMargin(0.18)
            canvs[-1].SetRightMargin(0.05)
            canvs[-1].SetBottomMargin(0.18)
            canvs[-1].SetTopMargin(0.05)
            canvs[-1].SetLogx(True)
            canvs[-1].SetLogy(True)
            hists[cat[0]].Draw('E1')

        if cat[0] == 'UU' : continue

        # canvas for asymmetry
        canvs.append( TCanvas( cat[0] + '_A' + var ) )
        canvs[-1].SetLeftMargin(0.18)
        canvs[-1].SetRightMargin(0.05)
        canvs[-1].SetBottomMargin(0.18)
        canvs[-1].SetTopMargin(0.05)
        if var == 'time' : canvs[-1].SetLogx(True)

        # create asymmetry graph
        asymVals = array( 'd', [ ] )
        asymErrs = array( 'd', [ ] )
        for bin in range( 1, hists['UU'].GetNbinsX() + 1 ) :
            val   = hists[cat[0]].GetBinContent(bin)
            err   = hists[cat[0]].GetBinError(bin)
            valUU = hists['UU'].GetBinContent(bin)
            errUU = hists['UU'].GetBinError(bin)
            asymVals.append( ( val - valUU ) / ( val + valUU ) )
            asymErrs.append( 0. if val + valUU <= 0. else 2. * sqrt( valUU**2 * err**2 + val**2 * errUU**2 ) / ( val + valUU )**2 )
        graph = Graph( len(binMeans), binMeans, asymVals, binErrsLow, binErrsHigh, asymErrs, asymErrs )
        hists[ cat[0] + '_A' + cat[0] ] = graph

        # asymmetry plot settings
        graph.GetXaxis().SetTitle( xTitles[var] )
        graph.GetYaxis().SetTitle( 'Asymmetry' )
        graph.GetXaxis().SetTitleOffset(1.15)
        graph.GetYaxis().SetTitleOffset(1.15)
        graph.GetXaxis().SetLimits( hists[cat[0]].GetXaxis().GetXmin(), hists[cat[0]].GetXaxis().GetXmax() )
        graph.SetMinimum(asymRange[0])
        graph.SetMaximum(asymRange[1])
        graph.SetMarkerStyle(kFullDotLarge)
        graph.SetMarkerSize(0.7)
        graph.SetLineWidth(3)

        # draw asymmetry
        graph.Draw('AP')

# print canvases
for it, canv in enumerate(canvs) :
    canv.Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len(canvs) - 1 else '' ) )
