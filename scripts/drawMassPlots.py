plotsFilePathIn = 'temp_plots.root'
plotsFilePathOut = 'temp_plots.ps'
plotLabelText = 'LHCb'

from ROOT import gStyle
from P2VV.Load import LHCbStyle
gStyle.SetLineStyleString( 5, ' 40 20 10 20'  )
gStyle.SetLineStyleString( 7, ' 40 20'        )
gStyle.SetLineStyleString( 9, ' 100 20'       )

if plotLabelText :
    from ROOT import TLatex
    plotLabel = TLatex()
    plotLabel.SetTextAlign(32)
    plotLabel.SetTextSize(0.072)

from ROOT import TFile, TPad, TCanvas
plotsFile = TFile.Open(plotsFilePathIn)
dumCanv = TCanvas('dummy')
dumCanv.Print( plotsFilePathOut + '[' )
for plotIt, name in enumerate( [ 'mass', 'mass_signal', 'mass_left', 'mass_right', 'mass_peakBkg' ] ) :
    plot = plotsFile.Get(name)
    residPlot = plotsFile.Get( name + '_resid' )
    canv = TCanvas(name)

    # draw mass plot
    canv.cd()
    plotName = name + '_plot1'
    plotPad = TPad( plotName, plotName, 0, 0.32, 1, 1 )
    if plotIt == 0 : plotPad.SetLogy(1)
    plotPad.SetNumber(1)
    plotPad.SetTopMargin(0.04)
    plotPad.SetBottomMargin(0.04)
    plotPad.Draw()
    canv.cd(1)
    plot.Draw()

    # draw residuals
    canv.cd()
    residName = name + '_resid1'
    residPad = TPad( residName, residName, 0, 0, 1, 0.32 )
    residPad.SetNumber(2)
    residPad.SetTopMargin(0.)
    residPad.SetBottomMargin(0.4)
    residPad.Draw()
    canv.cd(2)
    residPlot.Draw()

    if plotLabelText and plotIt < 2 :
        canv.cd()
        plotLabel.DrawLatexNDC( 0.87, 0.88, plotLabelText )

    # write canvas to file
    canv.Print(plotsFilePathOut)

dumCanv.Print( plotsFilePathOut + ']' )
plotsFile.Close()
