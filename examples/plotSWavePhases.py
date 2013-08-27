from math import sqrt
from P2VV.Load import LHCbStyle

plotsFilePath   = 'SWavePhases.ps'
deltaSAxisRange = ( None, None )
LHCbText1       = 'LHCb'
LHCbText2       = ''#'#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}'
drawLegend      = False
grayScale       = False

massBins = [ 990., 1008., 1016., 1020., 1024., 1032., 1050. ]

deltaSVals     = [ 1.31, 0.77, 0.51, -0.51, -0.46, -0.65 ]
deltaSLowErrs  = [ 0.50, 0.24, 0.36,  0.38,  0.26,  0.23 ]
deltaSHighErrs = [ 0.79, 0.39, 1.41,  0.26,  0.19,  0.19 ]

# theory values
theoryVals = None
#theoryVals = [ 3.0221, 2.8056, 2.0224, 0.8197, 0.2828, 0.1118 ]
#delS = ( sum(deltaSVals) - sum(theoryVals) ) / float( len(deltaSVals) )
#theoryVals = [ delS + val for val in theoryVals ]

from P2VV.Utilities.Plotting import plotSWaveBins
SWavePhaseCanv = plotSWaveBins( SAxisRange = deltaSAxisRange, LHCbTextLine1 = LHCbText1, LHCbTextLine2 = LHCbText2, MassBins = massBins
                               , SValues = deltaSVals, TheoryValues = theoryVals, SLowErrors = deltaSLowErrs, SHighErrors = deltaSHighErrs
                               , DrawLegend = drawLegend, GrayScale = grayScale
                              )
SWavePhaseCanv.Print(plotsFilePath)
