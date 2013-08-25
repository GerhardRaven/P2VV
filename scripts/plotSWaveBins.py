from math import sqrt
from P2VV.Load import LHCbStyle

plotsFilePath   = 'SWaveBins.ps'
LHCbText       = ''#'LHCb'
drawLegend      = True
grayScale       = False

#massBins = [ 990., 1008., 1016., 1020., 1024., 1032., 1050. ]
massBins = [ 990., 1008., 1020., 1032., 1050. ]

#fSVals     = [ 0.227, 0.067, 0.008, 0.016, 0.055, 0.168 ]
#fSLowErrs  = [ 0.076, 0.029, 0.009, 0.011, 0.026, 0.047 ]
#fSHighErrs = [ 0.083, 0.031, 0.015, 0.013, 0.028, 0.048 ]
fSVals     = [ 0.228, 0.026, 0.026, 0.169 ]
fSLowErrs  = [ 0.076, 0.014, 0.013, 0.047 ]
fSHighErrs = [ 0.083, 0.016, 0.014, 0.048 ]

#eventVals     = [ 135., 196.,  81., 134., 193., 358. ]
#eventLowErrs  = [  46.,  85.,  91.,  92.,  91., 100. ]
#eventHighErrs = [  50.,  91., 151., 109.,  98., 103. ]
eventVals     = [ 137., 340., 311., 361. ]
eventLowErrs  = [  46., 183., 155., 101. ]
eventHighErrs = [  50., 209., 167., 103. ]

#deltaSVals     = [ 1.31, 0.77, 0.51, -0.51, -0.46, -0.65 ]
#deltaSLowErrs  = [ 0.50, 0.24, 0.36,  0.38,  0.26,  0.23 ]
#deltaSHighErrs = [ 0.79, 0.39, 1.41,  0.26,  0.19,  0.19 ]
deltaSVals     = [ 1.33, 0.60, -0.53, -0.64 ]
deltaSLowErrs  = [ 0.51, 0.22,  0.26,  0.23 ]
deltaSHighErrs = [ 0.78, 0.36,  0.20,  0.19 ]

from P2VV.Utilities.Plotting import plotSWaveBins
SWavePhaseCanv = plotSWaveBins(  Mode = 'fractions', LHCbTextLine1 = LHCbText, GrayScale = grayScale
                               , MassBins = massBins, SValues = fSVals, SLowErrors = fSLowErrs, SHighErrors = fSHighErrs
                              )
SWavePhaseCanv.Print(plotsFilePath + '(')

SWavePhaseCanv = plotSWaveBins(  Mode = 'events', LHCbTextLine1 = LHCbText, GrayScale = grayScale, SAxisRange = ( 0., 80. )
                               , MassBins = massBins, SValues = eventVals, SLowErrors = eventLowErrs, SHighErrors = eventHighErrs
                              )
SWavePhaseCanv.Print(plotsFilePath)

SWavePhaseCanv = plotSWaveBins(  Mode = 'phases', LHCbTextLine1 = LHCbText, GrayScale = grayScale, DrawLegend = drawLegend
                               , MassBins = massBins, SValues = deltaSVals, SLowErrors = deltaSLowErrs, SHighErrors = deltaSHighErrs
                              )
SWavePhaseCanv.Print(plotsFilePath + ')')
