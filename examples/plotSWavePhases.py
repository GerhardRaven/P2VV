from math import sqrt

from ROOT import gStyle
from P2VVLoad import LHCbStyle
gStyle.SetEndErrorSize(4)

plotsFilePath   = 'SWavePhases.ps'
deltaSAxisRange = ( None, None ) # ( -1.5, 5.0 )
LHCbText1       = 'LHCb'
LHCbText2       = '' # '#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}'
drawLegend      = False

#massBins = [ 988., 1008., 1020., 1032., 1050. ]
massBins = [ 990., 1008., 1016., 1020., 1024., 1032., 1050. ]

# paper2012 values
deltaSVals     = [ 1.31, 0.77, 0.51, -0.51, -0.46, -0.65 ]
deltaSLowErrs  = [ 0.50, 0.24, 0.36,  0.38,  0.26,  0.23 ]
deltaSHighErrs = [ 0.79, 0.39, 1.41,  0.26,  0.19,  0.19 ]

# 6 bins, new analysis, new n-tuple, 1.0 1/fb
#deltaSVals     = [ 1.33767, 0.771849, 0.512522, -0.505443, -0.456398, -0.650681 ]
#deltaSLowErrs  = [ 0.513406, 0.23602, 0.309024, 0.353726, 0.260243, 0.223321 ]
#deltaSHighErrs = [ 0.756159, 0.378041, 1.37154, 0.21097,  0.185225, 0.184434 ]

# 4 bins, new analysis, new n-tuple, 1.0 1/fb
#deltaSVals     = [ 1.41716, 0.597858, -0.522053, -0.642156 ]
#deltaSLowErrs  = [ 0.524188, 0.194899, 0.242153, 0.224605 ]
#deltaSHighErrs = [ 0.661925, 0.35907,  0.165155, 0.179614 ]

# 4 bins, old analysis, new n-tuple, 1.0 1/fb
#deltaSVals     = [ 1.54757, 0.34823, -0.514583, -0.590583 ]
#deltaSLowErrs  = [ 0.447122, 0.108403, 0.355628, 0.223169 ]
#deltaSHighErrs = [ 0.454733, 0.136747, 0.1852,   0.194338 ]

# 4 bins, new analysis, new n-tuple, 0.37 1/fb
#deltaSVals     = [ 0.514488, 0.677262, -0.33752, -0.780554 ]
#deltaSLowErrs  = [ 0.512267, 0.36539, 0.542154, 2.13257  ]
#deltaSHighErrs = [ 0.84108,  2.32716, 0.187881, 0.280865 ]

# 4 bins, old analysis, new n-tuple, 0.37 1/fb
#deltaSVals     = [ 1.7303, 0.404123, -0.276128, -0.544511 ]
#deltaSLowErrs  = [ 1.02906,  0.167813, 0.254451, 0.298114 ]
#deltaSHighErrs = [ 0.863687, 0.240227, 0.192855, 0.251639 ]

# 4 bins, old analysis, old n-tuple, 0.37 1/fb
#deltaSVals     = [ 2.47014, 0.259854, -0.070761, -1.00204 ]
#deltaSLowErrs  = [ 0.742791, 0.135544, 0.215416, 0.595835 ]
#deltaSHighErrs = [ 0.47209,  0.156742, 0.19807,  0.344092 ]

# theory values
theoryVals = None
#theoryVals = [ 3.0221, 2.8056, 2.0224, 0.8197, 0.2828, 0.1118 ]
#delS = -0.0955 - 1.5589
#delS = ( sum(deltaSVals) - sum(theoryVals) ) / float( len(deltaSVals) )
#delS = -3.
#theoryVals = [ delS + val for val in theoryVals ]

from P2VVGeneralUtils import plotSWavePhases
SWavePhaseCanv = plotSWavePhases( DeltaSAxisRange = deltaSAxisRange, LHCbTextLine1 = LHCbText1, LHCbTextLine2 = LHCbText2
                                 , MassBins = massBins, DeltaSValues = deltaSVals, TheoryValues = theoryVals
                                 , DeltaSLowErrors = deltaSLowErrs, DeltaSHighErrors = deltaSHighErrs, DrawLegend = drawLegend
                                )
SWavePhaseCanv.Print(plotsFilePath)

