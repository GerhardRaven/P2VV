transAngles = False
momentsFiles = [
                  'hel_UB_UT_trueTime_BkgCat050_KK30_allOrders_Basis'
#                , 'hel_UB_UT_trueTime_BkgCat050_KK30_alt_allOrders_Basis'
#                , 'hel_UB_UT_trueTime_BkgCat050_KK30_PHSP_allOrders_Basis'
               ]
dataFile  = 'hel_UB_UT_trueTime_BkgCat050_KK30.root'
dataName  = 'DecayTree'
plotsFile = 'angularEfficiency.ps'

LHCbLabel = 'LHCb simulation'

numBins = ( 20, 20, 20 )

f10  = 462684. / ( 462684. + 112402. )
f11  = 112402. / ( 462684. + 112402. )
f20  = ( 462684. + 112402. ) / ( 462684. + 112402. + 119190. )
f21  = 119190.               / ( 462684. + 112402. + 119190. )
f2Sc = ( f10 * 3.5923835 + f11 * 3.5971249 ) / 3.5449077
addFactors = [  ( ) ]#, ( f10, f11 ), ( f20, f21 * f2Sc ) ]

from P2VV.Load import P2VVLibrary, LHCbStyle
from P2VV.RooFitWrappers import RooObject
from ROOT import TCanvas, gStyle
gStyle.SetEndErrorSize(4)

xLabels = ( 'cos#kern[0.1]{#theta_{K}}', 'cos#kern[0.1]{#theta_{#mu}}', '#varphi_{h} [rad]'  ) if not transAngles else\
          ( 'cos#kern[0.3]{#psi_{tr}}',  'cos#kern[0.3]{#theta_{tr}}',  '#varphi_{tr} [rad]' )
yLabels = [  (  '#varepsilon_{#Omega}(cos#kern[0.3]{#theta_{K}}, 0, 0) / #LT#varepsilon_{#Omega}#GT'
              , '#varepsilon_{#Omega}(0, cos#kern[0.3]{#theta_{#mu}}, 0) / #LT#varepsilon_{#Omega}#GT'
              , '#varepsilon_{#Omega}(0, 0, #varphi_{h}) / #LT#varepsilon_{#Omega}#GT'
             )
           , (  'Scaled acceptance integral'
              , 'Scaled acceptance integral'
              , 'Scaled acceptance integral'
             )
           #, (  '#int d_{}cos#theta_{#mu} d#varphi_{h} #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
           #   , '#int d_{}cos#theta_{K} d#varphi_{h} #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
           #   , '#int d_{}cos#theta_{K} dcos#theta_{#mu} #varepsilon_{#Omega}(#Omega) / (4 #LT#varepsilon_{#Omega}#GT)'
           #  ) if not transAngles else\
           #  (  '#int d_{}cos#kern[0.3]{#theta_{tr}} d#varphi_{tr} #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
           #   , '#int d_{}cos#kern[0.3]{#psi_{tr}} d#varphi_{tr} #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
           #   , '#int d_{}cos#kern[0.3]{#psi_{tr}} dcos#kern[0.3]{#theta_{tr}} #varepsilon_{#Omega}(#Omega) / (4 #LT#varepsilon_{#Omega}#GT)'
           #  )
          ]

ws = RooObject( workspace = 'workspace' ).ws()

if transAngles :
    from P2VV.Parameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

angles = [ angleFuncs.angles[ang] for ang in [ 'cpsi', 'ctheta', 'phi' ] ]

from P2VV.GeneralUtils import RealMomentsBuilder
from math import sqrt, pi
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(6) for YIndex0 in range(6)\
                                          for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#indices  = [ ( PIndex, YIndex0, 0 ) for PIndex in range(5) for YIndex0 in range(5) ]
moments = RealMomentsBuilder()
moments.appendPYList( angleFuncs.angles, indices )
for file, fac in zip( momentsFiles, addFactors ) :
    moments.read( file, AddMoments = fac )
moments.Print(  Scale = 1. / 2. / sqrt(pi), Names = 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0022'\
                                                    if transAngles else 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0040'
             )
moments.Print( Scale = 1. / 2. / sqrt(pi), MinSignificance = 2. )

momFuncTerms = moments.buildPDFTerms( CoefNamePrefix = 'transC_' if transAngles else 'helC_'
                                     , Names = 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0022'\
                                               if transAngles else 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0040' )
momFunc = momFuncTerms.buildAddition( 'efficiency' + ( 'Trans' if transAngles else 'Hel' ) )

# create efficiency functions with alternative angular values for slices
from ROOT import RooRealVar, RooCustomizer
angles1 = [ RooRealVar(ang._var) for ang in angles ]
angles2 = [ RooRealVar(ang._var) for ang in angles ]
momFuncCust1 = RooCustomizer( momFunc._var, '_1' )
momFuncCust2 = RooCustomizer( momFunc._var, '_2' )
for ang, ang1, ang2 in zip( angles, angles1, angles2 ) :
    momFuncCust1.replaceArg( ang._var, ang1 )
    momFuncCust2.replaceArg( ang._var, ang2 )
momFunc1 = momFuncCust1.build()
momFunc2 = momFuncCust2.build()

for ang, val in zip( angles,  (  0.,  0., 0. ) ) : ang.setVal(val)
for ang, val in zip( angles1, ( +1., +1., 0. ) ) : ang.setVal(val)
for ang, val in zip( angles2, ( -1., +1., 0. ) ) : ang.setVal(val)

print
momFunc.getVariables().Print('v')
print
momFunc1.getVariables().Print('v')
print
momFunc2.getVariables().Print('v')

# LHCb label
from ROOT import TPaveText
LHCbText = TPaveText( 0.31, 0.83, 0.67, 0.92, 'NDC' )
LHCbText.AddText(LHCbLabel)
LHCbText.SetShadowColor(0)
LHCbText.SetFillStyle(0)
LHCbText.SetBorderSize(0)
LHCbText.SetTextAlign(12)

# plot efficiency function slices
from P2VV.GeneralUtils import plot
from ROOT import kBlack, kBlue, kRed, kGreen, kFullDotLarge
canvs = [  TCanvas( 'cpsiFuncCanv',   'Angular Efficiency' )
         , TCanvas( 'cthetaFuncCanv', 'Angular Efficiency' )
         , TCanvas( 'phiFuncCanv',    'Angular Efficiency' )
        ]
for ( pad, angle, xTitle, yTitle, yScale )\
    in zip(  canvs
           , angles
           , xLabels
           , yLabels[0]
           , [ ( 0.84, 1.08 ), ( 0.96, 1.20 ), ( 0.86, 1.10 ) ]
          ) :
    pad.SetLeftMargin(0.28)
    plot(  pad, angle, None, momFunc
         #, addPDFs      = [ momFunc1, momFunc2 ]
         , xTitle       = xTitle
         , yTitle       = yTitle
         , yScale       = yScale
         , yTitleOffset = 1.0
         , frameOpts    = dict( Title = angle.GetTitle() )
         , pdfOpts      = dict( LineColor = kBlue, LineWidth = 4 )
         #, addPDFsOpts  = [ dict( LineColor = kRed, LineWidth = 4 ), dict( LineColor = kGreen + 2, LineWidth = 4 ) ]
        )
    LHCbText.Draw()

# plot efficiency function integrals
from ROOT import RooArgSet
canvs += [  TCanvas( 'cpsiIntCanv',   'Angular Efficiency' )
          , TCanvas( 'cthetaIntCanv', 'Angular Efficiency' )
          , TCanvas( 'phiIntCanv',    'Angular Efficiency' )
         ]
integrals = [  momFunc.createIntegral( RooArgSet( angles[1]._var, angles[2]._var ) )
             , momFunc.createIntegral( RooArgSet( angles[0]._var, angles[2]._var ) )
             , momFunc.createIntegral( RooArgSet( angles[0]._var, angles[1]._var ) )
            ]
for ( pad, func, angle, xTitle, yTitle, yScale, norm )\
    in zip(  canvs[ 3 : ]
           , integrals
           , angles
           , xLabels
           , yLabels[1]
           , [ ( 0.88, 1.12 ), ( 0.9328, 1.1872 ), ( 0.88, 1.12 ) ]
           , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
          ) :
    pad.SetLeftMargin(0.28)
    plot(  pad, angle, None, func
         #, addPDFs      = [ momFuncAdd ]
         , xTitle       = xTitle
         , yTitle       = yTitle
         , yScale       = yScale
         , yTitleOffset = 1.0
         , frameOpts    = dict( Title = angle.GetTitle() )
         , pdfOpts      = dict( LineColor = kBlue, LineWidth = 4, Normalization = norm )
         #, addPDFsOpts  = [ dict( LineColor = kRed ) ]
        )
    LHCbText.Draw()

if dataFile :
    # read data set
    from ROOT import TFile
    dataSetFile = TFile.Open(dataFile)
    data = dataSetFile.Get(dataName)
    dataSetFile.Close()

    # apply efficiency weights to events in data set
    from ROOT import RooDataSet
    data = RooDataSet( 'effWeightData', 'effWeightData', data.get(), Import = data, WeightVar = ( 'effWeight', True ) )

    # plot binned efficiency function integrals
    canvs += [  TCanvas( 'cpsiDataCanv',   'Angular Efficiency' )
              , TCanvas( 'cthetaDataCanv', 'Angular Efficiency' )
              , TCanvas( 'phiDataCanv',    'Angular Efficiency' )
             ]
    for ( pad, angle, xTitle, yTitle, yScale, norm, nBins )\
        in zip(  canvs[ 6 : ]
               , angles
               , xLabels
               , yLabels[1]
               , [ ( 0.88, 1.12 ), ( 0.9328, 1.1872 ), ( 0.88, 1.12 ) ]
               , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
               , numBins
              ) :
        pad.SetLeftMargin(0.28)
        plot(  pad, angle, data, None
             , xTitle       = xTitle
             , yTitle       = yTitle
             , yScale       = yScale
             , yTitleOffset = 1.0
             , frameOpts    = dict( Title = angle.GetTitle(), Bins = nBins )
             , dataOpts     = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.9, LineWidth = 3
                                   , Rescale = norm * float(nBins) / ( angle.getMax() - angle.getMin() ) )
            )
        LHCbText.Draw()

    # plot (binned) efficiency function integrals
    canvs += [  TCanvas( 'cpsiDataIntCanv',   'Angular Efficiency' )
              , TCanvas( 'cthetaDataIntCanv', 'Angular Efficiency' )
              , TCanvas( 'phiDataIntCanv',    'Angular Efficiency' )
             ]
    for ( pad, func, angle, xTitle, yTitle, yScale, norm, nBins )\
        in zip(  canvs[ 9 : ]
               , integrals
               , angles
               , xLabels
               , yLabels[1]
               , [ ( 0.88, 1.12 ), ( 0.9328, 1.1872 ), ( 0.88, 1.12 ) ]
               , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
               , numBins
              ) :
        pad.SetLeftMargin(0.28)
        plot(  pad, angle, data, func
             , xTitle       = xTitle
             , yTitle       = yTitle
             , yScale       = yScale
             , yTitleOffset = 1.0
             , frameOpts    = dict( Title = angle.GetTitle(), Bins = nBins )
             , dataOpts     = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.9, LineWidth = 3
                                   , Rescale = norm * float(nBins) / ( angle.getMax() - angle.getMin() ) )
             , pdfOpts      = dict( LineColor = kBlue, LineWidth = 4, Normalization = norm )
            )
        LHCbText.Draw()

for canvIt, canv in enumerate(canvs) : canv.Print( plotsFile + ( '(' if canvIt == 0 else ')' if canvIt == len(canvs) - 1 else '' ) )
