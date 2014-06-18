transAngles = False
momentsFiles = [ 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Basis_5thOrder' ]
nTupleFilePath  = 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_prodWeight.root'
nTupleName  = 'DecayTree'
plotsFile = 'angularEfficiency.pdf'

LHCbLabel = '' #'LHCb simulation'

numBins = ( 30, 30, 30 )
signMoms = [  ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 4, 0 )
            , ( 4, 0, 0 )#, ( 3, 0, 0 )
            #, ( 2, 2, 0 ), ( 2, 4, 0 )
            , ( 0, 2, 2 )#, ( 0, 3, -2), ( 0, 4, -2 )
            #, ( 4, 2, -1 ), ( 1, 5, 1 ), ( 1, 5, 3 ), ( 5, 4, -1 )
           ]

f10  = 462684. / ( 462684. + 112402. )
f11  = 112402. / ( 462684. + 112402. )
f20  = ( 462684. + 112402. ) / ( 462684. + 112402. + 119190. )
f21  = 119190.               / ( 462684. + 112402. + 119190. )
f2Sc = ( f10 * 3.5923835 + f11 * 3.5971249 ) / 3.5449077
addFactors = [  ( ) ]#, ( f10, f11 ), ( f20, f21 * f2Sc ) ]

from P2VV.Load import P2VVLibrary, LHCbStyle
from P2VV.RooFitWrappers import RooObject
from ROOT import TCanvas, gStyle
#gStyle.SetEndErrorSize(4)
gStyle.SetColorModelPS(1)
gStyle.SetPalette(1)
gStyle.SetNumberContours(50)

xLabels = ( 'cos#kern[0.1]{#theta_{K}}', 'cos#kern[0.1]{#theta_{#mu}}', '#varphi_{h} [rad]'  ) if not transAngles else\
          ( 'cos#kern[0.1]{#psi_{tr}}',  'cos#kern[0.1]{#theta_{tr}}',  '#varphi_{tr} [rad]' )
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
           , (  'Acceptance integral'
              , 'Acceptance integral'
              , 'Acceptance integral'
             )
          ]

ws = RooObject( workspace = 'workspace' ).ws()

if transAngles :
    from P2VV.Parameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
else :
    from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
    angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

angles = [ angleFuncs.angles[ang] for ang in [ 'cpsi', 'ctheta', 'phi' ] ]

from P2VV.Utilities.DataMoments import RealMomentsBuilder
from math import sqrt, pi
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(6) for YIndex0 in range(6)\
                                          for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#indices  = [ ( PIndex, YIndex0, 0 ) for PIndex in range(5) for YIndex0 in range(5) ]
#indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
#                                          for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#indices += [ ( 0, 4, 0 ), ( 0, 4, 2 ), ( 0, 4, 4 ) ]
moments = RealMomentsBuilder()
moments.appendPYList( angleFuncs.angles, indices )
for file, fac in zip( momentsFiles, addFactors ) :
    moments.read( file, AddMoments = fac )
moments.Print(  Scale = 1. / 2. / sqrt(pi)
              , Names = '|'.join( 'p2vvab_%d0%d%s%d' % ( moms[0], moms[1], 'm' if moms[2] < 0 else '', moms[2] ) for moms in signMoms )
             )
moments.Print( Scale = 1. / 2. / sqrt(pi), MinSignificance = 2.5 )

momFuncTerms = moments.buildPDFTerms( CoefNamePrefix = 'transC_' if transAngles else 'helC_'
                                     , Names = '|'.join( 'p2vvab_%d0%d%s%d' % ( moms[0], moms[1], 'm' if moms[2] < 0 else '', moms[2] )\
                                                        for moms in signMoms )
                                    )
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
if LHCbLabel :
    from ROOT import TPaveText
    LHCbText = TPaveText( 0.33, 0.81, 0.68, 0.89, 'BRNDC')
    LHCbText.AddText(LHCbLabel)
    LHCbText.SetFillColor(0)
    LHCbText.SetTextAlign(12)
    LHCbText.SetTextSize(0.072)
    LHCbText.SetBorderSize(0)

# plot efficiency function slices
from P2VV.Utilities.Plotting import plot
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
           , [ ( 0.84, 1.08 ), ( 0.96, 1.20 ), ( 0.86, 1.10 ) ] if not transAngles\
             else [ ( 1.01, 1.25 ), ( 0.96, 1.20 ), ( 0.97, 1.21 ) ]
          ) :
    pad.SetLeftMargin(0.28)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.18)
    pad.SetTopMargin(0.05)
    plot(  pad, angle, None, momFunc
         #, addPDFs      = [ momFunc1, momFunc2 ]
         , xTitle       = xTitle
         , yTitle       = yTitle
         , yScale       = yScale
         , xTitleOffset = 1.1
         , yTitleOffset = 1.0
         , frameOpts    = dict( Title = angle.GetTitle() )
         , pdfOpts      = dict( LineColor = kBlue, LineWidth = 3 )
         #, addPDFsOpts  = [ dict( LineColor = kRed, LineWidth = 3 ), dict( LineColor = kGreen + 2, LineWidth = 3 ) ]
        )
    if LHCbLabel : LHCbText.Draw()

# plot efficiency function integrals
from ROOT import RooArgSet
canvs += [  TCanvas( 'cpsiIntCanv',   'Angular Efficiency' )
          , TCanvas( 'cthetaIntCanv', 'Angular Efficiency' )
          , TCanvas( 'phiIntCanv',    'Angular Efficiency' )
         ]
integrals = [  momFunc.createIntegral( RooArgSet( angles[1]._var, angles[2]._var ) )
             , momFunc.createIntegral( RooArgSet( angles[0]._var, angles[2]._var ) )
             , momFunc.createIntegral( RooArgSet( angles[0]._var, angles[1]._var ) )
             , momFunc.createIntegral( RooArgSet( angles[0]._var                 ) )
             , momFunc.createIntegral( RooArgSet( angles[1]._var                 ) )
             , momFunc.createIntegral( RooArgSet( angles[2]._var                 ) )
            ]
for ( pad, func, angle, xTitle, yTitle, yScale, norm )\
    in zip(  canvs[ 3 : ]
           , integrals[ : 3 ]
           , angles
           , xLabels
           , yLabels[1]
           , [ ( 0.89, 1.11 ), ( 0.9345, 1.1655 ), ( 0.89, 1.11 ) ] if not transAngles\
             else [ ( 0.88, 1.12 ), ( 0.88, 1.12 ), ( 0.8976, 1.1424 ) ]
           , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
          ) :
    pad.SetLeftMargin(0.28)
    pad.SetRightMargin(0.05)
    pad.SetBottomMargin(0.18)
    pad.SetTopMargin(0.05)
    plot(  pad, angle, None, func
         #, addPDFs      = [ momFuncAdd ]
         , xTitle       = xTitle
         , yTitle       = yTitle
         , yScale       = yScale
         , xTitleOffset = 1.1
         , yTitleOffset = 1.0
         , frameOpts    = dict( Title = angle.GetTitle() )
         , pdfOpts      = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm )
         #, addPDFsOpts  = [ dict( LineColor = kRed ) ]
        )
    if LHCbLabel : LHCbText.Draw()

if nTupleFilePath :
    # read n-tuple from file
    from ROOT import TFile
    nTupleFile = TFile.Open(nTupleFilePath)
    nTuple = nTupleFile.Get(nTupleName)

    # convert n-tuple to dataset
    from ROOT import RooDataSet
    obsSet = RooArgSet( ang._var for ang in angles )
    from P2VV.RooFitWrappers import RealVar
    evWeight = RealVar( Name = 'prodWeight', Value = 1. )
    obsSet.add(evWeight._var)
    data = RooDataSet( 'effWeightData', 'effWeightData', obsSet, Import = nTuple, WeightVar = ( 'prodWeight', True ) )
    nTupleFile.Close()

    # plot binned efficiency function integrals
    canvs += [  TCanvas( 'cpsiDataCanv',   'Angular Efficiency' )
              , TCanvas( 'cthetaDataCanv', 'Angular Efficiency' )
              , TCanvas( 'phiDataCanv',    'Angular Efficiency' )
             ]
    for ( pad, angle, xTitle, yTitle, yScale, nBins )\
        in zip(  canvs[ 6 : ]
               , angles
               , xLabels
               , yLabels[1]
               , [ ( 0.89, 1.11 ), ( 0.9345, 1.1655 ), ( 0.89, 1.11 ) ] if not transAngles\
                 else [ ( 0.88, 1.12 ), ( 0.88, 1.12 ), ( 0.8976, 1.1424 ) ]
               , numBins
              ) :
        pad.SetLeftMargin(0.28)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.18)
        pad.SetTopMargin(0.05)
        plot(  pad, angle, data, None
             , xTitle       = xTitle
             , yTitle       = yTitle
             , yScale       = yScale
             , xTitleOffset = 1.1
             , yTitleOffset = 1.0
             , frameOpts    = dict( Title = angle.GetTitle(), Bins = nBins )
             , dataOpts     = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.8, LineWidth = 3
                                   , Rescale = float(nBins) / 8. / pi )
            )
        if LHCbLabel : LHCbText.Draw()

    # plot (binned) efficiency function integrals
    canvs += [  TCanvas( 'cpsiDataIntCanv',   'Angular Efficiency' )
              , TCanvas( 'cthetaDataIntCanv', 'Angular Efficiency' )
              , TCanvas( 'phiDataIntCanv',    'Angular Efficiency' )
             ]
    for ( pad, func, angle, xTitle, yTitle, yScale, norm, nBins )\
        in zip(  canvs[ 9 : ]
               , integrals[ : 3 ]
               , angles
               , xLabels
               , yLabels[1]
               , [ ( 0.89, 1.11 ), ( 0.9345, 1.1655 ), ( 0.89, 1.11 ) ] if not transAngles\
                 else [ ( 0.88, 1.12 ), ( 0.88, 1.12 ), ( 0.8976, 1.1424 ) ]
               , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
               , numBins
              ) :
        pad.SetLeftMargin(0.28)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.18)
        pad.SetTopMargin(0.05)
        plot(  pad, angle, data, func
             , xTitle       = xTitle
             , yTitle       = yTitle
             , yScale       = yScale
             , xTitleOffset = 1.1
             , yTitleOffset = 1.0
             , frameOpts    = dict( Title = angle.GetTitle(), Bins = nBins )
             , dataOpts     = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.8, LineWidth = 3
                                   , Rescale = float(nBins) / 8. / pi )
             , pdfOpts      = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm )
            )
        if LHCbLabel : LHCbText.Draw()

# plot 2D efficiency function integrals
from ROOT import RooArgList, RooConstVar, RooProduct
oneOverTwo   = RooConstVar( 'oneOverTwo',   'oneOverTwo',   0.5      )
oneOverTwoPi = RooConstVar( 'oneOverTwoPi', 'oneOverTwoPi', 0.5 / pi )
integrals += [  RooProduct( 'cthetaPhiEffInt',  'cthetaPhiEffInt',  RooArgList( integrals[-3], oneOverTwo   ) )
              , RooProduct( 'cpsiPhiEffInt',    'cpsiPhiEffInt',    RooArgList( integrals[-2], oneOverTwo   ) )
              , RooProduct( 'cpsiCthetaEffInt', 'cpsiCthetaEffInt', RooArgList( integrals[-1], oneOverTwoPi ) )
             ]
integrals += [  integrals[-3].asTF( RooArgList( angles[1], angles[2] ) )
              , integrals[-2].asTF( RooArgList( angles[0], angles[2] ) )
              , integrals[-1].asTF( RooArgList( angles[0], angles[1] ) )
             ]
canvs += [  TCanvas( 'cthetaPhiIntCanv',  'Angular Efficiency' )
          , TCanvas( 'cpsiPhiIntCanv',    'Angular Efficiency' )
          , TCanvas( 'cpsiCthetaIntCanv', 'Angular Efficiency' )
         ]
for integ, canv, xTitle, yTitle, zTitle, zScale\
        in zip( integrals[ -3 : ], canvs[ -3 : ], [ xLabels[1], xLabels[0], xLabels[0] ], [ xLabels[2], xLabels[2], xLabels[1] ]
               , [ yLabels[2][0], yLabels[2][1], yLabels[2][2] ]
               , [ ( 0.94, 1.17 ), ( 0.94, 1.17 ), ( 0.94, 1.17 ) ]
              ) :
    canv.cd()
    canv.SetLeftMargin(0.17)
    canv.SetRightMargin(0.12)
    canv.SetBottomMargin(0.10)
    canv.SetTopMargin(0.05)

    integ.SetMinimum( zScale[0] )
    integ.SetMaximum( zScale[1] )
    integ.GetXaxis().SetTitle(xTitle)
    integ.GetYaxis().SetTitle(yTitle)
    integ.GetZaxis().SetTitle(zTitle)
    integ.GetXaxis().SetTitleOffset(1.4)
    integ.GetYaxis().SetTitleOffset(1.4)
    integ.GetZaxis().SetTitleOffset(1.2)

    integ.Draw('SURF4')

    canv.SetTheta(30)
    canv.SetPhi(55)
    canv.Update()

if nTupleFilePath :
    # plot 2D efficiency binned function integrals
    from ROOT import TH2D
    integrals += [  TH2D( 'cthetaPhiEffHist',  'Angular Efficiency', 10, -1., +1., 10, -pi, +pi )
                  , TH2D( 'cpsiPhiEffHist',    'Angular Efficiency', 10, -1., +1., 10, -pi, +pi )
                  , TH2D( 'cpsiCthetaEffHist', 'Angular Efficiency', 10, -1., +1., 10, -1., +1. )
                 ]
    data.fillHistogram( integrals[-3], RooArgList( angles[1], angles[2] ) )
    data.fillHistogram( integrals[-2], RooArgList( angles[0], angles[2] ) )
    data.fillHistogram( integrals[-1], RooArgList( angles[0], angles[1] ) )
    integrals[-3].Scale( 100. / 8. / pi )
    integrals[-2].Scale( 100. / 8. / pi )
    integrals[-1].Scale( 100. / 8. / pi )

    canvs += [  TCanvas( 'cthetaPhiDataCanv',  'Angular Efficiency' )
              , TCanvas( 'cpsiPhiDataCanv',    'Angular Efficiency' )
              , TCanvas( 'cpsiCthetaDataCanv', 'Angular Efficiency' )
             ]
    for hist, canv, xTitle, yTitle, zTitle, zScale\
            in zip( integrals[ -3 : ], canvs[ -3 : ], [ xLabels[1], xLabels[0], xLabels[0] ], [ xLabels[2], xLabels[2], xLabels[1] ]
               , [ yLabels[2][0], yLabels[2][1], yLabels[2][2] ]
               , [ ( 0.94, 1.17 ), ( 0.94, 1.17 ), ( 0.94, 1.17 ) ]
              ) :
        canv.cd()
        canv.SetLeftMargin(0.17)
        canv.SetRightMargin(0.12)
        canv.SetBottomMargin(0.10)
        canv.SetTopMargin(0.05)

        hist.SetMinimum( zScale[0] )
        hist.SetMaximum( zScale[1] )
        hist.SetXTitle(xTitle)
        hist.SetYTitle(yTitle)
        hist.SetZTitle(zTitle)
        hist.SetTitleOffset( 1.4, 'X' )
        hist.SetTitleOffset( 1.4, 'Y' )
        hist.SetTitleOffset( 1.2, 'Z' )

        hist.Draw('LEGO2')

        canv.SetTheta(30)
        canv.SetPhi(55)
        canv.Update()

for canvIt, canv in enumerate(canvs) : canv.Print( plotsFile + ( '(' if canvIt == 0 else ')' if canvIt == len(canvs) - 1 else '' ) )
