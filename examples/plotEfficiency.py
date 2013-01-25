helEffMomentsFile = 'hel_UB_UT_trueTime_BkgCat050_KK30_allOrders_Basis'
trEffMomentsFile  = ''
plotsFile = 'angularEfficiency.ps'

from P2VVLoad import P2VVLibrary, LHCbStyle
from RooFitWrappers import RooObject
from ROOT import TCanvas

canvs = [  (  TCanvas( 'ctkFuncCanv',   'Angular Efficiency' )
            , TCanvas( 'ctlFuncCanv',   'Angular Efficiency' )
            , TCanvas( 'phihFuncCanv',  'Angular Efficiency' )
            , TCanvas( 'ctkIntCanv',    'Angular Efficiency' )
            , TCanvas( 'ctlIntCanv',    'Angular Efficiency' )
            , TCanvas( 'phihIntCanv',   'Angular Efficiency' )
           )
         , (  TCanvas( 'cpsiFuncCanv',  'Angular Efficiency' )
            , TCanvas( 'cthFuncCanv',   'Angular Efficiency' )
            , TCanvas( 'phitrFuncCanv', 'Angular Efficiency' )
            , TCanvas( 'cpsiIntCanv',   'Angular Efficiency' )
            , TCanvas( 'cthIntCanv',    'Angular Efficiency' )
            , TCanvas( 'phitrIntCanv',  'Angular Efficiency' )
           )
        ]

xLabels = [  ( 'cos#theta_{K}', 'cos#theta_{#mu}', '#phi_{h} [rad]'  )
           , ( 'cos#psi_{tr}',  'cos#theta_{tr}',  '#phi_{tr} [rad]' )
          ]
yLabels = [  (  '#int d_{}cos#theta_{#mu} d#phi_{h} #varepsilon_{#Omega}(#Omega) / #LT#varepsilon_{#Omega}#GT'
              , '#int d_{}cos#theta_{K} d#phi_{h} #varepsilon_{#Omega}(#Omega) / #LT#varepsilon_{#Omega}#GT'
              , '#int d_{}cos#theta_{K} dcos#theta_{#mu} #varepsilon_{#Omega}(#Omega) / #LT#varepsilon_{#Omega}#GT'
             )
           , (  '#int d_{}cos#theta_{tr} d#phi_{tr} #varepsilon_{#Omega}(#Omega) / #LT#varepsilon_{#Omega}#GT'
              , '#int d_{}cos#psi_{tr} d#phi_{tr} #varepsilon_{#Omega}(#Omega) / #LT#varepsilon_{#Omega}#GT'
              , '#int d_{}cos#psi_{tr} dcos#theta_{tr} #varepsilon_{#Omega}(#Omega) / #LT#varepsilon_{#Omega}#GT'
             )
          ]

for transAngles, canvTuple in zip( [ False, True ], canvs ) :
    if ( not transAngles and helEffMomentsFile ) or ( transAngles and trEffMomentsFile ) :
        ws = RooObject(workspace = ( 'trans' if transAngles else 'hel' ) + 'Workspace').ws()

        if transAngles :
            from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
            angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
        else :
            from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
            angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

        angles = [ angleFuncs.angles[ang] for ang in [ 'cpsi', 'ctheta', 'phi' ] ]

        from P2VVGeneralUtils import RealMomentsBuilder
        from math import sqrt, pi
        indices  = [(PIndex, YIndex0, YIndex1) for PIndex in range(9) for YIndex0 in range(9) for YIndex1 in range(-YIndex0, YIndex0 + 1)]
        moments = RealMomentsBuilder()
        moments.appendPYList( angleFuncs.angles, indices )
        moments.read( trEffMomentsFile if transAngles else helEffMomentsFile )
        moments.Print(  Scale = 1. / 2. / sqrt(pi)
                      , Names = 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0022' if transAngles else 'p2vvab_0000|p2vvab_2000|p2vvab_0020'
                     )
        moments.Print( Scale = 1. / 2. / sqrt(pi), MinSignificance = 3. )

        momFuncTerms = moments.buildPDFTerms( CoefNamePrefix = 'transC_' if transAngles else 'helC_', MinSignificance = 3. )
                #, Names = 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0022' if transAngles else 'p2vvab_0000|p2vvab_2000|p2vvab_0020' )
        momFunc = momFuncTerms.buildAddition( 'efficiency' + ( 'Trans' if transAngles else 'Hel' ) )

        #momFuncTermsAdd = moments.buildPDFTerms( CoefNamePrefix = 'transAddC_' if transAngles else 'helAddC_' )
        #momFuncAdd = momFuncTermsAdd.buildAddition( 'efficiencyAdd' + ( 'Trans' if transAngles else 'Hel' ) )

        momFunc.getVariables().Print('v')

        from P2VVGeneralUtils import plot
        from ROOT import kBlue, kRed
        for ( pad, angle, xTitle, yScale )\
            in zip(  canvTuple[ : 3 ]
                   , angles
                   , xLabels[1] if transAngles else xLabels[0]
                   , [ ( 0.84, 1.08 ), ( 0.96, 1.20 ), ( 0.86, 1.10 ) ]
                  ) :
            pad.SetLeftMargin(0.28)
            plot(  pad, angle, None, momFunc
                 #, addPDFs      = [ momFuncAdd ]
                 , xTitle       = xTitle
                 , yTitle       = '#varepsilon_{#Omega}(#Omega) / #LT#varepsilon_{#Omega}#GT'
                 , yScale       = yScale
                 , yTitleOffset = 0.9
                 , frameOpts    = dict( Title = angle.GetTitle() )
                 , pdfOpts      = dict( LineColor = kBlue, LineWidth = 4 )
                 #, addPDFsOpts  = [ dict( LineColor = kRed ) ]
                )

        from ROOT import RooArgSet
        integrals = [  momFunc.createIntegral( RooArgSet( angles[1]._var, angles[2]._var ) )
                     , momFunc.createIntegral( RooArgSet( angles[0]._var, angles[2]._var ) )
                     , momFunc.createIntegral( RooArgSet( angles[0]._var, angles[1]._var ) )
                    ]
        for ( pad, func, angle, xTitle, yTitle, yScale )\
            in zip(  canvTuple[ 3 : ]
                   , integrals
                   , angles
                   , xLabels[1] if transAngles else xLabels[0]
                   , yLabels[1] if transAngles else yLabels[0]
                   , [ ( 11.0, 14.0 ), ( 11.9, 14.9 ), ( 3.5, 4.5 ) ]
                  ) :
            pad.SetLeftMargin(0.28)
            plot(  pad, angle, None, func
                 #, addPDFs      = [ momFuncAdd ]
                 , xTitle       = xTitle
                 , yTitle       = yTitle
                 , yScale       = yScale
                 , yTitleOffset = 1.0
                 , frameOpts    = dict( Title = angle.GetTitle() )
                 , pdfOpts      = dict( LineColor = kBlue, LineWidth = 4 )
                 #, addPDFsOpts  = [ dict( LineColor = kRed ) ]
                )

    for canvIt, canv in enumerate(canvTuple) : canv.Print( plotsFile + ( '(' if not transAngles and canvIt == 0 else\
                                                                         ')' if transAngles and canvIt == len(canvTuple) - 1 else '' ) )
