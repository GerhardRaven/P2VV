from P2VVLoad import P2VVLibrary, ROOTStyle
from RooFitWrappers import RooObject
from ROOT import TCanvas
canv = TCanvas( 'effCanv', 'Angular Efficiency' )
canv.Divide( 3, 2 )

for transAngles in [ False, True ] :
    ws = RooObject(workspace = ( 'trans' if transAngles else 'hel' ) + 'Workspace').ws()

    angleNames = [ 'cos(#psi_{tr})', 'cos(#theta_{tr})', '#phi_{tr}' ] if transAngles\
                 else [ 'cos(#theta_{K})', 'cos(#theta_{l})', '#phi_{hel}' ]
    angEffMomentsFile = 'effMomentsTransBasis' if transAngles else 'effMomentsHelBasis'

    if transAngles :
        from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
        angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
    else :
        from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

    angles = [ angleFuncs.angles[ang] for ang in [ 'cpsi', 'ctheta', 'phi' ] ]

    from P2VVGeneralUtils import RealMomentsBuilder
    from math import sqrt, pi
    indices  = [ (PIndex, YIndex0, YIndex1) for PIndex in range(10) for YIndex0 in range(10) for YIndex1 in range(-YIndex0, YIndex0 + 1) ]
    moments = RealMomentsBuilder()
    moments.appendPYList( angleFuncs.angles, indices )
    moments.read(angEffMomentsFile)
    moments.Print( Scale = 1. / 2. / sqrt(pi),
                   Names = 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0022' if transAngles else 'p2vvab_0000|p2vvab_2000|p2vvab_0020' )
    moments.Print( Scale = 1. / 2. / sqrt(pi), MinSignificance = 2. )

    momFuncTerms = moments.buildPDFTerms( CoefNamePrefix = 'transC_' if transAngles else 'helC_',
            Names = 'p2vvab_0000|p2vvab_2000|p2vvab_0020|p2vvab_0022' if transAngles else 'p2vvab_0000|p2vvab_2000|p2vvab_0020' )
    momFunc = momFuncTerms.buildAddition( 'efficiency' + ( 'Trans' if transAngles else 'Hel' ) )
    #momFunc.Print()

    momFuncTermsAdd = moments.buildPDFTerms( CoefNamePrefix = 'transAddC_' if transAngles else 'helAddC_', MinSignificance = 3. )
    momFuncAdd = momFuncTermsAdd.buildAddition( 'efficiencyAdd' + ( 'Trans' if transAngles else 'Hel' ) )
    #momFuncAdd.Print()

    from P2VVGeneralUtils import plot
    from ROOT import kRed
    for ( pad, angle, angName ) in zip( canv.pads( predicate = lambda pad : pad in ( [ 4, 5, 6 ] if transAngles else [ 1, 2, 3 ] ) )
                                       , angles, angleNames ) :
        plot(  pad, angle, None, momFunc
             , addPDFs = [ momFuncAdd ]
             , xTitle = angName, yTitle = '#varepsilon / <#varepsilon>'
             , yScale = ( 0.85, 1.15 ), yTitleOffset = 1.6
             , frameOpts = dict( Title = angle.GetTitle() )
             , addPDFsOpts = [ dict( LineColor = kRed ) ]
            )

canv.Print('angularEfficiency.ps')
