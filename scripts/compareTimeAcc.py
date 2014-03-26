runPeriod = 2011
constrName = 'runPeriod_p%d_poisson' % runPeriod
binBounds = [ 0.3, 0.46, 0.63, 0.83, 1.05, 1.32, 1.65, 2.07, 2.67, 3.7, 14. ]
HLT1UB = 'runPeriod_p%d_hlt1_excl_biased_notExclB' % runPeriod + '_bin_%03d'
HLT1ExclB = 'runPeriod_p%d_hlt1_excl_biased_exclB' % runPeriod + '_bin_%03d'
HLT2UB = 'runPeriod_p%d_hlt2_unbiased_UB' % runPeriod + '_bin_%03d'
HLT2B = 'runPeriod_p%d_hlt2_biased_B' % runPeriod + '_bin_%03d'
plotsFilePath = 'timeAcc/plots%d.ps' % runPeriod
HLT1AccRange = ( 0.,   0.25 )
HLT2AccRange = ( 0.92, 1.   )

from math import sqrt
import P2VV.RooFitWrappers
from P2VV.Load import P2VVLibrary
from ROOT import TFile
dataFiles = [ TFile.Open('timeAcc/timeAcc_constr_noPrescale.root'), TFile.Open('timeAcc/timeAcc_constr_prescale.root') ]
prescales = [ 1., 0.2 ]
combDataFile = TFile.Open('timeAcc/timeAcc_constr_allData.root')

def setParameters( setting, parSet, preScale, parSetIn ) :
    assert setting in [ 'default', 'HLT1', 'HLT2', 'both' ]
    parSet.setRealValue( HLT1UB % 1, 0.7 )
    parSet.find( HLT1UB % 1 ).setConstant( setting not in [ 'HLT1', 'both' ] )
    if preScale :
        parSet.setRealValue( HLT2UB % 1, preScale )
        parSet.find( HLT2UB % 1 ).setConstant(True)
    else :
        parSet.find( HLT2UB % 1 ).setConstant(False)

    for binIt in range( 1, len(binBounds) ) :
        if setting in [ 'HLT1', 'both' ] :
            parSet.setRealValue( HLT1ExclB % binIt, parSetIn.getRealValue( HLT1ExclB % binIt ) )
            parSet.find( HLT1ExclB % binIt ).setConstant(True)
        else :
            parSet.find( HLT1ExclB % binIt ).setConstant(False)
        if setting in [ 'HLT2', 'both' ] :
            parSet.setRealValue( HLT2B % binIt, parSetIn.getRealValue( HLT2B % binIt ) )
            parSet.find( HLT2B % binIt ).setConstant(True)
        else :
            parSet.find( HLT2B % binIt ).setConstant(False)

from P2VV.Load import LHCbStyle
from array import array
from ROOT import TGraphAsymmErrors, TLatex, TCanvas, kFullDotLarge, kBlue
graphs = [ ]
canvs = [ ]
timeValsArr = array( 'd', range( 1, len(binBounds) ) )
timeErrArr = array( 'd', [ 0. ] * len(timeValsArr) )
label = TLatex()
label.SetTextAlign(32)
label.SetTextSize(0.072)
def createGraph( accType, accValsArr, accErrLoArr, accErrHiArr, minMax = ( None, None ) ) :
    graph = TGraphAsymmErrors( len(timeValsArr), timeValsArr, accValsArr, timeErrArr, timeErrArr, accErrLoArr, accErrHiArr )
    if minMax[0] != None : graph.SetMinimum( minMax[0] )
    if minMax[1] != None : graph.SetMaximum( minMax[1] )
    graph.GetXaxis().SetLimits( timeValsArr[0] - 1., timeValsArr[-1] + 1. )
    graph.GetYaxis().SetNdivisions( 5, 5, 0 )
    graph.SetMarkerStyle(kFullDotLarge)
    graph.SetLineColor(kBlue)
    graph.SetMarkerColor(kBlue)
    graph.SetLineWidth(3)
    graph.SetMarkerSize(1.0)
    graph.GetXaxis().SetTitle('time bin')
    graph.GetYaxis().SetTitle( 'HLT1 acceptance ratio (a.u.)' if accType == 'HLT1' else 'HLT2-biased acceptance' )
    graph.GetXaxis().SetTitleOffset(1.1)
    graph.GetYaxis().SetTitleOffset(1.1)
    graphs.append(graph)

    canv = TCanvas( 'canv%d' % len(canvs) )
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)
    graph.Draw('ALP')
    label.DrawLatex( timeValsArr[-1], graph.GetMinimum() + 0.15 * ( graph.GetMaximum() - graph.GetMinimum() ), 'LHCb unofficial' )
    canvs.append(canv)

from ROOT import RooArgSet
print 'reading data file %s' % combDataFile.GetName()
combConstr = combDataFile.Get(constrName)
combConstrSet = RooArgSet(combConstr)
combPars = combConstr.getVariables()

from ROOT import RooRealVar, RooDataSet
xVar = RooRealVar( 'x', 'x', 0., 0., 1. )
xSet = RooArgSet(xVar)
dummyData = RooDataSet( 'dummyData', 'dummyData', xSet )
dummyData.add(xSet)

from ROOT import RooUniform
combSet = RooArgSet(xVar)
combSet.add(combPars)
combPdf = RooUniform( 'combPdf', 'combPdf', combSet )
print '-' * 120
print 'fitting combined acceptance'
setParameters( 'default', combPars, None, None )
fitResult = combPdf.fitTo( dummyData, Minimizer = 'Minuit2', Save = True, PrintLevel = -1, ExternalConstraints = combConstrSet
                          , Minos = combPars )
fitResult.Print('v')
fitResult.correlationMatrix().Print()
combConstrChiSq = -2. * combConstr.getLogVal()
print 'chi^2 combined constraint = %.2f^2' % ( cmp( combConstrChiSq, 0. ) * sqrt( abs(combConstrChiSq) ) )

HLT1ValsArr  = array( 'd', [ 0. ] * len(timeValsArr) )
HLT2ValsArr  = array( 'd', [ 0. ] * len(timeValsArr) )
HLT1ErrLoArr = array( 'd', [ 0. ] * len(timeValsArr) )
HLT1ErrHiArr = array( 'd', [ 0. ] * len(timeValsArr) )
HLT2ErrLoArr = array( 'd', [ 0. ] * len(timeValsArr) )
HLT2ErrHiArr = array( 'd', [ 0. ] * len(timeValsArr) )
for binIt in range( 1, len(binBounds) ) :
    HLT1ValsArr[ binIt - 1 ]  = combPars.find( HLT1ExclB % binIt ).getVal()
    HLT2ValsArr[ binIt - 1 ]  = combPars.find( HLT2B % binIt ).getVal()
    HLT1ErrLoArr[ binIt - 1 ] = -combPars.find( HLT1ExclB % binIt ).getErrorLo()
    HLT1ErrHiArr[ binIt - 1 ] = +combPars.find( HLT1ExclB % binIt ).getErrorHi()
    HLT2ErrLoArr[ binIt - 1 ] = -combPars.find( HLT2B % binIt ).getErrorLo()
    HLT2ErrHiArr[ binIt - 1 ] = +combPars.find( HLT2B % binIt ).getErrorHi()
createGraph( 'HLT1', HLT1ValsArr, HLT1ErrLoArr, HLT1ErrHiArr, HLT1AccRange )
createGraph( 'HLT2', HLT2ValsArr, HLT2ErrLoArr, HLT2ErrHiArr, HLT2AccRange )

print '-' * 120

for fileIt, dataFile in enumerate(dataFiles) :
    print 'reading data file %s' % dataFile.GetName()
    constr = dataFile.Get(constrName)
    constrSet = RooArgSet(constr)
    pars = constr.getVariables()
    
    parSet = RooArgSet(xVar)
    parSet.add(pars)
    pdf = RooUniform( 'pdf', 'pdf', parSet )
    print '-' * 120
    print 'fitting acceptance %d' % fileIt
    setParameters( 'default', pars, prescales[fileIt], None )
    fitResult = pdf.fitTo( dummyData, Minimizer = 'Minuit2', Save = True, PrintLevel = -1, ExternalConstraints = constrSet, Minos = pars )
    fitResult.Print('v')
    fitResult.correlationMatrix().Print()
    constrChiSq = -2. * constr.getLogVal()
    print 'chi^2 constraint %d = %.2f^2' % ( fileIt, cmp( constrChiSq, 0. ) * sqrt( abs(constrChiSq) ) )

    for binIt in range( 1, len(binBounds) ) :
        HLT1ValsArr[ binIt - 1 ]  = pars.find( HLT1ExclB % binIt ).getVal()
        HLT2ValsArr[ binIt - 1 ]  = pars.find( HLT2B % binIt ).getVal()
        HLT1ErrLoArr[ binIt - 1 ] = -pars.find( HLT1ExclB % binIt ).getErrorLo()
        HLT1ErrHiArr[ binIt - 1 ] = +pars.find( HLT1ExclB % binIt ).getErrorHi()
        HLT2ErrLoArr[ binIt - 1 ] = -pars.find( HLT2B % binIt ).getErrorLo()
        HLT2ErrHiArr[ binIt - 1 ] = +pars.find( HLT2B % binIt ).getErrorHi()
    createGraph( 'HLT1', HLT1ValsArr, HLT1ErrLoArr, HLT1ErrHiArr, HLT1AccRange )
    createGraph( 'HLT2', HLT2ValsArr, HLT2ErrLoArr, HLT2ErrHiArr, HLT2AccRange )

    print '-' * 120
    print 'varying acceptance %d' % fileIt
    for setting in [ 'HLT1', 'HLT2', 'both' ] :
        ndof = len(binBounds) - 1
        if setting == 'both' : ndof *= 2
        setParameters( setting, pars, prescales[fileIt], combPars )
        constrChiSqConstr = -2. * constr.getLogVal()
        setParameters( 'default', pars, prescales[fileIt], None )
        fitResult = pdf.fitTo( dummyData, Minimizer = 'Minuit2', Save = True, PrintLevel = -1, ExternalConstraints = constrSet )
        print 'chi^2 constraint %d = %.2f^2' % ( fileIt, cmp( constrChiSqConstr, 0. ) * sqrt( abs(constrChiSqConstr) ) )
        print 'chi^2 difference / ndof = (%.1f - %.1f) / %d = %.2f^2'\
              % ( constrChiSqConstr, constrChiSq, ndof, sqrt( ( constrChiSqConstr - constrChiSq ) / float(ndof) ) )
    print '-' * 120

canvs[0].Print( plotsFilePath + '[' )
for canvIt in range( 0, len(canvs), 2 ) : canvs[canvIt].Print(plotsFilePath)
for canvIt in range( 1, len(canvs), 2 ) : canvs[canvIt].Print(plotsFilePath)
canvs[-1].Print( plotsFilePath + ']' )
