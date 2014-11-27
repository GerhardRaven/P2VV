DTTPath = 'DTT_merged.root'
treeName = 'DecayTree'
plotsFilePath = 'DTTPlots.pdf'
dataFilePath = 'DTTDataSets.root'
readData = True
massRange = ( 5200., 5550. )
massPlotRange = ( 5300., 5440. ) #( 5200., 5550. ) #( 5300., 5440. )
binSize = 2.5

#cutKeys = [ 'none', 'trigger', 'vertex', 'KKPT', 'time', 'all' ]
cutKeys = [ 'trigger', 'KKPT', 'vertex', 'time', 'all' ]
cuts = dict(  muPID       = 'min(muplus_PIDmu,muminus_PIDmu) > 0.'
            , muTrack     = 'muplus_TRACK_CHI2NDOF < 4. && muminus_TRACK_CHI2NDOF < 4.'
            , muPT        = 'min(muplus_PT,muminus_PT) > 500.'
            , KPID        = 'Kplus_PIDK > 0. && Kminus_PIDK > 0.'
            , KTrack      = 'Kplus_TRACK_CHI2NDOF < 4. && Kminus_TRACK_CHI2NDOF < 4.'
            , mumuMass    = 'J_psi_1S_MM > 3030. && J_psi_1S_MM < 3150.'
            , KKPT        = 'phi_1020_PT > 1000.'
            , KKVertex    = 'phi_1020_ENDVERTEX_CHI2/phi_1020_ENDVERTEX_NDOF < 20.'
            , KKMass      = 'phi_1020_MM > 990. && phi_1020_MM < 1050.'
            , psiKKVertex = 'B_s0_LOKI_DTF_CHI2NDOF > 0. && B_s0_LOKI_DTF_VCHI2NDOF > 0. && B_s0_LOKI_DTF_VCHI2NDOF < 5.'
            , psiKKMass   = 'B_s0_M > %f && B_s0_M < %f' % massRange
            , IPChiSq     = 'B_s0_IPCHI2_OWNPV < 25.'
            , nextPV      = '(B_s0_MINIPCHI2NEXTBEST > 50. || B_s0_MINIPCHI2NEXTBEST < 0.)'
            , clones      = '(muplus_TRACK_CloneDist < 0. || muplus_TRACK_CloneDist > 5000.)'
                            + ' && (muminus_TRACK_CloneDist < 0. || muminus_TRACK_CloneDist > 5000.)'
                            + ' && (Kplus_TRACK_CloneDist < 0. || Kplus_TRACK_CloneDist > 5000.)'
                            + ' && (Kminus_TRACK_CloneDist < 0. || Kminus_TRACK_CloneDist > 5000.)'
            , decTime     = 'B_s0_LOKI_DTF_CTAU/0.299792458 > 0.3'
            , decTimeMax  = 'B_s0_LOKI_DTF_CTAU/0.299792458 < 14.'
            , decTimeRes  = 'B_s0_LOKI_DTF_CTAUERR/0.299792458 < 0.12'
            , HLT1        = '(B_s0_Hlt1DiMuonHighMassDecision_TOS == 1 || B_s0_Hlt1TrackMuonDecision_TOS == 1 || B_s0_Hlt1TrackAllL0Decision_TOS == 1)'
            , HLT2        = '(B_s0_Hlt2DiMuonJPsiDecision_TOS == 1 || B_s0_Hlt2DiMuonDetachedJPsiDecision_TOS == 1)'
           )
cutCombs = dict(  none    = [ 'psiKKMass' ]
                , trigger = [ 'psiKKMass', 'HLT1', 'HLT2' ]
                , KKPT    = [ 'psiKKMass', 'HLT1', 'HLT2', 'KKPT' ]
                , vertex  = [ 'psiKKMass', 'HLT1', 'HLT2', 'KKPT', 'psiKKVertex' ]
                , time    = [ 'psiKKMass', 'HLT1', 'HLT2', 'KKPT', 'psiKKVertex', 'decTime' ]
                , all     = [ name for name in cuts.keys() ]
               )

from ROOT import TFile
from P2VV.RooFitWrappers import RooObject, RealVar
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()
mass = RealVar( Name = 'mass', Title = 'm(J/psiKK)', MinMax = ( massRange[0], massRange[1] ), Observable = True )
if not readData :
    print 'reading file "%s"' % DTTPath
    DTTFile = TFile.Open(DTTPath)
    DTT = DTTFile.Get(treeName)
    
    print 'create mass datasets'
    from ROOT import RooDataSet, RooArgSet, RooRealVar, gDirectory
    from P2VV.RooFitWrappers import __dref__
    massSet = RooArgSet( __dref__(mass) )
    dataSets = { }
    import sys
    for key, comb in cutCombs.iteritems() :
        print '--> creating dataset for "%s" selection' % key
        DTT.Draw( '>>selList', ' && '.join( cuts[cutName] for cutName in comb ), 'entrylist' )
        selList = gDirectory.Get('selList')
        print '%d selected events out of %d (%.1f%%)'\
              % ( selList.GetN(), DTT.GetEntries(), float( selList.GetN() ) / float( DTT.GetEntries() ) * 100. )
        dataSets[key] = RooDataSet( 'data_%s' % key, 'data_%s' % key, massSet )
        evIt = -1
        evInd = selList.Next()
        while evInd >= 0 :
            evIt += 1
            if evIt % 100000 == 0 :
                print 'event %d' % evIt
                sys.stdout.flush()
            DTT.GetEntry(evInd)
            mass.setVal(DTT.B_s0_M)
            dataSets[key].add(massSet)
            evInd = selList.Next()
        print '%d events in dataset' % dataSets[key].numEntries()
    DTTFile.Close()

    from ROOT import TObject
    dataFile = TFile.Open( dataFilePath, 'RECREATE' )
    for dataSet in dataSets.itervalues() : dataFile.Append(dataSet)
    dataFile.Write( dataFilePath, TObject.kOverwrite )
    dataFile.Close()

else :
    dataFile = TFile.Open(dataFilePath)
    dataSets = dict( [ ( key, dataFile.Get( 'data_' + key ) ) for key in cutKeys ] )
    dataFile.Close()

print 'build mass PDF'
from P2VV.RooFitWrappers import Component
sigMassComps = Component( 'sigMass', [ ], Yield = ( 85.e3, 0., 2.e5 ) )
bkgMassComps = Component( 'bkgMass', [ ], Yield = (  1.e5, 0., 2.e5 ) )
yields = [ comp.getYield() for comp in [ sigMassComps, bkgMassComps ] ]

from P2VV.Parameterizations.MassPDFs import DoubleGauss_Signal_Mass as SignalBMass
sigMass = SignalBMass( mass = mass )
sigMassComps += sigMass.pdf()
from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass
bkgMass = BackgroundBMass( mass = mass )
bkgMassComps += bkgMass.pdf()
from P2VV.RooFitWrappers import buildPdf
massPdf = buildPdf( [ sigMassComps, bkgMassComps ], Observables = [ mass ], Name = 'JpsiKKMass' )

ws['N_sigMass'].setVal(123.e3)
ws['m_sig_mean'].setVal(5368.2)
ws['m_sig_sigma_1'].setVal(5.6)
ws['m_sig_sigma_2'].setVal(14.)
ws['m_sig_frac'].setVal(0.65)
ws['m_bkg_exp'].setVal(-0.0015)
massPdf.getVariables().Print('v')

print 'draw mass plots'
from P2VV.Load import LHCbStyle
from ROOT import gStyle, TCanvas, kFullDotLarge, kBlack, kBlue, kGreen, kRed, kYellow
gStyle.SetColorModelPS(1)
colours = dict( none = kBlack, trigger = kBlack, KKPT = kBlue, vertex = kGreen + 3, time = kYellow + 3, all = kRed )
nBins = int( ( massPlotRange[1] - massPlotRange[0] ) / binSize )
binWidth = ( massPlotRange[1] - massPlotRange[0] ) / float(nBins)
plots = dict( [ ( key, dataSets[key].get().find('mass').frame( massPlotRange[0], massPlotRange[1], nBins ) ) for key in cutKeys ] )

canv = TCanvas('canv')
canv.SetLeftMargin(0.18)
canv.SetRightMargin(0.05)
canv.SetBottomMargin(0.18)
canv.SetTopMargin(0.05)
for keyIt, key in enumerate(cutKeys) :
    print '--> %s' % key
    nEv = dataSets[key].sumEntries()
    for y in yields : y.setMax( 1.5 * nEv )
    yields[1].setVal( nEv - yields[0].getVal() )
    fitResult = massPdf.fitTo( dataSets[key], Save = True, NumCPU = 8, Optimize = 2, Timer = True, Minimizer = 'Minuit2', Strategy = 1
                              , Offset = True )
    fitResult.PrintSpecial( text = True )

    dataSets[key].plotOn( plots[key], MarkerStyle = kFullDotLarge, MarkerSize = 0.6, MarkerColor = colours[key]
                         , LineWidth = 2, LineColor = colours[key] )
    massPdf.plotOn( plots[key], LineWidth = 3, LineColor = colours[key] )

    plots[key].SetMinimum(0.)
    plots[key].GetYaxis().SetNdivisions( 5, 5, 0 )
    plots[key].SetXTitle('m(J/#psi K^{+}K^{-}) [MeV/c^{2}]')
    plots[key].SetYTitle( 'Candidates / (%.1f MeV/c^{2})' % binWidth )
    plots[key].SetLabelOffset( 0.03, 'x' )
    plots[key].SetLabelOffset( 0.01, 'y' )
    plots[key].SetTitleOffset( 1.2,  'x' )
    plots[key].SetTitleOffset( 1.2,  'y' )
    plots[key].Draw( '' + ( '' if keyIt == 0 else 'SAMES' ) )
canv.Print(plotsFilePath)

dataFile.Close()
