###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
doFit = True

parFileIn  = ''
parFileOut = ''

dataSetName = 'JpsiKK_sigSWeight'
dataSetFile = 'data/P2VVDataSets2011Reco12_6KKMassBins_2TagCats.root'

# fit options
fitOpts = dict(  NumCPU    = 8
               , Optimize  = 2
               , Minimizer = 'Minuit2'
               , Offset    = True
#               , Hesse     = False
               , Timer     = True
#               , Verbose   = True
              )
pdfConfig['fitOptions'] = fitOpts

fitRange      = ''
corrSFitErr   = 'sumWeight' # '' / 'sumWeight' / ( 0.887, [ 0.566, 0.863, 0.956, 0.948, 0.855, 0.662 ] ) / 'matrix'
randomParVals = ( ) # ( 1., 12345 )
MinosPars     = [#  'AparPhase'
                 #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
                 #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
                ]

# PDF options
pdfConfig['numTimeResBins']     = 40
pdfConfig['timeResType']        = 'eventNoMean'
pdfConfig['constrainTResScale'] = 'constrain'
pdfConfig['timeEffType']        = 'paper2012'
pdfConfig['constrainDeltaM']    = 'constrain'

pdfConfig['timeEffHistFiles'] = dict(  file      = 'data/Bs_HltPropertimeAcceptance_Data-20120816.root'
                                     , hlt1UB    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                     , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                    )

pdfConfig['anglesEffType'] = 'weights'
pdfConfig['angEffMomsFiles'] = 'data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

pdfConfig['SSTagging']        = True
pdfConfig['condTagging']      = True
pdfConfig['contEstWTag']      = True
pdfConfig['constrainTagging'] = 'constrain'

pdfConfig['paramKKMass']     = 'simultaneous'
pdfConfig['KKMassBinBounds'] = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ]
pdfConfig['CSPValues']       = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]
KKMassPars = pdfConfig['obsDict']['KKMass']
pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                  , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )

pdfConfig['lambdaCPParam'] = 'lambPhi'

from P2VV.Imports import extConstraintValues
extConstraintValues.setVal( 'DM',      ( 17.63, 0.11 ) )
extConstraintValues.setVal( 'P0OS',    (  0.392, 0.008, 0.392 ) )
extConstraintValues.setVal( 'DelP0OS', (  0.0110, 0.0034 ) )
extConstraintValues.setVal( 'P1OS',    (  1.000,  0.023  ) )
extConstraintValues.setVal( 'DelP1OS', (  0.000,  0.001  ) )
extConstraintValues.setVal( 'P0SS',    (  0.350, 0.017, 0.350 ) )
extConstraintValues.setVal( 'DelP0SS', ( -0.019, 0.005   ) )
extConstraintValues.setVal( 'P1SS',    (  1.00,  0.16    ) )
extConstraintValues.setVal( 'DelP1SS', (  0.00,  0.01    ) )

# plot options
numBins = ( 60, 30, 30, 30 )
angleNames = (  ( 'helcosthetaK', 'cos#kern[0.1]{#theta_{K}}'   )
              , ( 'helcosthetaL', 'cos#kern[0.1]{#theta_{#mu}}' )
              , ( 'helphi',       '#varphi_{h} [rad]'           )
             )


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
worksp = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
sigData = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
pdfConfig['signalData'] = sigData
pdfConfig['readFromWS'] = True

# build the PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if not 'Optimize' in fitOpts or fitOpts['Optimize'] < 2 :
    # unset cache-and-track
    for par in pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

# data set with weights corrected for background dilution: for phi_s fit only!
if corrSFitErr == 'sumWeight'\
        or ( type(corrSFitErr) != str and hasattr( corrSFitErr, '__iter__' ) and hasattr( corrSFitErr, '__getitem__' ) ) :
    from P2VV.Utilities.DataHandling import correctSWeights
    fitData = correctSWeights( sigData, 'N_cbkgMass_sw'
                              , 'KKMassCat' if pdfConfig['paramKKMass'] == 'simultaneous' else ''
                              , CorrectionFactors = None if corrSFitErr == 'sumWeight' else corrSFitErr )

else :
    fitData = sigData

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)


###########################################################################################################################################
## fit data ##
##############

# float/fix values of some parameters
for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
    if not pdfConfig['SSTagging'] :
        CEvenOdds.setConstant('avgCEven.*')
        CEvenOdds.setConstant( 'avgCOdd.*', True )
    else :
        for CEvenOdd in CEvenOdds :
            CEvenOdd.setConstant('avgCEven.*')
            CEvenOdd.setConstant( 'avgCOdd.*', True )

pdfBuild['tagCatsOS'].parameter('wTagDelP1OS').setVal(0.)
pdfBuild['tagCatsSS'].parameter('wTagDelP1SS').setVal(0.)
pdfBuild['tagCatsOS'].setConstant('wTagDelP1')
pdfBuild['tagCatsSS'].setConstant('wTagDelP1')

pdfBuild['amplitudes'].setConstant('C_SP')

if randomParVals :
    # give parameters random offsets
    import random
    print 'plotProjectionsSixKKbins: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
          % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
    random.seed( randomParVals[1] if randomParVals[1] else None )
    for par in pdfPars :
        if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

# print parameters
print 120 * '='
print 'plotProjectionsSixKKbins: fit data:'
fitData.Print()
print 'plotProjectionsSixKKbins: observables in PDF:'
pdfObs.Print('v')
print 'plotProjectionsSixKKbins: parameters in PDF:'
pdfPars.Print('v')

if doFit :
    # fit data
    print 120 * '='
    print 'plotProjectionsSixKKbins: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    RooMinPars = [ ]
    if MinosPars :
        print 'plotProjectionsSixKKbins: running Minos for parameters',
        for parName in MinosPars :
            RooMinPars.append( pdfPars.find(parName) )
            print '"%s"' % RooMinPars[-1],
        print

    fitResult = pdf.fitTo( fitData, SumW2Error = True if corrSFitErr == 'matrix' else False, Minos = RooMinPars, Save = True
                          , Range = fitRange, **fitOpts )

    # print parameter values
    from P2VV.Imports import parNames, parValues2011 as parValues
    print 'plotProjectionsSixKKbins: parameters:'
    fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    fitResult.covarianceMatrix().Print()
    fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

if parFileOut :
    # write parameters to file
    pdfConfig.getParametersFromPdf( pdf, fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut )


###########################################################################################################################################
## make plots ##
################

# import plotting tools
from P2VV.Load import LHCbStyle
from P2VV.Utilities.Plotting import plot, CPcomponentsPlotingToolkit
from ROOT import TCanvas, kRed, kGreen, kMagenta, kBlue, kSolid

#Initialaze the CP components ploting toolkit
CpPlotsKit = CPcomponentsPlotingToolkit(pdf,sigData)

#Get dictionary with all the pdfs in the KKmass bins
# {'bin_i', {'total'=... , 'even'=... , 'odd'=... , 'swave'=...}}
pdfsDict = CpPlotsKit.getCPcompPdfKKbins()

#Get some useful stuff ncessesary for looping
KKbins      = CpPlotsKit.getNumKKbins()  #Get nummber of KKmass bins 
binNames    = CpPlotsKit.getKKbinNames() #Get list of KKmass bin names
CPcomps     = CpPlotsKit.getCpCompNames()#Get list of names of the CP components
observables = [ pdfBuild['observables'][name] for name in [ 'time', 'cpsi', 'ctheta', 'phi' ] ]

#Set plot options      
markStyle = 8
markSize  = 0.5
lineWidth = 4
CpPlotsKit.setLineColors( dict(total = kBlue, even=kRed, odd=kGreen+3, swave=kMagenta+3) )
CpPlotsKit.setLineStyles( dict(total = kSolid, even=9   , odd=7       , swave=5         ) )
CpPlotsKit.setLineWidth(lineWidth)


#LHCbLabel
from ROOT import TPaveText, gStyle
lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.14,
                         0.87 - gStyle.GetPadTopMargin(),
                         gStyle.GetPadLeftMargin() + 0.20,
                         0.95 - gStyle.GetPadTopMargin(),
                         "BRNDC")
lhcbName.AddText("LHCb")
lhcbName.SetFillColor(0)
lhcbName.SetTextAlign(12)
lhcbName.SetBorderSize(0)


#Get pdf and AdPdf options in case of bining and No bining the projData 
PDFoptsBin   = CpPlotsKit.getPdfOptsSixKKbins(BinData=True)
PDFoptsNOBin =  CpPlotsKit.getPdfOptsSixKKbins(BinData=False)
addPDFotpsBin = CpPlotsKit.getAddPdfsOptsSixKKbins(BinData=True)
addPDFotpsNOBin = CpPlotsKit.getAddPdfsOptsSixKKbins(BinData=False)

#Plot and save  lifetime and angles in all the KK mass bins
for bin in binNames:
         print '\n\nP2VV - INFO: Plotting decay time and angular distributions of ' + bin
         dataSlice = sigData.reduce('KKMassCat==KKMassCat::' + bin)
         pdfSlice  = pdfsDict[bin]['total']
         binIdx = binNames.index(bin)
         for ( pad, obs, nBins, xTitle, yScaleRel, logY )\
                 in zip(  [ TCanvas(o.GetName()+bin) for o in observables ]
                        , observables
                        , numBins
                        , ( observables[0].GetTitle()+' [ps]', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                        , 2 *(( 1., 1.2 ),) + 2 *(( 1. , 1.2 ),) 
                        , ( True, ) + 3 * ( False, )
                       ) :           
             plot(  pad, obs, dataSlice, pdfSlice, xTitle=xTitle, yScaleRel=yScaleRel, logy=logY
                  , frameOpts   = dict( Bins = nBins, Name = bin + obs.GetName() + 'Histo'      )
                  , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize          )
                  , pdfOpts     = PDFoptsNOBin[bin] if obs==observables[0] else PDFoptsBin[bin]
                  , addPDFs     = [ pdfsDict[bin][c] for c in CPcomps ]
                  , addPDFsOpts = [addPDFotpsNOBin[binIdx][c] for c in CPcomps ] if obs==observables[0]\
                            else  [addPDFotpsBin[binIdx][c] for c in CPcomps ]
                   )
             lhcbName.Draw()
        ## print canvas to file
             fName =  bin + '_' + obs.GetName() + '_sFit.ps'
             pad.Print(fName)
                      
# Save all the plots in a root file as RooPlot objects.
from P2VV.Utilities.Plotting import _P2VVPlotStash as rooplots
from ROOT import TFile
plotsFile = TFile('RooPlots6x4.root','recreate')
for plot in rooplots: plot.Write()
plotsFile.Close()
