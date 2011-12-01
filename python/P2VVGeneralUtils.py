###########################################################################################################################################
## P2VVGeneralUtils: General P2VV utilities                                                                                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

###########################################################################################################################################
## Handling Data                                                                                                                         ##
###########################################################################################################################################
def readData( filePath, dataSetName, NTuple = False, observables = None, tagCat = '', initTag = '' ) :
  """reads data from file (RooDataSet or TTree(s))
  """

  if NTuple :
    from ROOT import RooDataSet, TChain

    # create data set from NTuple file(s)
    print "P2VV - INFO: readData: reading NTuple(s) '%s' from file(s) '%s'" % ( dataSetName, filePath )
    files = TChain(dataSetName)
    files.Add(filePath)
    data = RooDataSet( dataSetName, dataSetName, files, observables )

  else :
    from ROOT import TFile

    # get data set from file
    print "P2VV - INFO: readData: reading RooDataset '%s' from file '%s'" % ( dataSetName, filePath )
    file = TFile.Open( filePath, 'READ' )
    data = file.Get(dataSetName)
    file.Close()

  return data


def writeData( filePath, dataSetName, data, NTuple = False ) :
  """writes data to file (RooDataSet or TTree)
  """

  from ROOT import TFile

  print "P2VV - INFO: writeData: writing RooDataSet '%s' to file '%s'" % ( dataSetName, filePath )

  file = TFile.Open( filePath, 'RECREATE' )
  if NTuple : data.tree().Write(dataSetName)
  else : data.Write(dataSetName)
  file.Close()


###########################################################################################################################################
## Plots                                                                                                                                 ##
###########################################################################################################################################

# plot stash: keep the relevant objects alive by keeping a reference to them
global _P2VVPlotStash
_P2VVPlotStash = []

# plotting function
def plot(  canv, obs, data, pdf, components, frameOpts = [ ], dataOpts = [ ], pdfOpts = [ ], plotResidHist = False
         , logy = False, normalize = True, symmetrize = True, usebar = True ) :
    """makes a P2VV plot

    exAxismple usage:

    canvas = plot( canvas.cd(1), observable, data, pdf
                  , {  'psi'    : [ RooFit.LineColor(RooFit.kGreen), RooFit.LineStyle(kDashed) ]
                     , 'nonpsi' : [ RooFit.LineColor(RooFit.kBlue),  RooFit.LineStyle(kDashed) ]
                    }
                  , frameOpts = [ RooFit.Bins(30) ]
                  , dataOpts  = [ RooFit.MarkerSize(0.4), RooFit.XErrorSize(0) ]
                  , pdfOpts   = [ RooFit.LineWidth(2) ]
                 )
    """
    from ROOT import TLine, TPad

    # create frame for observable
    obsFrame = obs.frame(*frameOpts)  if frameOpts else obs.frame()
    xAxis = obsFrame.GetXaxis()
    _P2VVPlotStash.append(obsFrame)

    # plot data
    data.plotOn( obsFrame, Name = 'data', *dataOpts )

    # plot PDF
    if components :
        for comp, opts in components.iteritems() :
            drawOpts = list(opts) + list(pdfOpts)
            pdf.plotOn(obsFrame, Components = comp, *drawOpts )
    pdf.plotOn( obsFrame, Name = 'pdf', *pdfOpts )

    # draw data after drawing the PDF
    obsFrame.drawAfter( 'pdf', 'data' )

    #TODO: add chisq/nbins
    #chisq = obsFrame.chiSquare( 'pdf', 'data' )
    #nbins = obsFrame.GetNbinsX()

    if plotResidHist :
        # get residuals histogram
        residHist = obsFrame.residHist( 'data', 'pdf', normalize )
        residHist.GetXaxis().SetLimits( xAxis.GetXmin(), xAxis.GetXmax() )
        _P2VVPlotStash.append(residHist)

        # create residuals frame
        residFrame = obsFrame.emptyClone( obsFrame.GetName() + '_resid' )
        xAxis = residFrame.GetXaxis()
        _P2VVPlotStash.append(residFrame)

        # set minimum for observable's frame if there is a log scale for y
        if logy : obsFrame.SetMinimum(0.1)

        # set residual plot options
        #TODO: if normalize : plot residHist as a filled histogram with fillcolor blue...
        #      or, maybe, with the 'bar chart' options: 'bar' or 'b'
        for opt in dataOpts :
            if opt.opcode() == 'MarkerSize'  : residHist.SetMarkerSize(opt.getDouble(0))
            if opt.opcode() == 'MarkerStyle' : residHist.SetMarkerStyle(opt.getInt(0))
            if opt.opcode() == 'MarkerColor' : residHist.SetMarkerColor(opt.getInt(0))
            if opt.opcode() == 'Title'       : residFrame.SetTitle(opt.getString(0))

        # residFrame.addPlotable( residHist, 'p' if not usebar else 'b' )
        # zz.plotOn(f,RooFit.DrawOption('B0'), RooFit.DataError( RooAbsData.None ) )
        #residFrame.SetBarWidth(1.0)
        #residHist.SetDrawOption("B HIST")
        residFrame.addPlotable( residHist, 'p' )  # , 'B HIST' )
        #residFrame.setDrawOptions(residHist.GetName(),'B')

        if symmetrize :
            # symmetrize y-axis residuals histogram
            maxY = max( abs(residHist.getYAxisMin()), abs(residHist.getYAxisMax()) )
            residFrame.SetMaximum(maxY)
            residFrame.SetMinimum(-maxY)

        if normalize :
            if residHist.getYAxisMin() > -5 : residFrame.SetMinimum(-5)
            if residHist.getYAxisMax() <  5 : residFrame.SetMaximum(5)

        # add a line at y=0
        zeroLine = TLine( xAxis.GetXmin(), 0, xAxis.GetXmax(), 0 )
        from ROOT import kRed
        zeroLine.SetLineColor(kRed)
        residFrame.addObject(zeroLine)
        #TODO: improve (remove?) axis labels from residFrame, move up against the initial plot

        # draw observable frame
        canv.cd()
        obsName = obs.GetName() + '_plot1'
        obsPad = TPad( obsName, obsName, 0, 0.2, 1, 1 )
        _P2VVPlotStash.append(obsPad)
        if logy: obsPad.SetLogy(1)
        obsPad.SetNumber(1)
        obsPad.Draw()
        canv.cd(1)
        obsFrame.Draw()

        # draw residuals frame
        canv.cd()
        residName = obs.GetName() + '_resid1'
        residPad = TPad( residName, residName, 0, 0, 1, 0.2 )
        _P2VVPlotStash.append(residPad)
        residPad.SetNumber(2)
        residPad.Draw()
        canv.cd(2)
        residFrame.Draw()

    else :
        # draw observable frame
        canv.cd()
        if logy: canv.SetLogy(1)
        obsFrame.Draw()

    canv.Update()
    return canv

