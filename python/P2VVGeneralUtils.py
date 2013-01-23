###########################################################################################################################################
## P2VVGeneralUtils: General P2VV utilities                                                                                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################


import sys
def numCPU( Max = sys.maxint ) :
    try : # needs >= 2.6
        import multiprocessing
        return min(Max,multiprocessing.cpu_count())
    except :
        import os
        ncpu = os.sysconf('SC_NPROCESSORS_ONLN')
        return min(Max,max( ncpu, 1 ))


###########################################################################################################################################
## Handling Data                                                                                                                         ##
###########################################################################################################################################

def readNTuple( filePath, **kwargs ) :
    treeName = kwargs.pop( 'TreeName', 'DecayTree' )
    MCNTuple = kwargs.pop( 'MCNTuple', False       )

    # options
    if 'Cuts' in kwargs :
        cuts = kwargs.pop('Cuts')
    else :
        cuts = dict(  muPlusTrack  = 'muplus_track_chi2ndof < 4.'
                    , muMinusTrack = 'muminus_track_chi2ndof < 4.'
                    , KPlusTrack   = 'Kplus_track_chi2ndof < 4.'
                    , KMinusTrack  = 'Kminus_track_chi2ndof < 4.'
                    , mumuKKMass   = 'mass > 5200. && mass < 5550.'
                    , mumuMass     = 'mdau1 > 3030. && mdau1 < 3150.'
                    , KKMass       = 'mdau2 > 990. && mdau2 < 1050.'
                    , decayTime    = 'time > 0.3 && time < 14.'
                    , decayTimeRes = 'sigmat < 0.12'
                    , selection    = 'sel == 1'
                    , hlt1         = '(hlt1_unbiased_dec || hlt1_biased)'
                    , hlt2         = 'hlt2_biased'
                   )

        if MCNTuple :
            cuts['bkgcat'] = '(bkgcat == 0 || bkgcat == 50)'

    for cutName in kwargs.pop( 'CutsMin', [ ] ) :
        if cutName in cuts : del cuts[cutName]
    for cutName, cut in kwargs.pop( 'CutsPlus', { } ) :
        cuts[cutName] = cut

    if 'Branches' in kwargs :
        branches = kwargs.pop('Branches')

    else :
        branches = [  'runNumber'
                    , 'eventNumber'
                    , 'muplus_track_chi2ndof'
                    , 'muminus_track_chi2ndof'
                    , 'Kplus_track_chi2ndof'
                    , 'Kminus_track_chi2ndof'
                    , 'mass'
                    , 'mdau1'
                    , 'mdau2'
                    , 'time'
                    , 'sigmat'
                    , 'tagdecision'
                    , 'tagdecision_os'
                    , 'tagdecision_ss'
                    , 'tagomega'
                    , 'tagomega_os'
                    , 'tagomega_ss'
                    , 'sel'
                    , 'hlt1_unbiased_dec'
                    , 'hlt1_biased'
                    , 'hlt1_excl_biased_dec'
                    , 'hlt2_unbiased'
                    , 'hlt2_biased'
                    #, 'muplus_PE'
                    , 'muplus_PX'
                    , 'muplus_PY'
                    , 'muplus_PZ'
                    #, 'muminus_PE'
                    , 'muminus_PX'
                    , 'muminus_PY'
                    , 'muminus_PZ'
                    #, 'Kplus_PE'
                    , 'Kplus_PX'
                    , 'Kplus_PY'
                    , 'Kplus_PZ'
                    #, 'Kminus_PE'
                    , 'Kminus_PX'
                    , 'Kminus_PY'
                    , 'Kminus_PZ'
                    , 'helcosthetaK'
                    , 'helcosthetaL'
                    , 'helphi'
                    , 'bkgcat'
                   ]

        if MCNTuple :
            branches += [  'truetime'
                         #, 'muplus_TRUEP_E'
                         , 'muplus_TRUEP_X'
                         , 'muplus_TRUEP_Y'
                         , 'muplus_TRUEP_Z'
                         #, 'muminus_TRUEP_E'
                         , 'muminus_TRUEP_X'
                         , 'muminus_TRUEP_Y'
                         , 'muminus_TRUEP_Z'
                         #, 'Kplus_TRUEP_E'
                         , 'Kplus_TRUEP_X'
                         , 'Kplus_TRUEP_Y'
                         , 'Kplus_TRUEP_Z'
                         #, 'Kminus_TRUEP_E'
                         , 'Kminus_TRUEP_X'
                         , 'Kminus_TRUEP_Y'
                         , 'Kminus_TRUEP_Z'
                        ]

    from ROOT import gROOT, TFile

    # read n-tuple
    file = TFile.Open(filePath)
    tree = file.Get(treeName)
    gROOT.cd('PyROOT:/')

    tree.SetBranchStatus( '*', False )
    for branch in branches : tree.SetBranchStatus( branch, True )

    cutStr = ' && '.join( cut for cut in cuts.values() )

    print 'P2VV INFO: readNTuple: selection cuts data:\n   %s' % cutStr
    print 'P2VV INFO: readNTuple: numbers of events:'
    if 'selection' in cuts and 'hlt1' in cuts and 'hlt2' in cuts :
        print '   n-tuple: %d\n   triggers: %d\n   pre-selection: %d\n   triggers + pre-selection: %d\n   full selection: %d'\
               % (  tree.GetEntries()
                  , tree.GetEntries( cuts['hlt1'] + ' && ' + cuts['hlt2'] )
                  , tree.GetEntries( cuts['selection'] )
                  , tree.GetEntries( cuts['hlt1'] + ' && ' + cuts['hlt2'] + ' && ' + cuts['selection'] )
                  , tree.GetEntries(cutStr)
                 )
    else :
        print '   n-tuple: %d\n   full selection: %d' % (  tree.GetEntries(), tree.GetEntries(cutStr) )

    tree = tree.CopyTree(cutStr)
    file.Close()
    print 'P2VV INFO: readNTuple: read n-tuple with %d events after selection' % tree.GetEntries()
    return tree

def readData( filePath, dataSetName, NTuple = False, observables = None, **kwargs ) :
    """reads data from file (RooDataSet or TTree(s))
    """
    from ROOT import RooFit
    if observables : noNAN = ' && '.join( '( %s==%s )' % ( obs, obs ) for obs in observables )
    cuts = kwargs.pop( 'cuts', '' )

    if observables :
        print 'P2VV - INFO: readData: reading data for observables [ %s ]' % ', '.join( obs.GetName() for obs in observables )

    if NTuple :
      from ROOT import RooDataSet, TChain
      assert observables != None, 'P2VV - ERROR: readData: set of observables is required for reading an n-tuple'

      # create data set from NTuple file(s)
      print 'P2VV - INFO: readData: reading NTuple(s) "%s" from file(s) "%s"' % ( dataSetName, filePath )
      chain = TChain(dataSetName)
      status = chain.Add( filePath, -1 )
      if status == 0 : raise RuntimeError( 'P2VV - ERROR: could not locate tree "%s" in file "%s"' % ( dataSetName, filePath ) )

      if 'ntupleCuts' in kwargs :
          ntupleCuts = kwargs.pop( 'ntupleCuts', '' )
          print 'P2VV - INFO: readData: applying cuts on n-tuple: %s' % ntupleCuts
          ntuple = chain.CopyTree(ntupleCuts)
      else :
          ntuple = chain

      if cuts : print 'P2VV - INFO: readData: applying cuts on data set: %s' % cuts
      data = RooDataSet( dataSetName, dataSetName
                       , [ obs._var for obs in observables ]
                       , Import = ntuple
                       , Cut = noNAN + ' && ' + cuts if cuts else noNAN )
      ntuple.IsA().Destructor(ntuple)
      if chain : chain.IsA().Destructor(chain)

    else :
      from ROOT import TFile

      # get data set from file
      print 'P2VV - INFO: readData: reading RooDataset "%s" from file "%s"' % ( dataSetName, filePath )
      file = TFile.Open( filePath, 'READ' )
      assert file, 'P2VV - ERROR: readData: file "%s" could not be opened' % filePath

      if observables :
          if cuts : print 'P2VV - INFO: readData: applying cuts: %s' % cuts
          from ROOT import RooDataSet
          data = RooDataSet( dataSetName, dataSetName
                           , [ obs._var for obs in observables ]
                           , Import = file.Get(dataSetName)
                           , Cut = noNAN + ' && ' + cuts if cuts else noNAN 
                           )
      else :
          data = file.Get(dataSetName)

      file.Close()

    print 'P2VV - INFO: read dataset with %s entries' % data.numEntries()

    # import data set into current workspace
    from RooFitWrappers import RooObject
    wsData = RooObject().ws().put( data, **kwargs )
    data.IsA().Destructor(data)
    return wsData


def writeData( filePath, dataSetName, data, NTuple = False ) :
    """writes data to file (RooDataSet or TTree)
    """

    from ROOT import TFile

    print "P2VV - INFO: writeData: writing RooDataSet '%s' to file '%s'" % ( dataSetName, filePath )

    f = TFile.Open( filePath, 'RECREATE' )
    assert f, "P2VV - ERROR: writeData: file '%s' could not be opened" % filePath
    if NTuple : data.tree().Write(dataSetName)
    else : data.Write(dataSetName)
    f.Close()


def correctSWeights( dataSet, bkgWeightName, splitCatName, **kwargs ) :
    """correct sWeights in dataSet for background dilution
    """

    # check if background weight variable exists in data set
    bkgWeight = dataSet.get().find(bkgWeightName)
    assert bkgWeight, 'P2VV - ERROR: correctSWeights: unknown background weight: "%s"' % bkgWeightName

    if splitCatName :
        # get category that splits data sample
        splitCat = dataSet.get().find(splitCatName)
        assert splitCat, 'P2VV - ERROR: correctSWeights: unknown spit category: "%s"' % splitCat

        indexDict = { }
        for iter, catType in enumerate( splitCat ) : indexDict[ iter ] = catType.getVal()

    corrFactors = kwargs.pop( 'CorrectionFactors', [ ] )
    if not corrFactors :
        if splitCatName :
            # initialize sums for the weights and the weights squared per category
            sumWeights   = [ 0. ] * ( splitCat.numTypes() + 1 )
            sumSqWeights = [ 0. ] * ( splitCat.numTypes() + 1 )
            posDict = { }
            for iter, catType in enumerate( splitCat ) : posDict[ catType.getVal() ] = iter

        else :
            # initialize sums for the weights and the weights squared
            sumWeights   = [ 0. ]
            sumSqWeights = [ 0. ]

        # loop over events and get sums of weights and weights squared
        for varSet in dataSet :
            weight = dataSet.weight()
            sumWeights[0]   += dataSet.weight()
            sumSqWeights[0] += dataSet.weight()**2
            if splitCatName :
                sumWeights[ posDict[ varSet.getCatIndex(splitCatName) ] + 1 ]   += dataSet.weight()
                sumSqWeights[ posDict[ varSet.getCatIndex(splitCatName) ] + 1 ] += dataSet.weight()**2

        # get correction factors
        corrFactors = [ sum / sumSq for sum, sumSq in zip( sumWeights, sumSqWeights ) ]

    # add corrected weights to data set
    from ROOT import RooCorrectedSWeight
    if splitCatName :
        from ROOT import std
        corrFactorsVec = std.vector('Double_t')()
        print 'P2VV - INFO: correctSWeights: multiplying sWeights (-ln(L)) to correct for background dilution with factors (overall factor %.4f):'\
              % corrFactors[0]
        for iter, fac in enumerate( corrFactors[ 1 : ] ) :
            corrFactorsVec.push_back(fac)
            print '    %d: %.4f' % ( indexDict[iter], fac )

        weightVar = RooCorrectedSWeight( 'weightVar', 'weight variable', bkgWeight, splitCat, corrFactorsVec, True )

    else :
        print 'P2VV - INFO: correctSWeights: multiplying sWeights (-ln(L)) to correct for background dilution with a factor %.4f'\
              % corrFactors[0]
        weightVar = RooCorrectedSWeight( 'weightVar', 'weight variable', bkgWeight, corrFactors[0], True )

    from ROOT import RooDataSet
    dataSet.addColumn(weightVar)
    dataSet = RooDataSet( dataSet.GetName() + '_corrErrs', dataSet.GetTitle() + ' corrected errors', dataSet.get()
                         , Import = dataSet, WeightVar = ( 'weightVar', True ) )

    # import data set into current workspace
    from RooFitWrappers import RooObject
    return RooObject().ws().put( dataSet, **kwargs )


def addTaggingObservables( dataSet, iTagName, tagCatName, tagDecisionName, estimWTagName, tagCatBins ) :
    """add tagging observables to data set
    """

    # get observables from data set
    obsSet      = dataSet.get(0)
    tagDecision = obsSet.find(tagDecisionName)
    estimWTag   = obsSet.find(estimWTagName)

    # create initial state tag
    from ROOT import RooTagDecisionWrapper
    iTagWrapper = RooTagDecisionWrapper(iTagName, 'Tagging Category', tagDecision)

    # create tagging category
    from ROOT import RooThresholdCategory
    binOneThresh = tagCatBins[1][2]
    tagCatFormula = RooThresholdCategory( tagCatName, 'P2VV tagging category', estimWTag, tagCatBins[0][0], tagCatBins[0][1] )
    for cat in range( 1, len(tagCatBins) ) : tagCatFormula.addThreshold( tagCatBins[cat][2], tagCatBins[cat][0], tagCatBins[cat][1] )

    # create new columns in data set
    dataSet.addColumn(iTagWrapper)
    dataSet.addColumn(tagCatFormula)

    # check tagging columns
    for obsSet in dataSet :
        assert obsSet.getCatIndex(iTagName) == +1 or obsSet.getCatIndex(iTagName) == -1,\
                'P2VV - ERROR: addTaggingObservables: initial state flavour tag has value %+d' % obsSet.getCatIndex(iTagName)
        assert obsSet.getCatIndex(tagDecisionName) == 0 or obsSet.getCatIndex(iTagName) == obsSet.getCatIndex(tagDecisionName),\
                'P2VV - ERROR: addTaggingObservables: initial state flavour tag and tag decision have different values: %+d and %+d'\
                % ( obsSet.getCatIndex(iTagName), obsSet.getCatIndex(tagDecisionName) )
        assert ( obsSet.getCatIndex(tagDecisionName) == 0 and obsSet.getRealValue(estimWTagName) >= binOneThresh )\
                or ( obsSet.getCatIndex(tagDecisionName) != 0 and obsSet.getRealValue(estimWTagName) < binOneThresh ),\
                'P2VV - ERROR: addTaggingObservables: tag decision = %+d, while tagging category = %d'\
                % ( obsSet.getCatIndex(tagDecisionName), obsSet.getCatIndex(tagCatName) )
        assert ( obsSet.getCatIndex(tagDecisionName) == 0 and obsSet.getCatIndex(tagCatName) == 0 )\
                or ( obsSet.getCatIndex(tagDecisionName) != 0 and obsSet.getCatIndex(tagCatName) > 0 ),\
                'P2VV - ERROR: addTaggingObservables: tag decision = %+d, while tagging category = %d'\
                % ( obsSet.getCatIndex(tagDecisionName), obsSet.getCatIndex(tagCatName) )


def addTransversityAngles( dataSet, cpsiName, cthetaTrName, phiTrName, cthetaKName, cthetalName, phiName ) :
    """add transversity angles to data set
    """

    # get observables from data set
    obsSet  = dataSet.get(0)
    cthetaK = obsSet.find(cthetaKName)
    cthetal = obsSet.find(cthetalName)
    phi     = obsSet.find(phiName)

    # create transversity angle functions
    from ROOT import RooTransAngle
    cpsi     = RooTransAngle( cpsiName,     'Cosine of kaon polarization angle',  cthetaK             )
    cthetaTr = RooTransAngle( cthetaTrName, 'Cosine of transversity polar angle', cthetal, phi, False )
    phiTr    = RooTransAngle( phiTrName,    'Transversity azimuthal angle',       cthetal, phi, True  )

    # create new columns in data set
    dataSet.addColumn(cpsi)
    dataSet.addColumn(cthetaTr)
    dataSet.addColumn(phiTr)


###########################################################################################################################################
## Plots                                                                                                                                 ##
###########################################################################################################################################

# plot stash: keep the relevant objects alive by keeping a reference to them
global _P2VVPlotStash
_P2VVPlotStash = []

# plotting function
def plot(  canv, obs, data = None, pdf = None, addPDFs = [ ], components = None, xTitle = '', yTitle = '', xTitleOffset = None
           , yTitleOffset = None, yScale = ( None, None ), yScaleRel = ( None, None ), frameOpts = { }, dataOpts = { }, pdfOpts = { }
           , addPDFsOpts = [ { } ], plotResidHist = False, logy = False, logx = False, normalize = True, symmetrize = True, usebar = True
        ) :
    """makes a P2VV plot

    example usage:

    canvas = plot( canvas.cd(1), observable, data, pdf
                  , {  'psi'    : { 'LineColor' : RooFit.kGreen, 'LineStyle' : RooFit.kDashed }
                     , 'nonpsi' : { 'LineColor' : RooFit.kBlue,  'LineStyle' : RooFit.kDashed }
                    }
                  , xTitle = 'M (MeV/c)'
                  , frameOpts = { 'Title'      : 'B mass', 'Bins        : 30 }
                  , dataOpts  = { 'MarkerSize' : 0.4,      'XErrorSize' : 0  }
                  , pdfOpts   = { 'LineWidth'  : 2                           }
                 )
    """
    from ROOT import TLine, TPad

    # create frame for observable
    obsFrame = obs.frame(**frameOpts)  if frameOpts else obs.frame()
    xAxis = obsFrame.GetXaxis()
    yAxis = obsFrame.GetYaxis()
    _P2VVPlotStash.append(obsFrame)

    # plot data
    if data :
        rooPlot = data.plotOn( obsFrame, Name = 'data', **dataOpts )
        # Set negative bins to 0 if logy is requested
        if logy:
            minimum = 0.
            hist = rooPlot.getHist()
            from ROOT import Double
            x = Double(0.)
            y = Double(0.)
            for i in range(hist.GetN()):
                r = hist.GetPoint(i, x, y)
                if y < minimum:
                    minimum = y
            #hist.SetMinimum(minimum + 0.1)
            obsFrame.SetMinimum( max( ( minimum, 0.1 ) ) )

    # plot PDF
    if pdf :
        # define function that parces the 'Slice(s)' argument and plots the pdf
        def plotPDFWithSlices( pdf, frame, name, **pdfOpts ) :
            if 'Slice' in pdfOpts or 'Slices' in pdfOpts :
                # get 'Slice(s)' argument from plot options
                origSlices = pdfOpts.pop( 'Slices', [ ] )
                if 'Slice' in pdfOpts : origSlices += [ pdfOpts.pop('Slice') ]

                # parse 'Slices' argument
                slicesList = [ [ ] ]
                for slice in origSlices :
                    tempList = [ ]
                    for slices in slicesList : tempList += [ slices + [( slice[0], catType.strip() )] for catType in slice[1].split(',') ]
                    slicesList = tempList

                for num, slices in enumerate(slicesList) :
                    # plot pdf for all slices
                    if num == 0 and len(slicesList) == 1 :
                        opts = dict( Name = name, Slices = slices, **pdfOpts )
                    elif num == 0 :
                        opts = dict( Name = name + '0', Invisible = None, Slices = slices, **pdfOpts )
                    elif num == len(slicesList) - 1 :
                        opts = dict( Name = name, AddTo = (name + '%d' % (num - 1), 1., 1.), Slices = slices, **pdfOpts )
                    else :
                        opts = dict(  Name = name + '%d' % num, AddTo = (name + '%d' % (num - 1), 1., 1.), Invisible = None
                                    , Slices = slices, **pdfOpts
                                   )

                    pdf.plotOn( obsFrame, **opts )

            else :
                pdf.plotOn( obsFrame, Name = name, **pdfOpts )

        if components :
            # plot separate components of the pdf
            for num, comp in enumerate( components.keys() ) :
                drawOpts = components[comp].copy()
                for opt, optVal in pdfOpts.iteritems() :
                    if opt not in drawOpts : drawOpts[opt] = optVal
                plotPDFWithSlices( pdf, obsFrame, 'comp%d' % num, Components = comp, **drawOpts )

        # plot total pdf
        drawOpts = pdfOpts.copy()
        plotPDFWithSlices( pdf, obsFrame, 'pdf', **drawOpts )

        # draw data after drawing the PDF
        if data and 'Asymmetry' not in pdfOpts : obsFrame.drawAfter( 'pdf',  'data' )

    # plot additional PDFs
    if addPDFs :
        for num, addPDF in enumerate(addPDFs) :
            addPDF.plotOn( obsFrame, Name = 'addPDF' + str(num), **(addPDFsOpts[num]) )
            if data and 'Asymmetry' not in addPDFsOpts[num] : obsFrame.drawAfter( 'addPDF' + str(num), 'data' )

    #TODO: add chisq/nbins
    #chisq = obsFrame.chiSquare( 'pdf', 'data' )
    #nbins = obsFrame.GetNbinsX()

    # set y scale
    if yScale[0]    : obsFrame.SetMinimum(yScale[0])
    if yScale[1]    : obsFrame.SetMaximum(yScale[1])
    if yScaleRel[0] : obsFrame.SetMinimum( yScaleRel[0] * obsFrame.GetMinimum() )
    if yScaleRel[1] : obsFrame.SetMaximum( yScaleRel[1] * obsFrame.GetMaximum() )
    if logy and obsFrame.GetMinimum() <= 0 : obsFrame.SetMinimum(0.1)

    # set axis titles
    if xTitle : xAxis.SetTitle(xTitle)
    if yTitle : yAxis.SetTitle(yTitle)

    # set axis title offsets
    if yTitleOffset: yAxis.SetTitleOffset(yTitleOffset)
    if xTitleOffset: xAxis.SetTitleOffset(xTitleOffset)

    # get residuals histogram
    if plotResidHist and data and pdf :
        residHist = obsFrame.residHist( 'data', 'pdf', normalize )
        residHist.GetXaxis().SetLimits( xAxis.GetXmin(), xAxis.GetXmax() )
        _P2VVPlotStash.append(residHist)

        xAxis.SetLabelOffset(0.1)
        yAxis.SetTitleSize(0.08)
        yAxis.SetLabelSize(0.07)
        yAxis.SetTitleOffset(0.75)

        # create residuals frame
        residFrame = obsFrame.emptyClone( obsFrame.GetName() + '_resid' )
        #if 'Title' in frameOpts: residFrame.SetTitle(frameOpts['Title'])
        #residFrame.SetTitle('')
        xAxis = residFrame.GetXaxis()
        xAxis.SetLabelSize(0.15)
        xAxis.SetTitleSize(0.15)
        xAxis.SetLabelOffset(0.02)
        xAxis.SetTitleOffset(1.0)
        yAxis = residFrame.GetYaxis()
        yAxis.SetTitle('')
        yAxis.SetLabelSize(0.13)
        yAxis.SetLabelOffset(0.01)
        _P2VVPlotStash.append(residFrame)

        # set minimum for observable's frame if there is a log scale for y
        #if logy : obsFrame.SetMinimum(0.1)

        # set residual plot options
        #TODO: if normalize : plot residHist as a filled histogram with fillcolor blue...
        #      or, maybe, with the 'bar chart' options: 'bar' or 'b'
        if dataOpts :
            fun = {  'MarkerSize'  : lambda x : residHist.SetMarkerSize(x)
                   , 'MarkerStyle' : lambda x : residHist.SetMarkerStyle(x)
                   , 'MarkerColor' : lambda x : residHist.SetMarkerColor(x)
                   , 'LineWidth'   : lambda x : residHist.SetLineWidth(x)
                   , 'Title'       : lambda x : residFrame.SetTitle(x)
                  }
            for k, v in dataOpts.iteritems() : 
                if k in fun : fun[k](v)

        # residFrame.addPlotable( residHist, 'p' if not usebar else 'b' )
        # zz.plotOn(f,RooFit.DrawOption('B0'), RooFit.DataError( RooAbsData.None ) )
        #residFrame.SetBarWidth(1.0)
        #residHist.SetDrawOption("B HIST")
        residFrame.addPlotable( residHist, 'P' )  # , 'B HIST' )
        #residFrame.setDrawOptions(residHist.GetName(),'B')

        if symmetrize :
            # symmetrize y-axis residuals histogram
            maxY = max( abs(residHist.getYAxisMin()), abs(residHist.getYAxisMax()) )
            residFrame.SetMaximum(maxY)
            residFrame.SetMinimum(-maxY)

        if normalize :
            if residHist.getYAxisMin() > -5.5 : residFrame.SetMinimum(-5.5)
            if residHist.getYAxisMax() < +5.5 : residFrame.SetMaximum(+5.5)

        # add a line at y=0
        zeroLine = TLine( xAxis.GetXmin(), 0, xAxis.GetXmax(), 0 )
        from ROOT import kRed
        zeroLine.SetLineColor(kRed)
        residFrame.addObject(zeroLine)
        #TODO: improve (remove?) axis labels from residFrame, move up against the initial plot

        # draw observable frame
        canv.cd()
        obsName = obs.GetName() + '_plot1'
        obsPad = TPad( obsName, obsName, 0, 0.32, 1, 1 )
        _P2VVPlotStash.append(obsPad)
        if logy: obsPad.SetLogy(1)
        if logx: obsPad.SetLogx(1)
        obsPad.SetNumber(1)
        #obsPad.SetLeftMargin(0.12)
        obsPad.SetTopMargin(0.04)
        obsPad.SetBottomMargin(0.04)
        obsPad.Draw()
        canv.cd(1)
        if 'Title' in frameOpts and not frameOpts['Title'] : obsFrame.SetTitle('')
        obsFrame.Draw()

        # draw residuals frame
        canv.cd()
        residName = obs.GetName() + '_resid1'
        residPad = TPad( residName, residName, 0, 0, 1, 0.32 )
        if logx: residPad.SetLogx(1)
        _P2VVPlotStash.append(residPad)
        residPad.SetNumber(2)
        #residPad.SetLeftMargin(0.12)
        residPad.SetTopMargin(0.)
        residPad.SetBottomMargin(0.4)
        residPad.Draw()
        canv.cd(2)
        if 'Title' in frameOpts and not frameOpts['Title']:
            residFrame.SetTitle("")
        residFrame.Draw()

    else :
        # draw observable frame
        canv.cd()
        if logy: canv.SetLogy(1)
        if logx: canv.SetLogx(1)
        title = frameOpts.get("Title", "")
        if title:
            obsFrame.SetTitle(title)
        obsFrame.Draw()

    canv.Update()
    return canv


# function for plotting a PDF in KK mass bins
def getCondObsPlotsInKKbins(pdf, data, canv):
    from ROOT import TCanvas, gPad, TH1F
    from RooFitWrappers import Category
    import RooFitDecorators

    condObs = pdf.ConditionalObservables()
    nPads = len(condObs)
    nBins =  pdf.indexCat().numTypes()

    #canv = TCanvas('CondObsCanvInKKbins')
    if (nPads & 1 == 0): canv.Divide(nPads/2, 2)
    else: canv.Divide((nPads+1)/2, 2)
    pad = 1
    #Dictionary with all the histogrms, the form is:  { observable, {'bin_i',hist}  }
    Histos = dict( (obs.GetName(), dict( ('bin{0}'.format(bin), TH1F() ) for bin in xrange(nBins)  ) )  for obs in condObs )

    t = data.buildTree()
    entries = t.GetEntries()

    for obs in condObs:
        canv.cd(pad)
        obsName = obs.GetName()
        for KKbin in xrange(nBins):
            binName = 'bin{0}'.format(KKbin)
            sImpose = 'same' if not KKbin == 0 else ''
            if isinstance(obs, Category): #Is CondObs discreate?  (It makes difference in the way Categories are listed in the TTree)
                if (obsName=='iTagOS' or 'iTagSS'):Histos[obsName][binName].SetBins(3,-1,2)
                else: Histos[obsName][binName].SetBins(2,0,2)
                Histos[obsName][binName].SetDefaultSumw2(True)
                Histos[obsName][binName].SetName(obs.GetName() + 'bin{0}'.format(KKbin))
                Histos[obsName][binName].SetLineColor(KKbin + 1)
                Histos[obsName][binName].SetLineWidth(1)
                Histos[obsName][binName].SetAxisRange(0, 350,'Y')
                Histos[obsName][binName].SetXTitle(obsName)

                for event in xrange(entries):
                    t.GetEntry(event)
                    whichKKbin = t.__getattr__('KKMassCat_idx')
                    if (whichKKbin == KKbin ):
                        Histos[obsName][binName].Fill( t.__getattr__(obsName+'_idx'), t.__getattr__('weightVar') )
            else:
                if (obsName=='sigmat'): Histos[obs.GetName()][binName].SetBins(15,0,0.12)
                else:                   Histos[obsName][binName].SetBins(15,0.05,0.55)
                Histos[obsName][binName].SetDefaultSumw2(True)
                Histos[obsName][binName].SetName(obs.GetName() + 'bin{0}'.format(KKbin))
                Histos[obsName][binName].SetLineColor(KKbin + 1)
                Histos[obsName][binName].SetLineWidth(1)
                Histos[obsName][binName].SetAxisRange(.01, 1e3,'Y')
                Histos[obsName][binName].SetXTitle(obsName)

                for event in xrange(entries):
                    t.GetEntry(event)
                    whichKKbin = t.__getattr__('KKMassCat_idx')
                    if (whichKKbin == KKbin ):
                        Histos[obsName][binName].Fill( t.__getattr__(obsName), t.__getattr__('weightVar') )


            int1 = Histos[obsName][binName].GetSumOfWeights()
            if     KKbin==0 : int0 = Histos[obsName]['bin0'].GetSumOfWeights()
            if not KKbin==0 : Histos[obsName][binName].Scale(int0 / int1)

            Histos[obsName][binName].Draw(sImpose)

        if not isinstance(obs, Category): gPad.SetLogy()
        pad += 1
        #assert(False)
    #return  Histos
    if obsName=='iTagSS': assert(False)

    return canv


# class for plotting a PDF in KK mass bins
class CPcomponentsPlotingToolkit():
    def __init__(self, pdf, data):
        #Initializer builds the CP component pdfs
        #Create objects
        self._data = data
        self._tPdf = pdf
        self._CpCompPdfs = dict(total = self._tPdf)
        self._comps = ['even','odd','swave']
        self._condObservables = self._tPdf.ConditionalObservables()
        self._observables = self._tPdf.Observables() - self._condObservables
        self._CPnormFracs = {}
        self._lineColors = dict(total=1,even=4,odd=4,swave=2)
        self._lineStyles = dict(total=1,even=9,odd=3,swave=5)
        self._lineWidth  = 2

        # Check if KK mass binning feature is active
        from RooFitWrappers import SimultaneousPdf
        self._flagKKbin = isinstance( self._tPdf, SimultaneousPdf )
        if (self._flagKKbin):
            self._nKKbins = self._tPdf.indexCat().numTypes()
            self._binNames = [('bin{0}'.format(bin))for bin in xrange(self._nKKbins)]
            self._pdfsSuperDict = {}
            self._sliceNormFracs = {}

        #Start CP components projection of pdf
        # Create pdf paramters dictionary and a set for the original paramters.
        parsDict    = dict((p.GetName(),p) for p in self._tPdf.Parameters()  )
        originalSet = set([ self._tPdf.ws().function("AparMag2"), parsDict['AperpMag2'], parsDict['A0Mag2']])
        if self._flagKKbin: # In the case of KK binning there are several f_S fractions
            f_sSet = set([ parsDict['f_S_bin{0}'.format(b)] for b in xrange(self._nKKbins) ])
        else: f_sSet = set([parsDict['f_S']])
        originalSet.update(f_sSet)

        #Construct CP pdf components
        from RooFitWrappers import ConstVar, Customizer
        for Comp in self._comps:
            #Dictionary with the values  of the parameters to be replaced
            replacementDict = {}
            if (Comp == "even"):
                replacementDict.update(dict(A0Mag2    = parsDict['A0Mag2'].getVal(),
                                            AperpMag2 =            0               ,
                                            AparMag2  = 1-parsDict['A0Mag2'].getVal()-parsDict['AperpMag2'].getVal()
                                            ))
                replacementDict.update(dict( (f_S_i.GetName(), 0 ) for f_S_i in f_sSet)  )

            if (Comp == "odd"):
                replacementDict.update(dict(A0Mag2    =                  0             ,
                                            AperpMag2 = parsDict['AperpMag2'].getVal() ,
                                            AparMag2  =                  0
                                            ))
                replacementDict.update(dict( (f_S_i.GetName(), 0 ) for f_S_i in f_sSet)  )

            if (Comp == "swave"):
                replacementDict.update(dict(A0Mag2    = 0 ,
                                            AperpMag2 = 0 ,
                                            AparMag2  = 0
                                            ))
                replacementDict.update(dict( (f_S_i.GetName(), f_S_i.getVal() ) for f_S_i in f_sSet)  )

            replacementSet = set([ConstVar(Name=k + Comp, Value=v) for k,v in replacementDict.iteritems()])
            CPcompPDF = Customizer( Pdf = self._tPdf, OriginalArgs = originalSet, SubstituteArgs = replacementSet
                                   , ReplaceByName = True, ArgumentSuffix = Comp )

            if (Comp == "even") : self._CpCompPdfs.update(dict(even  = CPcompPDF ))
            if (Comp == "odd")  : self._CpCompPdfs.update(dict(odd   = CPcompPDF ))
            if (Comp == "swave"): self._CpCompPdfs.update(dict(swave = CPcompPDF ))

        if self._flagKKbin:
            #Split the Cp Components further into individual KK mass copmponents
            tPdfs = dict( (bin, self._CpCompPdfs['total'].getPdf(bin)) for bin in self._binNames)
            ePdfs = dict( (bin, self._CpCompPdfs['even' ].getPdf(bin)) for bin in self._binNames)
            oPdfs = dict( (bin, self._CpCompPdfs['odd'  ].getPdf(bin)) for bin in self._binNames)
            sPdfs = dict( (bin, self._CpCompPdfs['swave'].getPdf(bin)) for bin in self._binNames)
        self._pdfsSuperDict.update(dict( (bin, dict(total = tPdfs[bin],
                                                    even  = ePdfs[bin],
                                                    odd   = oPdfs[bin],
                                                    swave = sPdfs[bin] ) ) for bin in self._binNames ))
        #Sneek these lines that set a unit for the phi, it helps RooFit
        #  when creating the RooPlot. for some reason the unit of phi was not set up
        try: list(self._observables)[[a.GetName() for a in list(self._observables)].index('helphi')].setUnit('rad')
        except ValueError: print 'helphi observable not in angles lsit, Failed to set, rad, as a unit. '

    #End of the initialazation

    #Internal methods
    def calculateNormFracs(self, pdfDict):
        from ROOT import RooArgSet
        obs = RooArgSet(o._target_() for o in self._observables )
        totInt = pdfDict['total'].getNorm(obs)
        fEven  = pdfDict['even' ].getNorm(obs) / totInt
        fOdd   = pdfDict['odd'  ].getNorm(obs) / totInt
        fSwave = pdfDict['swave'].getNorm(obs) / totInt
        return dict(even=fEven, odd=fOdd, swave=fSwave)

    def calculateCPnormFracs(self):
        print 'P2VV - INFO: Calculating relative normaliation fractions of the CP components.'
        if not self._flagKKbin:
             print 'P2VV - INFO: Finished calculating relative normaliation fractions of the CP components.'
             return calculateNormFracs(self._CpCompPdfs)
        if self._flagKKbin:
            self._CPnormFracs = dict( (bin ,self.calculateNormFracs(dict(
                                 total = self._CpCompPdfs['total'].getPdf(bin),
                                 even  = self._CpCompPdfs['even'].getPdf(bin),
                                 odd   = self._CpCompPdfs['odd'].getPdf(bin),
                                 swave = self._CpCompPdfs['swave'].getPdf(bin)
                                 )))for bin in self._binNames )
            print 'P2VV - INFO: Finished calculating relative normaliation fractions of the CP components.'
            return self._CPnormFracs

    def calculateKKslicesNormFracs(self):
        table = self._data.table(self._tPdf.indexCat())
        total = float(self._data.sumEntries()) if   self._data.isWeighted() \
                                               else float(self_.data.numEntries())
        self._sliceNormFracs =  dict( (bin,table.get(bin)/total)for bin in self._binNames )
        return self._sliceNormFracs

    def getProJWdata(self,bin,Bins):
        #Helping internal function to aviodavoid dublicating code,
        #  Usefull in the case where you make 6x4 observable plots
        if bin:
            projData = self.binDataSet(Bins)
            from ROOT import RooRealVar
            projVars = []
            for pV in projData.get():
                if isinstance(pV,RooRealVar): projVars.append(pV)
        else :
            projVars = list(self._condObservables)
            if self._flagKKbin: projVars.append(self._tPdf.indexCat())
            projData = self._data.reduce(ArgSet=projVars)

        return dict(data=projData, vars=projVars)

    def binDataSet(self, nBins):
        if self._flagKKbin: projVars = list(self._condObservables) + [self._tPdf.indexCat()]
        else              : projVars = list(self._condObservables)

        from RooFitWrappers import Category
        from ROOT import RooArgSet, RooDataHist
        binnedVarsList = []
        #Bin only the continous observables
        for pV in list(self._condObservables):
            if    isinstance(pV,Category):pass
            else: binnedVarsList.append(pV)
        for pV in binnedVarsList: pV.setBins(nBins)

        binnedVars =  RooArgSet(self._tPdf.indexCat(), *binnedVarsList)
        return RooDataHist('RDH', 'RDH', binnedVars, self._data.reduce(RooArgSet(*projVars)))

    #Interface
    def getCPcompPdf(self):        return self._CpCompPdfs
    def getNumKKbins(self):        return self._nKKbins
    def getCPcompPdfKKbins(self):  return self._pdfsSuperDict
    def getKKbinNames(self):       return self._binNames
    def getCpCompNames(self):      return self._comps
    def getCPnormFracs(self):
        if not self._CPnormFracs: self.calculateCPnormFracs()
        return self._CPnormFracs

    def getKKslicesNormFracs(self):
        if not self._sliceNormFracs: self.calculateKKslicesNormFracs()
        return self._sliceNormFracs

    def getPdfOpts(self, BinData=True,bins=20):
        if   BinData: projDataSet=self.binDataSet(bins)
        else:
            projVars = list(self._condObservables) + [self._tPdf.indexCat()]
            projDataSet = self._data.reduce(ArgSet=projVars)
        return dict( LineWidth = self._lineWidth
                   , LineColor = self._lineColors['total']
                   , ProjWData = (projDataSet, False)
                     )

    def getAddPdfs(self):
        return [self._CpCompPdfs['even' ].getPdf(b)for b in self._binNames] +\
               [self._CpCompPdfs['odd'  ].getPdf(b)for b in self._binNames] +\
               [self._CpCompPdfs['swave'].getPdf(b)for b in self._binNames]

    def getAddPdfsOpts(self, BinData=True,bins=20):
        if not self._CPnormFracs:    self.calculateCPnormFracs()
        if not self._sliceNormFracs: self.calculateKKslicesNormFracs()
        if BinData:
            data     = self.binDataSet(bins)
            projVars = self.getProJWdata(BinData,bins)['vars']
        else:
            data     = self._data
            projVars = list(self._condObservables)
        opts = []
        for comp in self._comps:
            for bin in self._binNames:
                binInd = self._binNames.index(bin)
                addPdfOpt_i = dict( ProjWData     = (data.reduce('KKMassCat==KKMassCat::' + bin),False),
                                    Normalization =  self._CPnormFracs[bin][comp] * self._sliceNormFracs[bin] )
                if not binInd==self._nKKbins-1:addPdfOpt_i.update(dict( Invisible = ()                       ))
                if     binInd==self._nKKbins-1:addPdfOpt_i.update(dict( LineColor = self._lineColors[comp],
                                                                        LineStyle = self._lineStyles[comp],
                                                                        LineWidth = self._lineWidth          ))
                if not binInd==0: #odd   Components First index = len(self._binNames)
                                  #swave Components First index = 2 *( len(self._binNames)
                                  if comp=='even' : argAddTo = ( 'addPDF{0}'.format(binInd-1),1.,1.)
                                  if comp=='odd'  : argAddTo = ( 'addPDF{0}'.format(  len(self._binNames) + binInd-1 ,1.,1.) )
                                  if comp=='swave': argAddTo = ( 'addPDF{0}'.format(2*len(self._binNames) + binInd-1 ,1.,1.) )
                                  addPdfOpt_i.update(dict(AddTo = argAddTo))
                opts.append(addPdfOpt_i)
        return opts

    def getPdfOptsSixKKbins(self, BinData=True, bins=20):
        projecting = self.getProJWdata(BinData,bins)
        projData   = projecting['data']
        projVars   = projecting['vars']
        if not BinData:  projVars.remove( self._tPdf.indexCat() )
        opts = {}
        KKCat = 'KKMassCat==KKMassCat::'
        for b in self._binNames:
            opts.update( {b : dict(  ProjWData = (projData.reduce(KKCat+b).reduce(ArgSet=projVars), False)
                                   , LineWidth = self._lineWidth
                                   , LineStyle = self._lineStyles['total']
                                   , LineColor = self._lineColors['total'])
                          }  )
        return opts

    def getAddPdfsOptsSixKKbins(self,BinData=True,bins=20):
        if not self._CPnormFracs:    self.calculateCPnormFracs()
        projecting = self.getProJWdata(BinData,bins)
        projData   = projecting['data']
        projVars   = projecting['vars']
        if not BinData:  projVars.remove( self._tPdf.indexCat() )
        opts = []
        for bin in self._binNames:
            ith_binOpts = { }
            for comp in self._comps:
                opt = dict(  ProjWData     = (projData.reduce('KKMassCat==KKMassCat::'+bin).reduce(ArgSet=projVars),False)
                           , LineColor     =  self._lineColors[comp]
                           , LineStyle     =  self._lineStyles[comp]
                           , LineWidth     =  self._lineWidth
                           , Normalization =  self._CPnormFracs[bin][comp]
                             )
                ith_binOpts.update( {comp:opt} )
            opts.append( ith_binOpts  )
        return opts

    def setLineColors(self,colors): self._lineColors = colors
    def setLineStyles(self,styles): self._lineStyles = styles
    def setLineWidth(self, width ): self._lineWidth  = width


# function for plotting the S-wave phases versus the (binned) KK mass
def plotSWavePhases( **kwargs ) :
    yAxisRange  = kwargs.pop( 'DeltaSAxisRange', ( None, None )                         )
    KKMassLabel = kwargs.pop( 'KKMassLabel',     'm(K^{+}K^{-}) [MeV/c^{2}]'            )
    deltaSLabel = kwargs.pop( 'DeltaSLabel',     '#delta_{S} - #delta_{#perp}    [rad]' )
    plotTitle   = kwargs.pop( 'PlotTitle',       ''                                     )
    LHCbText1   = kwargs.pop( 'LHCbTextLine1',   'LHCb'                                 )
    LHCbText2   = kwargs.pop( 'LHCbTextLine2',   ''                                     )
    drawLegend  = kwargs.pop( 'DrawLegend',      False                                  )

    if any( key not in kwargs for key in [ 'DeltaSValues', 'DeltaSLowErrors', 'DeltaSHighErrors' ] ) :
        raise KeyError, 'P2VV - ERROR: plotSWavePhases: "DeltaSValues", "DeltaSLowErrors" and "DeltaSHighErrors" arguments are required'
    massBins       = kwargs.pop( 'MassBins', [ 988., 1008., 1032., 1050. ] )
    theoryVals     = kwargs.pop( 'TheoryValues', None )
    deltaSVals     = kwargs.pop( 'DeltaSValues'       )
    deltaSLowErrs  = kwargs.pop( 'DeltaSLowErrors'    )
    deltaSHighErrs = kwargs.pop( 'DeltaSHighErrors'   )

    if kwargs :
        raise KeyError, 'P2VV - ERROR: plotSWavePhases: unexpected keyword arguments: %s' % kwargs

    from array import array
    KKMass         = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMassErr      = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] )                for it in range( len(massBins) - 1 ) ] )

    KKMass1        = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + 0.35 + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMass1LowErr  = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + 0.35                for it in range( len(massBins) - 1 ) ] )
    KKMass1HighErr = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - 0.35                for it in range( len(massBins) - 1 ) ] )

    KKMass2        = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - 0.35 + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMass2LowErr  = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - 0.35                for it in range( len(massBins) - 1 ) ] )
    KKMass2HighErr = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + 0.35                for it in range( len(massBins) - 1 ) ] )

    if theoryVals :
        from ROOT import TGraphErrors
        theory    = array( 'd', theoryVals           )
        theoryErr = array( 'd', [ 0. ] * len(KKMass) )
        theoryGraph = TGraphErrors( len(KKMass), KKMass, theory, KKMassErr, theoryErr )
    else :
        theoryGraph = None

    from ROOT import TGraphAsymmErrors
    deltaS1        = array( 'd', deltaSVals      )
    deltaS1LowErr  = array( 'd', deltaSLowErrs   )
    deltaS1HighErr = array( 'd', deltaSHighErrs )
    deltaS1Graph = TGraphAsymmErrors( len(KKMass1), KKMass1, deltaS1\
                                     , KKMass1LowErr, KKMass1HighErr, deltaS1LowErr, deltaS1HighErr )

    from math import pi
    deltaS2        = array( 'd', [ pi - val for val in deltaSVals ] )
    deltaS2LowErr  = array( 'd', deltaSHighErrs )
    deltaS2HighErr = array( 'd', deltaSLowErrs  )
    deltaS2Graph = TGraphAsymmErrors( len(KKMass2), KKMass2, deltaS2\
                                     , KKMass2LowErr, KKMass2HighErr, deltaS2LowErr, deltaS2HighErr )

    delSMin0 = min( delS for delS in theoryVals ) if theoryVals else +1.e32
    delSMax0 = max( delS for delS in theoryVals ) if theoryVals else -1.e32
    delSMin1 = min( delS - delSErr for delS, delSErr in zip( deltaSVals, deltaSLowErrs  ) )
    delSMax1 = max( delS + delSErr for delS, delSErr in zip( deltaSVals, deltaSHighErrs ) )
    delSMin  = min( [ delSMin0, delSMin1, pi - delSMax1 ] )
    delSMax  = max( [ delSMax0, delSMax1, pi - delSMin1 ] )

    from ROOT import gStyle, kSolid
    gStyle.SetLineStyleString( 11, ' 30 15' )
    if theoryGraph : theoryGraph.SetLineStyle(11)
    deltaS1Graph.SetLineStyle(kSolid)
    deltaS2Graph.SetLineStyle(kSolid)

    from ROOT import kBlack, kRed, kBlue, kWhite
    if theoryGraph : theoryGraph.SetLineColor(kBlack)
    deltaS1Graph.SetLineColor(kRed)
    deltaS2Graph.SetLineColor(kBlue)

    if theoryGraph : theoryGraph.SetMarkerColor(kBlack)
    deltaS1Graph.SetMarkerColor(kRed)
    deltaS2Graph.SetMarkerColor(kBlue)

    if theoryGraph : theoryGraph.SetLineWidth(3)
    deltaS1Graph.SetLineWidth(4)
    deltaS2Graph.SetLineWidth(4)

    from ROOT import kFullCircle, kFullSquare
    if theoryGraph : theoryGraph.SetMarkerStyle(kFullCircle)
    deltaS1Graph.SetMarkerStyle(kFullCircle)
    deltaS2Graph.SetMarkerStyle(kFullSquare)
    if theoryGraph : theoryGraph.SetMarkerSize(0.7)
    deltaS1Graph.SetMarkerSize(1.6)
    deltaS2Graph.SetMarkerSize(1.5)

    if theoryGraph : theoryGraph.SetMinimum( yAxisRange[0] if yAxisRange[0] else delSMin - 0.10 * ( delSMax - delSMin ) )
    if theoryGraph : theoryGraph.SetMaximum( yAxisRange[1] if yAxisRange[1] else delSMax + 0.15 * ( delSMax - delSMin ) )
    deltaS1Graph.SetMinimum( yAxisRange[0] if yAxisRange[0] else delSMin - 0.10 * ( delSMax - delSMin ) )
    deltaS1Graph.SetMaximum( yAxisRange[1] if yAxisRange[1] else delSMax + 0.15 * ( delSMax - delSMin ) )

    if theoryGraph : theoryGraph.GetXaxis().SetTitle(KKMassLabel)
    if theoryGraph : theoryGraph.GetYaxis().SetTitle(deltaSLabel)
    deltaS1Graph.GetXaxis().SetTitle(KKMassLabel)
    deltaS1Graph.GetYaxis().SetTitle(deltaSLabel)

    if theoryGraph : theoryGraph.GetXaxis().SetTitleOffset(1.0)
    if theoryGraph : theoryGraph.GetYaxis().SetTitleOffset(0.7)
    deltaS1Graph.GetXaxis().SetTitleOffset(1.0)
    deltaS1Graph.GetYaxis().SetTitleOffset(0.7)

    if theoryGraph : theoryGraph.SetTitle(plotTitle)
    deltaS1Graph.SetTitle(plotTitle)

    _P2VVPlotStash.append(deltaS1Graph)
    _P2VVPlotStash.append(deltaS2Graph)
    if theoryGraph : _P2VVPlotStash.append(theoryGraph)

    if drawLegend :
        from ROOT import gStyle, TLegend
        leg = TLegend( 0.67, 0.46, 0.93, 0.66 )
        leg.SetTextFont( gStyle.GetTextFont() )
        leg.SetMargin(0.45)
        leg.AddEntry( deltaS1Graph, '#Delta#Gamma_{s} > 0', 'LPE' )
        leg.AddEntry( deltaS2Graph, '#Delta#Gamma_{s} < 0', 'LPE' )
        leg.SetBorderSize(1)
        leg.SetFillStyle(0)
        _P2VVPlotStash.append(leg)
    else :
        leg = None

    if LHCbText1 or LHCbText2 :
        from ROOT import TPaveText
        LHCbText = TPaveText( 0.15, 0.77 if LHCbText1 and LHCbText2 else 0.84, 0.51, 0.93, 'NDC' )
        if LHCbText1 : LHCbText.AddText(LHCbText1)
        if LHCbText2 : LHCbText.AddText(LHCbText2)
        LHCbText.SetShadowColor(0)
        LHCbText.SetFillStyle(0)
        LHCbText.SetBorderSize(0)
        LHCbText.SetTextAlign(12)
        _P2VVPlotStash.append(LHCbText)
    else :
        LHCbText = None

    from ROOT import TCanvas
    SWavePhaseCanv = TCanvas( 'SWavePhaseCanv', 'S-Wave Phases' )
    SWavePhaseCanv.SetLeftMargin(0.12)
    SWavePhaseCanv.SetRightMargin(0.04)
    SWavePhaseCanv.SetTopMargin(0.04)
    SWavePhaseCanv.SetBottomMargin(0.17)
    deltaS1Graph.Draw('AP')
    deltaS2Graph.Draw('P sames')
    if theoryGraph : theoryGraph.Draw('P sames')
    if leg :      leg.Draw()
    if LHCbText : LHCbText.Draw()

    return SWavePhaseCanv

def splot( pdf, sdata ) :
    # switch off all yields, except current one
    from contextlib import contextmanager
    @contextmanager
    def __select_component( i, yields ):
        orig = dict( (j,j.getVal()) for j in yields )
        [ j.setVal(0) for j in orig.iterkeys() if j!=i ]
        try     : yield
        finally : [ j.setVal(v) for (j,v) in orig.iteritems() ]
    from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
    canvas = TCanvas(pdf.GetName() + '_splot')
    obs = [ o for o in pdf.Observables() if hasattr(o,'frame') and o not in sdata.usedObservables() ]
    for (p,o) in zip( canvas.pads(len(obs)), obs ) :
        # select yields
        _yields = [ y for y in pdf.Parameters() if y.getAttribute('Yield') ]
        # loop over components
        for (pp,i) in zip( p.pads(1,len(_yields)), _yields ) :
            # switch off all yields, except current one
            with __select_component( i, _yields ) :
                # plot both weighed data and PDF
                # TODO: add the same color coding as above...
                c_name = i.GetName()[2:]
                c_opts = { 'signal'             : dict( LineColor = kGreen )
                         , 'psi_background'     : dict( LineColor = kRed )
                         , 'cmb_background'     : dict( LineColor = kBlue )
                         }
                from P2VVGeneralUtils import plot
                plot( pp, o, sdata.data( c_name ), pdf, pdfOpts = c_opts[c_name] if c_name in c_opts else {})
    return canvas

def Oldsplot( canv, sdata, pdf, frameOpts = dict(), dataOpts = dict() , pdfOpts = dict() ) :
    obs = [ o for o in pdf.Observables() if hasattr(o,'frame') and o not in sdata.usedObservables() ]
    for (p,o) in zip( canv.pads(len(obs)), obs ) :
        # snapshot yeilds
        _yields = dict( (p,p.getVal()) for p in pdf.Parameters() if p.getAttribute('Yield')  )
        # loop over components
        for (pp,i) in zip( p.pads(1,len(_yields)), _yields.iterkeys() ) :
            # switch off all yields, except current one
            for j in filter( lambda x: x!=i, _yields.iterkeys() ) : j.setVal(0)
            # and plot both weighed data and PDF
            from P2VVGeneralUtils import plot
            c_name = i.GetName()[2:]
            if 'Title' not in frameOpts : frameOpts['Title'] =  '%s : %s' % ( c_name , o.GetTitle() )
            plot( pp, o, sdata.data( c_name ), pdf
                , pdfOpts = pdfOpts[c_name] if c_name in pdfOpts else {}
                , frameOpts = frameOpts
                , dataOpts = dataOpts
                )
            # and put back the original value!
            for (i,v) in _yields.iteritems() : i.setVal(v)

#To help sort the pads in the plot function
class Sorter(object):
    def __init__(self, d):
        self.__d = d
    def __call__(self, o):
        if o in self.__d:  return self.__d[o]
        else:              return len(self.__d) + 1


###########################################################################################################################################
## (Efficiency) Moments                                                                                                                  ##
###########################################################################################################################################

class RealMomentsBuilder ( dict ) :
    # TODO:  implement reduce: clone self, selecting a subset of available moments...
    # TODO:                    support as kw: MinSignificance, Names
    # def reduce( self, **kwargs ) :
    #
    #

    def __init__( self, **kwargs )   :
        self._basisFuncNames   = [ ]
        self._basisFuncIndices = { }
        self._basisFuncs       = { }
        self._coefficients     = { }
        if 'Moments' in kwargs :
            for mom in kwargs.pop('Moments') : self.append( Moment = mom )
        if 'Moment' in kwargs :
            self.append( Moment = kwargs.pop(Moment) )
        if kwargs :
            raise RuntimeError( 'unknown keyword arguments %s' % kwargs.keys() )

    def basisFuncs(self)       : return self._basisFuncs.copy()
    def basisFuncNames(self)   : return self._basisFuncNames[ : ]
    def basisFuncIndices(self) : return self._basisFuncIndices.copy()
    def coefficients(self)     : return self._coefficients.copy()

    def _iterFuncAndCoef( self, **kwargs ) :
        # get minimum significance for coefficient
        minSignif = kwargs.pop( 'MinSignificance', float('-inf') )

        # get name requirement for coefficient
        names = kwargs.pop('Names',None)
        import re
        nameExpr = re.compile(names) if names else None

        # check if there are no arguments left
        assert not kwargs, 'extraneous kw args: %s' % kwargs

        # yield name, coefficient and function for selected moments
        for funcName in self._basisFuncNames :
            if self._coefficients[funcName][2] < minSignif : continue
            if nameExpr and not nameExpr.match(funcName)   : continue
            yield ( funcName, self._coefficients[funcName], self._basisFuncs[funcName] )

    def appendPYList( self, Angles, IndicesList, PDF = None, IntSet = None, NormSet = None ) :
        # build moments from list of indices
        if not PDF and not IntSet and not NormSet :
            # build moment
            for inds in IndicesList :
                self.append(  Angles = Angles, PIndex = inds[0], YIndex0 = inds[1], YIndex1 = inds[2]
                                  , Norm = float( 2 * inds[0] + 1 ) / 2. )
        elif PDF and IntSet != None and NormSet != None :
            # TODO: default for IntSet should be an empty set and a set of angles for NormSet,
            #       but we need a dataset to determine this ;-(
            #       maybe set to 'None' , and then defer to compute???
            #       Problem is, at that point the moments have already been build...
            # build efficiency moment
            for inds in IndicesList :
                self.append(  Angles = Angles, PIndex = inds[0], YIndex0 = inds[1], YIndex1 = inds[2]
                                  , Norm = float( 2 * inds[0] + 1 ) / 2., PDF = PDF, IntSet = IntSet, NormSet = NormSet )
        else :
            print 'P2VV - ERROR: RealMomentsBuilder.appendList: a PDF, an integration set and a normalisation set are required for efficiency moments'

    def append( self, **kwargs ) :
        momIndices = None
        if 'Moment' in kwargs :
            # get moment directly from arguments
            func = None
            moment = kwargs.pop('Moment')

        elif 'Function' in kwargs or all( arg in kwargs for arg in ( 'Angles', 'PIndex', 'YIndex0', 'YIndex1' ) ):
            # build moment with function from arguments
            if 'Function' in kwargs :
                # get function from arguments
                func = kwargs.pop('Function')
            else :
                # build basis function
                from RooFitWrappers import P2VVAngleBasis
                momIndices = ( kwargs.pop('PIndex'), kwargs.pop('YIndex0'), kwargs.pop('YIndex1') )
                func = P2VVAngleBasis( kwargs.pop('Angles'), ( momIndices[0], 0, momIndices[1], momIndices[2] ) , 1. )

            if not 'PDF' in kwargs and not 'IntSet' in kwargs and not 'NormSet' in kwargs :
                # build moment
                from RooFitWrappers import RealMoment
                moment = RealMoment( func, kwargs.pop('Norm',1.) )
            elif 'PDF' in kwargs and 'IntSet' in kwargs and 'NormSet' in kwargs :
                # build efficiency moment
                from RooFitWrappers import RealEffMoment
                moment = RealEffMoment( func, kwargs.pop('Norm',1.), kwargs.pop('PDF'), kwargs.pop('IntSet'), kwargs.pop('NormSet') )
            else :
                print 'P2VV - ERROR: RealMomentsBuilder.append: a PDF, an integration set and a normalisation set are required for an efficiency moment'
                moment = None

        else :
            print 'P2VV - ERROR: RealMomentsBuilder.append: did not find required arguments'
            moment = None

        # check for unexpected arguments
        if kwargs :
            print 'P2VV - ERROR: RealMomentsBuilder.append: unknown arguments:', kwargs
            moment = None

        if moment :
            # append moment
            momName = moment.GetName()
            self._basisFuncNames.append(momName)
            if momIndices : self._basisFuncIndices[momName] = momIndices
            else :          self._basisFuncIndices[momName] = None
            if func : assert func == moment.basisFunc()
            self._basisFuncs[momName] = moment.basisFunc()
            self[momName] = moment

    def compute( self, data ) :
        """computes moments of data set (wrapper for C++ computeRooRealMoments)

        Looping over data in python is quite a bit slower than in C++. Hence, we
        adapt the arguments and then defer to the C++ computeRooRealMoments.
        """
        from P2VVLoad import P2VVLibrary
        from ROOT import std, computeRooRealMoments
        momVec = std.vector('RooAbsRealMoment*')()
        for func in self._basisFuncNames : momVec.push_back( self[func]._var )
        computeRooRealMoments( data, momVec )

        for func in self._basisFuncNames :
            self._coefficients[func] = ( self[func].coefficient(), self[func].stdDev(), self[func].significance() )

    def Print( self, **kwargs ) :
        # get maximum length of basis function name
        maxLenName = 15
        for func in self._basisFuncNames : maxLenName = min( max( len(func), maxLenName ), 60 )

        # get name requirements
        names = kwargs.pop('Names',None)
        import re
        nameExpr = re.compile(names) if names else None

        # get minimum significance
        minSignif = kwargs.pop('MinSignificance', float('-inf') )

        # get scale factors
        scale  = kwargs.pop( 'Scale',  None )
        scales = kwargs.pop( 'Scales', ( scale, scale, 1. ) if scale != None else None )

        # print header
        print 'P2VV - INFO: RealMomentsBuilder.printMoments:'
        print '  name requirement: \'' + ( names if names else '' ) + '\''
        print '  minimum significance = %.1f' % minSignif
        print '  scales = ' + ( str(scales) if scales else '(1., 1., 1.)' )
        print '  ' + '-' * (45 + maxLenName)
        print ( '  {0:<%d}   {1:<12}   {2:<12}   {3:<12}' % maxLenName )\
                .format( 'basis function', 'coefficient', 'std. dev.', 'significance' )
        print '  ' + '-' * (45 + maxLenName)

        # print moments
        for func in self._basisFuncNames :
            if ( nameExpr and not nameExpr.match(func) ) or ( func in self._coefficients and self._coefficients[func][2] < minSignif ) :
                continue

            print ( '  {0:<%d}' % maxLenName ).format(func if len(func) <= maxLenName else '...' + func[3 - maxLenName : ]),
            if func in self._coefficients :
                coef = self._coefficients[func]
                if scales :
                    print '  {0:<+12.5g}   {1:<12.5g}   {2:<12.5g}'.format( coef[0] * scales[0], coef[1] * scales[1], coef[2] * scales[2] )
                else :
                    print '  {0:<+12.5g}   {1:<12.5g}   {2:<12.5g}'.format( coef[0],             coef[1],             coef[2]             )
            else : print

        print '  ' + '-' * (45 + maxLenName)

    def write( self, filePath = 'moments', **kwargs ) :
        # get file path and name
        filePath = filePath.strip()
        fileName = filePath.split('/')[-1]

        # open file
        try :
            momFile = open( filePath, 'w' )
        except :
            raise IOError( 'P2VV - ERROR: RealMomentsBuilder.writeMoments: unable to open file \"%s\"' % filePath )

        # get maximum length of basis function name
        maxLenName = 13
        for func in self._basisFuncNames : maxLenName = max( len(func), maxLenName )

        # get minimum significance
        minSignif = kwargs.pop('MinSignificance',float('-inf'))

        # get scale factors
        scale  = kwargs.pop( 'Scale',  None )
        scales = kwargs.pop( 'Scales', ( scale, scale, 1. ) if scale != None else None )

        # get name requirements
        import re
        names = kwargs.pop('Names', None)
        nameExpr = re.compile(names) if names else None

        # write moments to content string
        cont = '# %s: angular moments\n' % fileName\
             + '# name requirement: \'{0}\'\n'.format( names if names else '' )\
             + '# minimum significance = {0:.1f}\n'.format(minSignif)\
             + '# scales = {0}\n'.format( str(scales) if scales else '(1., 1., 1.)' )\
             + '#\n'\
             + '# ' + '-' * (49 + maxLenName) + '\n'\
             + ( '# {0:<%s}   {1:<14}   {2:<13}   {3:<13}\n' % maxLenName )\
                   .format( 'basis function', 'coefficient', 'std. dev.', 'significance' )\
             + '# ' + '-' * (49 + maxLenName) + '\n'

        numMoments = 0
        for func in self._basisFuncNames :
            if ( nameExpr and not nameExpr.match(func) ) or ( func in self._coefficients and self._coefficients[func][2] < minSignif ) :
                continue

            cont += ( '  {0:<%s}' % maxLenName ).format(func)
            if func in self._coefficients :
                coef = self._coefficients[func]
                if scales :
                    cont += '   {0:<+14.8g}   {1:<13.8g}   {2:<13.8g}\n'.format( coef[0]*scales[0], coef[1]*scales[1], coef[2]*scales[2] )
                else :
                    cont += '   {0:<+14.8g}   {1:<13.8g}   {2:<13.8g}\n'.format( coef[0],           coef[1],           coef[2]           )
                numMoments += 1

            else :
                cont += '\n'

        cont += '# ' + '-' * (49 + maxLenName) + '\n'

        # write content to file
        momFile.write(cont)
        momFile.close()

        print 'P2VV - INFO: RealMomentsBuilder.writeMoments: %d efficiency moment%s written to file \"%s\"'\
                % ( numMoments, '' if numMoments == 1 else 's', filePath )

    def read( self, filePath = 'moments', **kwargs ) :
        self._coefficients = { }

        # get file path
        filePath = filePath.strip()

        # open file
        try :
          momFile = open(filePath, 'r')
        except :
          raise IOError( 'P2VV - ERROR: RealMomentsBuilder.readMoments: unable to open file \"%s\"' % filePath )

        # get minimum significance
        minSignif = kwargs.pop('MinSignificance',float('-inf'))

        # get name requirements
        import re
        nameExpr = re.compile(kwargs.pop('Names')) if 'Names' in kwargs else None

        # loop over lines and read moments
        numMoments = 0
        while True :
            # read next line
            line = momFile.readline()
            if not line : break

            # check for empty or comment lines
            line = line.strip()
            if not line or line[0] == '#' : continue

            # check moment format
            line = line.split()
            if len(line) < 4 or line[0] not in self._basisFuncNames : continue
            try :
              coef   = float(line[1])
              stdDev = float(line[2])
              signif = float(line[3])
            except :
              continue

            # check significance and name
            if ( nameExpr and not nameExpr.match(line[0]) ) or signif < minSignif : continue

            # get moment
            self._coefficients[line[0]] = ( coef, stdDev, signif )
            numMoments += 1

        momFile.close()

        print 'P2VV - INFO: RealMomentsBuilder.readMoments: %d efficiency moment%s read from file \"%s\"'\
                % ( numMoments, '' if numMoments == 1 else 's', filePath )

    def buildPDFTerms( self, **kwargs ) :
        # TODO: decide whether coefficients are ConstVar or RealVar?? (add keyword for that! -- what MinMax to give if RealVar?? x times their error??)
        # TODO: verify we've got moments, and not EffMoments???
        # TODO: verify we've either been computed or read

        # get minimum significance
        minSignif = kwargs.pop( 'MinSignificance', float('-inf') )

        # get name requirements
        import re
        names = kwargs.pop( 'Names', None )

        # get scale factors
        scale  = kwargs.pop( 'Scale',  None )
        scales = kwargs.pop( 'Scales', ( scale, scale ) if scale != None else ( None, None ) )

        # get number of standard deviations for range of the PDF coefficients
        numStdDevs = kwargs.pop( 'RangeNumStdDevs', 5. )

        # get prefix for PDF coefficient names
        namePref = kwargs.pop( 'CoefNamePrefix', 'C_' )

        # loop over PDF terms
        keys     = []
        angFuncs = {}
        angCoefs = {}
        from RooFitWrappers import ConstVar, RealVar
        for ( name, coef, func ) in self._iterFuncAndCoef( MinSignificance = minSignif, Names = names ) :
            # construct the key and the function for the term
            keys.append( ( name, None ) )
            angFuncs[( name, None )] = ( func, None )

            # get the coefficient value and standard deviation
            coefVal = coef[0] * scales[0] if scales[0] else coef[0]
            coefErr = coef[1] * scales[1] if scales[1] else coef[1]

            # create the coefficient parameter for the PDF term
            if self._basisFuncIndices[name] == ( 0, 0, 0 ) :
                angCoefs[( name, None )] = ( ConstVar( Name = namePref + name, Value = coefVal ), None )
            else :
                coefMin = coefVal - numStdDevs * coefErr
                coefMax = coefVal + numStdDevs * coefErr
                angCoefs[( name, None )] = ( RealVar(  namePref + name, Value = coefVal, MinMax = ( coefMin, coefMax ) ), None )

        from P2VVParameterizations.AngularPDFs import Coefficients_AngularPdfTerms
        return Coefficients_AngularPdfTerms( Keys = keys, AngFunctions = angFuncs, AngCoefficients = angCoefs )

    def createPDF( self, **kwargs ) :
        # TODO: decide whether coefficients are ConstVar or RealVar?? (add keyword for that! -- what MinMax to give if RealVar??)
        #        maybe take MinMax = ( -5 * c[1], 5*c[1] ) ??? and make the 5 settable??
        # TODO: verify we've got moments, and not EffMoments???
        # TODO: verify we've either been computed or read
        scale = kwargs.pop('Scale', 1. )
        name = kwargs.pop('Name')
        ( names, coefs, funs ) = zip( *self._iterFuncAndCoef( MinSignificance = kwargs.pop( 'MinSignificance', float('-inf') )
                                                         , Names           = kwargs.pop( 'Names', None)
                                                         ) )
        from RooFitWrappers import ConstVar,RealSumPdf
        # TODO: renornmalize wrt. 0,0,0? require 0,0,0 to be present??
        return RealSumPdf( name
                         , functions = funs
                         , coefficients = ( ConstVar( Name = ('C_%3.6f'%c[0]).replace('-','m').replace('.','d')
                                                    , Value = c[0]*scale ) for c in coefs ) 
                         )

    def __mul__( self, pdf ) :
        from RooFitWrappers import Pdf
        if not isinstance(pdf, Pdf) : raise RuntimeError( 'trying to multiply a %s with %s ... this is not supported!'%(type(pdf),type(self) ) )
        return self.multiplyPDFWithEff( pdf )


    def multiplyPDFWithEff( self, pdf, **kwargs ) :

        def _createProduct( f1, f2, c, namePF ) :
            assert not f1.prod()
            assert not f2.prod()
            assert f1.c()!=0
            assert f2.c()!=0
            assert c!=0
            d1 = dict( (type(i),i) for i in f1.components() )
            d2 = dict( (type(i),i) for i in f2.components() )
            from ROOT import RooLegendre, RooSpHarmonic
            for i in [ RooLegendre, RooSpHarmonic ] :
                assert d1[i].getVariables().equals( d2[i].getVariables() )
            (cpsi,) = d1[RooLegendre].getVariables()
            (ctheta,phi) = d1[RooSpHarmonic].getVariables()
            from RooFitWrappers import P2VVAngleBasis
            # WARNING: cpsi,ctheta, phi are barebones PyROOT objects!!!
            return P2VVAngleBasis( {'cpsi':cpsi,'ctheta':ctheta,'phi':phi}
                                 , ( f1.i(),f1.j(),f1.l(),f1.m() )
                                 , f1.c()*f2.c()*c
                                 , ( f2.i(),f2.j(),f2.l(),f2.m() )
                                 , NamePostFix = namePF
                                 ) # build a wrapped object inside workspace
            
        # TODO: check that 'we' contain efficiency moments?
        # TODO: and that we've actually either 'read' or 'compute'-ed them??
        from ROOT import RooP2VVAngleBasis
        subst = dict()
        # TODO: do not use type to recognize, but name??
        from RooFitWrappers import Addition,EditPdf
        effName = kwargs.pop( 'EffName', 'eff' )
        for comp in filter( lambda x : type(x) is RooP2VVAngleBasis, pdf.getComponents() )  :
            subst[comp] = Addition( '%s_x_%s' % ( comp.GetName(), effName )
                                  , [ _createProduct( comp, f, c[0], effName ) for n,c,f in self._iterFuncAndCoef( Names = 'p2vvab.*' )  ]
                                  )
        return EditPdf( Name = kwargs.pop( 'Name', '%s_x_Eff' % pdf.GetName() ), Original = pdf, Rules = subst )


###########################################################################################################################################
## S-weights                                                                                                                             ##
###########################################################################################################################################

class SData( object ) :
    def __init__( self, **kwargs ) :
        # get input arguments
        def getKwArg( keyword, member, kwargs ) :
            if keyword in kwargs : setattr( self, '_' + member, kwargs.pop(keyword) )
            else : raise KeyError, 'P2VV - ERROR: SData.__init__(): key %s not found in input arguments' % keyword
        getKwArg( 'Name', 'name',      kwargs )
        getKwArg( 'Data', 'inputData', kwargs )
        getKwArg( 'Pdf',  'pdf',       kwargs )

        # initialize dictionary for weighted data sets per specie
        self._data = dict()

        # get yields and observables
        self._yields = [ par for par in self._pdf.Parameters() if par.getAttribute('Yield') ]
        self._observables = self._pdf.Observables()

        # calculate sWeights
        from ROOT import RooStats, RooArgList, RooSimultaneous
        if isinstance( self._pdf._var, RooSimultaneous ) and kwargs.pop( 'Simultaneous', True ) :
            # split data set in categories of the simultaneous PDF
            splitCat        = self._pdf.indexCat()
            splitCatIter    = splitCat.typeIterator()
            splitData       = self._inputData.split(splitCat)
            self._sPlots    = [ ]
            self._sDataSets = [ ]
            splitCatState   = splitCatIter.Next()
            sDataVars       = None
            from ROOT import RooFormulaVar
            while splitCatState :
                # calculate sWeights per category
                cat = splitCatState.GetName()
                data = splitData.FindObject(cat)

                origYieldVals = [ ( par.GetName(), par.getVal(), par.getError() ) for par in self._yields if cat in par.GetName() ]
                self._sPlots.append(  RooStats.SPlot( self._name + '_sData_' + cat, self._name + '_sData_' + cat
                                                     , data, self._pdf.getPdf(cat)
                                                     , RooArgList( par._var for par in self._yields if cat in par.GetName() ) )
                                   )
                self._sDataSets.append( self._sPlots[-1].GetSDataSet() )

                print 'P2VV - INFO: SData.__init__(): yields category %s:' % cat
                print '    original:',
                for vals in origYieldVals : print '%s = %.2f +/- %.2f  ' % vals,
                print '\n    new:     ',
                for par in self._yields :
                    if cat in par.GetName() : print '%s = %.2f +/- %.2f  ' % ( par.GetName(), par.getVal(), par.getError() ),
                print

                # add column for splitting category/categories (was removed when data set was split)
                # FIXME: How can we do this more generally? These are special cases and it might go wrong here...
                from ROOT import RooSuperCategory
                splitCat.setLabel(cat)
                __dref = lambda o : o._target_() if hasattr(o,'_target_') else o
                if isinstance( __dref(splitCat), RooSuperCategory ) :
                    for fundCat in splitCat.inputCatList() :
                        if not self._sDataSets[-1].get().find( fundCat.GetName() ) : self._sDataSets[-1].addColumn(fundCat)
                elif splitCat.isFundamental() and not self._sDataSets[-1].get().find( splitCat.GetName() ) :
                    self._sDataSets[-1].addColumn(splitCat)

                # add general sWeight and PDF value columns (it must be possible to simplify this...)
                # FIXME: in some cases "par.GetName().strip( '_' + cat )" goes wrong:
                # use "par.GetName()[ : par.GetName().find(cat) - 1 ]" instead
                # (case: 'N_bkgMass_notExclBiased'.strip('_notExclBiased') --> 'N_bkgM' ?!!!!!!)
                weightVars = [ (  RooFormulaVar( par.GetName()[ : par.GetName().find(cat) - 1 ] + '_sw', '', '@0'
                                                , RooArgList( self._sDataSets[-1].get().find( par.GetName() + '_sw' ) ) )
                                , RooFormulaVar( 'L_' + par.GetName()[ : par.GetName().find(cat) - 1 ], '', '@0'
                                                , RooArgList( self._sDataSets[-1].get().find( 'L_' + par.GetName() ) ) )
                               ) for par in self._yields if cat in par.GetName()
                             ]

                for weight, pdfVal in weightVars :
                    self._sDataSets[-1].addColumn(weight)
                    self._sDataSets[-1].addColumn(pdfVal)

                if not sDataVars :
                    # get set of variables in data
                    sDataVars = self._sDataSets[-1].get()
                    for par in self._yields :
                        if cat in par.GetName() :
                            sDataVars.remove( sDataVars.find( par.GetName() + '_sw' ) )
                            sDataVars.remove( sDataVars.find( 'L_' + par.GetName()  ) )

                splitCatState = splitCatIter.Next()

            # merge data sets from categories
            from ROOT import RooDataSet
            self._sData = RooDataSet( self._name + '_splotdata', self._name + '_splotdata', sDataVars )
            for data in self._sDataSets : self._sData.append(data)

        else :
            # calculate sWeights with full data set
            if isinstance( self._pdf._var, RooSimultaneous ) :
                print 'P2VV - WARNING: SData.__init__(): computing sWeights with a simultaneous PDF'
            self._sPlot = RooStats.SPlot( self._name + '_splotdata', self._name + '_splotdata', self._inputData, self._pdf._var
                                         , RooArgList( par._var for par in self._yields ) )
            self._sData = self._sPlot.GetSDataSet()

        # check keyword arguments
        if kwargs : raise KeyError, 'P2VV - ERROR: SData.__init__(): got unknown keywords %s for %s' % ( kwargs, type(self) )

    def usedObservables( self ) : return self._observables
    def components( self )      : return [ y.GetName()[2:] for y in self._yields ]
    def Pdf( self )             : return self._pdf

    def Yield( self, Component ) :
        yName = 'N_%s' % Component
        for y in self._yields :
            if y.GetName() == yName : return y.getVal()

        raise KeyError, 'P2VV - ERROR: SData.__init__(): unknown component %s' % Component

    def data( self, Component = None ) :
        if not Component : return self._sData

        if Component not in self._data :
            # check if component exists
            yName = 'N_%s' % Component
            if not any( yName in y.GetName() for y in self._yields ) :
                raise KeyError, 'P2VV - ERROR: SData.__init__(): unknown component: %s' % Component
            wName = '%s_sw' % yName
            if wName not in [ w.GetName() for w in self._sData.get() ] :
                raise KeyError, 'no weight in dataset for component %s' % Component

            # create weighted data set
            dName = '%s_weighted_%s' % ( self._sData.GetName(), Component )
            from ROOT import RooDataSet
            from ROOTDecorators import ROOTversion
            if ROOTversion[0] <= 5 and ROOTversion[1] <= 34 and ROOTversion[2] < 2 :
                self._data[Component] = RooDataSet( dName, dName, self._sData.get(), Import = self._sData, WeightVar = ( wName ) )
            else :
                self._data[Component] = RooDataSet( dName, dName, self._sData.get(), Import = self._sData, WeightVar = ( wName, True ) )

        return self._data[Component]


def createSData( **kwargs ) :
    # make sweighted dataset using J/psi phi mass
    Observables = kwargs.pop('Observables')
    Data        = kwargs.pop('Data')
    FitOpts     = kwargs.pop('FitOpts')
    Components  = kwargs.pop('Components')
    Name        = kwargs.pop('Name')
    from RooFitWrappers import buildPdf
    pdf = buildPdf( Components, Observables = Observables, Name= Name + '_splot_pdf')
    c_m = pdf.fitTo(Data,**FitOpts)
    c_state = dict( ( p, p.isConstant() ) for p in pdf.Parameters() )
    for p in pdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
    sdata =  SData(Pdf = pdf, Data = Data, Name = Name + '_sdata')
    for p,c in c_state.iteritems() : p.setConstant(c)
    return sdata


###########################################################################################################################################
## Miscellaneous                                                                                                                         ##
###########################################################################################################################################

def createProfile(name,data,pdf,npoints,param1,param1min,param1max,param2,param2min,param2max,NumCPU=8,Extend=True):
    print '**************************************************'
    print 'making profile for %s and %s'%(param1.GetName(),param2.GetName())
    print '**************************************************'

    nll = pdf.createNLL(data,RooFit.NumCPU(NumCPU),RooFit.Extended(Extend))
    profile = nll.createProfile(RooArgSet( param1,param2))
    return profile.createHistogram(name,         param1, RooFit.Binning(npoints,param1_min,param1_max)
                                  , RooFit.YVar( param2, RooFit.Binning(npoints,param2_min,param2_max))
                                  , RooFit.Scaling(False)
                                  )


# function for finding a splitted parameter in a simultaneous PDF
def getSplitPar( parName, stateName, parSet ) :
    from itertools import permutations
    stateName = stateName[ 1 : -1 ].split(';') if stateName.startswith('{') and stateName.endswith('}') else [ stateName ]
    if len(stateName) > 1 :
        fullNames = [ '%s_{%s}' % ( parName, ';'.join( comp for comp in perm ) )\
                     for perm in permutations( stateName, len(stateName) ) ]
    else :
        fullNames = [ '%s_%s' % ( parName, stateName[0] ) ]

    name = lambda par : par if type(par) == str else par.GetName()
    for par in parSet :
        if name(par) in fullNames : return par
    return None
