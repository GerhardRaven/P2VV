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

def readData( filePath, dataSetName, NTuple = False, observables = None, **kwargs ) :
    """reads data from file (RooDataSet or TTree(s))
    """
    from ROOT import RooFit
    if observables : noNAN = ' && '.join( '( %s==%s )' % ( obs, obs ) for obs in observables )
    cuts = kwargs.pop( 'cuts', '' )

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

        # initialize sums for the weights and the weights squared per category
        sumWeights   = [ 0. ] * ( splitCat.numTypes() + 1 )
        sumSqWeights = [ 0. ] * ( splitCat.numTypes() + 1 )
        indexDict = { }
        posDict   = { }
        for iter, catType in enumerate( splitCat ) :
            posDict[ catType.getVal() ] = iter
            indexDict[ iter ] = catType.getVal()

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

    # add corrected weights to data set
    from ROOT import RooCorrectedSWeight
    if splitCatName :
        from ROOT import std
        alphaVec = std.vector('Double_t')()
        print 'P2VV - INFO: correctSWeights: multiplying sWeights (-ln(L)) to correct for background dilution with factors (overall factor %.4f):'\
              % ( sumWeights[0] / sumSqWeights[0] )
        for iter, ( sumW, sumSqW ) in enumerate( zip( sumWeights[ 1 : ], sumSqWeights[ 1 : ] ) ) :
            alphaVec.push_back( sumW / sumSqW )
            print '    %d: %.4f' % ( indexDict[iter], sumW / sumSqW )

        weightVar = RooCorrectedSWeight( 'weightVar', 'weight variable', bkgWeight, splitCat, alphaVec, True )

    else :
        print 'P2VV - INFO: correctSWeights: multiplying sWeights (-ln(L)) to correct for background dilution with a factor %.4f'\
              % ( sumWeights[0] / sumSqWeights[0] )
        weightVar = RooCorrectedSWeight( 'weightVar', 'weight variable', bkgWeight, sumWeights[0] / sumSqWeights[0], True )

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




def getCPprojectionOFpdf(pdf, cpComponent):
    #cpComponent: EVEN, ODD, SWAVE
    #Get the amplitudes
    parsDict = dict( (par.GetName(),par ) for par in pdf.Parameters()  )
    originalSet = set([pdf.ws().function("AparMag2"), parsDict['AperpMag2'], parsDict['A0Mag2'], parsDict['f_S']])

    if (cpComponent == "EVEN"):
        A0Mag2 = parsDict['A0Mag2'].getVal()
        AperpMag2 = 0
        AparMag2 = 1 - A0Mag2 - parsDict['AperpMag2'].getVal()
        f_S = 0 
    if (cpComponent == "ODD"):
        A0Mag2 = 0
        AperpMag2 = parsDict['AperpMag2'].getVal()
        AparMag2 = 0
        f_S = parsDict['f_S'].getVal()
    if (cpComponent == "SWAVE"):
        A0Mag2 = 0
        AperpMag2 = 0
        AparMag2 = 0
        f_S = parsDict['f_S'].getVal()

    #Create replacement ConstVars 
    from RooFitWrappers import ConstVar, CustomizePdf
    repl_AparMag2  = ConstVar(Name = "AparMag2" + cpComponent , Value = AparMag2  )
    repl_AperpMag2 = ConstVar(Name = "AperpMag2"+ cpComponent , Value = AperpMag2 )
    repl_A0Mag2    = ConstVar(Name = "A0Mag2"   + cpComponent , Value = A0Mag2    )
    repl_f_S       = ConstVar(Name = "f_S"      + cpComponent , Value =    f_S    )
    replacementSet = set([ repl_AparMag2,repl_AperpMag2, repl_A0Mag2, repl_f_S])

    CPcompPDF = CustomizePdf( pdf, origSet = originalSet, replSet = replacementSet, replString =  cpComponent )

    return CPcompPDF


def getNormOverObservables(pdf):
    #Do not integrate over the conditional observables thsi is done by plotOn
    observables = pdf.Observables() - pdf.ConditionalObservables()
    
    from ROOT import RooArgSet
    return pdf.getNorm(RooArgSet(obs._target_() for obs in observables))

 



###########################################################################################################################################
## Plots                                                                                                                                 ##
###########################################################################################################################################

# plot stash: keep the relevant objects alive by keeping a reference to them
global _P2VVPlotStash
_P2VVPlotStash = []

# plotting function
def plot(  canv, obs, data = None, pdf = None, addPDFs = [ ], components = None, xTitle = '', yTitle = '', xTitleOffset = None
           , yTitleOffset = None, yScale = ( None, None ), frameOpts = { }, dataOpts = { }, pdfOpts = { }, addPDFsOpts = [ { } ]
           , plotResidHist = False, logy = False, logx = False, normalize = True, symmetrize = True, usebar = True
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
        minimum = 0.
        if logy:
            hist = rooPlot.getHist()
            from ROOT import Double
            x = Double(0.)
            y = Double(0.)
            for i in range(hist.GetN()):
                r = hist.GetPoint(i, x, y)
                if y > 0 and y < minimum:
                    minimum = y
                if y < 0.:
                    hist.SetPoint(i, float(x), 0.)
                    minimum = 0.
            #hist.SetMinimum(minimum + 0.1)
            #rooPlot.SetMinimum(minimum + 0.1)

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
    if yScale[0] : obsFrame.SetMinimum(yScale[0])
    if yScale[1] : obsFrame.SetMaximum(yScale[1])

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

        # create residuals frame
        residFrame = obsFrame.emptyClone( obsFrame.GetName() + '_resid' )
        if 'Title' in frameOpts: residFrame.SetTitle(frameOpts['Title'])
        xAxis = residFrame.GetXaxis()
        _P2VVPlotStash.append(residFrame)

        # set minimum for observable's frame if there is a log scale for y
        if logy : obsFrame.SetMinimum(0.1)

        # set residual plot options
        #TODO: if normalize : plot residHist as a filled histogram with fillcolor blue...
        #      or, maybe, with the 'bar chart' options: 'bar' or 'b'
        if dataOpts :
            fun = { 'MarkerSize'  : lambda x : residHist.SetMarkerSize(x) 
                  , 'MarkerStyle' : lambda x : residHist.SetMarkerStyle(x)
                  , 'MarkerColor' : lambda x : residHist.SetMarkerColor(x)
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
        if logx: obsPad.SetLogx(1)
        obsPad.SetNumber(1)
        obsPad.SetLeftMargin(0.12)
        obsPad.Draw()
        canv.cd(1)
        if 'Title' in frameOpts and not frameOpts['Title']:
            obsFrame.SetTitle("")
        obsFrame.Draw()

        # draw residuals frame
        canv.cd()
        residName = obs.GetName() + '_resid1'
        residPad = TPad( residName, residName, 0, 0, 1, 0.2 )
        if logx: residPad.SetLogx(1)
        _P2VVPlotStash.append(residPad)
        residPad.SetNumber(2)
        residPad.SetLeftMargin(0.12)
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

        # get name requirements
        import re
        names = kwargs.pop('Names', None)
        nameExpr = re.compile(names) if names else None

        # write moments to content string
        cont = '# %s: angular moments\n' % fileName\
             + '# name requirement: \'{0}\'\n'.format( names if names else '' )\
             + '# minimum significance = {0:.1f}\n'.format(minSignif)\
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
                cont += '   {0:<+14.8g}   {1:<13.8g}   {2:<13.8g}\n'.format(coef[0], coef[1], coef[2])
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

