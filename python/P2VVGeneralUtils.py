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

def ROOTversion() :
    from ROOT import gROOT
    versionInt = gROOT.GetVersionInt()
    versionMajor = versionInt/10000
    versionMinor = versionInt/100 - versionMajor*100
    versionPatch = versionInt%100 
    return (versionMajor, versionMinor, versionPatch)


###########################################################################################################################################
## Handling Data                                                                                                                         ##
###########################################################################################################################################

def readData( filePath, dataSetName, NTuple = False, observables = None ) :
    """reads data from file (RooDataSet or TTree(s))
    """

    if NTuple :
      from ROOT import RooDataSet, TChain
      assert observables != None, "P2VV - ERROR: readData: set of observables is required for reading an NTuple"

      # create data set from NTuple file(s)
      print "P2VV - INFO: readData: reading NTuple(s) '%s' from file(s) '%s'" % ( dataSetName, filePath )
      chain = TChain(dataSetName)
      chain.Add(filePath)
      # TODO: how to check that filePath is a proper path? Explicitly checking for a file
      #       doesn't work, as it may be a wildcard...
      #       Also, file might not contain the right tree....

      noNAN = ' && '.join( '( %s==%s )' % ( obs, obs ) for obs in observables )
      data = RooDataSet( dataSetName, dataSetName, chain, [ obs._var for obs in observables ], noNAN )

    else :
      from ROOT import TFile

      # get data set from file
      print "P2VV - INFO: readData: reading RooDataset '%s' from file '%s'" % ( dataSetName, filePath )
      f = TFile.Open( filePath, 'READ' )
      assert f, "P2VV - ERROR: readData: file '%s' could not be opened" % filePath
      data = f.Get(dataSetName)
      f.Close()

    print 'P2VV - INFO: read dataset with %s entries' % data.numEntries()

    return data


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


def addTaggingObservables( dataSet, estimWTag, tagCatBinBoundaries = None ) :
    """add tagging observables to data set
    """

    # create binning
    from array import array
    from ROOT import RooBinning
    if not tagCatBinBoundaries : tagCatBinBoundaries = array( 'd', [ 0.50001, 0.49999, 0.38, 0.31, 0.24, 0.17, 0. ] )
    bins = RooBinning( len( tagCatBinBoundaries - 1 ), tagCatBinBoundaries, 'estWTagBins' )
    estimWTag.setBinning( bins, 'estWTagBins' )

    # create tagging category column
    from ROOT import RooBinningCategory
    tagCatFormula = RooBinningCategory( 'tagCatP2VV', 'P2VV tagging category', estimWTag, 'estWTagBins' )
    tagCat = dataSet.addColumn(tagCatFormula)
    tagCat.createFundamental()


###########################################################################################################################################
## Plots                                                                                                                                 ##
###########################################################################################################################################

# plot stash: keep the relevant objects alive by keeping a reference to them
global _P2VVPlotStash
_P2VVPlotStash = []

# plotting function
def plot(  canv, obs, data, pdf, addPDFs = [ ], components = None, xTitle = '', frameOpts = { }, dataOpts = { }, pdfOpts = { }
         , addPDFsOpts = [ { } ], plotResidHist = False, logy = False, normalize = True, symmetrize = True, usebar = True ) :
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
    _P2VVPlotStash.append(obsFrame)

    # set x-axis title
    if xTitle : xAxis.SetTitle(xTitle)

    # plot data
    if data : data.plotOn( obsFrame, Name = 'data', **dataOpts )

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
                drawOpts = dict( list(components[comp].items()) + list(pdfOpts.items()) )
                plotPDFWithSlices( pdf, obsFrame, 'comp%d' % num, Components = comp, **drawOpts )

        # plot total pdf
        plotPDFWithSlices( pdf, obsFrame, 'pdf', **pdfOpts )

        # draw data after drawing the PDF
        if data : obsFrame.drawAfter( 'pdf',  'data' )

    # plot additional PDFs
    if addPDFs :
        for num, addPDF in enumerate(addPDFs) :
            addPDF.plotOn( obsFrame, Name = 'addPDF' + str(num), **(addPDFsOpts[num]) )
            if data : obsFrame.drawAfter( 'addPDF' + str(num), 'data' )

    #TODO: add chisq/nbins
    #chisq = obsFrame.chiSquare( 'pdf', 'data' )
    #nbins = obsFrame.GetNbinsX()

    # get residuals histogram
    if plotResidHist and data and pdf :
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
        assert not kwargs

        # yield name, coefficient and function for selected moments
        for funcName in self._basisFuncNames :
            if self._coefficients[funcName][2] < minSignif : continue
            if nameExpr and not nameExpr.match(funcName)   : continue
            yield ( funcName, self._coefficients[funcName], self._basisFuncs[funcName] )

    def appendPYList( self, Angles, IndicesList, PDF = None, NormSet = None ) :
        # build moments from list of indices
        if not PDF and not NormSet :
            # build moment
            for inds in IndicesList :
                self.append(  Angles = Angles, PIndex = inds[0], YIndex0 = inds[1], YIndex1 = inds[2]
                                  , Norm = float( 2 * inds[0] + 1 ) / 2. )
        elif PDF and NormSet :
            # TODO: default for NormSet should be all observables in PDF except Angles, but
            #       we need a dataset to determine this ;-( 
            #       maybe set to 'None' , and then defer to compute???
            #       Problem is, at that point the moments have already been build...
            # build efficiency moment
            for inds in IndicesList :
                self.append(  Angles = Angles, PIndex = inds[0], YIndex0 = inds[1], YIndex1 = inds[2]
                                  , Norm = float( 2 * inds[0] + 1 ) / 2., PDF = PDF, NormSet = NormSet )
        else :
            print 'P2VV - ERROR: RealMomentsBuilder.appendList: both a PDF and a normalisation set are required for efficiency moments'

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

            if not 'PDF' in kwargs and not 'NormSet' in kwargs :
                # build moment
                from RooFitWrappers import RealMoment
                moment = RealMoment( func, kwargs.pop('Norm',1.) )
            elif 'PDF' in kwargs and 'NormSet' in kwargs :
                # build efficiency moment
                from RooFitWrappers import RealEffMoment
                moment = RealEffMoment( func, kwargs.pop('Norm',1.), kwargs.pop('PDF'), kwargs.pop('NormSet') )
            else :
                print 'P2VV - ERROR: RealMomentsBuilder.append: both a PDF and a normalisation set are required for an efficiency moment'
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
            self._basisFuncs[momName] = func
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
        scale = kwargs.pop('Scales', None )

        # print header
        print 'P2VV - INFO: RealMomentsBuilder.printMoments:'
        print '  name requirement: \'' + ( names if names else '' ) + '\''
        print '  minimum significance = %.1f' % minSignif
        print '  scale = ' + ( str(scale) if scale else '(1., 1., 1.)' )
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
                if scale :
                    print '  {0:<+12.4g}   {1:<12.4g}   {2:<12.4g}'.format( coef[0] * scale[0], coef[1] * scale[1], coef[2] * scale[2] )
                else :
                    print '  {0:<+12.4g}   {1:<12.4g}   {2:<12.4g}'.format( coef[0],            coef[1],            coef[2] )
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
            print 'P2VV - ERROR: RealMomentsBuilder.writeMoments: unable to open file \"%s\"' % filePath
            return

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

        cont += '# ' + '-' * (49 + maxLenName) + '\n\n'

        # write content to file
        momFile.write(cont)
        momFile.close()

        print 'P2VV - INFO: MomentsBuilder.writeMoments: %d efficiency moment%s written to file \"%s\"'\
                % ( numMoments, '' if numMoments == 1 else 's', filePath )

    def read( self, filePath = 'moments', **kwargs ) :
        self._coefficients = { }

        # get file path
        filePath = filePath.strip()

        # open file
        try :
          momFile = open(filePath, 'r')
        except :
          print 'P2VV - ERROR: MomentsBuilder.readMoments: unable to open file \"%s\"' % filePath
          return

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

        print 'P2VV - INFO: MomentsBuilder.readMoments: %d efficiency moment%s read from file \"%s\"'\
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
        scales = kwargs.pop( 'Scales', ( None, None ) )

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
        # TODO: verify we've got moments, and not EffMoments???
        # TODO: verify we've either been computed or read
        ( name, coef, fun ) = zip( *self._iterFuncAndCoef() )
        from RooFitWrappers import ConstVar,RealSumPdf
        return RealSumPdf( kwargs.pop('Name'), functions = fun, coefficients = ( ConstVar( Name = 'C_%3.6f'%c, Value = c ) for c in coef[0] ) )

    def __mul__( self, pdf ) :
        from RooFitWrappers import Pdf
        if not isinstance(pdf, Pdf) : raise RuntimeError( 'trying to multiply a %s with %s ... this is not supported!'%(type(pdf),type(self) ) )
        return self.multiplyPDFWithEff( pdf )


    def multiplyPDFWithEff( self, pdf, **kwargs ) :

        def _createProduct( f1, f2, c ) :
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
                                 ) # build a wrapped object inside workspace
            
        # TODO: check that 'we' contain efficiency moments?
        # TODO: and that we've actually either 'read' or 'compute'-ed them??
        from ROOT import RooP2VVAngleBasis
        subst = dict()
        # TODO: do not use type to recognize, but name??
        from RooFitWrappers import Addition,EditPdf
        for comp in filter( lambda x : type(x) is RooP2VVAngleBasis, pdf.getComponents() )  :
            subst[comp] = Addition( '%s_x_eff' % ( comp.GetName() )
                                  , [ _createProduct( comp, f, c[0] ) for n,c,f in self._iterFuncAndCoef( Names = 'p2vvab.*' )  ] 
                                  )
        return EditPdf( Name = kwargs.pop( 'Name', '%s_x_Eff' % pdf.GetName() ), Original = pdf, Rules = subst )
