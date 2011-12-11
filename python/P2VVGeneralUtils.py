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
    assert file
    data = file.Get(dataSetName)
    file.Close()

  return data


def writeData( filePath, dataSetName, data, NTuple = False ) :
  """writes data to file (RooDataSet or TTree)
  """

  from ROOT import TFile

  print "P2VV - INFO: writeData: writing RooDataSet '%s' to file '%s'" % ( dataSetName, filePath )

  file = TFile.Open( filePath, 'RECREATE' )
  assert file
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
def plot(  canv, obs, data, pdf, components = None, xTitle = '', frameOpts = { }, dataOpts = { }, pdfOpts = { }, plotResidHist = False
         , logy = False, normalize = True, symmetrize = True, usebar = True ) :
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
        if components :
            for comp, opts in components.iteritems() :
                drawOpts = dict( list(opts.items()) + list(pdfOpts.items()) )
                pdf.plotOn(obsFrame, Components = comp, **drawOpts )
        pdf.plotOn( obsFrame, Name = 'pdf', **pdfOpts )

    # draw data after drawing the PDF
    if data and pdf : obsFrame.drawAfter( 'pdf', 'data' )

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
            for k, v in dataOpts.iteritems() : fun[k](v)

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
    def __init__(self,**kw)   :
        self._basisFuncNames = [ ]
        self._basisFuncs     = { }
        self._coefficients   = { }
        if 'Moments' in kw :
            for i in kw.pop('Moments') : self.append( Moment = i )
        if 'Moment' in kw :
            self.append( Moment = kw.pop(Moment) )
        if kw :
            raise RuntimeError( 'unknown keyword arguments %s' % kw.keys() )

    def basisFuncs(self)     : return self._basisFuncs.copy()
    def basisFuncNames(self) : return self._basisFuncNames[ : ]
    def coefficients(self)   : return self._coefficients.copy()

    def appendPYList( self, Angles, IndicesList, PDF = None, NormSet = None ) :
        # build moments from list of indices
        if not PDF and not NormSet :
            # build moment
            for inds in IndicesList :
                self.append(  Angles = Angles, PIndex = inds[0], YIndex0 = inds[1], YIndex1 = inds[2]
                                  , Norm = float( 2 * inds[0] + 1 ) / 2. )
        elif PDF and NormSet :
            # build efficiency moment
            for inds in IndicesList :
                self.append(  Angles = Angles, PIndex = inds[0], YIndex0 = inds[1], YIndex1 = inds[2]
                                  , Norm = float( 2 * inds[0] + 1 ) / 2., PDF = PDF, NormSet = NormSet )
        else :
            print 'P2VV - ERROR: RealMomentsBuilder.appendList: both a PDF and a normalisation set are required for efficiency moments'

    def append( self, **kwargs ) :
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
                func = P2VVAngleBasis(kwargs.pop('Angles'), kwargs.pop('PIndex'), 0, kwargs.pop('YIndex0'), kwargs.pop('YIndex1'), 1.)

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
        scale = kwargs.pop('Scale', None )

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
        from P2VVParameterizations.AngularPDFs import Coefficients_AngularPdfTerms
        return Coefficients_AngularPdfTerms()

    def _iterFuncAndCoef( self ) :
        for fn in self._basisFuncNames :
            yield (self._basisFuncs[fn], self._coefficients[fn])

    def multiplyPDFWithEff( self, pdf, **kwargs ) :
        pdfName = kwargs.pop('Name','%s_x_Eff'%pdf.GetName())
        # TODO: check that 'we' contain efficiency moments?
        # TODO: and that we've actually either 'read' or 'compute'-ed them??

        from ROOT import RooP2VVAngleBasis,RooAddition,RooArgSet
        ws = pdf.ws()
        subst = dict()
        for comp in pdf.getComponents() :
            if type(comp) is not RooP2VVAngleBasis : continue  # TODO: do not use type to recognize, but name??
            name  = '%s_%s_eff' % ( pdfName, comp.GetName() )
            effTerms = RooArgSet( comp.createProduct( f,c ) for f,c in self._iterFuncAndCoef() )
            subst[comp.GetName()] = ws.put( RooAddition( name, name, effTerms, True ) ).GetName() 
        # TODO: the returned object ought to be wrapped in a Pdf class...
        return ws.factory('EDIT::%s(%s,%s)' % ( pdfName, pdf.GetName(), ','.join( '%s=%s'%(v,k) for k,v in subst.iteritems() ) ) )
