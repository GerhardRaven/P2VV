###########################################################################################################################################
## Utilities.Studies: P2VV utilities for toys and systematics studies                                                                    ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                                                                            ##
##                                                                                                                                       ##
###########################################################################################################################################

# class to generate parameter variations based on measured values and covariance matrix
class parameterGen(object) :
    def __init__( self, Values, Covariances, ParNames = [ ], StartID = 10**9, ParFile = '' ) :
        # get number of parameters
        self._nPars = len(Values)
        if self._nPars < 1 : return

        # read parameter file
        self._filePars = [ ]
        self.readParameters(ParFile)

        # set parameter names
        if ParNames :
            assert len(ParNames) == self._nPars
            self._parNames = ParNames
        else :
           self._parNames = [ 'par%02d' % parIt for parIt in range(self._nPars) ]

        # parse covariance matrix
        from ROOT import TMatrixD
        if hasattr( Covariances, 'ClassName' ) and Covariances.ClassName() in [ 'TMatrixT<double>', 'TMatrixTSym<double>' ] :
            assert Covariances.GetNrows() == self._nPars and Covariances.GetNcols() == self._nPars\
                   , 'P2VV - ERROR: parameterGen: size of covariance matrix does not match number of parameters'
            covMat = Covariances
        else :
            assert len(Covariances) == self._nPars and len( Covariances[0] ) == self._nPars\
                   , 'P2VV - ERROR: parameterGen: size of covariance matrix does not match number of parameters'
            covMat = TMatrixD( self._nPars, self._nPars )
            for it0, covs in enumerate(Covariances) :
                for it1, cov in enumerate(covs) :
                    covMat[it0][it1] = cov

        # transform parameters to an uncorrelated set
        from ROOT import TVectorD
        vals = TVectorD(self._nPars)
        ranges = [ ]
        for it, val in enumerate(Values) :
            if hasattr( val, '__iter__' ) :
                assert len(val) == 3
                vals[it] = val[0]
                ranges.append( val[ 1 : ] )
            else :
                vals[it] = val
        assert not ranges or len(ranges) == self._nPars
        diagVars = TVectorD(self._nPars)
        invTrans = covMat.EigenVectors(diagVars)
        trans = TMatrixD(invTrans).Invert()
        diagVals = trans * vals

        # set members needed for parameter generation
        from math import sqrt
        self._diagVals = [       diagVals[it]   for it in range(self._nPars) ]
        self._diagErrs = [ sqrt( diagVars[it] ) for it in range(self._nPars) ]
        self._genDiagVals = TVectorD(self._nPars)
        self._genVals = [ ]
        self._invTrans = invTrans
        self._ranges = ranges
        self._nGen = 0
        self._nOutRange = 0
        self._sumVals = [ 0. ] * self._nPars
        self._sumSqVals = [ 0. ] * self._nPars

        # initialize random number generator
        assert type(StartID) == int and StartID > 0 and StartID < 2**32\
               , 'P2VV - ERROR: parameterGen: StartID is required to be between 0 and 2^32 (got %d)' % StartID
        from ROOT import TRandom3
        self._parSetID = StartID
        self._randGen = TRandom3(StartID)

    def setParSetID( self, ID ) :
        if self._nPars < 1 : return
        self._parSetID = ID
        self._randGen.SetSeed(ID)

    def parSetID(self) :
        if self._nPars < 1 : return None
        return self._parSetID

    def readParameters( self, File ) :
        if File :
            from P2VV.Parameterizations.FullPDFs import PdfConfiguration
            self._filePars.append( PdfConfiguration() )
            self._filePars[-1].readParametersFromFile( filePath = File )
            return len(self._filePars) - 1

    def writeParameters( self, File, TransFuncs = { }, Index = -1 ) :
        if self._nPars < 1 or not self._genVals or not self._filePars or Index >= len(self._filePars) or Index < -len(self._filePars) :
            print 'P2VV - ERROR: parameterGen.writeParameters(): no parameters to write to file'
            return

        # set parameters for output file
        funcs = TransFuncs if TransFuncs else dict( [ ( name, val ) for name, val in zip( self._parNames, self._genVals ) ] )
        pars = self._filePars[Index].parameters()
        parsFound = True
        for name, func in funcs.iteritems() :
            if name not in pars :
                parsFound = False
                continue
            pars[name] = tuple( [ func(self._genVals) if callable(func) else func ] + [ v for v in pars[name][ 1 : -1 ] ] + [ False ] )

        if not parsFound and self._parsFoundWarn :
            print 'P2VV - WARNING: parameterGen.generateParSet(): not all generated parameters found in parameter set for output file'
            self._parsFoundWarn = False

        # write parameter file
        return self._filePars[Index].writeParametersToFile( filePath = File, Verbose = False )

    def generateParSet(self) :
        if self._nPars < 1 : return [ ]

        while True :
            # generate values for uncorrelated parameters
            for parIt, ( val, err ) in enumerate( zip( self._diagVals, self._diagErrs ) ) :
                self._genDiagVals[parIt] = self._randGen.Gaus( val, err )

            # transform generated values to original correlated-parameter values
            genVals = self._invTrans * self._genDiagVals

            # get generated values and check if they are in range
            self._genVals = [ ]
            for parIt in range(self._nPars) :
                val = genVals[parIt]
                if self._ranges and ( val < self._ranges[parIt][0] or val > self._ranges[parIt][1] ) : break
                self._genVals.append(val)
            if len(self._genVals) == self._nPars : break
            self._nOutRange += 1

        # update generation statistics
        self._nGen += 1
        self._parSetID += 1
        for parIt, val in enumerate(self._genVals) :
            self._sumVals[parIt]   += val
            self._sumSqVals[parIt] += val**2

        return self._genVals

    def genVals(self) :
        if self._nPars < 1 : return [ ]
        return self._genVals

    def nGen(self) :
        if self._nPars < 1 : return 0
        return self._nGen

    def nOutRange(self) :
        if self._nPars < 1 : return 0
        return self._nOutRange

    def stats(self) :
        if self._nPars < 1 : return [ ]
        from math import sqrt
        return [ ( name, sumVals / float(self._nGen), sqrt( sumSqVals / float(self._nGen) - ( sumVals / float(self._nGen) )**2 ) )\
                 for name, sumVals, sumSqVals in zip( self._parNames, self._sumVals, self._sumSqVals ) ]

    def printStats(self) :
        print 'P2VV - INFO: parameterGen.printStats(): generation statistics:'
        if self._nPars < 1 :
            print '  no parameters to generate'
            return

        print '  generated %d parameter sets' % self._nGen
        print '  %d sets generated outside range' % self._nOutRange
        nameLen = min( 30, max( len(name) for name in self._parNames ) )
        from math import log10, ceil
        for name, val, err in self.stats() :
            prec = max( 0, 3 - int( ceil( log10(err) ) ) )
            print ( '  {0:<%ds}  {1:<+10.%df} +/- {2:<10.%df}' % ( nameLen, prec, prec ) ).format( name, val, err )

# class to analyse fit results
class fitResultsAnalysis(object) :
    def __init__( self, ParNames, AnaFiles = [ ], RefFile = '' ) :
        assert ParNames and hasattr( ParNames, '__iter__' ) and hasattr( ParNames, '__len__' ) and len(ParNames) > 0, 'P2VV'\
               'P2VV - ERROR: fitResultsAnalysis: no parameter names specified'
        self._parNames = ParNames
        from P2VV.Parameterizations.FullPDFs import PdfConfiguration
        self._refPars = PdfConfiguration()
        self._anaPars = PdfConfiguration()
        self._refParVals = [ None ] * len(self._parNames)
        self._anaParVals = { }
        self._parHists = { }
        self._nFiles = 0
        self.readRefFile(RefFile)
        self.readAnaFiles(AnaFiles)

    def readRefFile( self, file ) :
        self._refParVals = [ None ] * len(self._parNames)
        if file :
            self._refPars.readParametersFromFile( filePath = file )
            parVals = self._refPars.parameters()
            for parIt, name in enumerate(self._parNames) :
                assert name in parVals\
                       , 'P2VV - ERROR: fitResultsAnalysis.readRefFile(): parameter "%s" not found in file "%s"' % ( name, file )
                self._refParVals[parIt] = ( parVals[name][0], parVals[name][1] )

    def readAnaFiles( self, files ) :
        self._anaParVals = dict( [ ( name, [ ] ) for name in self._parNames ] )
        self._nFiles = 0
        if files :
            fileStats = { }
            for file in files :
                fitStatus = self._anaPars.readParametersFromFile( filePath = file, Verbose = False )
                if fitStatus[0] :
                    print 'P2VV - WARNING: fitResultsAnalysis.readAnaFiles(): fit status %d for file "%s"' % ( fitStatus[0], file )
                if not fitStatus[0] in fileStats :
                    fileStats[ fitStatus[0] ] = 1
                else :
                    fileStats[ fitStatus[0] ] += 1

                parVals = self._anaPars.parameters()
                for name, vals in self._anaParVals.iteritems() :
                    assert name in parVals\
                           , 'P2VV - ERROR: fitResultsAnalysis.readAnaFiles(): parameter "%s" not found in file "%s"' % ( name, file )
                    vals.append( parVals[name][ : 2 ] )

            self._nFiles = sum( stats for stats in fileStats.itervalues() )
            print 'P2VV - INFO: fitResultsAnalysis.readAnaFiles(): read %d parameter files for analysis (fit status: %s)'\
                  % ( self._nFiles, ', '.join( '%d: %d' % ( stat, count ) for stat, count in fileStats.iteritems() ) )

    def processResults( self, **kwargs ) :
        if not self._anaParVals or self._nFiles == 0 :
            print 'P2VV - WARNING: fitResultsAnalysis.processResults(): no parameter values available for analysis'
            return

        histsFile = kwargs.pop( 'HistsFile', '' ).split('.')[0]
        anaParSetts = kwargs.pop( 'ParSettings', { } )
        toyMode = kwargs.pop( 'ToyMode', False )

        self._nParVals = [ [ 0, 0, 0 ] for name in self._parNames ]
        self._maxPulls = [ 0. for name in self._parNames ]
        self._parVals = [ [ ] for name in self._parNames ]
        self._parSums = [ [ 0., 0. ] for name in self._parNames ]
        self._errSums = [ 0. for name in self._parNames ]
        self._pullVals = [ [ ] for name in self._parNames ]
        self._pullSums = [ [ 0., 0., 0., 0., 0., 0. ] for name in self._parNames ]
        if histsFile :
            from P2VV.Load import LHCbStyle
            from ROOT import TH1D, TLine, TCanvas, kFullDotLarge, kBlue
            dCanv = TCanvas('dummy')
            dCanv.Print( histsFile + '_pars.pdf[' )
            dCanv.Print( histsFile + '_pulls.pdf[' )

            def drawHist( parName, parHist, fileName, doFit = False, refVal = None ) :
                print 'P2VV - INFO: fitResultsAnalysis.processResults(): drawing "%s" histogram' % parName
                parSetts = anaParSetts[parName] if parName in anaParSetts else { }
                self._parHists[parName] = parHist
                parHist.GetXaxis().SetNdivisions( 7, 5, 0 )
                parHist.SetTitleOffset( 1.15, 'x' )
                parHist.SetTitleOffset( 1.05, 'y' )
                parHist.SetXTitle( parSetts['name'] if 'name' in parSetts else parName )
                parHist.SetYTitle('Entries / %.2g' % parHist.GetBinWidth(1) )
                parHist.SetMarkerStyle(kFullDotLarge)
                parHist.SetMarkerSize(0.7)
                parHist.SetLineWidth(2)
                if doFit :
                    parHist.Fit('gaus')
                    fitFunc = parHist.GetFunction('gaus')
                    fitFunc.SetNpx(1000)
                    fitFunc.SetLineWidth(3)
                    fitFunc.SetLineColor(kBlue)

                histMax = -1.
                for bin in range( 1, parHist.GetNbinsX() + 1 ) :
                    val = parHist.GetBinContent(bin) + parHist.GetBinErrorUp(bin)
                    if val > histMax : histMax = val
                assert histMax > 0.
                histMax *= 1.05
                parHist.SetMinimum(0.)
                parHist.SetMaximum(histMax)

                canv = TCanvas(parName)
                canv.SetLeftMargin(0.19)
                canv.SetRightMargin(0.05)
                canv.SetBottomMargin(0.21)
                canv.SetTopMargin(0.05)
                parHist.Draw('E1')
                if refVal != None :
                    line = TLine()
                    line.DrawLine( refVal, 0., refVal, histMax )
                canv.Print(fileName)

        for parIt, name in enumerate(self._parNames) :
            isPhase = name in anaParSetts and 'isPhase' in anaParSetts[name] and anaParSetts[name]['isPhase']
            from math import floor, pi
            for vals in self._anaParVals[name] :
                self._nParVals[parIt][0] += 1
                val = vals[0]
                if self._refParVals[parIt] != None :
                    diff = val - self._refParVals[parIt][0]
                    if isPhase and abs(diff) >= pi :
                        nShift = floor( 0.5 * ( diff / pi + 1. ) )
                        val -= nShift * 2. * pi
                        diff -= nShift * 2. * pi
                    pull = diff / ( vals[1] if toyMode else self._refParVals[parIt][1] )
                    self._nParVals[parIt][ 2 if pull < 0. else 1 ] += 1
                    if abs(pull) > self._maxPulls[parIt] : self._maxPulls[parIt] = abs(pull)
                    self._pullVals[parIt].append(pull)
                    self._pullSums[parIt][ 1 if pull < 0. else 0 ] += pull
                    self._pullSums[parIt][ 3 if pull < 0. else 2 ] += pull**2
                    self._pullSums[parIt][ 5 if pull < 0. else 4 ] += pull**4
                self._parVals[parIt].append(val)
                self._parSums[parIt][0] += val
                self._parSums[parIt][1] += val**2
                self._errSums[parIt] += vals[1]

            if name in self._parHists and self._parHists[name] :
                self._parHists[name].Delete()
            if histsFile :
                parMinMax = anaParSetts[name]['parMinMax'] if name in anaParSetts and 'parMinMax' in anaParSetts[name] else None
                parHistBins = anaParSetts[name]['parHistBins'] if name in anaParSetts and 'parHistBins' in anaParSetts[name]\
                              else int( float(self._nParVals[parIt][0]) / 250. ) if self._nParVals[parIt][0] > 1000 else 10
                plotShift = anaParSetts[name]['plotShift'] if name in anaParSetts and 'plotShift' in anaParSetts[name] else 0.
                parHistMin = parMinMax[0] if parMinMax else min( self._parVals[parIt] ) + plotShift
                parHistMax = parMinMax[1] if parMinMax else max( self._parVals[parIt] ) + plotShift
                parHistRange = parHistMax - parHistMin
                if not parMinMax and parHistRange > 0. :
                    parHistMin = parHistMin - 0.01 * parHistRange
                    parHistMax = parHistMax + 0.01 * parHistRange
                parHist = TH1D( name + '_par', name, parHistBins, parHistMin, parHistMax )
                for val in self._parVals[parIt] : parHist.Fill( val + plotShift )
                doParFit = name in anaParSetts and 'doParFit' in anaParSetts[name] and anaParSetts[name]['doParFit']
                drawHist( name, parHist, histsFile + '_pars.pdf', doParFit
                         , self._refParVals[parIt][0] + plotShift if self._refParVals[parIt] != None else None )

                pullMinMax = anaParSetts[name]['pullMinMax'] if name in anaParSetts and 'pullMinMax' in anaParSetts[name] else None
                pullHistBins = anaParSetts[name]['pullHistBins'] if name in anaParSetts and 'pullHistBins' in anaParSetts[name]\
                               else int( float(self._nParVals[parIt][0]) / 250. ) if self._nParVals[parIt][0] > 1000 else 10
                pullHistMin = pullMinMax[0] if pullMinMax else min( self._pullVals[parIt] )
                pullHistMax = pullMinMax[1] if pullMinMax else max( self._pullVals[parIt] )
                pullHistRange = pullHistMax - pullHistMin
                if not pullMinMax and pullHistRange > 0. :
                    pullHistMin = pullHistMin - 0.01 * pullHistRange
                    pullHistMax = pullHistMax + 0.01 * pullHistRange
                pullHist = TH1D( name + '_pull', name, pullHistBins, pullHistMin, pullHistMax )
                for val in self._pullVals[parIt] : pullHist.Fill(val)
                doPullFit = name in anaParSetts and 'doPullFit' in anaParSetts[name] and anaParSetts[name]['doPullFit']
                drawHist(name, pullHist, histsFile + '_pulls.pdf', doPullFit, 0. )

        if histsFile :
            dCanv.Print( histsFile + '_pars.pdf]' )
            dCanv.Print( histsFile + '_pulls.pdf]' )

        for valCounts in self._nParVals :
            assert valCounts[0] == self._nFiles and ( valCounts[1] + valCounts[2] == self._nFiles or self._refParVals[0] == None )

        from math import sqrt, log10, ceil
        nameLen = min( 30, max( len(name) for name in self._parNames ) )
        sepStr = ' ' + '-' * ( nameLen + ( 167 if self._refParVals[0] != None else 32 ) )
        print 'P2VV - INFO: fitResultsAnalysis.processResults(): parameter statistics for %d files:' % self._nFiles
        print sepStr
        print ( '  {0:<%ds}   {1:<22s}' % nameLen ).format( 'name', 'measured' ),
        if self._refParVals[0] != None :
            print '   {0:<18s}   {1:<20s}   {2:<19s}   {3:<19s}   {4:<19s}   {5:<6s}   {6:<6s}   {7:<8s}'\
                  .format( 'reference', 'pull mean', 'mean error', 'pull RMS', 'width error', '+RMS', '-RMS', 'max pull' )
        else :
            print
        print sepStr
        for name, nParVals, maxPull, parSums, errSums, pullSums, refVal\
                in zip( self._parNames, self._nParVals, self._maxPulls, self._parSums, self._errSums, self._pullSums, self._refParVals ) :
            meanVal = parSums[0] / float( nParVals[0] )
            meanSqVal = parSums[1] / float( nParVals[0] )
            stdDev = sqrt( meanSqVal - meanVal**2 )
            precDev = max( 0, 3 - int( ceil( log10(stdDev) ) ) )
            prec = max( 0, 3 - int( ceil( log10( refVal[1] ) ) ) ) if refVal != None else precDev
            print ( '  {0:<%ds}   {1:<+8.%df}   {2:<11.%df}' % ( nameLen, prec, precDev ) ).format( name, meanVal, stdDev ),
            if refVal != None :
                refErr = errSums / float( nParVals[0] ) if toyMode else refVal[1]
                pullMean = ( pullSums[0] + pullSums[1] ) / float( nParVals[0] )
                pullVar = ( pullSums[2] + pullSums[3] ) / float( nParVals[0] )
                pullVarVar = ( pullSums[4] + pullSums[5] ) / float( nParVals[0] )
                meanErr = sqrt( ( pullVar - pullMean**2 ) / float( nParVals[0] ) )
                width = sqrt(pullVar)
                widthErr = 0.5 / sqrt(pullVar) * sqrt( ( pullVarVar - pullVar**2 ) / float( nParVals[0] ) )
                widthPlus = sqrt( pullSums[2] / float( nParVals[1] ) ) if nParVals[1] > 0 else 0.
                widthMin  = sqrt( pullSums[3] / float( nParVals[2] ) ) if nParVals[2] > 0 else 0.
                print ( '   {0:<+8.%df}   {1:<7.%df}   {2:<+6.3f} / {3:<+11.%df}   {4:<5.3f} / {5:<11.%df}   {6:<5.3f} / {7:<11.%df}   {8:<5.3f} / {9:<11.%df}   {10:<6.3f}   {11:<6.3f}   {12:<8.3f}'\
                        % ( prec, prec, precDev, precDev, precDev, precDev ) )\
                        .format( refVal[0], refErr, pullMean, pullMean * refErr, meanErr, meanErr * refErr
                                , width, width * refErr, widthErr, widthErr * refErr, widthPlus, widthMin, maxPull )
            else :
                print
        print sepStr

# distributions of conditional observables
condBins = dict(  KKMass = [  1008.0, 1012.0, 1013.7, 1014.7, 1015.5, 1016.0, 1016.4, 1016.8, 1017.2, 1017.5, 1017.7, 1017.9, 1018.1
                            , 1018.3, 1018.5, 1018.7, 1018.9, 1019.1, 1019.3, 1019.5, 1019.7, 1019.8, 1019.9, 1020.0, 1020.1, 1020.2
                            , 1020.3, 1020.5, 1020.7, 1020.9, 1021.2, 1021.5, 1021.8, 1022.1, 1022.4, 1022.7, 1023.1, 1023.5, 1024.0
                            , 1024.5, 1025.2, 1026.0, 1027.0, 1028.2, 1029.9, 1032.0, 1034.8, 1038.7, 1043.6
                           ]
                , sigmat = [  0.0167, 0.0184, 0.0196, 0.0206, 0.0214, 0.0222, 0.0228, 0.0235, 0.0240
                            , 0.0246, 0.0251, 0.0256, 0.0261, 0.0265, 0.0270, 0.0275, 0.0279, 0.0283, 0.0288
                            , 0.0292, 0.0297, 0.0301, 0.0305, 0.0310, 0.0314, 0.0319, 0.0323, 0.0328, 0.0332
                            , 0.0337, 0.0342, 0.0346, 0.0352, 0.0357, 0.0362, 0.0367, 0.0373, 0.0380, 0.0386
                            , 0.0394, 0.0401, 0.0410, 0.0420, 0.0431, 0.0444, 0.0461, 0.0482, 0.0512, 0.0571
                           ]
                , etaOS  = [ 0.336, 0.400, 0.450, 0.4999999 ]
                , etaSS  = [ 0.397, 0.437, 0.456, 0.469, 0.481, 0.492, 0.4999999 ]
               )
condDists = dict(  runPeriod = dict(  sig  = { 2011 : 0.310319, 2012 : 0.689681 }
                                    , cbkg = { 2011 : 0.281569, 2012 : 0.718431 }
                                   )
                 , HLT1Cat   = dict(  sig       = { 0 : 0.818302, 1 : 0.181698 }
                                    , cbkg      = { 0 : 0.805423, 1 : 0.194577 }
                                    , runPeriod = dict(  sig  = {  2011 : { 0 : 0.843546, 1 : 0.156454 }
                                                                 , 2012 : { 0 : 0.806944, 1 : 0.193056 }
                                                                }
                                                       , cbkg = {  2011 : { 0 : 0.856037, 1 : 0.143963 }
                                                                 , 2012 : { 0 : 0.785586, 1 : 0.214414 }
                                                                }
                                                      )
                                   )
                 , KKMass    = dict(  sig  = [  0.025814, 0.020853, 0.020693, 0.020291, 0.023258, 0.019357, 0.018538, 0.021672
                                              , 0.026782, 0.022449, 0.017277, 0.018419, 0.018845, 0.020021, 0.022168, 0.022035
                                              , 0.022217, 0.023525, 0.023312, 0.023698, 0.024279, 0.011650, 0.012220, 0.012185
                                              , 0.011229, 0.012069, 0.010915, 0.022498, 0.022111, 0.020329, 0.028616, 0.026382
                                              , 0.024512, 0.021714, 0.019688, 0.017839, 0.021438, 0.017724, 0.018935, 0.017509
                                              , 0.019177, 0.019780, 0.018714, 0.018264, 0.019433, 0.018454, 0.020260, 0.020182
                                              , 0.020574, 0.020097
                                             ]
                                    , cbkg = [  0.180933, 0.055585, 0.026909, 0.016060, 0.014051, 0.008910, 0.007640, 0.007740
                                              , 0.008299, 0.006528, 0.004168, 0.004494, 0.004574, 0.004746, 0.004987, 0.004478
                                              , 0.004782, 0.005055, 0.005344, 0.004877, 0.004845, 0.002650, 0.002628, 0.002381
                                              , 0.002346, 0.002329, 0.002657, 0.005422, 0.004872, 0.004865, 0.007424, 0.007199
                                              , 0.006828, 0.006684, 0.006643, 0.006871, 0.008626, 0.008385, 0.010162, 0.009698
                                              , 0.013198, 0.014997, 0.018620, 0.022847, 0.030789, 0.040792, 0.052295, 0.076255
                                              , 0.099805, 0.136729
                                             ]
                                   )
                 , sigmat    = dict(  sig  = [  0.019796, 0.019495, 0.019904, 0.020579, 0.018815, 0.021543, 0.018069, 0.023126
                                              , 0.017016, 0.022413, 0.019565, 0.019905, 0.021382, 0.016410, 0.021358, 0.022354
                                              , 0.018107, 0.018131, 0.022602, 0.017596, 0.022001, 0.019149, 0.018955, 0.022684
                                              , 0.018234, 0.023206, 0.017898, 0.022820, 0.017114, 0.021308, 0.021153, 0.016184
                                              , 0.023660, 0.020557, 0.018899, 0.018535, 0.019962, 0.022308, 0.017698, 0.022120
                                              , 0.018305, 0.021128, 0.020204, 0.019412, 0.019030, 0.020753, 0.019572, 0.019835
                                              , 0.020160, 0.018988
                                             ]
                                    , cbkg = [  0.010222, 0.013613, 0.016119, 0.018777, 0.018120, 0.020468, 0.018103, 0.023072
                                              , 0.017782, 0.022384, 0.019973, 0.021010, 0.021865, 0.017964, 0.021964, 0.023093
                                              , 0.017952, 0.018280, 0.022960, 0.018680, 0.023386, 0.018765, 0.018531, 0.023194
                                              , 0.017584, 0.021278, 0.018691, 0.022198, 0.016584, 0.020672, 0.019400, 0.015299
                                              , 0.023210, 0.017422, 0.017133, 0.016945, 0.020554, 0.021354, 0.017000, 0.021909
                                              , 0.017957, 0.021104, 0.021434, 0.021143, 0.020627, 0.024306, 0.022801, 0.023320
                                              , 0.025919, 0.027882
                                             ]
                                   )
                 , etaTag    = dict(  sig  = [  [ 0.006587, 0.006158, 0.005884, 0.005483, 0.005284, 0.005387, 0.005798, 0.036124 ]
                                              , [ 0.006215, 0.005937, 0.006096, 0.005390, 0.005717, 0.005780, 0.006070, 0.036398 ]
                                              , [ 0.006366, 0.005806, 0.006637, 0.005780, 0.006117, 0.005834, 0.005890, 0.034395 ]
                                              , [ 0.006626, 0.006453, 0.006281, 0.006433, 0.006117, 0.006635, 0.006568, 0.031342 ]
                                              , [ 0.050068, 0.051382, 0.050736, 0.051315, 0.049790, 0.053596, 0.052917, 0.332608 ]
                                             ]
                                    , cbkg = [  [ 0.003327, 0.003606, 0.004014, 0.004860, 0.005355, 0.005540, 0.004862, 0.037073 ]
                                              , [ 0.003360, 0.004326, 0.004820, 0.006033, 0.006819, 0.007343, 0.005648, 0.046975 ]
                                              , [ 0.003701, 0.005168, 0.005987, 0.007426, 0.008269, 0.007753, 0.007018, 0.057499 ]
                                              , [ 0.004464, 0.006245, 0.007493, 0.009373, 0.009916, 0.010178, 0.008936, 0.075733 ]
                                              , [ 0.021447, 0.029485, 0.033599, 0.042014, 0.048202, 0.049706, 0.044735, 0.341693 ]
                                             ]
                                   )
                 , OSSSTagEq = dict(  sig  = (0.500531, 0.499469)
                                    , cbkg = (0.499288, 0.500712)
                                    , etaTag = dict(  sig  = [#  [  (0.518531, 0.481469), (0.511729, 0.488271), (0.532317, 0.467683)
                                                              #   , (0.501278, 0.498722), (0.537036, 0.462964), (0.505679, 0.494321)
                                                              #   , (0.503547, 0.496453), (0.505935, 0.494065)
                                                              #  ]
                                                              #, [  (0.501955, 0.498045), (0.533521, 0.466479), (0.565687, 0.434313)
                                                              #   , (0.504880, 0.495120), (0.502906, 0.497094), (0.501393, 0.498607)
                                                              #   , (0.508728, 0.491272), (0.498901, 0.501099)
                                                              #  ]
                                                              #, [  (0.500270, 0.499730), (0.516512, 0.483488), (0.495842, 0.504158)
                                                              #   , (0.467002, 0.532998), (0.501260, 0.498740), (0.502659, 0.497341)
                                                              #   , (0.484945, 0.515055), (0.498325, 0.501675)
                                                              #  ]
                                                              #, [  (0.496065, 0.503935), (0.530474, 0.469526), (0.470131, 0.529869)
                                                              #   , (0.496751, 0.503249), (0.499674, 0.500326), (0.489182, 0.510818)
                                                              #   , (0.505240, 0.494760), (0.498645, 0.501355)
                                                              #  ]
                                                              #, [  (0.477089, 0.522911), (0.507735, 0.492265), (0.493996, 0.506004)
                                                              #   , (0.498700, 0.501300), (0.506971, 0.493029), (0.506913, 0.493087)
                                                              #   , (0.498733, 0.501267), (0.499510, 0.500490)
                                                              #  ]
                                                             ]
                                                    , cbkg = [  [  (0.432368, 0.567632), (0.450853, 0.549147), (0.502949, 0.497051)
                                                                 , (0.477038, 0.522962), (0.488115, 0.511885), (0.460788, 0.539212)
                                                                 , (0.542491, 0.457509), (0.497941, 0.502059)
                                                                ]
                                                              , [  (0.431933, 0.568067), (0.448352, 0.551648), (0.462688, 0.537312)
                                                                 , (0.520658, 0.479342), (0.473874, 0.526126), (0.514179, 0.485821)
                                                                 , (0.523592, 0.476408), (0.499034, 0.500966)
                                                                ]
                                                              , [  (0.443730, 0.556270), (0.493874, 0.506126), (0.496520, 0.503480)
                                                                 , (0.478078, 0.521922), (0.485809, 0.514191), (0.513988, 0.486012)
                                                                 , (0.479254, 0.520746), (0.501633, 0.498367)
                                                                ]
                                                              , [  (0.519854, 0.480146), (0.466892, 0.533108), (0.512752, 0.487248)
                                                                 , (0.503124, 0.496876), (0.481143, 0.518857), (0.491513, 0.508487)
                                                                 , (0.470787, 0.529213), (0.507020, 0.492980)
                                                                ]
                                                              , [  (0.517407, 0.482593), (0.514572, 0.485428), (0.478907, 0.521093)
                                                                 , (0.502737, 0.497263), (0.498347, 0.501653), (0.495888, 0.504112)
                                                                 , (0.509286, 0.490714), (0.501677, 0.498323)
                                                                ]
                                                             ]
                                                   )
                                   )
                )

def findBin( obs, val ) :
    bin = 0
    for bound in condBins[obs] :
        if val < bound : break
        bin += 1
    return bin

def condObsDist( obsNames, obsVals, component ) :
    if obsNames[0] == 'OSSSTagEq' and obsNames[1] == 'etaTag' and component == 'sig' :
        sign = 1. if obsVals[0][0] == obsVals[0][1] else -1.
        return 0.5 * ( 1. + sign * ( 1. - 2. * obsVals[1][0] ) * ( 1. - 2. * obsVals[1][1] ) )

    dists = condDists[ obsNames[0] ]
    obsList = [ ( obsNames[0], obsVals[0] ) ]
    for name, val in zip( obsNames[ 1 : ], obsVals[ 1 : ] ) :
        dists = dists[name]
        obsList.append( ( name, val ) )

    dists = dists[component]
    for obs in reversed(obsList) :
        bin = ( findBin( 'etaOS', obs[1][0] ), findBin( 'etaSS', obs[1][1] ) ) if obs[0] == 'etaTag'\
              else findBin( obs[0], obs[1] ) if obs[0] in [ 'KKMass', 'sigmat' ]\
              else int( obs[1][0] != obs[1][1] ) if obs[0] == 'OSSSTagEq'\
              else obs[1]
        dists = dists[bin[0]][bin[1]] if type(bin) == tuple else dists[bin]

    return dists
