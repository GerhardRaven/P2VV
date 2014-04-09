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
                    vals.append( parVals[name][0] )

            self._nFiles = sum( stats for stats in fileStats.itervalues() )
            print 'P2VV - INFO: fitResultsAnalysis.readAnaFiles(): read %d parameter files for analysis (fit status: %s)'\
                  % ( self._nFiles, ', '.join( '%d: %d' % ( stat, count ) for stat, count in fileStats.iteritems() ) )

    def processResults( self, histsFile = '' ) :
        if not self._anaParVals :
            print 'P2VV - WARNING: fitResultsAnalysis.processResults(): no parameter values available for analysis'
            return

        self._parSums = [ [ sum( self._anaParVals[name] ), 0., 0., 0. ] for name in self._parNames ]
        self._nParVals = [ [ len( self._anaParVals[name] ), 0., 0. ] for name in self._parNames ]
        for parIt, name in enumerate(self._parNames) :
            if name in self._parHists and self._parHists[name] :
                self._parHists[name].Delete()
            if histsFile :
                histBins = int( float(self._nFiles) / 100. ) if self._nFiles > 1000 else 10
                histMin = min( self._anaParVals[name] )
                histMax = max( self._anaParVals[name] )
                histRange = histMax - histMin
                if histRange > 0. :
                    histMin = histMin - 0.01 * histRange
                    histMax = histMax + 0.01 * histRange
                from P2VV.Load import LHCbStyle
                from ROOT import TH1D, kFullDotLarge
                self._parHists[name] = TH1D( name, name, histBins, histMin, histMax )
                self._parHists[name].SetXTitle(name)
                self._parHists[name].SetYTitle('Entries / %.2g' % self._parHists[name].GetBinWidth(1) )
                self._parHists[name].SetMarkerStyle(kFullDotLarge)
                self._parHists[name].SetMarkerSize(0.6)
            for val in self._anaParVals[name] :
                self._parSums[parIt][1] += val**2
                if self._refParVals[parIt] != None :
                    valDiff = val - self._refParVals[parIt][0]
                    self._parSums[parIt][ 3 if valDiff < 0. else 2 ] += valDiff**2
                    self._nParVals[parIt][ 2 if valDiff < 0. else 1 ] += 1
                if histsFile :
                    self._parHists[name].Fill(val)
        for valCounts in self._nParVals :
            assert valCounts[0] == self._nFiles and ( valCounts[1] + valCounts[2] == self._nFiles or self._refParVals[0] == None )

        from math import sqrt, log10, ceil
        nameLen = min( 30, max( len(name) for name in self._parNames ) )
        sepStr = '  ' + '-' * ( nameLen + ( 105 if self._refParVals[0] != None else 31 ) )
        print 'P2VV - INFO: fitResultsAnalysis.processResults(): parameter statistics for %d files:' % self._nFiles
        print sepStr
        print ( '  {0:<%ds}   {1:<8s}   {2:<11s}' % nameLen ).format( 'name', 'mean', 'std dev' ),
        if self._refParVals[0] != None :
            print '   {0:<8s}   {1:<7s}   {2:<11s}   {3:<19s}   {4:<8s}   {5:<8s}'\
                  .format( 'value', 'uncert', 'mean dev', 'dev rel / abs', '+dev rel', '-dev rel' )
        else :
            print
        print sepStr
        for name, parSums, nParVals, refVal in zip( self._parNames, self._parSums, self._nParVals, self._refParVals ) :
            meanVal = parSums[0] / float( nParVals[0] )
            meanSqVal = parSums[1] / float( nParVals[0] )
            stdDev = sqrt( meanSqVal - meanVal**2 )
            precDev = max( 0, 3 - int( ceil( log10(stdDev) ) ) )
            prec = max( 0, 2 - int( ceil( log10( refVal[1] ) ) ) ) if refVal != None else precDev
            print ( '  {0:<%ds}   {1:<+8.%df}   {2:<11.%df}' % ( nameLen, prec, precDev ) ).format( name, meanVal, stdDev ),
            if refVal != None :
                dev = sqrt( ( parSums[2] + parSums[3] ) / float( nParVals[0] ) ) / refVal[1]
                devPlus = sqrt( parSums[2] / float( nParVals[1] ) ) / refVal[1]
                devMin  = sqrt( parSums[3] / float( nParVals[2] ) ) / refVal[1]
                print ( '   {0:<+8.%df}   {1:<7.%df}   {2:<+11.%df}   {3:<5.3f} / {4:<11.%df}   {5:<8.3f}   {6:<8.3f}'\
                        % ( prec, prec, precDev, precDev ) )\
                        .format( refVal[0], refVal[1], meanVal - refVal[0], dev, dev * refVal[1], devPlus, devMin )
            else :
                print
        print sepStr

        if histsFile :
            from ROOT import TCanvas
            dCanv = TCanvas('dummy')
            dCanv.Print( histsFile + '[' )
            for name in self._parNames :
                canv = TCanvas(name)
                self._parHists[name].Draw('E1')
                canv.Print(histsFile)
            dCanv.Print( histsFile + ']' )
