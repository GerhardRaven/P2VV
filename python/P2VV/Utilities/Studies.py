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

    def setParSetID( self, id ) :
        if self._nPars < 1 : return
        self._parSetID = id
        self._randGen.SetSeed(id)

    def parSetID(self) :
        if self._nPars < 1 : return None
        return self._parSetID

    def readParameters( self, file ) :
        self._filePars = None
        if file :
            from P2VV.Parameterizations.FullPDFs import PdfConfiguration
            self._filePars = PdfConfiguration()
            self._filePars.readParametersFromFile( filePath = file )

    def writeParameters( self, file ) :
        if self._filePars == None :
            print 'P2VV - ERROR: parameterGen.writeParameters(): no parameters to write to file'
            return
        self._filePars.writeParametersToFile( filePath = file, Verbose = False )

    def generateParSet(self) :
        if self._nPars < 1 : return [ ]

        while True :
            # generate values for uncorrelated parameters
            for parIt, ( val, err ) in enumerate( zip( self._diagVals, self._diagErrs ) ) :
                self._genDiagVals[parIt] = self._randGen.Gaus( val, err )

            # transform generated values to original correlated-parameter values
            genVals = self._invTrans * self._genDiagVals

            # get generated values and check if they are in range
            genValsRet = [ ]
            for parIt in range(self._nPars) :
                val = genVals[parIt]
                if self._ranges and ( val < self._ranges[parIt][0] or val > self._ranges[parIt][1] ) : break
                genValsRet.append(val)
            if len(genValsRet) == self._nPars : break
            self._nOutRange += 1

        # update generation statistics
        self._nGen += 1
        self._parSetID += 1
        for parIt, val in enumerate(genValsRet) :
            self._sumVals[parIt]   += val
            self._sumSqVals[parIt] += val**2

        if self._filePars != None :
            # set parameters for output file
            pars = self._filePars.parameters()
            parsFound = True
            for name, val in zip( self._parNames, genValsRet ) :
                if name in pars :
                    pars[name] = tuple( [ val ] + [ v for v in pars[name][ 1 : -1 ] ] + [ False ] )
                else :
                    parsFound = False
            if not parsFound and self._parsFoundWarn :
                print 'P2VV - WARNING: parameterGen.generateParSet(): not all generated parameters found in parameter set for output file'
                self._parsFoundWarn = False

        return genValsRet

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
        self._anaParVals = [ ]
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
        self._anaParVals = [ ]
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

                self._anaParVals.append( [ ] )
                parVals = self._anaPars.parameters()
                for name in self._parNames :
                    assert name in parVals\
                           , 'P2VV - ERROR: fitResultsAnalysis.readAnaFiles(): parameter "%s" not found in file "%s"' % ( name, file )
                    self._anaParVals[-1].append( parVals[name][0] )

            self._nFiles = sum( stats for stats in fileStats.itervalues() )
            print 'P2VV - INFO: fitResultsAnalysis.readAnaFiles(): read %d parameter files for analysis (fit status: %s)'\
                  % ( self._nFiles, ', '.join( '%d: %d' % ( stat, count ) for stat, count in fileStats.iteritems() ) )

    def processResults( self, histsFile = '' ) :
        if not self._anaParVals :
            print 'P2VV - WARNING: fitResultsAnalysis.processResults(): no parameter values available for analysis'
            return

        self._parSums = [ [ 0., 0., 0., 0. ] for par in self._parNames ]
        self._nParVals = [ [ 0., 0., 0. ] for par in self._parNames ]
        for vals in self._anaParVals :
            for valIt, val in enumerate(vals) :
                self._parSums[valIt][0] += val
                self._parSums[valIt][1] += val**2
                self._nParVals[valIt][0] += 1
                if self._refParVals[valIt] != None :
                    valDiff = val - self._refParVals[valIt][0]
                    self._parSums[valIt][ 3 if valDiff < 0. else 2 ] += valDiff**2
                    self._nParVals[valIt][ 2 if valDiff < 0. else 1 ] += 1
        for valCounts in self._nParVals :
            assert valCounts[0] == self._nFiles and ( valCounts[1] + valCounts[2] == self._nFiles or self._refParVals[0] == None )

        from math import sqrt, log10, ceil
        nameLen = min( 30, max( len(name) for name in self._parNames ) )
        sepStr = '  ' + '-' * ( nameLen + ( 97 if self._refParVals[0] != None else 31 ) )
        print 'P2VV - INFO: fitResultsAnalysis.processResults(): parameter statistics for %d files:' % self._nFiles
        print sepStr
        print ( '  {0:<%ds}   {1:<8s}   {2:<11s}' % nameLen ).format( 'name', 'mean', 'std. dev.' ),
        if self._refParVals[0] != None :
            print '   {0:<8s}   {1:<7s}   {2:<11s}   {3:<10s}   {4:<11s}   {5:<11s}'\
                  .format( 'value', 'uncert.', 'mean dev.', 'dev. (sig)', '+dev. (sig)', '-dev. (sig)' )
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
                print ( '   {0:<+8.%df}   {1:<7.%df}   {2:<+11.%df}   {3:<10.3f}   {4:<11.3f}   {5:<11.3f}' % ( prec, prec, precDev ) )\
                      .format( refVal[0], refVal[1], meanVal - refVal[0], dev, devPlus, devMin )
            else :
                print
        print sepStr
