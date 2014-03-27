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
