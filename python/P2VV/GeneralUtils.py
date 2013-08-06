###########################################################################################################################################
## P2VVGeneralUtils: General P2VV utilities                                                                                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

# clever switch construct from http://code.activestate.com/recipes/410692/
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args:
            self.fall = True
            return True
        else:
            return False



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
    noNAN = ( ' && '.join( '( %s==%s )' % ( obs, obs ) for obs in observables ) ) if hasattr( observables, '__iter__' ) else ''
    cuts = kwargs.pop( 'cuts', '' )
    tmp_file = None
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
          import tempfile
          import os
          from ROOT import TFile, gFile
          orig_file = gFile
          d = None
          if 'TMP' in os.environ:
              d = os.environ['TMP']
          elif 'TMPDIR' in os.environ:
              d = os.environ['TMPDIR']
          elif os.access(os.path.dirname(filePath), os.W_OK):
              d = os.path.dirname(filePath)
          else:
              d = '/tmp'
          fd, temp_name = tempfile.mkstemp(suffix = '.root', dir = d)
          os.close(fd)
          os.remove(temp_name)
          tmp_file = TFile.Open(temp_name, 'recreate')
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
      if cuts : print 'P2VV - INFO: readData: applying cuts: %s' % cuts

      # loop over category states
      states = tuple( [ [ ( cat[0], ind ) for ind in cat[1] ] for cat in kwargs.pop( 'Categories', [ ( '', [ '' ] ) ] ) ] )
      from itertools import product
      for it, state in enumerate( product(*states) ) :
          # get data set
          dsName = '_'.join( str(catSt[0]) + str(catSt[1]) for catSt in state )
          dsName = dataSetName + ( ( '_' + dsName ) if dsName else '' )
          dataSet = file.Get(dsName)
          assert dataSet, 'P2VV - ERROR: data set "%s" not found' % dsName
          if it == 0 :
              if observables :
                  from ROOT import RooDataSet
                  data = RooDataSet( dataSetName, dataSetName
                                   , [ obs._var for obs in observables ]
                                   , Import = dataSet
                                   , Cut = noNAN + ' && ' + cuts if cuts else noNAN
                                   )
              else :
                  data = dataSet

          else :
              data.append(dataSet)

      file.Close()

    print 'P2VV - INFO: read dataset with %s entries (%.1f weighted)' % ( data.numEntries(), data.sumEntries() )

    # import data set into current workspace
    from P2VV.RooFitWrappers import RooObject
    wsData = RooObject().ws().put( data, **kwargs )
    data.IsA().Destructor(data)
    if tmp_file:
        tmp_file.Close()
        os.remove(tmp_file.GetName())
        if orig_file: orig_file.cd()
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
        corrFactors = (  sumWeights[0] / sumSqWeights[0]
                       , [ sum / sumSq for sum, sumSq in zip( sumWeights[ 1 : ], sumSqWeights[ 1 : ] ) ]
                      )

    # add corrected weights to data set
    from P2VV.Load import P2VVLibrary
    from ROOT import RooCorrectedSWeight
    if splitCatName :
        from ROOT import std
        corrFactorsVec = std.vector('Double_t')()
        print 'P2VV - INFO: correctSWeights: multiplying sWeights (-ln(L)) to correct for background dilution with factors (overall factor %.4f):'\
              % corrFactors[0]
        for iter, fac in enumerate( corrFactors[1] ) :
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

    importIntoWS = kwargs.pop( 'ImportIntoWS', True )
    if importIntoWS :
        # import data set into current workspace
        from P2VV.RooFitWrappers import RooObject
        return RooObject().ws().put( dataSet, **kwargs )
    else :
        return dataSet


def addTaggingObservables( dataSet, iTagName, tagCatName, tagDecisionName, estimWTagName, tagCatBins ) :
    """add tagging observables to data set
    """

    # get observables from data set
    obsSet      = dataSet.get(0)
    tagDecision = obsSet.find(tagDecisionName)
    estimWTag   = obsSet.find(estimWTagName)

    # create initial state tag
    from P2VV.Load import P2VVLibrary
    from ROOT import RooTagDecisionWrapper
    iTagWrapper = RooTagDecisionWrapper(iTagName, 'Flavour tag', tagDecision)

    # create tagging category
    from ROOT import RooThresholdCategory
    tagCatFormula = RooThresholdCategory( tagCatName, 'P2VV tagging category', estimWTag, tagCatBins[0][0], tagCatBins[0][1] )
    for cat in range( 1, len(tagCatBins) ) : tagCatFormula.addThreshold( tagCatBins[cat][2], tagCatBins[cat][0], tagCatBins[cat][1] )

    # add tagging category binning to estimated wrong-tag probability variable
    from array import array
    from ROOT import RooBinning
    binBounds = array( 'd', [ 0. ] + [ tagCatBins[ len(tagCatBins) - it ][2] for it in range( 1, len(tagCatBins) + 1 ) ] )
    tagCatBinning = RooBinning( len(binBounds) - 1, binBounds, 'tagCats' )
    estimWTag.setBinning( tagCatBinning, 'tagCats' )

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
        assert ( obsSet.getCatIndex(tagDecisionName) == 0 and obsSet.getRealValue(estimWTagName) >= binBounds[-2] )\
                or ( obsSet.getCatIndex(tagDecisionName) != 0 and obsSet.getRealValue(estimWTagName) < binBounds[-2] ),\
                'P2VV - ERROR: addTaggingObservables: tag decision = %+d, while estimated wrong-tag probability = %.10f (threshold = %.10f)'\
                % ( obsSet.getCatIndex(tagDecisionName), obsSet.getRealValue(estimWTagName), binBounds[-2] )
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


def printEventYields( **kwargs ) :
    # get arguments
    parSet     = kwargs.pop( 'ParameterSet',        [ ] )
    yieldNames = kwargs.pop( 'YieldNames',          [ ] )
    splitCats  = kwargs.pop( 'SplittingCategories', [ ] )

    assert parSet,     'P2VV - ERROR: printEventYields: no parameter set with yield variables found in arguments ("ParameterSet")'
    assert yieldNames, 'P2VV - ERROR: printEventYields: no yield names found in arguments ("YieldNames")'

    if splitCats : splitCats = list( set( cat for cat in splitCats ) )

    # variables for looping over splitting category states
    iters = { }
    inds  = { }
    labs  = { }
    for cat in splitCats :
        iters[cat] = 0
        inds[cat]  = [ ]
        labs[cat]  = [ ]
        catIter = cat.typeIterator()
        catState = catIter.Next()
        while catState :
            inds[cat].append( catState.getVal() )
            labs[cat].append( catState.GetName() )
            catState = catIter.Next()

    # print yields and fractions with error from S/(S+B) fraction only (no Poisson error for number of events!)
    print
    print '-'.join( dashes for dashes in [ ' ' * 4 + '-' * 22, '-' * 8 ] + [ '-' * num for name in yieldNames for num in ( 23, 19 ) ] )
    print ' '.join( ( '{0:^%ds}' % width ).format(title) for ( title, width ) in [ ( '', 26 ), ( ' total', 8 ) ]\
                   + [ ( name if num == 23 else 'f(' + name + ')', num ) for name in yieldNames for num in ( 23, 19 ) ] )
    print '-'.join( dashes for dashes in [ ' ' * 4 + '-' * 22, '-' * 8 ] + [ '-' * num for name in yieldNames for num in ( 23, 19 ) ] )

    cont = True
    while cont :
        stateName = ';'.join( labs[cat][ iters[cat] ] for cat in splitCats )
        yields = [ getSplitPar( name, ( '{%s}' % stateName ) if stateName else '', parSet ) for name in yieldNames ]

        from math import sqrt
        nEv        = [ yieldVar.getVal()   for yieldVar in yields ]
        nEvErr     = [ yieldVar.getError() for yieldVar in yields ]
        nEvTot     = sum(nEv)
        frac       = [ num / nEvTot if nEvTot > 0. else 0. for num in nEv ]
        nEvCorrErr = [ sqrt( numErr**2 - num**2 / nEvTot ) for num, numErr in zip( nEv, nEvErr ) ]
        fracErr    = [ err / nEvTot if nEvTot > 0. else 0. for err in nEvCorrErr ]

        print '     {0:>20s}   {1:>6.0f}  '.format( stateName, nEvTot )\
              + ' '.join( ' {0:>9.2f} +/- {1:>7.2f}   {2:>6.4f} +/- {3:>6.4f} '.format( num, err, fr, frErr )\
                         for ( num, err, fr, frErr ) in zip( nEv, nEvCorrErr, frac, fracErr ) )

        if not splitCats : break

        iters[ splitCats[-1] ] += 1
        for catIt in range( len(splitCats) ) :
            if iters[ splitCats[ -catIt - 1 ] ] >= splitCats[ -catIt - 1 ].numTypes() :
                if catIt == len(splitCats) - 1 :
                    cont = False
                else :
                    iters[ splitCats[ -catIt - 1 ] ] = 0
                    iters[ splitCats[ -catIt - 2 ] ] +=1
            else :
                continue

    print '-'.join( dashes for dashes in [ ' ' * 4 + '-' * 22, '-' * 8 ] + [ '-' * num for name in yieldNames for num in ( 23, 19 ) ] )
    print


def printEventYieldsData( **kwargs ) :
    # get arguments
    fullDataSet    = kwargs.pop( 'FullDataSet',         None )
    weightDataSets = kwargs.pop( 'WeightedDataSets',    [ ]  )
    dataSetNames   = kwargs.pop( 'DataSetNames',        [ ]  )
    splitCats      = kwargs.pop( 'SplittingCategories', [ ]  )

    assert fullDataSet,    'P2VV - ERROR: printEventYieldsData: no data set found in arguments ("FullDataSet")'
    assert weightDataSets, 'P2VV - ERROR: printEventYieldsData: no weighted data sets found in arguments ("WeightedDataSets")'
    if not dataSetNames : dataSetNames = [ 'data set %d' % it for it, dataSet in enumerate(weightDataSets) ]

    if splitCats : splitCats = list( set( cat for cat in splitCats ) )

    # print total numbers of events
    from math import sqrt
    nEvTot = fullDataSet.sumEntries()
    nEv    = [ dataSet.sumEntries() for dataSet in weightDataSets ]
    frac   = [ num / nEvTot       if nEvTot > 0. else 0. for num in nEv ]
    signif = [ num / sqrt(nEvTot) if nEvTot > 0. else 0. for num in nEv ]

    print
    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 31 ] + [ '-' * 30 for dataSet in weightDataSets ] )
    print ' ' * 35 + '|' + '|'.join( ' {0:^28} '.format(name) for name in dataSetNames )
    print ' ' * 27 + '  total |' + '|'.join( ' {0:^9s}   {1:^6s}   {2:^9s} '\
          .format( 'N_%d' % it, 'f_%d' % it, u'N_%d/\u221AN'.encode('utf-8') % it ) for it, dataSet in enumerate(weightDataSets) )
    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 31 ] + [ '-' * 30 for dataSet in weightDataSets ] )
    print ' ' *  4 + ' {0:>20s}   {1:>6.0f} |'.format( 'total', nEvTot ) + '|'.join( ' {0:>9.2f}   {1:>6.4f}   {2:>7.3f} '\
          .format( num, fr, sig ) for ( num, fr, sig ) in zip( nEv, frac, signif ) )
    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 31 ] + [ '-' * 30 for dataSet in weightDataSets ] )

    if not splitCats :
        print
        return

    # print numbers of events per splitting category
    iters = { }
    inds  = { }
    labs  = { }
    for cat in splitCats :
        iters[cat] = 0
        inds[cat]  = [ ]
        labs[cat]  = [ ]
        catIter = cat.typeIterator()
        catState = catIter.Next()
        while catState :
            inds[cat].append( catState.getVal() )
            labs[cat].append( catState.GetName() )

            cut    = '!(%s-%d)' % ( cat.GetName(), inds[cat][-1] )
            nEvTot = fullDataSet.sumEntries(cut)
            nEv    = [ dataSet.sumEntries(cut) for dataSet in weightDataSets ]
            frac   = [ num / nEvTot       if nEvTot > 0. else 0. for num in nEv ]
            signif = [ num / sqrt(nEvTot) if nEvTot > 0. else 0. for num in nEv ]
            print ' ' *  4 + ' {0:>20s}   {1:>6.0f} |'.format( labs[cat][-1], nEvTot ) + '|'.join( ' {0:>9.2f}   {1:>6.4f}   {2:>7.3f} '\
                  .format( num, fr, sig ) for ( num, fr, sig ) in zip( nEv, frac, signif ) )

            catState = catIter.Next()

        print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 31 ] + [ '-' * 30 for dataSet in weightDataSets ] )

    if len(splitCats) < 2 :
        print
        return

    # print numbers of events for each combination of splitting categories
    cont = True
    while cont :
        stateName = ';'.join( labs[cat][ iters[cat] ] for cat in splitCats )
        cut       = '&&'.join( '!(%s-%d)' % ( cat.GetName(), inds[cat][ iters[cat] ] ) for cat in splitCats )
        nEvTot    = fullDataSet.sumEntries(cut)
        nEv       = [ dataSet.sumEntries(cut) for dataSet in weightDataSets ]
        frac      = [ num / nEvTot       if nEvTot > 0. else 0. for num in nEv ]
        signif    = [ num / sqrt(nEvTot) if nEvTot > 0. else 0. for num in nEv ]
        print ' ' *  4 + ' {0:>20s}   {1:>6.0f} |'.format( stateName, nEvTot ) + '|'.join( ' {0:>9.2f}   {1:>6.4f}   {2:>7.3f} '\
              .format( num, fr, sig ) for ( num, fr, sig ) in zip( nEv, frac, signif ) )

        iters[ splitCats[-1] ] += 1
        for catIt in range( len(splitCats) ) :
            if iters[ splitCats[ -catIt - 1 ] ] >= splitCats[ -catIt - 1 ].numTypes() :
                if catIt == len(splitCats) - 1 :
                    cont = False
                else :
                    iters[ splitCats[ -catIt - 1 ] ] = 0
                    iters[ splitCats[ -catIt - 2 ] ] +=1
            else :
                continue

    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 31 ] + [ '-' * 30 for dataSet in weightDataSets ] )
    print


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

    frames = []
    frames.append(obsFrame)

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
        if data and 'Asymmetry' not in pdfOpts : obsFrame.drawAfter( 'data', 'pdf' )

    # plot additional PDFs
    if addPDFs :
        for num, addPDF in enumerate(addPDFs) :
            addPDF.plotOn( obsFrame, Name = 'addPDF' + str(num), **(addPDFsOpts[num]) )
            if data and 'Asymmetry' not in addPDFsOpts[num] : obsFrame.drawAfter( 'data', 'addPDF' + str(num) )

    #TODO: add chisq/nbins
    #chisq = obsFrame.chiSquare( 'pdf', 'data' )
    #nbins = obsFrame.GetNbinsX()

    # set y scale
    if yScale[0]    != None : obsFrame.SetMinimum(yScale[0])
    if yScale[1]    != None : obsFrame.SetMaximum(yScale[1])
    if yScaleRel[0] != None : obsFrame.SetMinimum( yScaleRel[0] * obsFrame.GetMinimum() )
    if yScaleRel[1] != None : obsFrame.SetMaximum( yScaleRel[1] * obsFrame.GetMaximum() )
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
        #yAxis.SetTitleSize(0.10)
        #yAxis.SetLabelSize(0.08)
        yAxis.SetTitleOffset( 0.7 * yAxis.GetTitleOffset() )

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
        frames.append(residFrame)
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
        residFrame.addPlotable( residHist, 'P' if not type(plotResidHist) == str else plotResidHist )
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
    return frames


# function for plotting a PDF in KK mass bins
def getCondObsPlotsInKKbins(pdf, data, canv):
    from ROOT import TCanvas, gPad, TH1F
    from P2VV.RooFitWrappers import Category

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
        from P2VV.RooFitWrappers import SimultaneousPdf
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
        from P2VV.RooFitWrappers import ConstVar, Customizer
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

        from P2VV.RooFitWrappers import Category
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


# function for plotting the S-wave parameters versus the (binned) KK mass
def plotSWaveBins( **kwargs ) :
    mode = kwargs.pop( 'Mode', 'phases' )
    if mode not in [ 'phases', 'fractions', 'events' ] :
        raise KeyError, 'P2VV - ERROR: plotSWaveBins: possible plot modes: "phases", "fractions", "events"'
    if any( key not in kwargs for key in [ 'SValues', 'SLowErrors', 'SHighErrors' ] ) :
        raise KeyError, 'P2VV - ERROR: plotSWaveBins: "SValues", "SLowErrors" and "SHighErrors" arguments are required'

    defRange  = ( None, None ) if mode == 'phases' else ( 0., None )
    defSLabel = '#delta_{S} - #delta_{#perp}    [rad]' if mode == 'phases' else\
                'F_{S}' if mode == 'fractions' else\
                'N_{S} / (MeV/c^{2})'
    yAxisRange  = kwargs.pop( 'SAxisRange',    defRange                      )
    KKMassLabel = kwargs.pop( 'KKMassLabel',   'm(K^{+}K^{-}) [MeV/c^{2}]'   )
    SLabel      = kwargs.pop( 'SLabel',        defSLabel                     )
    plotTitle   = kwargs.pop( 'PlotTitle',     ''                            )
    LHCbText1   = kwargs.pop( 'LHCbTextLine1', 'LHCb'                        )
    LHCbText2   = kwargs.pop( 'LHCbTextLine2', ''                            )
    drawLegend  = kwargs.pop( 'DrawLegend',    False                         )
    massBins    = kwargs.pop( 'MassBins',      [ 988., 1008., 1032., 1050. ] )
    theoryVals  = kwargs.pop( 'TheoryValues',  None                          )
    SVals       = kwargs.pop( 'SValues'                                      )
    SLowErrs    = kwargs.pop( 'SLowErrors'                                   )
    SHighErrs   = kwargs.pop( 'SHighErrors'                                  )
    gray        = kwargs.pop( 'GrayScale',     False                         )

    if kwargs :
        raise KeyError, 'P2VV - ERROR: plotSWaveBins: unexpected keyword arguments: %s' % kwargs

    from array import array
    KKMass         = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMassErr      = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] )                for it in range( len(massBins) - 1 ) ] )

    offs = 0.35 if mode == 'phases' else 0.
    KKMass1        = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + offs + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMass1LowErr  = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + offs                for it in range( len(massBins) - 1 ) ] )
    KKMass1HighErr = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - offs                for it in range( len(massBins) - 1 ) ] )

    KKMass2        = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - offs + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMass2LowErr  = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - offs                for it in range( len(massBins) - 1 ) ] )
    KKMass2HighErr = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + offs                for it in range( len(massBins) - 1 ) ] )

    from ROOT import TGraphAsymmErrors
    SGraphs = [ ]

    S1        = array( 'd', [ SVals[it]     / (massBins[it+1] - massBins[it]) for it in range(len(massBins)-1) ] if mode == 'events' else SVals     )
    S1LowErr  = array( 'd', [ SLowErrs[it]  / (massBins[it+1] - massBins[it]) for it in range(len(massBins)-1) ] if mode == 'events' else SLowErrs  )
    S1HighErr = array( 'd', [ SHighErrs[it] / (massBins[it+1] - massBins[it]) for it in range(len(massBins)-1) ] if mode == 'events' else SHighErrs )
    SGraphs.append( TGraphAsymmErrors( len(KKMass1), KKMass1, S1, KKMass1LowErr, KKMass1HighErr, S1LowErr, S1HighErr ) )

    if mode == 'phases' :
        from math import pi
        S2        = array( 'd', [ pi - val for val in S1 ] )
        S2LowErr  = array( 'd', S1HighErr )
        S2HighErr = array( 'd', S1LowErr  )
        SGraphs.append( TGraphAsymmErrors( len(KKMass2), KKMass2, S2, KKMass2LowErr, KKMass2HighErr, S2LowErr, S2HighErr ) )

    theory = None
    if theoryVals :
        from ROOT import TGraphErrors
        theory    = array( 'd', theoryVals           )
        theoryErr = array( 'd', [ 0. ] * len(KKMass) )
        SGraphs.append( TGraphErrors( len(KKMass), KKMass, theory, KKMassErr, theoryErr ) )

    SMin0 = min( val for val in theory ) if theory else +1.e32
    SMax0 = max( val for val in theory ) if theory else -1.e32
    SMin1 = min( val - err for val, err in zip( S1, S1LowErr  ) )
    SMax1 = max( val + err for val, err in zip( S1, S1HighErr ) )
    SMin  = min( [ SMin0, SMin1, ( pi - SMax1 ) if mode == 'phases' else +1.e32 ] )
    SMax  = max( [ SMax0, SMax1, ( pi - SMin1 ) if mode == 'phases' else -1.e32 ] )

    from ROOT import gStyle, kSolid
    gStyle.SetEndErrorSize(3)
    gStyle.SetLineStyleString( 11, ' 30 15' )
    SGraphs[0].SetLineStyle(kSolid)
    if len(SGraphs) > 1 : SGraphs[1].SetLineStyle(kSolid)
    if len(SGraphs) > 2 : SGraphs[2].SetLineStyle(11)

    from ROOT import kBlack, kBlue, kRed
    SGraphs[0].SetLineColor( kBlue + 2 )
    if len(SGraphs) > 1 : SGraphs[1].SetLineColor( kRed - 4 )
    if len(SGraphs) > 2 : SGraphs[2].SetLineColor(kBlack)

    SGraphs[0].SetMarkerColor( kBlue + 2 )
    if len(SGraphs) > 1 : SGraphs[1].SetMarkerColor( kRed - 4 )
    if len(SGraphs) > 2 : SGraphs[2].SetMarkerColor(kBlack)

    SGraphs[0].SetLineWidth(3)
    if len(SGraphs) > 1 : SGraphs[1].SetLineWidth(3)
    if len(SGraphs) > 2 : SGraphs[2].SetLineWidth(3)

    from ROOT import kFullCircle, kFullTriangleDown
    SGraphs[0].SetMarkerStyle(kFullCircle)
    if len(SGraphs) > 1 : SGraphs[1].SetMarkerStyle(kFullTriangleDown)
    if len(SGraphs) > 2 : SGraphs[2].SetMarkerStyle(kFullCircle)
    SGraphs[0].SetMarkerSize(1.2)
    if len(SGraphs) > 1 : SGraphs[1].SetMarkerSize(1.4)
    if len(SGraphs) > 2 : SGraphs[2].SetMarkerSize(0.7)

    for graph in SGraphs :
        graph.SetMinimum( yAxisRange[0] if yAxisRange[0] != None else SMin - 0.10 * ( SMax - SMin ) )
        graph.SetMaximum( yAxisRange[1] if yAxisRange[1] != None else SMax + 0.15 * ( SMax - SMin ) )

        graph.GetXaxis().SetTitle(KKMassLabel)
        graph.GetYaxis().SetTitle(SLabel)

        graph.GetXaxis().SetTitleOffset(1.1)
        graph.GetYaxis().SetTitleOffset( 0.8 if mode == 'phases' else 1.1 if mode == 'fractions' else 0.9 )

        graph.SetTitle(plotTitle)

        _P2VVPlotStash.append(graph)

    if mode == 'phases' and drawLegend :
        from ROOT import gStyle, TLegend
        leg = TLegend( 0.65, 0.46, 0.91, 0.66 )
        leg.SetTextFont( gStyle.GetTextFont() )
        leg.SetTextSize(0.072)
        leg.SetMargin(0.45)
        leg.AddEntry( SGraphs[0], '#Delta#Gamma_{s} > 0', 'LPE' )
        leg.AddEntry( SGraphs[1], '#Delta#Gamma_{s} < 0', 'LPE' )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        _P2VVPlotStash.append(leg)
    else :
        leg = None

    if LHCbText1 or LHCbText2 :
        from ROOT import TPaveText
        LHCbText = TPaveText( 0.24, 0.73 if LHCbText1 and LHCbText2 else 0.81, 0.70 if LHCbText1 and LHCbText2 else 0.37, 0.89, 'BRNDC' )
        if LHCbText1 : LHCbText.AddText(LHCbText1)
        if LHCbText2 : LHCbText.AddText(LHCbText2)
        LHCbText.SetShadowColor(0)
        LHCbText.SetFillStyle(0)
        LHCbText.SetBorderSize(0)
        LHCbText.SetTextAlign(12)
        LHCbText.SetTextSize(0.072)
        _P2VVPlotStash.append(LHCbText)
    else :
        LHCbText = None

    from ROOT import TCanvas
    canv = TCanvas( 'SWavePhaseCanv' if mode == 'phases' else 'SWaveFracCanv' if mode == 'fractions' else 'SWaveEventsCanv', 'S-Wave' )
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)
    if gray :
        canv.SetGrayscale()
    if len(SGraphs) > 1 :
        SGraphs[1].Draw('AP')
        SGraphs[0].Draw('P sames')
    else :
        SGraphs[0].Draw('AP')
    if len(SGraphs) > 2 :
        SGraphs[2].Draw('P sames')
    if leg :
        leg.Draw()
    if LHCbText :
        LHCbText.Draw()

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
                from P2VV.GeneralUtils import plot
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
            from P2VV.GeneralUtils import plot
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
def angularMomentIndices(label,angleFuncs) :
        from P2VV.Parameterizations.AngularFunctions import JpsiphiTransversityAngles,  JpsiphiHelicityAngles
        transAngles = { JpsiphiTransversityAngles : True, JpsiphiHelicityAngles : False  }[ type(angleFuncs) ]
        for case in switch(label):
            if case('weights') :
                return [ ( 0, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ), ( 2, 0, 0 ), ( 0, 2, 1 ), ( 0, 2, -1 ), ( 0, 2, -2 )
                       , ( 1, 0, 0 ), ( 1, 2, 1 ), ( 1, 2, -1 ) ]
            if case('basis012') :
                return [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
                             for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
            if case('basis012Plus') :
                return [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
                              for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ] + [ ( 0, 4, 0 ) ]
            if case( 'basis012Thetal' ) :
                return   [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
                               for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ] \
                       + [ ( 0, YIndex0, 0 ) for YIndex0 in range( 3, 5 ) ]
            if case('basis012ThetalPhi') :
                return  [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
                               for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ] \
                      + [ ( 0, YIndex0, YIndex1 ) for YIndex0 in range( 3, 5 ) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
            if case('basis0123') :
                return [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(4) for YIndex0 in range(4)\
                              for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
            if case('basis01234') :
                return  [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(5) for YIndex0 in range(5)\
                              for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
            if case('basisSig3') :
                return [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if not transAngles\
                              else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
            if case('basisSig4') :
                if transAngles :
                    raise RuntimeError('P2VV - ERROR: angularMomentIndices: not a valid angular efficiency configuration with transversity angles: %s'\
                                       % multiplyByAngEff)
                return [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 4, 0 ) ]
            raise RuntimeError('P2VV - ERROR: angularMomentIndices: not a valid angular efficiency configuration: %s'\
                                   % label)


def getMomentFuncArgs( funcName, kwargs ) :
    # new moments or append to existing moments?
    kwargs['addFacs'] = kwargs.pop( 'AddMoments', ( ) )

    # process all moments?
    kwargs['procAll'] = kwargs.pop( 'ProcessAll', False )

    # moments and correlations
    funcNames    = kwargs.pop( 'BasisFuncNames', [ ] )
    moments      = kwargs.pop( 'Moments',        { } )
    correlations = kwargs.pop( 'Correlations',   { } )
    if ( kwargs['addFacs'] or not kwargs['procAll'] ) and not funcNames :
        print 'P2VV - ERROR: %s(): no basis function names specified' % ( funcName if funcName else 'getMomentFuncArgs' )
        return False
    for name in funcNames :
        if moments and name not in moments :
            print 'P2VV - ERROR: %s(): could not find moment of function "%s"' % ( funcName if funcName else 'getMomentFuncArgs', name )
            return False
        if correlations and name not in correlations :
            print 'P2VV - ERROR: %s(): could not find correlations for function "%s"'\
                  % ( funcName if funcName else 'getMomentFuncArgs', name )
            return False
    kwargs.update( dict( funcNames = funcNames, moments = moments, correlations = correlations ) )

    # maximum length of basis function name
    maxName = 0
    for func in funcNames : maxName = max( len(func), maxName )
    kwargs['maxLenName'] = max( 15, min( [ kwargs.pop( 'MaxLenName', 20 ), maxName ] ) )

    # name requirements
    kwargs['names'] = kwargs.pop( 'Names', None )
    import re
    kwargs['nameExpr'] = re.compile( kwargs['names'] ) if kwargs['names'] else None

    # minimum significance
    kwargs['minSignif'] = kwargs.pop( 'MinSignificance', float('-inf') )

    # scale factors
    kwargs['scale']  = kwargs.pop( 'Scale',  1. )
    kwargs['scales'] = kwargs.pop( 'Scales', ( kwargs['scale'], kwargs['scale'], 1. ) )

    return True


def printMoments( **kwargs ) :
    # parse arguments
    assert getMomentFuncArgs( 'printMoments', kwargs ), 'P2VV - ERROR: printMoments(): unable to parse function arguments'

    # print header
    print 'P2VV - INFO: printMoments():'
    print '  name requirement: \'' + ( kwargs['names'] if kwargs['names'] else '' ) + '\''
    print '  minimum significance = %.1f' % kwargs['minSignif']
    print '  scales = ({0:.5g}, {1:.5g}, {2:.5g})'.format( kwargs['scales'][0], kwargs['scales'][1], kwargs['scales'][2] )
    print '  ' + '-' * ( 45 + kwargs['maxLenName'] )
    print ( '  {0:<%d}   {1:>11}   {2:>9}   {3:>12}' % kwargs['maxLenName'] )\
            .format( 'basis function', 'coefficient', 'std. dev.', 'significance' )
    print '  ' + '-' * ( 45 + kwargs['maxLenName'] )

    # print moments
    for func in kwargs['funcNames'] :
        if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func) )\
                or ( func in kwargs['moments'] and kwargs['moments'][func][2] < kwargs['minSignif'] ) :
            continue

        print ( '  {0:<%d}' % kwargs['maxLenName'] )\
                .format( func if len(func) <= kwargs['maxLenName'] else '...' + func[ 3 - kwargs['maxLenName'] : ] ),
        if func not in kwargs['moments'] :
            print
            continue

        coef = kwargs['moments'][func]
        print '  {0:>+11.5f}   {1:>9.5f}   {2:>12.1f}'\
              .format( coef[0] * kwargs['scales'][0], coef[1] * kwargs['scales'][1], coef[2] * kwargs['scales'][2] )

    print '  ' + '-' * ( 45 + kwargs['maxLenName'] ) + '\n'

    if not kwargs['correlations'] : return

    # print correlation matrix
    print '  correlation matrix:'
    print ' ' * ( 2 + kwargs['maxLenName'] ),
    for func2 in kwargs['funcNames'] :
        if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func2) )\
                or ( func2 in kwargs['moments'] and kwargs['moments'][func2][2] < kwargs['minSignif'] ) :
            continue
        print ( '  {0:>11}' ).format( func2 if len(func2) <= 11 else '.' + func2[ 1 - 11 : ] ),
    print

    for it1, func1 in enumerate(kwargs['funcNames']) :
        if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func1) )\
                or ( func1 in kwargs['moments'] and kwargs['moments'][func1][2] < kwargs['minSignif'] ) :
            continue
        print ( '  {0:<%d}' % kwargs['maxLenName'] )\
                .format( func1 if len(func1) <= kwargs['maxLenName'] else '...' + func1[ 3 - kwargs['maxLenName'] : ] ),
        if func1 not in kwargs['correlations'] :
            print
            continue

        for it2, func2 in enumerate(kwargs['funcNames']) :
            if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func2) )\
                    or ( func2 in kwargs['moments'] and kwargs['moments'][func2][2] < kwargs['minSignif'] ) :
                continue
            if it2 >= it1 :
                print '  {0:>+11.3f}'.format( kwargs['correlations'][func1][func2] ),
            else :
                print ' ' * 13,
        print

    print


def writeMoments( filePath = 'moments', **kwargs ) :
    # parse arguments
    assert getMomentFuncArgs( 'writeMoments', kwargs ), 'P2VV - ERROR: writeMoments(): unable to parse function arguments'

    # get file path and name
    filePath = filePath.strip()
    fileName = filePath.split('/')[-1]

    # open file
    try :
        momFile = open( filePath, 'w' )
    except :
        raise IOError( 'P2VV - ERROR: writeMoments(): unable to open file \"%s\"' % filePath )

    # write moments to content string
    cont = '# %s: angular moments\n' % fileName\
         + '# name requirement: \'{0}\'\n'.format( kwargs['names'] if kwargs['names'] else '' )\
         + '# minimum significance = {0:.1f}\n'.format(kwargs['minSignif'])\
         + '# scales = ({0:.5g}, {1:.5g}, {2:.5g})\n'.format( kwargs['scales'][0], kwargs['scales'][1], kwargs['scales'][2] )\
         + '#\n'\
         + '# ' + '-' * ( 49 + kwargs['maxLenName'] ) + '\n'\
         + ( '# {0:<%s}   {1:<14}   {2:<13}   {3:<13}\n' % kwargs['maxLenName'] )\
               .format( 'basis function', 'coefficient', 'std. dev.', 'significance' )\
         + '# ' + '-' * ( 49 + kwargs['maxLenName'] ) + '\n'

    numMoments = 0
    for func in kwargs['funcNames'] :
        if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func) )\
                or ( func in kwargs['moments'] and kwargs['moments'][func][2] < kwargs['minSignif'] ) :
            continue

        cont += ( '  {0:<%s}' % kwargs['maxLenName'] ).format(func)
        if func in kwargs['moments'] :
            coef = kwargs['moments'][func]
            cont += '   {0:<+14.8g}   {1:<13.8g}   {2:<13.8g}\n'\
                    .format( coef[0] * kwargs['scales'][0], coef[1] * kwargs['scales'][1], coef[2] * kwargs['scales'][2] )
            numMoments += 1

        else :
            cont += '\n'

    cont += '# ' + '-' * (49 + kwargs['maxLenName']) + '\n'

    if kwargs['correlations'] :
        # write correlation matrix to content string
        cont += '#\n# correlation matrix:\n'
        cont += ' ' * ( 2 + kwargs['maxLenName'] )
        for func2 in kwargs['funcNames'] :
            if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func2) )\
                    or ( func2 in kwargs['moments'] and kwargs['moments'][func2][2] < kwargs['minSignif'] ) :
                continue
            cont += ( '  {0:>12}' ).format( func2 if len(func2) <= 12 else '.' + func2[ 1 - 12 : ] )
        cont += '\n'

        for it1, func1 in enumerate(kwargs['funcNames']) :
            if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func1) )\
                    or ( func1 in kwargs['moments'] and kwargs['moments'][func1][2] < kwargs['minSignif'] ) :
                continue
            cont += ( '  {0:<%d}' % kwargs['maxLenName'] ).format(func1)
            if func1 not in kwargs['correlations'] :
                cont += '\n'
                continue

            for it2, func2 in enumerate(kwargs['funcNames']) :
                if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(func2) )\
                        or ( func2 in kwargs['moments'] and kwargs['moments'][func2][2] < kwargs['minSignif'] ) :
                    continue
                cont += '  {0:>+12.6g}'.format( kwargs['correlations'][func1][func2] )
            cont += '\n'

    # write content to file
    momFile.write(cont)
    momFile.close()

    print 'P2VV - INFO: writeMoments(): %d efficiency moment%s written to file \"%s\"'\
            % ( numMoments, '' if numMoments == 1 else 's', filePath )


def readMoments( filePath = 'moments', **kwargs ) :
    # parse arguments
    assert getMomentFuncArgs( 'readMoments', kwargs ), 'P2VV - ERROR: readMoments(): unable to parse function arguments'

    # reset moments and correlations dictionaries
    if not kwargs['addFacs'] :
        if kwargs['procAll'] : del kwargs['funcNames'][ : ]
        kwargs['moments'].clear()
        kwargs['correlations'].clear()

    # get file path
    filePath = filePath.strip()

    # open file
    try :
      momFile = open(filePath, 'r')
    except :
      raise IOError( 'P2VV - ERROR: readMoments(): unable to open file \"%s\"' % filePath )

    # loop over lines and read moments
    from math import sqrt
    numMoments = 0
    corrMatrix = False
    funcNamesMoms = [ ]
    funcNamesCorrs = [ ]
    correlations = [ ]
    while True :
        # read next line
        line = momFile.readline()
        if not line : break

        # check for empty or comment lines and start of correlation matrix
        line = line.strip()
        if 'correlation' in line : corrMatrix = True
        if not line or line[0] == '#' : continue
        line = line.split()

        if not corrMatrix :
            # check moment format
            if len(line) < 4 or ( not kwargs['addFacs'] and not kwargs['procAll'] and line[0] not in kwargs['funcNames'] ) : continue
            try :
              coef   = float(line[1])
              stdDev = float(line[2])
              signif = float(line[3])
            except :
              continue

            # check significance and name
            if ( kwargs['nameExpr'] and not kwargs['nameExpr'].match(line[0]) ) or signif < kwargs['minSignif'] : continue

            # get moment
            if kwargs['addFacs'] :
                assert line[0] in kwargs['moments'], 'P2VV - ERROR: readMoments(): moment %s does not exist yet' % line[0]
                origCoef, origStdDev = kwargs['moments'][line[0]][ : 2 ]
                totCoef = kwargs['addFacs'][0] * origCoef + kwargs['addFacs'][1] * coef
                totStdDev = sqrt( kwargs['addFacs'][0]**2 * origStdDev**2 + kwargs['addFacs'][1]**2 * stdDev**2 )
                kwargs['moments'][line[0]] = ( totCoef, totStdDev, abs(totCoef) / totStdDev )
            else :
                kwargs['moments'][line[0]] = ( coef, stdDev, signif )
            numMoments += 1
            if kwargs['procAll'] : kwargs['funcNames'].append(line[0])
            funcNamesMoms.append(line[0])

        else :
            # read correlation matrix row
            if len(line) < 1 + numMoments or type(line[0]) != str : continue
            try :
                corrs = [ float(corr) for corr in line[ 1 : ] ]
            except :
                continue
            funcNamesCorrs.append(line[0])
            correlations.append(corrs)

    momFile.close()

    # process correlation matrix
    if corrMatrix :
        # check which moments to include in correlation matrix
        skipInds = [ ]
        for it, func in enumerate(funcNamesCorrs) :
            if func not in funcNamesMoms : skipInds.append(it)
        assert len(funcNamesCorrs) - len(skipInds) == numMoments\
               , 'P2VV - ERROR: readMoments(): number of rows in correlation matrix is not consistent with number of moments read'

        if kwargs['addFacs'] :
            print 'P2VV - WARNING: readMoments(): appending correlations matrices not implemented yet'

        else :
            # construct correlation matrix
            for func1, corrs in zip( funcNamesCorrs, correlations ) :
                if func1 not in funcNamesMoms : continue
                kwargs['correlations'][func1] = { }
                for it, func2 in enumerate(funcNamesCorrs) :
                    if func2 not in funcNamesMoms : continue
                    kwargs['correlations'][func1][func2] = corrs[it]

                assert len(kwargs['correlations'][func1]) == numMoments\
                     , 'P2VV - ERROR: readMoments(): number of columns in correlation matrix is not consistent with number of moments read'

            assert len(kwargs['correlations']) == numMoments\
                   , 'P2VV - ERROR: readMoments(): number of rows in correlation matrix is not consistent with number of moments read'

    print 'P2VV - INFO: readMoments(): %d efficiency moment%s read from file \"%s\"'\
          % ( numMoments, '' if numMoments == 1 else 's', filePath )


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
        self._correlations     = { }
        if 'Moments' in kwargs :
            for mom in kwargs.pop('Moments') : self.append( Moment = mom )
        if 'Moment' in kwargs :
            self.append( Moment = kwargs.pop(Moment) )
        self._covMatrixSize = 0
        if kwargs :
            raise RuntimeError( 'unknown keyword arguments %s' % kwargs.keys() )

    def basisFuncs(self)       : return self._basisFuncs.copy()
    def basisFuncNames(self)   : return self._basisFuncNames[ : ]
    def basisFuncIndices(self) : return self._basisFuncIndices.copy()
    def coefficients(self)     : return self._coefficients.copy()
    def correlations(self)     : return dict( [ ( key, val.copy() ) for key, val in self._correlations.iteritems() ] )

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
                from P2VV.RooFitWrappers import P2VVAngleBasis
                momIndices = ( kwargs.pop('PIndex'), kwargs.pop('YIndex0'), kwargs.pop('YIndex1') )
                func = P2VVAngleBasis( kwargs.pop('Angles'), ( momIndices[0], 0, momIndices[1], momIndices[2] ) , 1. )

            if not 'PDF' in kwargs and not 'IntSet' in kwargs and not 'NormSet' in kwargs :
                # build moment
                from P2VV.RooFitWrappers import RealMoment
                moment = RealMoment( Name = func.GetName(), BasisFunc = func, Norm = kwargs.pop( 'Norm', 1. ) )
            elif 'PDF' in kwargs and 'IntSet' in kwargs and 'NormSet' in kwargs :
                # build efficiency moment
                from P2VV.RooFitWrappers import RealEffMoment
                pdf = kwargs.pop('PDF')
                moment = RealEffMoment( Name = func.GetName(), BasisFunc = func, Norm = kwargs.pop( 'Norm', 1. ), PDF = pdf
                                       , IntSet = kwargs.pop('IntSet'), NormSet = kwargs.pop('NormSet') )
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

    def initCovariances(self) :
        self._covMatrixSize = len(self._basisFuncNames)
        for it, name1 in enumerate(self._basisFuncNames) :
            self[name1].clearMoments()
            moments = [ self[name2] for name2 in self._basisFuncNames[ it + 1 : ] ]
            self[name1].appendMoments(moments)

    def compute( self, data ) :
        """computes moments of data set (wrapper for C++ computeRooRealMoments)

        Looping over data in python is quite a bit slower than in C++. Hence, we
        adapt the arguments and then defer to the C++ computeRooRealMoments.
        """
        from P2VV.Load import P2VVLibrary
        from ROOT import std, computeRooRealMoments
        momVec = std.vector('RooRealMoment*')()
        for func in self._basisFuncNames : momVec.push_back( self[func]._var )
        computeRooRealMoments( data, momVec )

        for it1, func1 in enumerate(self._basisFuncNames) :
            self._coefficients[func1] = ( self[func1].coefficient(), self[func1].stdDev(), self[func1].significance() )
            if self._covMatrixSize == len(self._basisFuncNames) :
                corrs = { }
                for it2, func2 in enumerate(self._basisFuncNames) :
                    corrs[func2] = self[func2].correlation( it1 - it2 ) if it2 < it1 else self[func1].correlation( it2 - it1)
                self._correlations[func1] = corrs

    def Print( self, **kwargs ) :
        printMoments( BasisFuncNames = self._basisFuncNames, Moments = self._coefficients, Correlations = self._correlations, **kwargs )

    def write( self, filePath = 'moments', **kwargs ) :
        writeMoments( filePath, BasisFuncNames = self._basisFuncNames, Moments = self._coefficients, Correlations = self._correlations,
                     **kwargs )

    def read( self, filePath = 'moments', **kwargs ) :
        readMoments( filePath, BasisFuncNames = self._basisFuncNames, Moments = self._coefficients, Correlations = self._correlations,
                    **kwargs )
        self._covMatrixSize = len(self._correlations)

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
        from P2VV.RooFitWrappers import ConstVar, RealVar
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

        from P2VV.Parameterizations.AngularPDFs import Coefficients_AngularPdfTerms
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
        from P2VV.RooFitWrappers import ConstVar,RealSumPdf
        # TODO: renornmalize wrt. 0,0,0? require 0,0,0 to be present??
        return RealSumPdf( name
                         , functions = funs
                         , coefficients = ( ConstVar( Name = ('C_%3.6f'%c[0]).replace('-','m').replace('.','d')
                                                    , Value = c[0]*scale ) for c in coefs )
                         )

    def __mul__( self, pdf ) :
        from P2VV.RooFitWrappers import Pdf
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
            from P2VV.RooFitWrappers import P2VVAngleBasis
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
        from P2VV.RooFitWrappers import Addition,EditPdf
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
            splitData       = self._inputData.split(splitCat)
            self._sPlots    = [ ]
            self._sDataSets = [ ]
            sDataVars       = None
            from ROOT import RooFormulaVar
            for splitCatState in splitCat:
                # calculate sWeights per category
                cat = splitCatState.GetName()
                data = splitData.FindObject(cat)

                origYieldVals = [ ( par.GetName(), par.getVal(), par.getError() ) for par in self._yields if par.GetName().endswith(cat) ]
                self._sPlots.append(  RooStats.SPlot( self._name + '_sData_' + cat, self._name + '_sData_' + cat
                                                     , data, self._pdf.getPdf(cat)
                                                     , RooArgList( par._var for par in self._yields if par.GetName().endswith(cat) ) )
                                   )
                self._sDataSets.append( self._sPlots[-1].GetSDataSet() )

                print 'P2VV - INFO: SData.__init__(): yields category %s:' % cat
                print '    original:',
                for vals in origYieldVals : print '%s = %.2f +/- %.2f  ' % vals,
                print '\n    new:     ',
                for par in self._yields :
                    if par.GetName().endswith(cat) : print '%s = %.2f +/- %.2f  ' % ( par.GetName(), par.getVal(), par.getError() ),
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
                               ) for par in self._yields if par.GetName().endswith(cat)
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
            from P2VV.ROOTDecorators import ROOTversion
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
    from P2VV.RooFitWrappers import buildPdf
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
        fullNames = [ ( '%s_%s' % ( parName, stateName[0] ) ) if stateName[0] else parName ]

    name = lambda par : par if type(par) == str else par.GetName()
    for par in parSet :
        if name(par) in fullNames : return par
    return None

def make_binning(data, var, n_bins):
    tmp = data.get().find(var.GetName())
    values = []
    for i in range(data.numEntries()):
        data.get(i)
        values.append((tmp.getVal(), data.weight()))
    
    s = sum(e[1] for e in values)
    d = s / float(n_bins)
    
    from operator import itemgetter
    values = sorted(values, key = itemgetter(0))
    
    bounds = [var.getMin()]
    total = 0
    
    for v, w in values:
        total += w
        if total >= d:
            total = 0
            bounds.append(v)
    bounds.append(var.getMax())
    return bounds
