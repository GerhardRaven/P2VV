###########################################################################################################################################
## Utilities.DataHandling: P2VV utilities to read, write and manipulate data sets                                                        ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##   RA,  Roel Aaij,          Nikhef,      r.aaij@nikhef.nl                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
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

    dataSetArgs = { }
    if 'WeightVar' in kwargs : dataSetArgs['WeightVar'] = kwargs.pop('WeightVar')
    if NTuple :
      from ROOT import RooDataSet, TFile
      assert observables != None, 'P2VV - ERROR: readData: set of observables is required for reading an n-tuple'

      # create data set from NTuple file(s)
      print 'P2VV - INFO: readData: reading NTuple "%s" from file "%s"' % ( dataSetName, filePath )
      ntupleFile = TFile.Open(filePath)
      ntupleOrig = ntupleFile.Get(dataSetName)
      if not ntupleOrig : raise RuntimeError( 'P2VV - ERROR: could not locate tree "%s" in file "%s"' % ( dataSetName, filePath ) )

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
          ntuple = ntupleOrig.CopyTree(ntupleCuts)
      else :
          ntuple = ntupleOrig

      if cuts : print 'P2VV - INFO: readData: applying cuts on data set: %s' % cuts
      data = RooDataSet(  dataSetName, dataSetName
                       , [ obs._var for obs in observables ]
                       , Import = ntuple
                       , Cut = noNAN + ' && ' + cuts if cuts else noNAN
                       , **dataSetArgs
                       )

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
                                   , **dataSetArgs
                                   )
              else :
                  data = dataSet
          else :
              data.append(dataSet)

      file.Close()

    print 'P2VV - INFO: read dataset with %s entries (%.1f weighted)' % ( data.numEntries(), data.sumEntries() )

    # import data set into current workspace
    from P2VV.RooFitWrappers import RooObject
    importIntoWS = kwargs.pop( 'ImportIntoWS', True )
    rData = None
    if importIntoWS :
        # import data set into current workspace
        from P2VV.RooFitWrappers import RooObject
        wsData = RooObject().ws().put( data, **kwargs )
        rData = wsData
    else :
        rData = data

    if tmp_file:
        tmp_file.Close()
        os.remove(tmp_file.GetName())
        if orig_file: orig_file.cd()
    return rData


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


def correctWeights( dataSet, splitCatNames, **kwargs ) :
    """correct weights in dataSet for background dilution
    """

    if splitCatNames :
        # get category that splits data sample
        if type(splitCatNames) == str : splitCatNames = [ splitCatNames ]
        if len(splitCatNames) == 1 :
            splitCat = dataSet.get().find( splitCatNames[0] )
            assert splitCat, 'P2VV - ERROR: correctWeights: unknown split category: "%s"' % splitCatNames[0]

        else :
            from ROOT import RooArgSet
            catSet = RooArgSet()
            for catName in splitCatNames :
                cat = dataSet.get().find(catName)
                assert cat, 'P2VV - ERROR: correctWeights: unknown split category: "%s"' % catName
                catSet.add(cat)

            from ROOT import RooSuperCategory
            splitCat = RooSuperCategory( 'splitCat', 'splitCat', catSet)

        indexDict = { }
        for iter, catType in enumerate( splitCat ) : indexDict[ iter ] = ( catType.GetName(), catType.getVal() )

    # get variable that contains original event weights
    origWeightVar = None
    origWeightName = kwargs.pop( 'WeightName', '' )
    if origWeightName :
        origWeightVar = dataSet.get().find(origWeightName)
        from ROOT import RooAbsReal
        assert origWeightVar and isinstance( origWeightVar, RooAbsReal )\
               , 'P2VV - ERROR: correctWeights: unknown weight variable: "%s"' % origWeightName

    corrFactors = kwargs.pop( 'CorrectionFactors', [ ] )
    if not corrFactors :
        if splitCatNames :
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
            weight = origWeightVar.getVal() if origWeightVar else dataSet.weight()
            sumWeights[0]   += weight
            sumSqWeights[0] += weight**2
            if splitCatNames :
                sumWeights[   posDict[ splitCat.getIndex() ] + 1 ] += weight
                sumSqWeights[ posDict[ splitCat.getIndex() ] + 1 ] += weight**2

        # get correction factors
        corrFactors = (  ( sumWeights[0] / sumSqWeights[0], sumWeights[0] )
                       , [ ( sum / sumSq, sum ) for sum, sumSq in zip( sumWeights[ 1 : ], sumSqWeights[ 1 : ] ) ]
                      )

    # add corrected weights to data set
    from P2VV.Load import P2VVLibrary
    from ROOT import RooCorrectedWeight
    if splitCatNames :
        from ROOT import std
        corrFactorsVec = std.vector('Double_t')()
        print 'P2VV - INFO: correctWeights: multiplying sWeights (-ln(L)) to correct for background dilution with factors (overall factor %.4f for %.1f events):'\
              % corrFactors[0]
        maxLenStateName = max( len( ind[0] ) for ind in indexDict.values() )
        for iter, fac in enumerate( corrFactors[1] ) :
            corrFactorsVec.push_back( fac[0] )
            print ( '    {0:%ds} : {1:.4f} (events: {2:.1f})' % maxLenStateName ).format( indexDict[iter][0], fac[0], fac[1] )

        if origWeightVar :
            weightVar = RooCorrectedWeight( 'weightVar', 'weight variable', origWeightVar, splitCat, corrFactorsVec )
        else :
            weightVar = RooCorrectedWeight( 'weightVar', 'weight variable', dataSet, splitCat, corrFactorsVec )

    else :
        print 'P2VV - INFO: correctWeights: multiplying sWeights (-ln(L)) to correct for background dilution with a factor %.4f for %.1f events'\
              % corrFactors[0]
        if origWeightVar :
            weightVar = RooCorrectedWeight( 'weightVar', 'weight variable', origWeightVar, corrFactors[0][0] )
        else :
            weightVar = RooCorrectedWeight( 'weightVar', 'weight variable', dataSet, corrFactors[0][0] )

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
        if len(binBounds) < 3 : continue
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


def addGlobalTagCat( dataSet, tagCatOSName, tagCatSSName ) :
    tagCatOS = dataSet.get().find(tagCatOSName)
    tagCatSS = dataSet.get().find(tagCatSSName)
    assert tagCatOS and tagCatSS\
           , 'P2VV - ERROR: addGlobalTagCat(): one or more tagging categories not found in data set ("%s" and "%s")'\
             % ( tagCatOSName, tagCatSSName )

    from ROOT import RooMappedCategory, RooMultiCategory, RooArgSet
    tagCatMult = RooMultiCategory( 'tagCatMult', 'tagCatMult', RooArgSet( tagCatOS, tagCatSS ) )
    tagCatP2VV = RooMappedCategory( 'tagCatP2VV', 'tagCatP2VV', tagCatMult, 'Tagged', 1 )

    tagCatOS.setIndex(0)
    tagCatSS.setIndex(0)
    tagCatP2VV.map( tagCatMult.getLabel(), 'Untagged', 0 )

    tagCatP2VV = dataSet.addColumn(tagCatP2VV)

    for obsSet in dataSet :
        assert (obsSet.getCatIndex('tagCatP2VV') == 0 and obsSet.getCatIndex(tagCatOSName) == 0 and obsSet.getCatIndex(tagCatSSName) == 0)\
            or (obsSet.getCatIndex('tagCatP2VV') == 1 and (obsSet.getCatIndex(tagCatOSName) > 0 or obsSet.getCatIndex(tagCatSSName) > 0))\
            , 'P2VV - ERROR: addGlobalTagCat(): "tagCatP2VV" index (%d) does not correspond to "%s" and "%s" indices (%d and %d)'\
              % ( obsSet.getCatIndex('tagCatP2VV'), tagCatOSName, tagCatSSName
                 , obsSet.getCatIndex(tagCatOSName), obsSet.getCatIndex(tagCatSSName) )


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


def addTrackMomenta( dataSet ):
    """ add K+, K-, mu+ and mu- momentum magnitudes to data set
    """

    # get observables from data set
    from ROOT import RooArgList, RooFormulaVar
    KplusList   = RooArgList('KplusList')
    KminusList  = RooArgList('KminusList')
    muplusList  = RooArgList('muplusList')
    muminusList = RooArgList('muminusList')

    for suf in [ 'X','Y','Z' ] :
        KplusList.add( dataSet.get().find('Kplus_P'  + suf) )
        KminusList.add( dataSet.get().find('Kminus_P'  + suf) )
        muplusList.add( dataSet.get().find('muplus_P' + suf) )
        muminusList.add( dataSet.get().find('muminus_P' + suf) )

    #  create formulas
    KplusMomFunc   = RooFormulaVar( 'Kplus_P',  'Kplus_P',   'TMath::Sqrt(Kplus_PX**2 + Kplus_PY**2 + Kplus_PZ**2)',       KplusList   )
    KminusMomFunc  = RooFormulaVar( 'Kminus_P', 'Kminus_P',  'TMath::Sqrt(Kminus_PX**2 + Kminus_PY**2 + Kminus_PZ**2)',    KminusList  )
    muplusMomFunc  = RooFormulaVar( 'muplus_P', 'muplus_P',  'TMath::Sqrt(muplus_PX**2 + muplus_PY**2 + muplus_PZ**2)',    muplusList  )
    muminusMomFunc = RooFormulaVar( 'muminus_P','muminus_P', 'TMath::Sqrt(muminus_PX**2 + muminus_PY**2 + muminus_PZ**2)', muminusList )

    # create new columns in data set
    dataSet.addColumn(KplusMomFunc)
    dataSet.addColumn(KminusMomFunc)
    dataSet.addColumn(muplusMomFunc)
    dataSet.addColumn(muminusMomFunc)


def printEventYields( **kwargs ) :
    # get arguments
    parSet     = kwargs.pop( 'ParameterSet',        [ ] )
    yieldNames = kwargs.pop( 'YieldNames',          [ ] )
    splitCats  = kwargs.pop( 'SplittingCategories', [ ] )

    assert parSet,     'P2VV - ERROR: printEventYields: no parameter set with yield variables found in arguments ("ParameterSet")'
    assert yieldNames, 'P2VV - ERROR: printEventYields: no yield names found in arguments ("YieldNames")'

    if splitCats :
        splitCats = sorted( list( set( cat for cat in splitCats ) ), cmp = lambda a, b : cmp( a.GetName(), b.GetName() ) )

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
    print '-'.join( dashes for dashes in [ ' ' * 4 + '-' * 26, '-' * 8 ] + [ '-' * num for name in yieldNames for num in ( 23, 19 ) ] )
    print ' '.join( ( '{0:^%ds}' % width ).format(title) for ( title, width ) in [ ( '', 30 ), ( ' total', 8 ) ]\
                   + [ ( name if num == 23 else 'f(' + name + ')', num ) for name in yieldNames for num in ( 23, 19 ) ] )
    print '-'.join( dashes for dashes in [ ' ' * 4 + '-' * 26, '-' * 8 ] + [ '-' * num for name in yieldNames for num in ( 23, 19 ) ] )

    from P2VV.Utilities.General import getSplitPar
    cont = True
    while cont :
        stateName = ';'.join( labs[cat][ iters[cat] ] for cat in splitCats )
        yields = [ getSplitPar( name, ( '{%s}' % stateName ) if stateName else '', parSet ) for name in yieldNames ]

        from math import sqrt
        nEv        = [ yieldVar.getVal()   for yieldVar in yields ]
        nEvErr     = [ yieldVar.getError() for yieldVar in yields ]
        nEvTot     = sum(nEv)
        frac       = [ num / nEvTot if nEvTot > 0. else 0. for num in nEv ]
        nEvCorrErr = [ sqrt( ( numErr**2 - num**2 / nEvTot ) if numErr**2 > num**2 / nEvTot else 0. ) for num, numErr in zip(nEv, nEvErr) ]
        fracErr    = [ err / nEvTot if nEvTot > 0. else 0. for err in nEvCorrErr ]

        print '     {0:>24s}   {1:>6.0f}  '.format( stateName, nEvTot )\
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

    print '-'.join( dashes for dashes in [ ' ' * 4 + '-' * 26, '-' * 8 ] + [ '-' * num for name in yieldNames for num in ( 23, 19 ) ] )
    print


def printEventYieldsData( **kwargs ) :
    # get arguments
    weightDataSets = kwargs.pop( 'WeightedDataSets',    [ ]  )
    dataSetNames   = kwargs.pop( 'DataSetNames',        [ ]  )
    splitCats      = kwargs.pop( 'SplittingCategories', [ ]  )

    assert weightDataSets, 'P2VV - ERROR: printEventYieldsData: no weighted data sets found in arguments ("WeightedDataSets")'
    if not dataSetNames : dataSetNames = [ 'data set %d' % it for it, dataSet in enumerate(weightDataSets) ]

    if splitCats :
        splitCats = sorted( list( set( cat for cat in splitCats ) ), cmp = lambda a, b : cmp( a.GetName(), b.GetName() ) )

    # print total numbers of events
    from math import sqrt
    nEv    = [ dataSet.sumEntries() for dataSet in weightDataSets ]
    nEvTot = sum(nEv)
    frac   = [ num / nEvTot       if nEvTot > 0. else 0. for num in nEv ]
    signif = [ num / sqrt(nEvTot) if nEvTot > 0. else 0. for num in nEv ]

    print
    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 41 ] + [ '-' * 30 for dataSet in weightDataSets ] )
    print ' ' * 45 + '|' + '|'.join( ' {0:^28} '.format(name) for name in dataSetNames )
    print ' ' * 37 + '  total |' + '|'.join( ' {0:^9s}   {1:^6s}   {2:^9s} '\
          .format( 'N_%d' % it, 'f_%d' % it, u'N_%d/\u221AN'.encode('utf-8') % it ) for it, dataSet in enumerate(weightDataSets) )
    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 41 ] + [ '-' * 30 for dataSet in weightDataSets ] )
    print ' ' *  4 + ' {0:>30s}   {1:>6.0f} |'.format( 'total', nEvTot ) + '|'.join( ' {0:>9.2f}   {1:>6.4f}   {2:>7.3f} '\
          .format( num, fr, sig ) for ( num, fr, sig ) in zip( nEv, frac, signif ) )
    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 41 ] + [ '-' * 30 for dataSet in weightDataSets ] )

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
            nEv    = [ dataSet.sumEntries(cut) for dataSet in weightDataSets ]
            nEvTot = sum(nEv)
            frac   = [ num / nEvTot       if nEvTot > 0. else 0. for num in nEv ]
            signif = [ num / sqrt(nEvTot) if nEvTot > 0. else 0. for num in nEv ]
            print ' ' *  4 + ' {0:>30s}   {1:>6.0f} |'.format( labs[cat][-1], nEvTot ) + '|'.join( ' {0:>9.2f}   {1:>6.4f}   {2:>7.3f} '\
                  .format( num, fr, sig ) for ( num, fr, sig ) in zip( nEv, frac, signif ) )

            catState = catIter.Next()

        print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 41 ] + [ '-' * 30 for dataSet in weightDataSets ] )

    if len(splitCats) < 2 :
        print
        return

    # print numbers of events for each combination of splitting categories
    cont = True
    while cont :
        stateName = ';'.join( labs[cat][ iters[cat] ] for cat in splitCats )
        cut       = '&&'.join( '!(%s-%d)' % ( cat.GetName(), inds[cat][ iters[cat] ] ) for cat in splitCats )
        nEv       = [ dataSet.sumEntries(cut) for dataSet in weightDataSets ]
        nEvTot    = sum(nEv)
        frac      = [ num / nEvTot       if nEvTot > 0. else 0. for num in nEv ]
        signif    = [ num / sqrt(nEvTot) if nEvTot > 0. else 0. for num in nEv ]
        print ' ' *  4 + ' {0:>30s}   {1:>6.0f} |'.format( stateName, nEvTot ) + '|'.join( ' {0:>9.2f}   {1:>6.4f}   {2:>7.3f} '\
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

    print ' ' *  4 + '|'.join( dashes for dashes in [ '-' * 41 ] + [ '-' * 30 for dataSet in weightDataSets ] )
    print
