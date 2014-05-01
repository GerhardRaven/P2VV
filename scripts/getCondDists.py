import P2VV.RooFitWrappers
from ROOT import TFile
dataFile = TFile.Open('/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root')
#dataFile = TFile.Open('/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_cbkg_I2Mass_6KKMassBins_2TagCats_HLT2B_20140421.root')
data = dataFile.Get('JpsiKK_sigSWeight')
#data = dataFile.Get('JpsiKK_cbkgSWeight')
dataObs = data.get()
dataFile.Close()

condObsNames = [ 'runPeriod', 'HLT1Cat', 'KKMassCat', 'KKMass', 'sigmat', 'etaOS', 'etaSS', 'tagOS', 'tagSS' ]
condObs = dict(  runPeriod = dataObs.find('runPeriod')
               , HLT1Cat   = dataObs.find('hlt1_excl_biased')
               , KKMassCat = dataObs.find('KKMassCat')
               , KKMass    = dataObs.find('mdau2')
               , sigmat    = dataObs.find('sigmat')
               , etaOS     = dataObs.find('tagomega_os_cb')
               , etaSS     = dataObs.find('tagomega_ss_nn')
               , tagOS     = dataObs.find('iTagOS')
               , tagSS     = dataObs.find('iTagSS')
              )
condObsBins = dict(  sigmat = [  0.0167, 0.0184, 0.0196, 0.0206, 0.0214, 0.0222, 0.0228, 0.0235, 0.0240
                               , 0.0246, 0.0251, 0.0256, 0.0261, 0.0265, 0.0270, 0.0275, 0.0279, 0.0283, 0.0288
                               , 0.0292, 0.0297, 0.0301, 0.0305, 0.0310, 0.0314, 0.0319, 0.0323, 0.0328, 0.0332
                               , 0.0337, 0.0342, 0.0346, 0.0352, 0.0357, 0.0362, 0.0367, 0.0373, 0.0380, 0.0386
                               , 0.0394, 0.0401, 0.0410, 0.0420, 0.0431, 0.0444, 0.0461, 0.0482, 0.0512, 0.0571
                              ]
		   , KKMass = [  1008.0, 1012.0, 1013.7, 1014.7, 1015.5, 1016.0, 1016.4, 1016.8, 1017.2, 1017.5, 1017.7, 1017.9, 1018.1
                               , 1018.3, 1018.5, 1018.7, 1018.9, 1019.1, 1019.3, 1019.5, 1019.7, 1019.8, 1019.9, 1020.0, 1020.1, 1020.2
                               , 1020.3, 1020.5, 1020.7, 1020.9, 1021.2, 1021.5, 1021.8, 1022.1, 1022.4, 1022.7, 1023.1, 1023.5, 1024.0
                               , 1024.5, 1025.2, 1026.0, 1027.0, 1028.2, 1029.9, 1032.0, 1034.8, 1038.7, 1043.6
                              ]
                   , etaOS  = [ 0.336, 0.400, 0.450, 0.4999999 ]
                   , etaSS  = [ 0.397, 0.437, 0.456, 0.469, 0.481, 0.492, 0.4999999 ]
                  )

sumWTot = 0.
sumWSqTot = 0.
dists = dict( [ ( name, [ [0., 0.] for it in range( len(condObsBins[name]) + 1 if name in condObsBins else condObs[name].numTypes() ) ] )\
                for name in condObsNames if name[ -2 : ] ] )
dists['eta'] = [ [ ( [ 0., 0. ], ( [ 0., 0. ], [ 0., 0. ], [ 0., 0. ], [ 0., 0. ] ) )\
                 for it1 in range( len( condObsBins['etaSS'] ) + 1 ) ] for it0 in range( len( condObsBins['etaOS'] ) + 1 ) ]
for ev in data :
    #if condObs['runPeriod'].getIndex() == 2012 : continue

    etaOSInd = 0
    etaSSInd = 0
    tagOSInd = 0
    tagSSInd = 0
    evW = data.weight()
    evWSq = data.weight()**2
    sumWTot += evW
    sumWSqTot += evWSq
    for name in condObsNames :
        ind = 0
        if name in condObsBins :
            for thresh in condObsBins[name] :
                if condObs[name].getVal() < thresh : break
                ind += 1
            if   name == 'etaOS' : etaOSInd = ind
            elif name == 'etaSS' : etaSSInd = ind
        else :
            for state in condObs[name] :
                if condObs[name].getIndex() == state.getVal() : break
                ind += 1
                if   name == 'tagOS' : tagOSInd = ind
                elif name == 'tagSS' : tagSSInd = ind

        dists[name][ind][0] += evW
        dists[name][ind][1] += evWSq

    dists['eta'][etaOSInd][etaSSInd][0][0] += evW
    dists['eta'][etaOSInd][etaSSInd][0][1] += evWSq
    dists['eta'][etaOSInd][etaSSInd][1][ 2 * ( tagOSInd != tagSSInd ) + tagOSInd ][0] += evW
    dists['eta'][etaOSInd][etaSSInd][1][ 2 * ( tagOSInd != tagSSInd ) + tagOSInd ][1] += evWSq

from math import sqrt, log10, ceil
for name in condObsNames :
    print '%s: [%s]' % ( name, ', '.join( '%.6f' % ( sumW / sumWTot ) for [ sumW, sumWSq ] in dists[name] ) )
    for it, [ sumW, sumWSq ] in enumerate( dists[name] ) :
        frac = sumW / sumWTot
        fracErr = sqrt( ( sumWTot - sumW )**2 * sumWSq + sumW**2 * ( sumWSqTot - sumWSq ) ) / sumWTot**2
        prec = max( 0, 1 - int( ceil( log10(fracErr) ) ) ) if fracErr > 0. else 0
        print ( '{0:>2d}: {1:<6.%df} +/- {2:<6.%df}' % ( prec, prec ) ).format( it, frac, fracErr )
    print

print 'eta:'
for OSBin, SSVals in enumerate( dists['eta'] ) :
    print '%d: [%s]' % ( OSBin, ', '.join( '%.6f' % ( vals[0][0] / sumWTot ) for vals in SSVals ) )
for OSBin, SSVals in enumerate( dists['eta'] ) :
    print '%d: [%s]' % ( OSBin, ', '.join( '(%.6f, (%.6f, %.6f))'\
          % ( vals[0][0] / sumWTot, ( vals[1][0][0] + vals[1][1][0] ) / vals[0][0], ( vals[1][2][0] + vals[1][3][0] ) / vals[0][0] )\
          for vals in SSVals ) )
for OSBin, SSVals in enumerate( dists['eta'] ) :
    for SSBin, vals in enumerate(SSVals) :
        sumWCat   = vals[0][0]
        sumWSqCat = vals[0][1]
        frac = sumWCat / sumWTot
        fracErr = sqrt( ( sumWTot - sumWCat )**2 * sumWSqCat + sumWCat**2 * ( sumWSqTot - sumWSqCat ) ) / sumWTot**2
        prec = max( 0, 1 - int( ceil( log10(fracErr) ) ) ) if fracErr > 0. else 0

        assert abs( sumWCat - sum( v[0] for v in vals[1] ) ) < 1.e-4
        fracs = [ v[0] / sumWCat for v in vals[1] ]
        fracErrs = [ sqrt( ( sumWCat - v[0] )**2 * v[1] + v[0]**2 * ( sumWSqCat - v[1] ) ) / sumWCat**2 for v in vals[1] ]
        precs = [ ( max( 0, 1 - int( ceil( log10(err) ) ) ) if err > 0. else 0 ) for err in fracErrs ]

        fracEq = ( vals[1][0][0] + vals[1][1][0] ) / sumWCat
        fracEqErr = sqrt( ( vals[1][2][0] + vals[1][3][0] )**2 * ( vals[1][0][1] + vals[1][1][1] )\
                         + ( vals[1][0][0] + vals[1][1][0] )**2 * ( vals[1][2][1] + vals[1][3][1] ) ) / sumWCat**2
        precEq = max( 0, 1 - int( ceil( log10(fracEqErr) ) ) ) if fracEqErr > 0. else 0

        print ( '{0:d}-{1:d}: {2:<6.%df} +/- {3:<6.%df}   eq: {4:<6.%df} +/- {5:<6.%df}   {6}' % ( prec, prec, precEq, precEq ) )\
                .format( OSBin, SSBin, frac, fracErr, fracEq, fracEqErr, '   '.join( ('{0:<6.%df} +/- {1:<6.%df}' % ( p, p ) )\
                .format( f, e ) for f, e, p in zip( fracs, fracErrs, precs ) ) )
