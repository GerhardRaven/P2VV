###########################################################################################################################################
## set script parameters ##
###########################

floatVars = [  ( 'tagomega_os_cb', 10 )
             , ( 'tagomega_ss_nn', 10 )
             , ( 'sigmat',         10 )
             , ( 'helcosthetaK',   33 )
             , ( 'helcosthetaL',   33 )
             , ( 'helphi',         33 )
            ]
catVars = [  ( 'runPeriod',            2011 )
           , ( 'hlt1_excl_biased_dec', 0    )
           , ( 'KKMassCat',            0    )
           , ( 'tagCatP2VVOS',         1    )
           , ( 'tagCatP2VVSS',         1    )
           , ( 'iTagOS',               1    )
           , ( 'iTagSS',               1    )
          ]

import sys
if len(sys.argv) > 1 :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--floats', nargs = '+' )
    parser.add_argument( '--categories', nargs = '+' )

    args = parser.parse_args()
    floats = args.floats
    categories = args.categories
    floatVars = [ ]
    catVars = [ ]
    for varVal in floats     : floatVars.append( ( varVal.split(':')[0], int( varVal.split(':')[1] ) )  )
    for varVal in categories : catVars.append( ( varVal.split(':')[0], int( varVal.split(':')[1] ) )  )

print 'floats:'
print floatVars
print 'categories:'
print catVars

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
parFileIn   = 'genValues.par' # '20112012Reco14DataPlotValues_6KKMassBins.par'
dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetName = 'JpsiKK_sigSWeight'
dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2DiegoMass_6KKMassBins_2TagCats.root'

# PDF options
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file']\
        = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file']\
        = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
pdfConfig['anglesEffType'] = 'basisSig4'
pdfConfig['angEffMomsFiles'] = 'genAngAcc' # dataPath + 'Sim08_hel_UB_UT_trueTime_BkgCat050_KK30_Basis'

pdfConfig['lambdaCPParam'] = 'lambPhi'


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
pdfConfig['signalData'] = dataSet
pdfConfig['readFromWS'] = True

# build the PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(dataSet)
pdfPars = pdf.getParameters(dataSet)


###########################################################################################################################################
## get maximum of angular PDF ##
################################

# get observables
scanObsSet = [ ws[var] for var in [ 'helcosthetaK', 'helcosthetaL', 'helphi' ] ]
intObsSet  = [ ws['time'] ]

# get conditional observables
from P2VV.RooFitWrappers import SimultaneousPdf
condObsSet = [ obs for obs in pdf.ConditionalObservables() ]
if isinstance( pdf, SimultaneousPdf ) :
    if pdf.indexCat().isFundamental() :
        condObsSet += [ pdf.indexCat() ]
    else :
        condObsSet += [ cat for cat in pdf.indexCat().inputCatList() ]

# split conditional observables into categories and floats
from ROOT import RooCategory
from P2VV.RooFitWrappers import Category
condCatSet   = [ ]
condFloatSet = [ ]
for obs in condObsSet :
    if isinstance( obs, Category ) or isinstance( obs, RooCategory ) :
        condCatSet.append(obs)
    else :
        condFloatSet.append(obs)

# set conditional-category values
catVarsDict = dict(catVars)
for obs in condCatSet :
    assert obs.GetName() in catVarsDict
    obs.setIndex( catVarsDict.pop( obs.GetName() ) )
assert not catVarsDict

# get scan variables
obsSet = scanObsSet + condFloatSet
scanSet = [ ]
for varVal in floatVars :
    scanVar = None
    for var in scanObsSet + condFloatSet :
        if var.GetName() == varVal[0] :
            scanVar = var
            obsSet.remove(var)
            break
    assert var
    scanSet.append(var)
assert not obsSet

from P2VV.RooFitWrappers import getFuncMaxVal
maxVal = getFuncMaxVal( Function = pdf, IntegrationSet = intObsSet, NormalizationSet = intObsSet + scanObsSet
                   , ScanSet = scanSet, NumScanPoints = dict(floatVars) )

print '\nPDF:'
pdf.Print()

print '\nscan set:'
for obs in condFloatSet + scanObsSet : obs.Print()

print '\nconditional categories:'
for obs in condCatSet : obs.Print()

print '\nintegration set:'
for obs in intObsSet : obs.Print()

print '\nnormalization set:'
for obs in intObsSet + scanObsSet : obs.Print()

print '\nmaximum PDF value = %.3g' % maxVal
