#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-o', '--rewSteps',    default = 'Bmom_mkk_phys_KKmom' )
parser.add_argument( '-s', '--sample',      default = '20112012'            )
parser.add_argument( '-t', '--trTimeCut',   default = False                 )               
parser.add_argument( '-n', '--numIters',    default = 8,     type=int       )
parser.add_argument( '-N', '--oneIter',     default = 0,     type=int       )
parser.add_argument( '-c', '--numCpu',      default = 8,     type=int       )
parser.add_argument( '-i', '--initPhPars',  default = None                  )
parser.add_argument( '-p', '--makePlots',   default = False, action = 'store_true'  )
parser.add_argument( '-f', '--fit',         default = True,  action = 'store_false' )
parser.add_argument( '-r', '--paralRew',    default = True , action = 'store_false' )
parser.add_argument( '-w', '--writeData',   default = False, action = 'store_true'  )
parser.add_argument( '-m', '--sevdaImpmnt', default = True,  action = 'store_false' )
parser.add_argument( '-b', '--Bmom2DRew',   default = False, action = 'store_true'  )
parser.add_argument( '-D', '--BmomMkk2D',   default = True,  action = 'store_false' )
parser.add_argument( '-e', '--eqStatBins',  default = False, action = 'store_true'  )               
options = parser.parse_args()

# paths and paramteres
path               = '/project/bfys/vsyropou/data/Bs2JpsiPhi/nominalFitResults/'
numberOfIterations = options.numIters
oneIterationScript = 'python /project/bfys/vsyropou/repos/p2vv/scripts/iterativeAngAcc.py'
fittingScript      = 'python /project/bfys/vsyropou/repos/p2vv/scripts/vsFit.py'
fitData            = '/project/bfys/jleerdam/data/Bs2Jpsiphi/angEff/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
parameterEstimates = options.initPhPars if options.initPhPars else path + 'corrAngAccDEC/20112012Reco14DataFitValues_6KKMassBins.par' 
                                                                 # path + 'uncorrAngAccTOS/20112012Reco14DataFitValues_6KKMassBins.par'

# file names 
correctedAngAccBaseName = 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_weights_'
parameterEstimatesName  = lambda n, u: '20112012Reco14DataFitValues_6KKMassBins.par'.replace('.par','_%s_unbl.par'%n) if u else\
                                       '20112012Reco14DataFitValues_6KKMassBins.par'.replace('.par','_%s.par'%n) 

# set up subproceses options
import subprocess, shlex, select, sys
processes  = []
parallelReweighting = options.paralRew
if '2012' in options.sample and not '2011' in options.sample: parallelReweighting = False

combMomOpt     = ' -cTrue'
writeOpt       = ' -wTrue' if options.writeData else ' -wFalse'
plotOpt        = ' -pTrue' if options.makePlots else ' -pFalse'
Bmom2DRewOpt   = options.Bmom2DRew
BmomMkk2DRew   = options.BmomMkk2D
equalStatBins  = options.eqStatBins
KKmomWghtsOnly = options.sevdaImpmnt
rewSpetOpt     = options.rewSteps 
trTimeCut      = options.trTimeCut # cut on true time to calculate ang. acc. in bins of true time.
finalIterOpts  = ' ' + plotOpt + writeOpt

rewOpts = lambda s, n: '-n%i -s%s -d%s -o%s -b%s -e%s -m%s -D%s -t%s'%( n, s, parameterEstimatesName(n-1,True), rewSpetOpt, Bmom2DRewOpt, equalStatBins, KKmomWghtsOnly, BmomMkk2DRew, trTimeCut ) if n!=1 else \
                       '-n%i -s%s -d%s -o%s -b%s -e%s -m%s -D%s -t%s'%( 1, s, parameterEstimates.replace('.par','_unbl.par'), rewSpetOpt, Bmom2DRewOpt, equalStatBins, KKmomWghtsOnly, BmomMkk2DRew, trTimeCut )
fitOpts = lambda n:  '-d%s -a%s -i%s -o%s -c%s'%( fitData, (correctedAngAccBaseName + str(n)), parameterEstimatesName(n-1,False), parameterEstimatesName(n,False), options.numCpu ) if n!=1 else \
                     '-d%s -a%s -i%s -o%s -c%s'%( fitData, (correctedAngAccBaseName + str(n)), parameterEstimates, parameterEstimatesName(n,False), options.numCpu )

# info printing dictionaries
rewOptsLegend = {'-c' : 'Combine eff. moments      ',
                 '-w' : 'Write weightied mc to file',
                 '-p' : 'Plot after reweighting    ',
                 '-n' : 'Iteration index           ',
                 '-s' : 'Monte Carlo dataset       ',
                 '-d' : 'Input physics parameters  ',
                 '-o' : 'Reweighting steps         ',
                 '-b' : '2D B(p,p_T) reweighting   ',
                 '-e' : 'Equal statistics binning  ',
                 '-m' : 'Reweight with w_pkk only  ',
                 '-D' : 'Reweight (B_p,mkk) in 2D  ',
                 '-t' : 'TrueTime cut              '
                 }
fitOptsLegend = {'-d' : 'Fiting dataset          ', 
                 '-a' : 'Input angular acceptance', 
                 '-i' : 'Initial parametes values', 
                 '-o' : 'Fitted parameters values',
                 '-c' : 'Number of cpus          '
                 }

# info prining function
def _info( s, n, opts, what, indent=False ):
    if 'rew' in what:
        indnt = '    ' if indent else ''
        print indnt + 'Iteration number %s. Reweighting mc%s with options:'%(n,s)
        print indnt + 50*'-'
        for o in opts[2:]: 
            for entry in rewOptsLegend.keys():
                if entry in o:
                    print indnt + ' ', rewOptsLegend[entry] + ' (' + o.partition(entry)[1] + '): ', o.partition(entry)[2]
        print
    elif 'fit' in what:
        print '\nIteration number %s. Running 3fb Fit script with options:'%n
        print 50*'-'
        for o in opts[2:]:
            for entry in fitOptsLegend.keys():
                if entry in o:
                    print ' ' + fitOptsLegend[entry] + ' (' + o.partition(entry)[1] + '): ', o.partition(entry)[2]
        print 

###############################
# begin iterative prcedure ####
###############################
# either run 1 to n iterations or run a certain iteration N
whichIterations = range(1, numberOfIterations + 1) if options.numIters and not options.oneIter else \
                  range(options.oneIter, options.oneIter+1)

print 'P2VV - INFO: Begin Iteartive procedure, running %s iteration(s)'%len(whichIterations)
for itNum in whichIterations:

    # set script options
    rew11_Opts = rewOpts( 2011, itNum) + finalIterOpts if itNum==numberOfIterations else rewOpts( 2011, itNum)
    rew12_Opts = rewOpts( 2012, itNum) + finalIterOpts if itNum==numberOfIterations else rewOpts( 2012, itNum)

    # reweight 2011 mc
    if '2011' in options.sample:
        cmd_11 = shlex.split( oneIterationScript + ' ' + rew11_Opts )
        if parallelReweighting:
            print '\nP2VV - INFO: Reweighting of mc2011 and mc2012 samples will run in parallel.'
            _info( '2011', itNum, cmd_11, 'rew', indent=parallelReweighting )
            sys.stdout.flush()
            rew_11 = subprocess.Popen(cmd_11, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
            processes += [ [rew_11, True] ] # append parallel processes to save the output
        else: 
            rew_11 = subprocess.call(cmd_11, stdin = None, stdout = None, stderr = subprocess.STDOUT)
            print
            _info( '2011', itNum, cmd_11, 'rew' )
            sys.stdout.flush()

    if '2012' in options.sample:
        # reweight 2012 mc
        cmd_12 = shlex.split( oneIterationScript + ' ' + rew12_Opts + combMomOpt )
        _info( '2012', itNum, cmd_12, 'rew', indent=parallelReweighting )
        sys.stdout.flush()
        subprocess.call(cmd_12, stdin = None, stdout = None, stderr = subprocess.STDOUT)
    
    # perform 3fb fit
    if options.fit:
        cmdfit = shlex.split( fittingScript + ' ' + fitOpts(itNum) )
        _info( '', itNum, cmdfit, 'fit' )
        sys.stdout.flush()
        subprocess.call(cmdfit, stdin = None, stdout = None, stderr = subprocess.STDOUT)

# print the terminal output into a file
logs = [open('logRew2011_%d' % i, 'w') for i in range(len(processes))]
while any(e[1] for e in processes):
    for i, (p, d) in enumerate(processes):
        if not d:
            continue
        readable = select.select([p.stdout], [], [])
        for line in readable[0][0]:
            if not line.strip(): continue
            logs[i].write(line)
        logs[i].flush()
        if p.poll() != None:
            processes[i][1] = False
            logs[i].close()
