#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-n', '--numIters',   dest='numIters',  default = 8,     type=int,       help='number of iterations')
parser.add_option('-N', '--oneIter',    dest='oneIter',   default = 0,     type=int,       help='run a specific iteration')
parser.add_option('-c', '--numCpu',     dest='numCpu',    default = 8,     type=int,       help='fit with specified number of cpus')
parser.add_option('-p', '--makePlots',  dest='makePlots', default = 'False',               help='switch on/off plotting')
parser.add_option('-s', '--sample',     dest='sample',    default = '20112012',            help='reweight 11/12 or 11+12')
parser.add_option('-f', '--fit',        dest='fit',       default = 'True',                help='switch on/off fitting')
parser.add_option('-o', '--rewSteps',   dest='rewSteps',  default = 'Bmom_mkk_phys_KKmom', help='reweghting steps order')
parser.add_option('-r', '--paralRew',   dest='paralRew',  default = 'True',                help='switch on/off plotting')
parser.add_option('-m', '--sevdaImpmnt',dest='sevdaImpmnt',default = 'True',               help='use only w_pkk weights to calcllate eff. oments')
parser.add_option('-w', '--writeData',  dest='writeData', default = 'False',               help='save mc datasets to file')
parser.add_option('-b', '--Bmom2DRew',  dest='Bmom2DRew', default = 'False',               help='2 dimentional Bmom reweighting switch')
parser.add_option('-D', '--BmomMkk2D',  dest='BmomMkk2D', default = 'True',                help='2 dimentional (Bmom,mkk) reweighting switch')
parser.add_option('-e', '--eqStatBins', dest='eqStatBins',default = 'False',               help='2 dimentional Bmom reweighting switch')
(options, args) = parser.parse_args()

# paths and paramteres
numberOfIterations = options.numIters
oneIterationScript = 'python /project/bfys/vsyropou/Repository/p2vv/scripts/iterativeAngAcc.py'
fittingScript      = 'python /project/bfys/vsyropou/Repository/p2vv/scripts/vsFit.py'
fitData            = '/project/bfys/jleerdam/data/Bs2Jpsiphi/angEff/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
parameterEstimates = '/project/bfys/vsyropou/data/Bs2JsiPhi/nominalFitResults/20112012Reco14DataFitValues_6KKMassBins.par'

# file names 
correctedAngAccBaseName = 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_weights_'
parameterEstimatesName  = lambda n, u: '20112012Reco14DataFitValues_6KKMassBins.par'.replace('.par','_%s_unbl.par'%n) if u else\
                                       '20112012Reco14DataFitValues_6KKMassBins.par'.replace('.par','_%s.par'%n) 

# set up subproceses options
import subprocess, shlex, select, sys
processes  = []
parallelReweighting = True if 'True' in options.paralRew else False
if '2012' in options.sample and not '2011' in options.sample: parallelReweighting = False

combMomOpt     = ' -cTrue'
writeOpt       = ' -wTrue' if 'True' in options.writeData else ' -wFalse'
plotOpt        = ' -pTrue' if 'True' in options.makePlots else ' -pFalse'
Bmom2DRewOpt   = True if 'True' in options.Bmom2DRew else False
BmomMkk2DRew   = True if 'True' in options.BmomMkk2D else False
equalStatBins  = True if 'True' in options.eqStatBins else False
KKmomWghtsOnly = True if 'True' in options.sevdaImpmnt else False
rewSpetOpt     = options.rewSteps 
finalIterOpts  = ' ' + plotOpt + writeOpt

rewOpts = lambda s, n: '-n%i -s%s -d%s -o%s -b%s -e%s -m%s -D%s'%( n, s, parameterEstimatesName(n-1,True), rewSpetOpt, Bmom2DRewOpt, equalStatBins, KKmomWghtsOnly, BmomMkk2DRew ) if n!=1 else \
                       '-n%i -s%s -d%s -o%s -b%s -e%s -m%s -D%s'%( 1, s, parameterEstimates.replace('.par','_unbl.par'), rewSpetOpt, Bmom2DRewOpt, equalStatBins, KKmomWghtsOnly, BmomMkk2DRew )
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
                 '-D' : 'Reweight (B_p,mkk) in 2D  '
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

print 'P2VV - INFO: Begin Iteartive procedure, %s iteration(s)'%len(whichIterations)
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
    if 'True' in options.fit:
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
