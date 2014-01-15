#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-n', '--numIters',  dest='numIters',  default=8,     type=int,         help='number of iterations'  )
parser.add_option('-p', '--makePlots', dest='makePlots', default='False',                   help='switch on/off plotting')
parser.add_option('-o', '--rewSteps',  dest='rewSteps',  default = 'Bmom_mkk_phys_KKmom', help='reweghting steps order')
parser.add_option('-r', '--paralRew',  dest='paralRew',  default='True',                    help='switch on/off plotting')
(options, args) = parser.parse_args()

# paths and paramteres
numberOfIterations = options.numIters
oneIterationScript = 'python /project/bfys/vsyropou/Repository/p2vv/scripts/iterativeAngAcc.py'
fittingScript      = 'python /project/bfys/vsyropou/Repository/p2vv/scripts/vsFit.py'
fitData            = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
parameterEstimates = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/nominalFitResults/20112012Reco14DataFitValues_6KKMassBins.par'
startingAngAccFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'

correctedAngAccBaseName    = 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm_'
parameterEstimatesName = lambda n: '20112012Reco14DataFitValues_6KKMassBins.par'.replace('.par','_%s.par'%n)

# set up subproceses options
import subprocess, shlex, select
processes  = []
parallelReweighting = True if 'True' in options.paralRew else False

combMomOpt    = ' -cTrue'
writeOpt      = ' -wTrue'
plotOpt       = ' -pTrue' if 'True' in options.makePlots else ' -pFalse'
rewSpetOpt    = options.rewSteps 
finalIterOpts = ' ' + plotOpt + writeOpt

rewOpts = lambda s, n: '-n%i -s%s -d%s -o%s -fFalse'%( n, s, parameterEstimatesName(n-1), rewSpetOpt ) if n!=1 else \
                       '-n%i -s%s -d%s -o%s -fFalse'%( 1, s, parameterEstimates, rewSpetOpt )
fitOpts = lambda n:  '-d%s -a%s -i%s -o%s'%( fitData, (correctedAngAccBaseName + str(n)).replace('Phys_norm','weights'), parameterEstimatesName(n-1), parameterEstimatesName(n) ) if n!=1 else \
                     '-d%s -a%s -i%s -o%s'%( fitData, (correctedAngAccBaseName + str(n)).replace('Phys_norm','weights'), parameterEstimates, parameterEstimatesName(n) )

rewOptsLegend = {'-c' : 'Combine eff. moments      ',
                 '-w' : 'Write weightied mc to file',
                 '-p' : 'Plot after reweighting    ',
                 '-n' : 'Iteration index           ',
                 '-s' : 'Monte Carlo dataset       ',
                 '-d' : 'Input physics parameters  ',
                 '-o' : 'Reweighting steps         ',
                 '-f' : 'Fit after reweighting     '
                 }
fitOptsLegend = {'-d' : 'Fiting dataset          ', 
                 '-a' : 'Input angular acceptance', 
                 '-i' : 'Initial parametes values', 
                 '-o' : 'Fitted parameters values', 
                 }

# info printing function
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
        print 'Iteration number %s. Running 3fb Fit script with options:'%n
        print 50*'-'
        for o in opts[2:]:
            for entry in fitOptsLegend.keys():
                if entry in o:
                    print ' ' + fitOptsLegend[entry] + ' (' + o.partition(entry)[1] + '): ', o.partition(entry)[2]
        print 

###############################
# begin iterative prcedure ####
###############################
print 'P2VV - INFO: Begin Iteartive procedure'
for itNum in range(1, numberOfIterations + 1):

    # set script options
    rew11_Opts = rewOpts( 2011, itNum) + finalIterOpts if itNum==numberOfIterations else rewOpts( 2011, itNum)
    rew12_Opts = rewOpts( 2012, itNum) + finalIterOpts if itNum==numberOfIterations else rewOpts( 2012, itNum)

    # reweight 2011 mc
    cmd_11 = shlex.split( oneIterationScript + ' ' + rew11_Opts )
    if parallelReweighting:
        print 'P2VV - INFO: Reweighting of mc2011 and mc2012 samples will run in parallel.'
        _info( '2011', itNum, cmd_11, 'rew', indent=parallelReweighting )
        rew_11 = subprocess.Popen(cmd_11, stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    else: 
        rew_11 = subprocess.call(cmd_11, stdin = None, stdout = None, stderr = subprocess.STDOUT)
        _info( '2011', itNum, cmd_11, 'rew' )

    # reweight 2012 mc
    cmd_12 = shlex.split( oneIterationScript + ' ' + rew12_Opts + combMomOpt )
    _info( '2012', itNum, cmd_12, 'rew', indent=parallelReweighting )
    rew_12 = subprocess.call(cmd_12, stdin = None, stdout = None, stderr = subprocess.STDOUT)
    
    # perform 3fb fit
    cmdfit = shlex.split( fittingScript + ' ' + fitOpts(itNum) )
    _info( '', itNum, cmdfit, 'fit' )
    fit3fb = subprocess.call(cmdfit, stdin = None, stdout = None, stderr = subprocess.STDOUT)

logs = [open('log_%d' % i, 'w') for i in range(len(processes))]
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
