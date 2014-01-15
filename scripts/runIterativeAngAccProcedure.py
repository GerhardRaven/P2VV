#!/usr/bin/env python
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-n', '--numIters',  dest='numIters',  default=8,     type=int, help='number of iterations'  )
parser.add_option('-p', '--makePlots', dest='makePlots', default=False,           help='switch on/off plotting')
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
rewOpts = lambda s, n: '-n%i -s%s -d%s -fFalse'%( n, s, parameterEstimatesName(n-1) ) if n!=1 else \
                       '-n%i -s%s -d%s -fFalse'%( 1, s, parameterEstimates )
fitOpts = lambda n:  '-d%s -a%s -i%s -o%s'%( fitData, correctedAngAccBaseName + str(n), parameterEstimatesName(n-1), parameterEstimatesName(n) ) if n!=1 else \
                     '-d%s -a%s -i%s -o%s'%( fitData, correctedAngAccBaseName + str(n), parameterEstimates, parameterEstimatesName(n) )
combMomOpt = ' -cTrue'
writeOpt   = ' -wTrue'
plotOpt    = ' -pTrue' if options.makePlots else ''

processes  = []

def _info( s, n, opts, what ):
    if 'rew' in what:
        print 'Iteration number %s. Reweighting mc%s with options:'%(s,n)
        for o in opts[2:]: print o
    elif 'fit' in what:
        print 'Iteration number %s. Running 3fb Fit script with options:'%n
        for o in opts[2:]: print o

# _info( '', itNum, cmd_, 'rew' )
# _info( '', itNum, cmdfit, 'fit' )

###############################
# begin iterative prcedure ####
###############################
print 'P2VV - INFO: Begin Iteartive procedure'
for itNum in range(1, numberOfIterations):
    cmd_11 = shlex.split( oneIterationScript + ' ' + rewOpts( 2011, itNum) )
    # print 'Iteration number %s. Reweighting mc2011 with options:'%itNum 
    # for o in cmd_11[2:]: print o
    _info( '2011', itNum, cmd_11, 'rew' )
    rew_11 = subprocess.call(cmd_11, stdin = None, stdout = None, stderr = subprocess.STDOUT)

    cmd_12 = shlex.split( oneIterationScript + ' ' + rewOpts( 2012, itNum) + combMomOpt )
    _info( '2012', itNum, cmd_12, 'rew' )
    # print 'Iteration number %s. Reweighting mc2012 with options:'%itNum 
    # for o in cmd_12[2:]: print o
    rew_12 = subprocess.call(cmd_12, stdin = None, stdout = None, stderr = subprocess.STDOUT)

    cmdfit = shlex.split( fittingScript + ' ' + fitOpts(itNum) )
    # print 'Iteration number %s. Running 3fb Fit script with options:'%itNum 
    # for o in cmdfit[2:]: print o
    _info( '', itNum, cmd_fit, 'fit' )
    fit3fb = subprocess.call(cmdfit, stdin = None, stdout = None, stderr = subprocess.STDOUT)

cmd_11 = shlex.split( oneIterationScript + ' ' + rewOpts( 2011, numberOfIterations) + plotOpt + writeOpt )
#print 'Iteration number %s. Reweighting mc2011 with options:'%numberOfIterations
#for o in cmd_11[2:]: print o
_info( '2011', itNum, cmd_11, 'rew' )
rew_11 = subprocess.call(cmd_11, stdin = None, stdout = None, stderr = subprocess.STDOUT)

cmd_12 = shlex.split( oneIterationScript + ' ' + rewOpts( 2012, numberOfIterations) + plotOpt + combMomOpt + writeOpt )
# print 'Iteration number %s. Reweighting mc2012 with options:'%numberOfIterations
# for o in cmd_12[2:]: print o
_info( '2012', itNum, cmd_12, 'rew' )
rew_12 = subprocess.call(cmd_12, stdin = None, stdout = None, stderr = subprocess.STDOUT)

cmdfit = shlex.split( fittingScript + ' ' + fitOpts(numberOfIterations) )
# print 'Iteration number %s. Running 3fb Fit script with options:'%numberOfIterations 
# for o in cmdfit[2:]: print o
_info( '', itNum, cmdfit, 'fit' )
fit3fb = subprocess.call(cmdfit, stdin = None, stdout = None, stderr = subprocess.STDOUT)


assert False




testScrpt = 'python /project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/_junk/subproc/sleep.py'

cmd = shlex.split( oneIterationScript + ' -o%s -c%s'%('True', 'True')  )
p1 = subprocess.call(cmd, stdin = None, stdout = None, stderr = subprocess.STDOUT)
#p1.stdin.close()
#p1.wait()
processes.append([p1, True])


# cmd = shlex.split( testScrpt + ' -s%d'%10 )
# p2 = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE,
#                      stderr = subprocess.STDOUT)
# p2.stdin.close()
# processes.append([p2, True])


# for n in range( 1, numberOfIterations+1 ):
#     cmd = shlex.split( testScrpt + ' %d'%5 )
#     p = subprocess.Popen(cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE,
#                          stderr = subprocess.STDOUT)
#     p.stdin.close()
#     processes.append([p, True])

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
