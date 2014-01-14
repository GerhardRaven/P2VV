#!/usr/bin/env python
import subprocess
import shlex
#from time import sleep
import select

# configuration
numberOfIterations = 4
oneIterationScript = 'python /project/bfys/vsyropou/Repository/p2vv/scripts/iterativeAngAcc.py'
fittingScript      = 'python /project/bfys/vsyropou/Repository/p2vv/scripts/vsFit.py'
fiData             = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
parameterEstimates = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/nominalFitResults/20112012Reco14DataFitValues_6KKMassBins.par'
startingAngAccFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
correctedAngAccBaseName = 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_norm_'


testScrpt = 'python /project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/_junk/subproc/sleep.py'

# begin
print 'begin'
processes = []
opts = {
    1 : dict( rew11 = '-n%i -s%s -a%s -d%s -f%s'%(1,'2011',startingAngAccFile,parameterEstimates,'False'),
              rew12 = '-n%i -s%s -a%s -d%s -f%s -c%s'%(1,'2011',startingAngAccFile,parameterEstimates,'False','True'),
              dofit = '-d%s -a%s -i%s -o%s'%(fidData,correctedAngAccBaseName + '1',parameterEstimates, '20112012Reco14DataFitValues_6KKMassBins_1.par' )
              )
    2 : dict(
        )
    }

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
