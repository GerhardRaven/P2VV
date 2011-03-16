from ROOT import *
import os.path

################
### Settings ###
################

RDSnameIn = 'profile.root'
RDSnameOut = 'MergedRDS'
nameinRDS = 'dlogLLDataSet'

job = 57
npoints = 40

####################################
### open the first RDS, clone it ###
####################################
tfile = TFile('/data/bfys/dveijk/gangadir/workspace/dveijk/LocalXML/'+str(job)+'/'+str(0)+'/output/'+RDSnameIn)
RDS = tfile.Get(nameinRDS)
tfile.Close()
RDSMerge = RDS.Clone()
    
for j in range(1,npoints):
    filename = '/data/bfys/dveijk/gangadir/workspace/dveijk/LocalXML/'+str(job)+'/'+str(j)+'/output/'+RDSnameIn
    if os.path.exists(filename):
        tfile = TFile(filename)
        RDS = tfile.Get(nameinRDS)
        tfile.Close()
        RDSMerge.append(RDS)
    else:
        print "File %s does not exist"%(filename)
        continue

RDSMerge.Print()
    
writefile = TFile.Open(RDSnameOut+'.root',"RECREATE")
RDSMerge.Write(nameinRDS)
writefile.Close()
