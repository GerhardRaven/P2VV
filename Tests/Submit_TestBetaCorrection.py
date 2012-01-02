from math import pi

g = Executable()
g.exe = 'python'
g.args = []
g.env = {"ROOTSYS":"/project/bfys/jleerdam/local/root/5.30.02",
         "PYTHONPATH":"/user/dveijk/LHCb/P2VV/FreshStart/p2vv/python:/project/bfys/jleerdam/local/root/5.30.02/lib:/cvmfs/lhcb.cern.ch/lib/lhcb/LBSCRIPTS/LBSCRIPTS_v6r6p3/InstallArea/python",
         "LD_LIBRARY_PATH":"/user/dveijk/LHCb/P2VV/FreshStart/p2vv/lib:/project/bfys/jleerdam/local/root/5.30.02/lib:/cvmfs/lhcb.cern.ch/lib/lcg/external/gcc/4.3.5/x86_64-slc5/lib64:/cvmfs/lhcb.cern.ch/lib/lhcb/COMPAT/COMPAT_v1r8/CompatSys/x86_64-slc5-gcc43-opt/lib:/cvmfs/lhcb.cern.ch/lib/lcg/external/Python/2.6.5/x86_64-slc5-gcc43-opt/lib",
         "PATH":"/project/bfys/jleerdam/local/root/5.30.02/bin:/cvmfs/lhcb.cern.ch/lib/lhcb/LBSCRIPTS/LBSCRIPTS_v6r6p3/InstallArea/scripts:/cvmfs/lhcb.cern.ch/lib/lcg/external/Python/2.6.5/x86_64-slc5-gcc43-opt/bin:/cvmfs/lhcb.cern.ch/lib/lcg/external/gcc/4.3.5/x86_64-slc5/bin:/cvmfs/lhcb.cern.ch/lib/lhcb/COMPAT/COMPAT_v1r8/CompatSys/x86_64-slc5-gcc43-opt/bin:/cvmfs/lhcb.cern.ch/lib/contrib/CMT/v1r20p20090520/Linux-x86_64:/usr/lib64/qt-3.3/bin:/usr/kerberos/bin:/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin"
         }
#g.project = 'DaVinci'
#g.version = 'v2r2'
#g.user_release_area = '/user/dveijk/LHCb/P2VV/Generate'


###############################
# Make loop for arguments
###############################
seedmax = 2
list = []

for j in range(0,seedmax):
    list.append(['TestBetaCorrection.py','-s',str(j+1)])

print 'list =',list
J = Job(name = 'TestBetaCorrection_Merge', application = g, backend=PBS())
s = ArgSplitter(args=list)
J.splitter = s
J.outputdata = []
J.outputsandbox = ['*.root']
J.inputsandbox = [ File (
    name = '/user/dveijk/LHCb/P2VV/FreshStart/p2vv/Tests/TestBetaCorrection.py',
    subdir = '.'
    ),]
J.merger = CustomMerger(
    files = ['Toy.root'],
    ignorefailed = False,
    overwrite = False,
    module = File(
    name = '/project/bfys/raaij/cmtuser/Ganga_v505r9/Utils/MergeDataSets.py' ,
    subdir = '.'
    )
    )
J.submit()
