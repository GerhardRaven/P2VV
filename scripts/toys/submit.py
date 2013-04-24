import os, subprocess
import bz2, select
import sys, atexit

# Save the current working directory and change back to it on exit
location = os.path.realpath(os.curdir)
atexit.register(os.chdir, location)

# change directory to the git repository
os.chdir('/bfys/raaij/p2vv/code')

# Open bz2 file
cmd = 'git archive --format=tar master | bzip2 > snapshot.tar.bz2'
op = subprocess.Popen(cmd, shell = True)

r = op.poll()
if r == None:
    r = op.wait()

if r:
    print "Error from git command %d" % r
    print op.stderr.read()
    sys.exit(r)

# Create the job
j = Job()
j.application = Executable()
j.application.exe = 'python'

# Create the right environment
env = {}
root_location = '/project/bfys/raaij/sw/root-5.34-patches'
env_vars = {'PATH' : 'bin',
            'PYTHONPATH' : 'lib',
            'LD_LIBRARY_PATH' : 'lib',
            'ROOTSYS' : ''}
for var, d in env_vars.iteritems():
    val = os.environ[var]
    s = val.split(':')
    found = False
    for k in s:
        if k.startswith(root_location):
            found = True
            break
    if not d:
        env[var] = root_location
    elif found:
        env[var] = val
    else:
        env[var] = os.path.join(root_location, d) + ':' + val
j.application.env = env

# Add the inputsandbox
j.inputsandbox = ['/project/bfys/raaij/p2vv/code/python/ToyMCUtils.py',
                  '/project/bfys/raaij/p2vv/code/toys/SFit_toy.py',
                  '/project/bfys/raaij/p2vv/code/toys/CFit_unbiased.pars',
                  '/project/bfys/raaij/p2vv/code/lib/libP2VV.so',
                  '/project/bfys/raaij/p2vv/code/snapshot.tar.bz2']

# Add the outputsandbox
j.outputsandbox = ['*.root']

# The merger
j.merger = CustomMerger(
    files = ['toy.root'],
    module = '/project/bfys/raaij/cmtuser/Ganga_v505r9/Utils/MergeDataSets.py'
    )

# Add the splitter
args = ['SFit_toy.py', '--ncpu=1', '-n',
        '10', '--nevents=20000', '-s', 'snapshot.tar.bz2']
j.splitter = GenericSplitter(
    attribute = 'application.args',
    values = [args for i in range(100)]
    )
j.name = 'SFit_toys'

# backend
j.backend = PBS(queue = 'stbcq')

# change back to original location
os.chdir(location)