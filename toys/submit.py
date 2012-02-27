import os, subprocess
import bz2, select
import sys, atexit

# Save the current working directory and change back to it on exit
location = os.path.realpath(os.curdir)
atexit.register(os.chdir, location)

# change directory to the git repository
os.chdir('/project/bfys/raaij/p2vv/code')

# Create subprocess and 
cmd = ['git', 'archive', '--format=tar', 'HEAD']
# Run the lines script to get the available Hlt lines
p = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                     stderr = subprocess.PIPE)

# Open bz2 file
snapshot_file = bz2.BZ2File(os.path.join(location, 'snapshot.tar.bz2'), 'w')

dy = '0.05'

# Read the git output and put it into the bz2 file
data = None
while True:
    ready = select.select([p.stdout], [], [])
    if ready[0]:
        data = ready[0][0].read(1024)
    else:
        data = None
    if data:
        snapshot_file.write(data)
    else:
        break

# Check the result of the git command
r = p.poll()
if r:
    print "Error from git command %d" % r
    print p.stderr.read()
    sys.exit(r)

# Create the job
j = Job()
j.application = Executable()
j.application.exe = 'python'

# Create the right environment
env = {}
root_location = '/project/bfys/jleerdam/local/root/5.30.05'
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
                  '/project/bfys/raaij/p2vv/code/toys/toyMC_BinnedPdf.py',
                  '/project/bfys/raaij/p2vv/code/lib/libP2VV.so',
                  os.path.join(location, 'snapshot.tar.bz2')]

# Add the outputsandbox
j.outputsandbox = ['*.root']

# The merger
j.merger = CustomMerger(
    files = ['toy.root'],
    module = '/project/bfys/raaij/cmtuser/Ganga_v505r9/Utils/MergeDataSets.py'
    )

# Add the splitter
args = ['toyMC_BinnedPdf.py', '--ncpu=1', '-n',
        '50', '-s', 'snapshot.tar.bz2', dy]
j.splitter = GenericSplitter(
    attribute = 'application.args',
    values = [args for i in range(100)]
    )
j.name = 'BinnedPdf_toy_vbins_%s' % dy.replace('.', '_')

# backend
j.backend = PBS(queue = 'stbcq')

# change back to original location
os.chdir(location)
