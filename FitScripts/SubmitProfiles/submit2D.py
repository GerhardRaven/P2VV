import os

acceptance = 0

if acceptance:
    fitname = 'cfit_acc'
else:
    fitname = 'cfit_noacc'

scriptname = 'Profile2D.py'
varname1 = 'dGamma'
varname2 = 'phi_s'

npoints = 40

#for xstep in range(0,20):
for xstep in range(20,npoints):
    # Create the job
    j = Job()
    j.application = Executable()
    j.application.exe = 'python'

    # Create the right environment
    env = {}
    root_location = '/project/bfys/jleerdam/local/root/5.32.00-patches_P2VV'
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
    j.inputsandbox = ['/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/%s'%(scriptname)
                      ,'/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/%s_result'%(fitname)
                      ,'/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/%sparams.txt'%(fitname)
                      ]

    # Add the outputsandbox
    j.outputsandbox = ['*.txt']


    # Add the splitter
    list = []
    for ystep in range(0,npoints):
        list.append([scriptname, '-a', acceptance, '-n', str(npoints), '-x', str(xstep), '-y',str(ystep),'-v',varname1,'-w',varname2,])

    s = ArgSplitter(args=list)
    j.splitter = s

    j.name = '2DProfiles_inversion_%s_%s'%(str(xstep),fitname)

    # backend
    #j.backend = PBS(queue = 'stbcq')
    #j.backend = PBS(queue = 'qlong')
    j.backend = PBS()

    j.submit()
